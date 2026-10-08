"""Unit tests for variant effect annotation (app.services.annotation)."""

import pytest

from app.services.alignment import Alignment, global_align
from app.services.annotation import annotate_variants
from app.services.variants import Variant, call_variants


def annotate(reference: str, query: str):
    """Align with the real aligner, then annotate."""
    aln = global_align(reference, query)
    return annotate_variants(
        reference, aln, call_variants(aln.reference_aligned, aln.query_aligned)
    )


def annotate_aligned(ref_aln: str, query_aln: str):
    """Annotate a hand-written alignment (to control indel placement exactly)."""
    aln = Alignment(ref_aln, query_aln, score=0)
    return annotate_variants(ref_aln.replace("-", ""), aln, call_variants(ref_aln, query_aln))


def only(results):
    assert len(results) == 1, results
    return results[0]


# --------------------------------------------------------------------------- SNVs

#        codon:  1   2   3   4   5
REF = "ATGGCCAAATGGTAA"  # M A K W *


def test_synonymous():
    a = only(annotate(REF, "ATGGCGAAATGGTAA"))  # GCC > GCG, both Ala
    assert (a.effect, a.codon, a.ref_codon, a.alt_codon) == ("synonymous", 2, "GCC", "GCG")
    assert a.protein_change == "p.A2="


def test_missense():
    a = only(annotate(REF, "ATGGCCGAATGGTAA"))  # AAA (K) > GAA (E)
    assert (a.effect, a.ref_aa, a.alt_aa, a.protein_change) == ("missense", "K", "E", "p.K3E")


def test_nonsense():
    a = only(annotate(REF, "ATGGCCAAATGATAA"))  # TGG (W) > TGA (stop)
    assert (a.effect, a.protein_change) == ("nonsense", "p.W4*")


def test_stop_lost():
    a = only(annotate(REF, "ATGGCCAAATGGTAC"))  # TAA > TAC (Tyr)
    assert (a.effect, a.protein_change) == ("stop_lost", "p.*5Yext*?")


def test_stop_to_other_stop_is_synonymous():
    a = only(annotate(REF, "ATGGCCAAATGGTAG"))
    assert (a.effect, a.protein_change) == ("synonymous", "p.*5=")


def test_start_codon_change_is_flagged():
    a = only(annotate(REF, "ACGGCCAAATGGTAA"))
    assert a.effect == "missense"
    assert any("start" in note for note in a.notes)


def test_snvs_in_same_codon_are_evaluated_together():
    # AGC (Ser) > TCC (Ser): each SNV alone would be missense (TGC Cys, ACC Thr),
    # but the combined codon is synonymous.
    ref = "ATGAGCTAA"
    results = annotate(ref, "ATGTCCTAA")
    assert [r.variant.position for r in results] == [4, 5]
    for r in results:
        assert (r.effect, r.alt_codon) == ("synonymous", "TCC")
        assert any("same codon" in note for note in r.notes)


def test_snv_in_codon_with_n_is_undetermined():
    a = only(annotate(REF, "ATGGCCANTTGGTAA"))  # codon 3 = ANT
    assert a.variant.position == 9
    assert (a.effect, a.alt_codon, a.protein_change) == ("undetermined", "ANT", None)
    assert any("N" in note for note in a.notes)


def test_snv_in_codon_hit_by_deletion_is_undetermined():
    #            codon 3 = "-AT": base 7 deleted, base 9 A>T
    results = annotate_aligned("ATGGCCAAAGGG", "ATGGCC-ATGGG")
    snv = next(r for r in results if r.variant.type == "SNV")
    assert snv.effect == "undetermined"
    assert any("insertion or deletion" in note for note in snv.notes)


def test_snv_in_codon_interrupted_by_insertion_is_undetermined():
    results = annotate_aligned("ATGGC-CAAA", "ATGGCTCTAA")  # ins after 5, SNV at 7 (A>T)
    snv = next(r for r in results if r.variant.type == "SNV")
    assert snv.variant.position == 7
    assert snv.effect != "undetermined"  # codon 3 (7-9) is untouched by the insertion
    other = annotate_aligned("ATGGC-CAAA", "ATGGTTCAAA")  # SNV at 5 + ins after 5
    snv = next(r for r in other if r.variant.type == "SNV")
    assert snv.effect == "undetermined"


def test_snv_outside_last_complete_codon():
    a = only(annotate("ATGGCCA", "ATGGCCG"))
    assert (a.effect, a.codon) == ("undetermined", None)


# --------------------------------------------------------------------------- indels


def test_codon_aligned_inframe_deletion():
    a = only(annotate(REF, "ATGGCCTGGTAA"))  # AAA (K3) deleted
    assert (a.effect, a.codon, a.protein_change) == ("inframe_deletion", 3, "p.K3del")


def test_inframe_deletion_spanning_codons_is_delins():
    a = only(annotate_aligned("ATGGCCAAATGGTAA", "ATGG---AATGGTAA"))
    # ATG G|CC A|AA -> ATG GAA: M A K W * becomes M E W *
    assert (a.effect, a.protein_change) == ("inframe_deletion", "p.A2_K3delinsE")


def test_inframe_insertion_between_codons():
    a = only(annotate_aligned("ATGGCC---AAATGGTAA", "ATGGCCCTGAAATGGTAA"))
    assert (a.effect, a.codon, a.protein_change) == ("inframe_insertion", 3, "p.A2_K3insL")


def test_inframe_insertion_creating_stop_is_flagged():
    a = only(annotate_aligned("ATG---GCCAAATGGTAA", "ATGTAAGCCAAATGGTAA"))
    assert a.effect == "inframe_insertion"
    assert a.protein_change == "p.M1_A2ins*"
    assert any("premature stop" in note for note in a.notes)


def test_frameshift_deletion():
    a = only(annotate_aligned("ATGGCCAAATGGTAA", "ATGG-CAAATGGTAA"))
    # ATG GCA AAT GGT AA -> M A N G (no stop in the new frame)
    assert a.effect == "frameshift"
    assert a.codon == 2
    assert a.protein_change == "p.K3Nfs*?"


def test_frameshift_with_new_stop_counts_residues():
    # ref ATG GCC AAA TGG TAA GTG A  -> M A K W *
    # del A at 9: ATG GCC AAT GGT AAG TGA -> M A N G K *
    # first change K3N; counting N as 1, the new stop is residue 4.
    ref = "ATGGCCAAATGGTAAGTGA"
    a = only(annotate_aligned(ref, "ATGGCCAA-TGGTAAGTGA"))
    assert (a.ref_aa, a.alt_aa, a.protein_change) == ("K", "N", "p.K3Nfs*4")


def test_frameshift_creating_immediate_stop():
    ref = "ATGTCAAGGTAA"  # M S R *
    a = only(annotate_aligned(ref, "ATGT-AAGGTAA"))  # ATG TAA GGT AA
    assert (a.effect, a.protein_change) == ("frameshift", "p.S2*")


def test_frameshift_insertion():
    a = only(annotate_aligned("ATGGCC-AAATGGTAA", "ATGGCCTAAATGGTAA"))
    assert a.effect == "frameshift"
    assert a.protein_change == "p.K3*"


def test_indel_after_stop_codon_does_not_change_protein():
    a = only(annotate_aligned("ATGTAAGGGCCC", "ATGTAAGGG-CC"))
    assert (a.effect, a.protein_change) == ("frameshift", "p.=")


# --------------------------------------------------------------------------- context notes


def test_variant_downstream_of_frameshift_is_flagged():
    ref = "ATGGCCAAATGGCGTTAA"
    results = annotate_aligned(ref, "ATGG-CAAATGGCATTAA")  # del at 5, SNV at 14
    assert [r.variant.type for r in results] == ["deletion", "SNV"]
    assert any("Downstream of the frameshift at position 5" in n for n in results[1].notes)


def test_compensating_frameshift_restores_frame():
    # -1 deletion at 5, then +1 insertion after 12: the frame is restored.
    results = annotate_aligned("ATGGCCAAATGG-CGTTAA", "ATGG-CAAATGGACGTTAA")
    assert [r.variant.type for r in results] == ["deletion", "insertion"]
    assert [r.effect for r in results] == ["frameshift", "frameshift"]
    assert any("Restores the reading frame" in n for n in results[1].notes)


def test_variant_downstream_of_premature_stop_is_flagged():
    ref = "ATGGCCTGGAAACGTTAA"
    results = annotate(ref, "ATGGCCTGAAAACGCTAA")  # W3* then CGT>CGC
    assert results[0].effect == "nonsense"
    assert any("premature stop at codon 3" in n for n in results[1].notes)


@pytest.mark.parametrize("vtype", ["insertion", "deletion"])
def test_indel_beyond_last_codon(vtype):
    ref = "ATGGCCA"
    variant = (
        Variant(7, "-", "T", "insertion")
        if vtype == "insertion"
        else Variant(7, "A", "-", "deletion")
    )
    query_aln = "ATGGCCAT" if vtype == "insertion" else "ATGGCC-"
    ref_aln = "ATGGCCA-" if vtype == "insertion" else ref
    a = only(annotate_variants(ref, Alignment(ref_aln, query_aln, 0), [variant]))
    assert a.codon is None
    assert a.effect == "frameshift"
