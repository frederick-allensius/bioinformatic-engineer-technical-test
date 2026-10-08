"""Unit tests for variant calling from aligned strings (app.services.variants)."""

import pytest

from app.services.variants import Variant, call_variants


def calls(ref_aln: str, query_aln: str) -> list[tuple]:
    return [(v.position, v.ref, v.alt, v.type) for v in call_variants(ref_aln, query_aln)]


def test_no_variants():
    assert calls("ACGT", "ACGT") == []


def test_snv_position_is_one_based():
    assert calls("ACGT", "AGGT") == [(2, "C", "G", "SNV")]


def test_n_in_query_is_not_an_snv():
    assert calls("ACGT", "ANGN") == []


def test_n_in_reference_is_not_an_snv():
    assert calls("ANGT", "ACGT") == []


def test_adjacent_snvs_are_reported_separately():
    assert calls("ACGT", "AGCT") == [(2, "C", "G", "SNV"), (3, "G", "C", "SNV")]


def test_deletion_groups_consecutive_bases_and_uses_first_deleted_position():
    assert calls("ACGTAC", "A---AC") == [(2, "CGT", "-", "deletion")]


def test_deletion_at_start():
    assert calls("ACGT", "--GT") == [(1, "AC", "-", "deletion")]


def test_deletion_at_end():
    assert calls("ACGT", "AC--") == [(3, "GT", "-", "deletion")]


def test_insertion_position_is_preceding_reference_base():
    assert calls("AC--GT", "ACTTGT") == [(2, "-", "TT", "insertion")]


def test_insertion_before_first_base_is_position_zero():
    assert calls("--ACGT", "GGACGT") == [(0, "-", "GG", "insertion")]


def test_insertion_after_last_base():
    assert calls("ACGT--", "ACGTAA") == [(4, "-", "AA", "insertion")]


def test_separate_gap_runs_are_separate_events():
    assert calls("ACGTA", "A-G-A") == [(2, "C", "-", "deletion"), (4, "T", "-", "deletion")]


def test_inserted_n_is_kept_in_alt():
    assert calls("AC-GT", "ACNGT") == [(2, "-", "N", "insertion")]


def test_mixed_events_are_sorted_by_position():
    #           pos: 1 2 3 . . 4 5 6 7 8
    ref_aln = "ACG--TACGA"
    qry_aln = "AGGCC-ACTA"
    assert calls(ref_aln, qry_aln) == [
        (2, "C", "G", "SNV"),
        (3, "-", "CC", "insertion"),
        (4, "T", "-", "deletion"),
        (7, "G", "T", "SNV"),
    ]


def test_deletion_then_insertion_at_same_locus_keeps_alignment_order():
    # Deletion of base 2 followed directly by an insertion after base 2.
    assert calls("AC-G", "A-TG") == [(2, "C", "-", "deletion"), (2, "-", "T", "insertion")]


def test_coordinates_account_for_earlier_indels():
    # After a 2 bp insertion and a 1 bp deletion, the SNV is still at reference base 6.
    ref_aln = "A--CGTACG"
    qry_aln = "ATTC-TAGG"
    assert calls(ref_aln, qry_aln) == [
        (1, "-", "TT", "insertion"),
        (3, "G", "-", "deletion"),
        (6, "C", "G", "SNV"),
    ]


def test_unequal_lengths_rejected():
    with pytest.raises(ValueError):
        call_variants("ACG", "AC")


def test_variant_length():
    assert Variant(5, "-", "ATC", "insertion").length == 3
    assert Variant(5, "AC", "-", "deletion").length == 2
    assert Variant(5, "A", "C", "SNV").length == 1
