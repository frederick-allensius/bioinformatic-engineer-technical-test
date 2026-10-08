"""Unit tests for the Needleman-Wunsch/Gotoh implementation (app.services.alignment).

Besides hand-made cases, two independent oracles are used:

* brute-force enumeration of *every* global alignment for tiny inputs, which checks
  both optimality and the "fewest gap openings" tie-break;
* Biopython's PairwiseAligner (test-only dependency) for larger random inputs.
"""

from functools import cache

import pytest
from hypothesis import given, settings
from hypothesis import strategies as st

from app.services.alignment import (
    DEFAULT_SCORING,
    ScoringScheme,
    global_align,
    identity,
    score_alignment,
)

dna = st.text(alphabet="ACGT", min_size=1, max_size=60)
dna_with_n = st.text(alphabet="ACGTN", min_size=1, max_size=60)
tiny_dna = st.text(alphabet="ACGN", min_size=0, max_size=5)


def gap_runs(aligned_ref: str, aligned_query: str) -> int:
    """Number of separate gap runs (gap openings) in an alignment."""
    runs, prev = 0, None
    for r, q in zip(aligned_ref, aligned_query, strict=True):
        state = "D" if q == "-" else "I" if r == "-" else "M"
        if state != "M" and state != prev:
            runs += 1
        prev = state
    return runs


def all_alignments(a: str, b: str) -> list[tuple[str, str]]:
    @cache
    def rec(i: int, j: int) -> tuple[tuple[str, str], ...]:
        if i == 0 and j == 0:
            return (("", ""),)
        out = []
        if i and j:
            out += [(x + a[i - 1], y + b[j - 1]) for x, y in rec(i - 1, j - 1)]
        if i:
            out += [(x + a[i - 1], y + "-") for x, y in rec(i - 1, j)]
        if j:
            out += [(x + "-", y + b[j - 1]) for x, y in rec(i, j - 1)]
        return tuple(out)

    return list(rec(len(a), len(b)))


# --------------------------------------------------------------------------- basics


def test_identical_sequences():
    aln = global_align("ACGTACGT", "ACGTACGT")
    assert aln.reference_aligned == aln.query_aligned == "ACGTACGT"
    assert aln.score == 8
    assert aln.identity == 1.0


def test_single_mismatch_prefers_substitution_over_gaps():
    aln = global_align("ACGTACGT", "ACGAACGT")
    assert aln.query_aligned == "ACGAACGT"
    assert aln.score == 7 - 1


def test_n_scores_zero_and_is_not_identical():
    aln = global_align("ACGT", "ANGN")
    assert aln.query_aligned == "ANGN"
    assert aln.score == 2  # A/A + G/G; both N columns score 0
    assert aln.identity == 0.5


def test_n_against_n_is_not_a_match():
    assert identity("ANG", "ANG") == pytest.approx(2 / 3)
    assert score_alignment("ANG", "ANG") == 2


def test_identity_counts_gap_columns_in_denominator():
    assert identity("ACGT", "AC-T") == 0.75


def test_score_alignment_rejects_gap_gap_column():
    with pytest.raises(ValueError):
        score_alignment("A-", "A-")


def test_long_n_run_stays_on_the_diagonal():
    ref = "ACGTACGTACGTACGTACGT"
    query = "ACGT" + "N" * 12 + "ACGT"
    aln = global_align(ref, query)
    assert "-" not in aln.reference_aligned + aln.query_aligned
    assert aln.score == 8


def test_custom_scoring_scheme():
    scoring = ScoringScheme(match=2, mismatch=-3, gap=-1, n=0)
    aln = global_align("AAA", "AGA", scoring)
    # One mismatch costs -3, but two gap columns cost only -2.
    assert aln.score == 2 + 2 - 2
    assert "-" in aln.query_aligned


def test_empty_inputs():
    assert global_align("", "").score == 0
    aln = global_align("ACG", "")
    assert (aln.reference_aligned, aln.query_aligned) == ("ACG", "---")
    aln = global_align("", "AC")
    assert (aln.reference_aligned, aln.query_aligned) == ("--", "AC")


# --------------------------------------------------------------------------- indels


def test_deletion_is_one_contiguous_gap():
    # ref GGAACCCTG vs query GGCCTG: "GG--C-CTG" scores the same as "GG---CCTG"
    # under a linear penalty; the tie-break must keep the 3 bp deletion together.
    aln = global_align("GGAACCCTG", "GGCCTG")
    assert aln.query_aligned == "GG---CCTG"
    assert gap_runs(aln.reference_aligned, aln.query_aligned) == 1


def test_insertion_is_one_contiguous_gap():
    aln = global_align("GGCCTG", "GGAACCCTG")
    assert aln.reference_aligned == "GG---CCTG"


@pytest.mark.parametrize(
    ("ref", "query", "expected_query_aligned"),
    [
        ("GAAAT", "GAAT", "G-AAT"),  # homopolymer deletion -> leftmost A
        ("GCACACAT", "GCACAT", "G--CACAT"),  # dinucleotide repeat -> leftmost "CA"
    ],
)
def test_deletions_in_repeats_are_left_aligned(ref, query, expected_query_aligned):
    assert global_align(ref, query).query_aligned == expected_query_aligned


def test_insertion_in_homopolymer_is_left_aligned():
    aln = global_align("GAAT", "GAAAT")
    assert aln.reference_aligned == "G-AAT"


def test_terminal_gaps_are_allowed():
    assert global_align("ACGTACGT", "ACGTAC").query_aligned == "ACGTAC--"
    assert global_align("ACGT", "TTACGT").reference_aligned == "--ACGT"


def test_alignment_is_deterministic():
    results = {
        (a.reference_aligned, a.query_aligned)
        for a in (global_align("ACGTTGCA", "AGTTTGA") for _ in range(5))
    }
    assert len(results) == 1


# --------------------------------------------------------------------------- properties


@given(ref=dna_with_n, query=dna_with_n)
def test_alignment_preserves_both_sequences(ref, query):
    aln = global_align(ref, query)
    assert aln.reference_aligned.replace("-", "") == ref
    assert aln.query_aligned.replace("-", "") == query
    columns = zip(aln.reference_aligned, aln.query_aligned, strict=True)
    assert all(not (r == q == "-") for r, q in columns)


@given(ref=dna_with_n, query=dna_with_n)
def test_reported_score_matches_columns(ref, query):
    aln = global_align(ref, query)
    assert aln.score == score_alignment(aln.reference_aligned, aln.query_aligned)
    assert 0.0 <= aln.identity <= 1.0


@pytest.mark.slow
@settings(max_examples=300, deadline=None)
@given(ref=tiny_dna, query=tiny_dna)
def test_optimal_and_fewest_gap_openings_vs_brute_force(ref, query):
    candidates = all_alignments(ref, query)
    best_score = max(score_alignment(r, q) for r, q in candidates)
    fewest_runs = min(gap_runs(r, q) for r, q in candidates if score_alignment(r, q) == best_score)

    aln = global_align(ref, query)
    assert aln.score == best_score
    assert gap_runs(aln.reference_aligned, aln.query_aligned) == fewest_runs


@pytest.mark.slow
@settings(max_examples=150, deadline=None)
@given(ref=dna_with_n, query=dna_with_n)
def test_score_matches_biopython(ref, query):
    align = pytest.importorskip("Bio.Align")
    aligner = align.PairwiseAligner(
        mode="global",
        match_score=DEFAULT_SCORING.match,
        mismatch_score=DEFAULT_SCORING.mismatch,
        open_gap_score=DEFAULT_SCORING.gap,
        extend_gap_score=DEFAULT_SCORING.gap,
    )
    aligner.wildcard = "N"  # Biopython scores wildcard columns as 0, like us
    assert global_align(ref, query).score == aligner.score(ref, query)
