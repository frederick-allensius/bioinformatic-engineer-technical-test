"""Unit tests for translation (app.services.translation)."""

import pytest
from hypothesis import given, settings
from hypothesis import strategies as st

from app.services.translation import (
    CODON_TABLE,
    STOP_CODONS,
    InvalidFrameError,
    codons,
    translate,
    translate_codon,
    translate_full,
)

# --------------------------------------------------------------------------- codon table


def test_codon_table_has_all_64_codons():
    assert len(CODON_TABLE) == 64


def test_stop_codons():
    assert {"TAA", "TAG", "TGA"} == STOP_CODONS


@pytest.mark.parametrize(
    ("codon", "aa"),
    [("ATG", "M"), ("TGG", "W"), ("TTT", "F"), ("GGC", "G"), ("AAA", "K"), ("TGT", "C")],
)
def test_known_codons(codon, aa):
    assert CODON_TABLE[codon] == aa


def test_codon_table_matches_biopython_standard_table():
    codon_data = pytest.importorskip("Bio.Data.CodonTable")
    standard = codon_data.unambiguous_dna_by_id[1]
    expected = dict(standard.forward_table) | dict.fromkeys(standard.stop_codons, "*")
    assert expected == CODON_TABLE


def test_codon_with_n_is_x_even_when_inferable():
    # GCN is alanine whatever N is, but the rules require X.
    assert translate_codon("GCN") == "X"
    assert translate_codon("NNN") == "X"


def test_translate_codon_rejects_wrong_length():
    with pytest.raises(ValueError):
        translate_codon("AT")


# --------------------------------------------------------------------------- frames


def test_codons_per_frame():
    seq = "ATGGCCTAA"
    assert codons(seq, 1) == ["ATG", "GCC", "TAA"]
    assert codons(seq, 2) == ["TGG", "CCT"]
    assert codons(seq, 3) == ["GGC", "CTA"]


@pytest.mark.parametrize("frame", [0, 4, -1, 10])
def test_invalid_frame(frame):
    with pytest.raises(InvalidFrameError):
        translate("ATG", frame)


def test_each_frame_translates_from_its_offset():
    seq = "AATGGCCAAAGG"
    assert translate(seq, 1).protein == "NGQR"  # AAT GGC CAA AGG
    assert translate(seq, 2).protein == "MAK"  # ATG GCC AAA (GG ignored)
    assert translate(seq, 3).protein == "WPK"  # TGG CCA AAG (G ignored)


# --------------------------------------------------------------------------- stops


def test_stops_before_first_stop_codon():
    result = translate("ATGGCCTAAGGG")
    assert result.protein == "MA"
    assert result.has_stop_codon is True
    assert result.length == 2


@pytest.mark.parametrize("stop", ["TAA", "TAG", "TGA"])
def test_each_stop_codon(stop):
    result = translate("ATG" + stop + "ATG")
    assert (result.protein, result.has_stop_codon) == ("M", True)


def test_stop_as_first_codon_gives_empty_protein_with_stop():
    result = translate("TAAATG")
    assert (result.protein, result.length, result.has_stop_codon) == ("", 0, True)


def test_no_stop_codon():
    result = translate("ATGGCC")
    assert (result.protein, result.has_stop_codon) == ("MA", False)


def test_stop_out_of_frame_is_ignored():
    # TAA appears at offset 1, but frame 1 reads ATA AGC.
    assert translate("ATAAGC", 1).protein == "IS"


def test_ambiguous_stop_is_x_not_stop():
    # TNA could be TAA (stop) or TCA (Ser): not an *unambiguous* stop.
    result = translate("ATGTNAGCC")
    assert (result.protein, result.has_stop_codon) == ("MXA", False)


def test_n_codons_are_counted_in_length():
    result = translate("ATGNNNGCN")
    assert result.protein == "MXX"
    assert result.length == 3


# --------------------------------------------------------------------------- incomplete codons


@pytest.mark.parametrize(
    ("seq", "frame", "protein"),
    [
        ("ATGGC", 1, "M"),  # 2 trailing bases ignored
        ("ATGG", 1, "M"),  # 1 trailing base ignored
        ("AATGGC", 2, "M"),  # ATG + GC ignored
    ],
)
def test_trailing_partial_codon_is_ignored(seq, frame, protein):
    assert translate(seq, frame).protein == protein


@pytest.mark.parametrize(("seq", "frame"), [("A", 1), ("AT", 1), ("ATG", 2), ("ATGG", 3)])
def test_no_complete_codon(seq, frame):
    result = translate(seq, frame)
    assert (result.protein, result.length, result.has_stop_codon) == ("", 0, False)


def test_translate_full_keeps_stops():
    assert translate_full("ATGTAAGGG") == "M*G"


# --------------------------------------------------------------------------- vs Biopython


@pytest.mark.slow
@settings(max_examples=200, deadline=None)
@given(seq=st.text(alphabet="ACGT", min_size=0, max_size=120), frame=st.sampled_from([1, 2, 3]))
def test_matches_biopython(seq, frame):
    seq_module = pytest.importorskip("Bio.Seq")
    framed = seq[frame - 1 :]
    framed = framed[: len(framed) - len(framed) % 3]
    expected = str(seq_module.Seq(framed).translate(to_stop=True)) if framed else ""
    assert translate(seq, frame).protein == expected
