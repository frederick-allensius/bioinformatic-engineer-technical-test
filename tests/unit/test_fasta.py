"""Unit tests for FASTA parsing / validation (app.services.fasta)."""

import pytest

from app.services.fasta import (
    FastaValidationError,
    decode_fasta_bytes,
    load_fasta_file,
    parse_fasta,
)


def error_code(text: str, **kwargs) -> str:
    with pytest.raises(FastaValidationError) as exc_info:
        parse_fasta(text, **kwargs)
    return exc_info.value.code


# --------------------------------------------------------------------------- happy paths


def test_single_line_record():
    record = parse_fasta(">seq1\nACGT\n")
    assert record.sequence_id == "seq1"
    assert record.sequence == "ACGT"


def test_multiline_sequence_is_concatenated():
    record = parse_fasta(">seq1\nACGT\nTTGG\nCA\n")
    assert record.sequence == "ACGTTTGGCA"


def test_blank_lines_are_ignored_everywhere():
    record = parse_fasta("\n\n>seq1\n\nACGT\n\n\nTT\n\n")
    assert record.sequence == "ACGTTT"


def test_lowercase_and_mixed_case_are_uppercased():
    assert parse_fasta(">s\nacgtn\nAcGtN").sequence == "ACGTNACGTN"


def test_windows_and_old_mac_line_endings():
    assert parse_fasta(">s\r\nACGT\r\nGG\r\n").sequence == "ACGTGG"
    assert parse_fasta(">s\rACGT\rGG\r").sequence == "ACGTGG"


def test_trailing_whitespace_on_lines_is_stripped():
    # "Surrounding" whitespace (e.g. editor artefacts) is removed per line.
    assert parse_fasta(">s  \n  ACGT \t\nGG  ").sequence == "ACGTGG"


def test_no_trailing_newline():
    assert parse_fasta(">s\nACGT").sequence == "ACGT"


def test_sequence_id_is_first_token_of_header():
    record = parse_fasta(">REF_001 synthetic CDS, 300 bp\nACGT")
    assert record.sequence_id == "REF_001"
    assert record.description == "REF_001 synthetic CDS, 300 bp"


def test_whitespace_between_gt_and_identifier_is_tolerated():
    assert parse_fasta(">   patient_x  extra words\nACGT").sequence_id == "patient_x"


def test_tab_delimited_header():
    assert parse_fasta(">id\tdescription\nACGT").sequence_id == "id"


def test_n_is_valid_even_when_sequence_is_all_n():
    assert parse_fasta(">s\n" + "N" * 1000).sequence == "N" * 1000


def test_length_limit_is_inclusive():
    assert len(parse_fasta(">s\n" + "A" * 1000)) == 1000


def test_length_is_counted_after_removing_line_breaks():
    # 1000 nt wrapped at 60 columns is ~1016 characters of text but still valid.
    seq = "ACGTN" * 200
    wrapped = "\n".join(seq[i : i + 60] for i in range(0, len(seq), 60))
    assert len(parse_fasta(">s\n" + wrapped)) == 1000


def test_max_length_none_disables_limit():
    assert len(parse_fasta(">s\n" + "A" * 5000, max_length=None)) == 5000


def test_custom_max_length():
    assert error_code(">s\nACGTA", max_length=4) == "sequence_too_long"


# --------------------------------------------------------------------------- rejections


@pytest.mark.parametrize("text", ["", "   ", "\n\n", " \r\n \t "])
def test_empty_file(text):
    assert error_code(text) == "empty_file"


def test_missing_header():
    assert error_code("ACGTACGT\n") == "malformed_fasta"


def test_sequence_before_header():
    assert error_code("ACGT\n>s\nACGT") == "malformed_fasta"


@pytest.mark.parametrize("header", [">", ">   ", ">\t"])
def test_empty_header(header):
    assert error_code(f"{header}\nACGT") == "empty_header"


def test_header_only_has_empty_sequence():
    assert error_code(">seq1\n") == "empty_sequence"
    assert error_code(">seq1\n\n\n") == "empty_sequence"


def test_multiple_records():
    assert error_code(">a\nACGT\n>b\nGGCC\n") == "multiple_records"


def test_multiple_records_reports_count_and_line():
    with pytest.raises(FastaValidationError) as exc_info:
        parse_fasta(">a\nAC\n>b\nGG\n>c\nTT")
    assert "found 3" in exc_info.value.message
    assert "line 3" in exc_info.value.message


def test_two_headers_without_sequence_is_multiple_records():
    assert error_code(">a\n>b\nACGT") == "multiple_records"


@pytest.mark.parametrize("bad", ["U", "R", "Y", "X", "-", "*", "1", ".", ">", "é"])
def test_invalid_characters(bad):
    assert error_code(f">s\nACG{bad}T") == "invalid_character"


def test_invalid_character_message_points_at_location():
    with pytest.raises(FastaValidationError) as exc_info:
        parse_fasta(">s\nACGT\nACXT")
    assert "'X'" in exc_info.value.message
    assert "line 3, column 3" in exc_info.value.message


@pytest.mark.parametrize(
    "line", ["ACG T", "ACG\tT", "AC\N{NO-BREAK SPACE}GT", "AC\x0bGT", "AC\x0cGT"]
)
def test_embedded_whitespace(line):
    assert error_code(f">s\n{line}") == "embedded_whitespace"


def test_too_long():
    assert error_code(">s\n" + "A" * 1001) == "sequence_too_long"


def test_too_long_counts_across_lines():
    assert error_code(">s\n" + "A" * 600 + "\n" + "C" * 401) == "sequence_too_long"


# --------------------------------------------------------------------------- bytes / files


def test_decode_strips_utf8_bom():
    assert parse_fasta(decode_fasta_bytes(b"\xef\xbb\xbf>s\nACGT")).sequence_id == "s"


def test_decode_rejects_binary():
    with pytest.raises(FastaValidationError) as exc_info:
        decode_fasta_bytes(b"\xff\xfe\x00\x81binary")
    assert exc_info.value.code == "invalid_encoding"


def test_load_reference_file(reference_path):
    record = load_fasta_file(reference_path)
    assert record.sequence_id == "REF_001"
    assert len(record) == 300
    assert set(record.sequence) <= set("ACGT")
