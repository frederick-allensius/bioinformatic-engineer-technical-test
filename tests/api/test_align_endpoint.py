"""End-to-end tests for POST /api/align."""

import pytest

from app.core.config import Settings
from app.main import create_app
from tests.helpers import apply_variants

ENDPOINT = "/api/align"


def snv(pos, ref, alt):
    return {"position": pos, "ref": ref, "alt": alt, "type": "SNV"}


def deletion(pos, ref):
    return {"position": pos, "ref": ref, "alt": "-", "type": "deletion"}


def insertion(pos, alt):
    return {"position": pos, "ref": "-", "alt": alt, "type": "insertion"}


# Expected calls for the supplied samples, verified by hand against the reference
# (see test_expected_calls_are_consistent_with_reference below).
EXPECTED = {
    "patient_01": {"score": 300.0, "identity": 1.0, "variants": []},
    "patient_02": {
        "score": 294.0,
        "identity": 0.99,
        "variants": [snv(52, "C", "A"), snv(177, "A", "T"), snv(245, "T", "G")],
    },
    "patient_03": {
        # 296 matches - 4 gap columns * 2 = 288; identity 296 / 300 columns
        "score": 288.0,
        "identity": 0.9867,
        "variants": [deletion(106, "AAC"), deletion(208, "A")],
    },
    "patient_04": {
        # 170 N columns score 0 and are never SNVs; 130 matches remain
        "score": 130.0,
        "identity": 0.4333,
        "variants": [],
    },
    "patient_05": {
        # 300 matches - 4 gap columns * 2 = 292; identity 300 / 304 columns
        "score": 292.0,
        "identity": 0.9868,
        "variants": [insertion(18, "ATC"), insertion(230, "T")],
    },
}


# --------------------------------------------------------------------------- supplied samples


def test_patient_02_reproduces_the_example_response(post_sample):
    response = post_sample(ENDPOINT, "patient_02")
    assert response.status_code == 200
    body = response.json()
    assert body["sequence_id"] == "patient_02"
    assert body["reference_id"] == "REF_001"
    assert body["alignment"]["score"] == 294.0
    assert body["alignment"]["identity"] == 0.99
    assert body["alignment"]["reference_aligned"].startswith("ATGAGTTCG")
    assert body["alignment"]["query_aligned"].startswith("ATGAGTTCG")
    assert body["variants"] == EXPECTED["patient_02"]["variants"]


@pytest.mark.parametrize("name", sorted(EXPECTED))
def test_supplied_samples(post_sample, name):
    response = post_sample(ENDPOINT, name)
    assert response.status_code == 200
    body = response.json()
    expected = EXPECTED[name]
    assert body["sequence_id"] == name
    assert body["alignment"]["score"] == expected["score"]
    assert body["alignment"]["identity"] == expected["identity"]
    assert body["variants"] == expected["variants"]


@pytest.mark.parametrize("name", sorted(EXPECTED))
def test_aligned_strings_are_complete(post_sample, samples_dir, reference_seq, name):
    body = post_sample(ENDPOINT, name).json()["alignment"]
    query = "".join((samples_dir / f"{name}.fasta").read_text().splitlines()[1:])
    assert len(body["reference_aligned"]) == len(body["query_aligned"])
    assert body["reference_aligned"].replace("-", "") == reference_seq
    assert body["query_aligned"].replace("-", "") == query


@pytest.mark.parametrize("name", sorted(EXPECTED))
def test_variants_rebuild_the_query(post_sample, samples_dir, reference_seq, name):
    """Applying the reported variants to the reference must give back the query.

    N positions are not variants, so they are compared as wildcards.
    """
    variants = post_sample(ENDPOINT, name).json()["variants"]
    query = "".join((samples_dir / f"{name}.fasta").read_text().splitlines()[1:])
    rebuilt = apply_variants(reference_seq, variants)
    assert len(rebuilt) == len(query)
    assert all(q in ("N", r) for r, q in zip(rebuilt, query, strict=True))


def test_expected_calls_are_consistent_with_reference(reference_seq):
    """Sanity-check the hand-written expectations against the reference itself."""
    assert reference_seq[52 - 1] == "C"
    assert reference_seq[177 - 1] == "A"
    assert reference_seq[245 - 1] == "T"
    assert reference_seq[106 - 1 : 106 - 1 + 3] == "AAC"
    assert reference_seq[208 - 1] == "A"
    # Insertion at 18 = after reference base 18.
    assert reference_seq[:18] == "ATGAGTTCGTTCGTGCCT"


def test_variants_are_sorted_by_position(post_sample):
    for name in EXPECTED:
        positions = [v["position"] for v in post_sample(ENDPOINT, name).json()["variants"]]
        assert positions == sorted(positions)


# --------------------------------------------------------------------------- synthetic inputs


def test_identical_to_reference_returns_empty_variant_list(post_fasta, reference_seq):
    body = post_fasta(ENDPOINT, f">same\n{reference_seq}").json()
    assert body["variants"] == []
    assert body["alignment"]["identity"] == 1.0


def test_lowercase_input_gives_same_result(post_sample, post_fasta, samples_dir):
    text = (samples_dir / "patient_02.fasta").read_text()
    header, *seq_lines = text.splitlines()
    lower = header + "\n" + "\n".join(line.lower() for line in seq_lines)
    assert post_fasta(ENDPOINT, lower).json() == post_sample(ENDPOINT, "patient_02").json()


def test_multiline_crlf_and_blank_lines(post_fasta, reference_seq):
    chunks = [reference_seq[i : i + 70] for i in range(0, 300, 70)]
    text = "\r\n\r\n>wrapped some description\r\n" + "\r\n\r\n".join(chunks) + "\r\n"
    body = post_fasta(ENDPOINT, text).json()
    assert body["sequence_id"] == "wrapped"
    assert body["variants"] == []


def test_deletion_at_sequence_start(post_fasta, reference_seq):
    # ATG deleted; "TGA" at 2-4 would also fit, the left-most placement is reported.
    body = post_fasta(ENDPOINT, ">q\n" + reference_seq[3:]).json()
    assert body["variants"] == [deletion(1, "ATG")]


def test_insertion_before_first_base(post_fasta, reference_seq):
    body = post_fasta(ENDPOINT, ">q\nGG" + reference_seq).json()
    assert body["variants"] == [insertion(0, "GG")]


def test_insertion_after_last_base(post_fasta, reference_seq):
    body = post_fasta(ENDPOINT, ">q\n" + reference_seq + "TT").json()
    assert body["variants"] == [insertion(300, "TT")]


def test_multi_base_deletion_is_grouped(post_fasta, reference_seq):
    # Remove 6 bases (positions 151-156) in one block.
    query = reference_seq[:150] + reference_seq[156:]
    body = post_fasta(ENDPOINT, ">q\n" + query).json()
    assert len(body["variants"]) == 1
    event = body["variants"][0]
    assert event["type"] == "deletion"
    assert len(event["ref"]) == 6
    assert apply_variants(reference_seq, body["variants"]) == query


def test_mixed_snv_insertion_and_deletion(post_fasta, reference_seq):
    query = (
        reference_seq[:9]
        + "G"  # SNV at position 10 (T>G)
        + reference_seq[10:100]
        + "CCC"  # insertion after base 100
        + reference_seq[100:199]
        + reference_seq[201:]  # deletion of bases 200-201
    )
    assert reference_seq[9] == "T"
    body = post_fasta(ENDPOINT, ">q\n" + query).json()
    types = [v["type"] for v in body["variants"]]
    assert types == ["SNV", "insertion", "deletion"]
    assert body["variants"][0] == snv(10, "T", "G")
    assert apply_variants(reference_seq, body["variants"]) == query


def test_n_at_variant_site_is_not_an_snv(post_sample, post_fasta, samples_dir):
    text = (samples_dir / "patient_02.fasta").read_text()
    header, *lines = text.splitlines()
    seq = "".join(lines)
    masked = seq[:51] + "N" + seq[52:]  # position 52 carries the C>A SNV
    body = post_fasta(ENDPOINT, f"{header}\n{masked}").json()
    assert [v["position"] for v in body["variants"]] == [177, 245]


def test_all_n_query_is_valid(post_fasta):
    body = post_fasta(ENDPOINT, ">unknown\n" + "N" * 300)
    assert body.status_code == 200
    assert body.json()["variants"] == []
    assert body.json()["alignment"]["identity"] == 0.0


def test_max_length_query_is_accepted(post_fasta, reference_seq):
    response = post_fasta(ENDPOINT, ">long\n" + (reference_seq * 4)[:1000])
    assert response.status_code == 200


def test_single_base_query(post_fasta):
    response = post_fasta(ENDPOINT, ">tiny\nA")
    assert response.status_code == 200
    assert response.json()["alignment"]["query_aligned"].replace("-", "") == "A"


# --------------------------------------------------------------------------- errors


def assert_error(response, status: int, code: str):
    assert response.status_code == status, response.text
    body = response.json()
    assert body["error"]["code"] == code
    assert body["error"]["message"]
    return body


def test_missing_file_returns_400(client):
    body = assert_error(client.post(ENDPOINT), 400, "missing_file")
    assert body["error"]["field"] == "file"


def test_file_sent_as_text_field(client):
    assert_error(client.post(ENDPOINT, data={"file": "ACGT"}), 422, "invalid_request")


def test_wrong_field_name_is_missing_file(client):
    response = client.post(ENDPOINT, files={"fasta": ("q.fasta", b">q\nACGT")})
    assert_error(response, 400, "missing_file")


@pytest.mark.parametrize(
    ("content", "code"),
    [
        ("", "empty_file"),
        ("ACGTACGT\n", "malformed_fasta"),
        (">\nACGT", "empty_header"),
        (">only_header\n", "empty_sequence"),
        (">a\nACGT\n>b\nACGT\n", "multiple_records"),
        (">q\nACGU", "invalid_character"),
        (">q\nACGT-ACGT", "invalid_character"),
        (">q\nACG TACGT", "embedded_whitespace"),
        (">q\n" + "A" * 1001, "sequence_too_long"),
    ],
)
def test_invalid_fasta_returns_422(post_fasta, content, code):
    body = assert_error(post_fasta(ENDPOINT, content), 422, code)
    assert body["error"]["field"] == "file"


def test_non_utf8_file(post_fasta):
    assert_error(post_fasta(ENDPOINT, b">q\n\xff\xfeAC"), 422, "invalid_encoding")


def test_oversized_upload_is_rejected_before_parsing():
    from fastapi.testclient import TestClient

    small = TestClient(create_app(Settings(max_upload_bytes=64)))
    response = small.post(ENDPOINT, files={"file": ("big.fasta", b">q\n" + b"A" * 100)})
    assert_error(response, 422, "file_too_large")


def test_get_is_not_allowed(client):
    assert_error(client.get(ENDPOINT), 405, "method_not_allowed")
