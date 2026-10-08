"""End-to-end tests for POST /api/translate."""

import pytest

ENDPOINT = "/api/translate"

EXAMPLE_PROTEIN = (
    "MSSFVPLLYVSRVGTSNRLFSLLTKFTWVSVTRLRNPAEFAACVRHPYLAHSSQSHEHNLALVITSGQSTPDSKGEQNWVS"
)
# Full reference CDS (= patient_01) translated in frame 1: 99 aa then TGA at codon 100.
REFERENCE_PROTEIN = (
    "MSSFVPLLYVSRVGTSNRLFSLLTKFTWVSVTRLRNPAEFAACVRHPYLAHSSQSHEHKLALVITSGQSTPDSKGEQNWVSLRLAPSE"
    "LESQPISDQCQ"
)


def translate_body(response, status=200):
    assert response.status_code == status, response.text
    return response.json()


# --------------------------------------------------------------------------- supplied samples


def test_patient_02_reproduces_the_example_response(post_sample):
    body = translate_body(post_sample(ENDPOINT, "patient_02", frame=1))
    assert body == {
        "sequence_id": "patient_02",
        "frame": 1,
        "protein": EXAMPLE_PROTEIN,
        "length": 81,
        "has_stop_codon": True,
    }


def test_frame_defaults_to_1(post_sample):
    body = translate_body(post_sample(ENDPOINT, "patient_02"))
    assert body["frame"] == 1
    assert body["protein"] == EXAMPLE_PROTEIN


def test_reference_like_sample_translates_full_cds(post_sample):
    body = translate_body(post_sample(ENDPOINT, "patient_01"))
    assert body["protein"] == REFERENCE_PROTEIN
    assert body["length"] == 99
    assert body["has_stop_codon"] is True


def test_n_rich_sample_uses_x(post_sample):
    body = translate_body(post_sample(ENDPOINT, "patient_04"))
    assert body["length"] == 99  # same CDS length as the reference: N never stops translation
    assert body["protein"].startswith("MSSFVPX")
    assert body["protein"].count("X") > 50
    assert "*" not in body["protein"]


@pytest.mark.parametrize("name", ["patient_01", "patient_02", "patient_03", "patient_05"])
@pytest.mark.parametrize("frame", [1, 2, 3])
def test_samples_match_biopython_in_every_frame(post_sample, samples_dir, name, frame):
    seq_module = pytest.importorskip("Bio.Seq")
    seq = "".join((samples_dir / f"{name}.fasta").read_text().splitlines()[1:])
    framed = seq[frame - 1 :]
    framed = framed[: len(framed) - len(framed) % 3]
    full = str(seq_module.Seq(framed).translate())

    body = translate_body(post_sample(ENDPOINT, name, frame=frame))
    assert body["frame"] == frame
    assert body["protein"] == full.split("*")[0]
    assert body["has_stop_codon"] == ("*" in full)
    assert body["length"] == len(body["protein"])


# --------------------------------------------------------------------------- frames


@pytest.mark.parametrize(
    ("frame", "protein"),
    [
        (1, "NGQR"),  # AAT GGC CAA AGG
        (2, "MAK"),  # ATG GCC AAA + GG ignored
        (3, "WPK"),  # TGG CCA AAG + G ignored
    ],
)
def test_each_reading_frame(post_fasta, frame, protein):
    body = translate_body(post_fasta(ENDPOINT, ">q\nAATGGCCAAAGG", frame=frame))
    assert body["frame"] == frame
    assert body["protein"] == protein
    assert body["length"] == len(protein)
    assert body["has_stop_codon"] is False


def test_frame_with_surrounding_whitespace_is_accepted(post_fasta):
    assert translate_body(post_fasta(ENDPOINT, ">q\nAATGGC", frame=" 2 "))["frame"] == 2


@pytest.mark.parametrize("frame", ["0", "4", "-1", "1.5", "abc", "one", " ", "1,2"])
def test_invalid_frame_returns_422(post_fasta, frame):
    body = translate_body(post_fasta(ENDPOINT, ">q\nATGGCC", frame=frame), status=422)
    assert body["error"]["code"] == "invalid_frame"
    assert body["error"]["field"] == "frame"


# --------------------------------------------------------------------------- codon rules


def test_stops_at_first_stop_codon(post_fasta):
    body = translate_body(post_fasta(ENDPOINT, ">q\nATGGCCTAGGGCTAA"))
    assert (body["protein"], body["length"], body["has_stop_codon"]) == ("MA", 2, True)


def test_stop_codon_first_gives_empty_protein(post_fasta):
    body = translate_body(post_fasta(ENDPOINT, ">q\nTGAATGGCC"))
    assert (body["protein"], body["length"], body["has_stop_codon"]) == ("", 0, True)


def test_no_stop_codon(post_fasta):
    body = translate_body(post_fasta(ENDPOINT, ">q\nATGGCCGGG"))
    assert (body["protein"], body["has_stop_codon"]) == ("MAG", False)


def test_n_codons_become_x_and_count_towards_length(post_fasta):
    # GCN is Ala for any N, but must still be X; NNN is X.
    body = translate_body(post_fasta(ENDPOINT, ">q\nATGGCNNNNTGG"))
    assert (body["protein"], body["length"]) == ("MXXW", 4)


def test_ambiguous_stop_does_not_stop(post_fasta):
    # TNA may or may not be a stop codon, so it is X and translation continues.
    body = translate_body(post_fasta(ENDPOINT, ">q\nATGTNAGCC"))
    assert (body["protein"], body["has_stop_codon"]) == ("MXA", False)


@pytest.mark.parametrize(("seq", "protein"), [("ATGGCCG", "MA"), ("ATGGCCGG", "MA")])
def test_incomplete_final_codon_is_ignored(post_fasta, seq, protein):
    body = translate_body(post_fasta(ENDPOINT, f">q\n{seq}"))
    assert (body["protein"], body["has_stop_codon"]) == (protein, False)


@pytest.mark.parametrize(
    ("seq", "frame"), [("A", 1), ("AT", 1), ("ATG", 2), ("ATG", 3), ("ATGC", 3)]
)
def test_no_complete_codon(post_fasta, seq, frame):
    body = translate_body(post_fasta(ENDPOINT, f">q\n{seq}", frame=frame))
    assert body == {
        "sequence_id": "q",
        "frame": frame,
        "protein": "",
        "length": 0,
        "has_stop_codon": False,
    }


def test_lowercase_input(post_fasta):
    upper = translate_body(post_fasta(ENDPOINT, ">q\nATGGCCNNNTAA"))
    lower = translate_body(post_fasta(ENDPOINT, ">q\natggccnnntaa"))
    assert lower == upper
    assert lower["protein"] == "MAX"


def test_multiline_input(post_fasta):
    body = translate_body(post_fasta(ENDPOINT, ">q desc\nATGG\n\nCCT\nAA\n"))
    assert body["sequence_id"] == "q"
    assert (body["protein"], body["has_stop_codon"]) == ("MA", True)


# --------------------------------------------------------------------------- validation


def test_length_limit(post_fasta):
    assert translate_body(post_fasta(ENDPOINT, ">q\n" + "ATG" * 333 + "A"))["length"] == 333
    body = translate_body(post_fasta(ENDPOINT, ">q\n" + "A" * 1001), status=422)
    assert body["error"]["code"] == "sequence_too_long"


def test_missing_file_returns_400(client):
    response = client.post(ENDPOINT, data={"frame": "1"})
    assert response.status_code == 400
    assert response.json()["error"]["code"] == "missing_file"


@pytest.mark.parametrize(
    ("content", "code"),
    [
        (">a\nATG\n>b\nATG", "multiple_records"),
        (">q\nATGXYZ", "invalid_character"),
        (">q\nATG GCC", "embedded_whitespace"),
        ("ATGGCC", "malformed_fasta"),
        (">q\n", "empty_sequence"),
    ],
)
def test_invalid_fasta_returns_422(post_fasta, content, code):
    body = translate_body(post_fasta(ENDPOINT, content), status=422)
    assert body["error"]["code"] == code
