"""End-to-end tests for the bonus POST /api/analyze endpoint."""

import pytest

ENDPOINT = "/api/analyze"


def effects(body):
    return [(v["position"], v["type"], v["effect"], v["protein_change"]) for v in body["variants"]]


def test_patient_02_snv_effects(post_sample):
    body = post_sample(ENDPOINT, "patient_02").json()
    assert effects(body) == [
        (52, "SNV", "synonymous", "p.R18="),  # CGA > AGA, both Arg
        (177, "SNV", "missense", "p.K59N"),  # AAA > AAT
        (245, "SNV", "nonsense", "p.L82*"),  # TTA > TGA
    ]
    nonsense = body["variants"][2]
    assert (nonsense["codon"], nonsense["ref_codon"], nonsense["alt_codon"]) == (82, "TTA", "TGA")
    # The nonsense call agrees with the translation: 81 residues, then the stop.
    assert body["translation"]["length"] == 81
    assert body["translation"]["has_stop_codon"] is True


def test_patient_03_deletions(post_sample):
    body = post_sample(ENDPOINT, "patient_03").json()
    assert effects(body) == [
        (106, "deletion", "inframe_deletion", "p.N36del"),
        (208, "deletion", "frameshift", "p.T70Rfs*13"),
    ]


def test_patient_05_insertions(post_sample):
    body = post_sample(ENDPOINT, "patient_05").json()
    assert effects(body) == [
        (18, "insertion", "inframe_insertion", "p.P6_L7insI"),
        (230, "insertion", "frameshift", "p.Q77Hfs*19"),
    ]


@pytest.mark.parametrize("name", ["patient_01", "patient_04"])
def test_samples_without_variants(post_sample, name):
    body = post_sample(ENDPOINT, name).json()
    assert body["variants"] == []
    assert set(body["effect_counts"].values()) == {0}


def test_effect_counts(post_sample):
    counts = post_sample(ENDPOINT, "patient_02").json()["effect_counts"]
    assert counts["synonymous"] == counts["missense"] == counts["nonsense"] == 1
    assert sum(counts.values()) == 3


def test_reference_translation_is_included(post_sample):
    body = post_sample(ENDPOINT, "patient_02").json()
    assert body["reference_translation"]["length"] == 99
    assert body["reference_translation"]["protein"].startswith("MSSFVP")


def test_analyze_is_consistent_with_align_and_translate(post_sample):
    for name in ["patient_02", "patient_03", "patient_05"]:
        analyzed = post_sample(ENDPOINT, name).json()
        aligned = post_sample("/api/align", name).json()
        translated = post_sample("/api/translate", name, frame=1).json()

        assert analyzed["alignment"] == aligned["alignment"]
        core = [{k: v[k] for k in ("position", "ref", "alt", "type")} for v in analyzed["variants"]]
        assert core == aligned["variants"]
        assert analyzed["translation"]["protein"] == translated["protein"]


def test_n_codon_snv_is_undetermined(post_fasta, reference_seq):
    # Codon 82 (244-246) = TTA. Make it TGN: SNV at 245 plus an N at 246.
    query = reference_seq[:244] + "GN" + reference_seq[246:]
    body = post_fasta(ENDPOINT, ">q\n" + query).json()
    (variant,) = body["variants"]
    assert (variant["position"], variant["effect"]) == (245, "undetermined")
    assert variant["alt_codon"] == "TGN"


def test_invalid_input_uses_shared_validation(post_fasta, client):
    response = post_fasta(ENDPOINT, ">a\nACGT\n>b\nACGT")
    assert response.status_code == 422
    assert response.json()["error"]["code"] == "multiple_records"
    assert client.post(ENDPOINT).status_code == 400
