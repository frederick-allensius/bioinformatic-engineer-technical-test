# DNA Sequence Analysis API

[![CI](https://github.com/frederick-allensius/bioinformatic-engineer-technical-test/actions/workflows/ci.yml/badge.svg)](https://github.com/frederick-allensius/bioinformatic-engineer-technical-test/actions/workflows/ci.yml)

A REST API built with FastAPI that analyses a DNA sequence against the reference `REF_001` (300 bp):

- `POST /api/align`: global alignment to the reference and variant detection
- `POST /api/translate`: translation to protein in frame 1, 2 or 3
- `POST /api/analyze` (bonus): alignment + translation + the effect of each variant
- Web interface at `/` (bonus), plus a Dockerfile (bonus)

The alignment (Needleman-Wunsch) and the translation are implemented without bioinformatics libraries. Biopython is only used in the tests, to double-check the results.

## Installation

Requires Python 3.11 or newer.

```bash
git clone https://github.com/frederick-allensius/bioinformatic-engineer-technical-test.git
cd bioinformatic-engineer-technical-test
python -m venv .venv
source .venv/bin/activate        # Windows: .venv\Scripts\activate
pip install -r requirements.txt
```

## Running the API

```bash
uvicorn app.main:app --reload
```

- API: <http://localhost:8000>
- Interactive docs: <http://localhost:8000/docs>
- Web interface: <http://localhost:8000>

With Docker:

```bash
docker build -t dna-analysis-api .
docker run --rm -p 8000:8000 dna-analysis-api
```

## Usage

Both endpoints take `multipart/form-data` with the FASTA file in the field `file`.

### Align

```bash
curl -X POST http://localhost:8000/api/align -F "file=@samples/patient_02.fasta"
```

```json
{
  "sequence_id": "patient_02",
  "reference_id": "REF_001",
  "alignment": {
    "score": 294.0,
    "identity": 0.99,
    "reference_aligned": "ATGAGTTCG...",
    "query_aligned": "ATGAGTTCG..."
  },
  "variants": [
    {"position": 52, "ref": "C", "alt": "A", "type": "SNV"},
    {"position": 177, "ref": "A", "alt": "T", "type": "SNV"},
    {"position": 245, "ref": "T", "alt": "G", "type": "SNV"}
  ]
}
```

(The aligned strings are shortened here; the API returns them in full.)

### Translate

`frame` is optional and defaults to 1.

```bash
curl -X POST http://localhost:8000/api/translate -F "file=@samples/patient_02.fasta" -F "frame=1"
```

```json
{
  "sequence_id": "patient_02",
  "frame": 1,
  "protein": "MSSFVPLLYVSRVGTSNRLFSLLTKFTWVSVTRLRNPAEFAACVRHPYLAHSSQSHEHNLALVITSGQSTPDSKGEQNWVS",
  "length": 81,
  "has_stop_codon": true
}
```

### Analyze (bonus)

```bash
curl -X POST http://localhost:8000/api/analyze -F "file=@samples/patient_02.fasta"
```

Returns the alignment, the translation and, for each variant, its effect, for example `{"position": 245, "type": "SNV", "effect": "nonsense", "protein_change": "p.L82*", ...}`.

### Errors

Invalid requests return a JSON error:

```json
{"error": {"code": "multiple_records", "message": "Expected exactly one FASTA record but found 2 (additional header on line 3).", "field": "file"}}
```

- **400**: no file was uploaded.
- **422**: the file is not valid FASTA (no header, more than one record, invalid characters, whitespace inside the sequence, empty, or longer than 1,000 nt), or `frame` is not 1, 2 or 3.

## Alignment

Global alignment with the Needleman-Wunsch algorithm (`app/services/alignment.py`).

| | Score |
|---|---|
| Match | +1 |
| Mismatch | −1 |
| Gap | −2 per gap character |
| Any position with N | 0 |

N is an unknown base, so it is neither rewarded nor penalised.

**Ties:** when several alignments have the same best score, the API prefers the one with the fewest separate gaps, so a 3-base deletion is reported as one event instead of being split up. If there is still a tie, gaps are placed as far left as possible. The result is therefore always the same for the same input.

## Identity and variants

```
identity = identical A/C/G/T positions / all alignment columns (including gaps and N)
```

N/N does not count as identical.

| Type | Position (1-based, on the reference) | ref | alt |
|---|---|---|---|
| SNV | the changed base | reference base | new base |
| deletion | first deleted base | deleted bases | `-` |
| insertion | base before the insertion (0 = before base 1) | `-` | inserted bases |

- An SNV is only called between two different A/C/G/T bases, so N is never an SNV.
- Consecutive inserted or deleted bases are grouped into one variant.
- Variants are sorted by position. If none are found, `variants` is `[]`.

## Translation

- Standard genetic code. Frame 1, 2 or 3 starts at nucleotide 1, 2 or 3.
- Only complete codons are translated; 1–2 leftover bases at the end are ignored.
- Any codon containing N becomes `X`.
- Translation stops before the first stop codon (TAA, TAG, TGA). The stop is not included in `protein`, and `has_stop_codon` is `true`.
- If the frame has no complete codon, the result is `protein: ""`, `length: 0`, `has_stop_codon: false`.

## Variant effects (bonus)

`/api/analyze` classifies each variant using reference frame 1: `synonymous`, `missense`, `nonsense`, `stop_lost`, `frameshift`, `inframe_insertion` or `inframe_deletion`.

- If the codon of an SNV contains N, the effect is `undetermined`.
- SNVs in the same codon are evaluated together, as one codon change.
- Each insertion/deletion is classified by its length (a multiple of 3 is in-frame, otherwise frameshift).

## Web interface (bonus)

Upload or paste a FASTA file (or pick a sample), choose Align / Translate / Analyze, and view the results and any error messages.

![Web interface](docs/ui-align.png)

## Running the tests

```bash
pytest
```

The tests cover the five sample files (exact variants and positions), invalid and multiple-record FASTA, lowercase input, N-containing codons, each reading frame, stop codons, incomplete final codons and the length limit. The alignment and translation results are also checked against Biopython.

```
$ pytest
============================= test session starts =============================
platform win32 -- Python 3.12.5, pytest-9.1.1, pluggy-1.6.0
rootdir: D:\1. Work\task\dna-analysis-api
configfile: pyproject.toml
testpaths: tests
plugins: anyio-4.15.1, hypothesis-6.168.5, cov-7.1.0
collected 258 items

tests\api\test_align_endpoint.py ....................................... [ 15%]
......                                                                   [ 17%]
tests\api\test_analyze_endpoint.py ..........                            [ 21%]
tests\api\test_error_handlers.py ..                                      [ 22%]
tests\api\test_translate_endpoint.py ................................... [ 35%]
..............                                                           [ 41%]
tests\api\test_web_ui.py .....                                           [ 43%]
tests\unit\test_alignment.py ....................                        [ 50%]
tests\unit\test_annotation.py .........................                  [ 60%]
tests\unit\test_fasta.py ............................................... [ 78%]
..                                                                       [ 79%]
tests\unit\test_translation.py ...................................       [ 93%]
tests\unit\test_variants.py ..................                           [100%]

============================ 258 passed in 15.74s =============================
```

The tests also run automatically on GitHub Actions for every push.

## Project structure

```
app/
  routers/    API endpoints
  schemas/    request/response models
  services/   logic: FASTA parsing, alignment, variant calling, translation, annotation
  static/     web interface
data/         reference sequence (REF_001)
samples/      the five patient FASTA files
tests/        unit and API tests
```

## Assumptions

- All sequences are the forward strand of the same complete fragment, so global alignment is used.
- Spaces at the start or end of a line are ignored; spaces inside the sequence are rejected.
- If `frame` is empty or not sent, frame 1 is used.

## Limitations

- The alignment is written in plain Python, which is fast enough for sequences up to 1,000 nt but would be slow for much longer ones.
- Gaps use a simple per-base penalty; a gap-open + gap-extend penalty is more realistic biologically.
- When an insertion or deletion falls inside a region of Ns, its exact position is arbitrary.
- Apart from SNVs in the same codon, each variant's effect is evaluated on its own. For example, two frameshifts that cancel each other out are both still reported as `frameshift`.
