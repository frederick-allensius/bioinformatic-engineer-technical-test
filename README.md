# DNA Analysis API

REST API for DNA sequence analysis against a fixed reference (`REF_001`, a 300 bp synthetic CDS):

1. **`POST /api/align`**: global alignment to the reference plus variant detection (SNVs, insertions, deletions).
2. **`POST /api/translate`**: translation to protein in reading frame 1, 2 or 3.

> Work in progress. Full documentation follows once the endpoints are implemented.

## Repository layout

```
data/       reference sequence loaded by the server at start-up
samples/    the five supplied patient FASTA files (sample inputs)
```
