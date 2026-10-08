"""Variant calling from a reference/query global alignment.

Conventions (all positions are 1-based reference coordinates):

* **SNV**: position of the substituted base. Only called between two *different*
  unambiguous bases (A/C/G/T); a column involving N is never an SNV.
* **deletion**: position of the first deleted reference base; ``ref`` is the deleted
  sequence and ``alt`` is ``"-"``.
* **insertion**: position of the reference base *preceding* the insertion (``0`` means
  before base 1); ``ref`` is ``"-"`` and ``alt`` is the inserted sequence.
* Consecutive gap columns of the same kind are grouped into one event.
* Events are returned sorted by position; ties keep alignment order, so a
  deletion/SNV at position p comes before an insertion after base p.
"""

from __future__ import annotations

from dataclasses import dataclass
from typing import Literal

from app.services.alignment import GAP, UNAMBIGUOUS

VariantType = Literal["SNV", "insertion", "deletion"]


@dataclass(frozen=True, slots=True)
class Variant:
    position: int
    ref: str
    alt: str
    type: VariantType

    @property
    def length(self) -> int:
        """Number of reference (deletion) or query (insertion/SNV) bases involved."""
        return len(self.alt) if self.type == "insertion" else len(self.ref)


def call_variants(reference_aligned: str, query_aligned: str) -> list[Variant]:
    if len(reference_aligned) != len(query_aligned):
        raise ValueError("aligned sequences must have equal length")

    variants: list[Variant] = []
    ref_pos = 0  # number of reference bases consumed so far
    col, n_cols = 0, len(reference_aligned)

    while col < n_cols:
        r, q = reference_aligned[col], query_aligned[col]

        if r == GAP:  # insertion run: gap in the reference
            start = col
            while col < n_cols and reference_aligned[col] == GAP:
                col += 1
            variants.append(Variant(ref_pos, GAP, query_aligned[start:col], "insertion"))
            continue

        if q == GAP:  # deletion run: gap in the query
            start, first_deleted = col, ref_pos + 1
            while col < n_cols and query_aligned[col] == GAP:
                col += 1
            ref_pos += col - start
            variants.append(Variant(first_deleted, reference_aligned[start:col], GAP, "deletion"))
            continue

        ref_pos += 1
        if r != q and r in UNAMBIGUOUS and q in UNAMBIGUOUS:
            variants.append(Variant(ref_pos, r, q, "SNV"))
        col += 1

    # Alignment order is already position order; the stable sort makes it explicit.
    variants.sort(key=lambda v: v.position)
    return variants
