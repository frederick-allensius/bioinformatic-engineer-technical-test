"""Use-case level functions that combine the lower-level services.

Routers call these so they stay thin: no bioinformatics logic lives in the API layer.
"""

from __future__ import annotations

from dataclasses import dataclass

from app.services.alignment import DEFAULT_SCORING, Alignment, ScoringScheme, global_align
from app.services.variants import Variant, call_variants


@dataclass(frozen=True, slots=True)
class AlignmentReport:
    alignment: Alignment
    variants: list[Variant]


def align_to_reference(
    reference: str, query: str, scoring: ScoringScheme = DEFAULT_SCORING
) -> AlignmentReport:
    alignment = global_align(reference, query, scoring)
    variants = call_variants(alignment.reference_aligned, alignment.query_aligned)
    return AlignmentReport(alignment=alignment, variants=variants)
