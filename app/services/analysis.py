"""Use-case level functions that combine the lower-level services.

Routers call these so they stay thin: no bioinformatics logic lives in the API layer.
"""

from __future__ import annotations

from collections import Counter
from dataclasses import dataclass

from app.services.alignment import DEFAULT_SCORING, Alignment, ScoringScheme, global_align
from app.services.annotation import EFFECTS, AnnotatedVariant, Effect, annotate_variants
from app.services.translation import Translation, translate
from app.services.variants import Variant, call_variants

# Variant effects are always annotated in reference frame 1 (the CDS starts at base 1).
ANNOTATION_FRAME = 1


@dataclass(frozen=True, slots=True)
class AlignmentReport:
    alignment: Alignment
    variants: list[Variant]


@dataclass(frozen=True, slots=True)
class AnalysisReport:
    alignment: Alignment
    variants: list[AnnotatedVariant]
    translation: Translation
    reference_translation: Translation

    @property
    def effect_counts(self) -> dict[Effect, int]:
        counts = Counter(v.effect for v in self.variants)
        return {effect: counts.get(effect, 0) for effect in EFFECTS}


def align_to_reference(
    reference: str, query: str, scoring: ScoringScheme = DEFAULT_SCORING
) -> AlignmentReport:
    alignment = global_align(reference, query, scoring)
    variants = call_variants(alignment.reference_aligned, alignment.query_aligned)
    return AlignmentReport(alignment=alignment, variants=variants)


def analyze(reference: str, query: str, scoring: ScoringScheme = DEFAULT_SCORING) -> AnalysisReport:
    """Align, call variants, annotate their effects and translate the query (frame 1)."""
    report = align_to_reference(reference, query, scoring)
    return AnalysisReport(
        alignment=report.alignment,
        variants=annotate_variants(reference, report.alignment, report.variants),
        translation=translate(query, ANNOTATION_FRAME),
        reference_translation=translate(reference, ANNOTATION_FRAME),
    )
