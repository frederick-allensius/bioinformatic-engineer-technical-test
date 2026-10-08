from fastapi import APIRouter

from app.routers.align import to_alignment_result
from app.routers.dependencies import FastaUpload, Reference
from app.schemas.analyze import AnalyzeResponse, AnnotatedVariantOut, TranslationOut
from app.schemas.common import ERROR_RESPONSES
from app.services.analysis import analyze as analyze_sequence
from app.services.translation import Translation

router = APIRouter(prefix="/api", tags=["analysis"])


def _translation_out(t: Translation) -> TranslationOut:
    return TranslationOut(
        frame=t.frame, protein=t.protein, length=t.length, has_stop_codon=t.has_stop_codon
    )


@router.post(
    "/analyze",
    response_model=AnalyzeResponse,
    responses=ERROR_RESPONSES,
    summary="Align, translate and annotate variant effects (reference frame 1)",
)
def analyze(record: FastaUpload, reference: Reference) -> AnalyzeResponse:
    """Combined endpoint: alignment + variants (as in `/api/align`), each variant
    annotated as synonymous, missense, nonsense, stop_lost, frameshift,
    in-frame insertion/deletion or undetermined, plus the frame-1 translation.
    """
    report = analyze_sequence(reference.sequence, record.sequence)
    return AnalyzeResponse(
        sequence_id=record.sequence_id,
        reference_id=reference.sequence_id,
        alignment=to_alignment_result(report.alignment),
        variants=[
            AnnotatedVariantOut(
                position=a.variant.position,
                ref=a.variant.ref,
                alt=a.variant.alt,
                type=a.variant.type,
                effect=a.effect,
                codon=a.codon,
                ref_codon=a.ref_codon,
                alt_codon=a.alt_codon,
                ref_aa=a.ref_aa,
                alt_aa=a.alt_aa,
                protein_change=a.protein_change,
                notes=list(a.notes),
            )
            for a in report.variants
        ],
        translation=_translation_out(report.translation),
        reference_translation=_translation_out(report.reference_translation),
        effect_counts=report.effect_counts,
    )
