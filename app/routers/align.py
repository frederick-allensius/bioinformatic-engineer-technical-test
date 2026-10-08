from fastapi import APIRouter

from app.routers.dependencies import FastaUpload, Reference
from app.schemas.align import AlignmentResult, AlignResponse, VariantOut
from app.schemas.common import ERROR_RESPONSES
from app.services.alignment import Alignment
from app.services.analysis import align_to_reference

router = APIRouter(prefix="/api", tags=["alignment"])


def to_alignment_result(alignment: Alignment) -> AlignmentResult:
    return AlignmentResult(
        score=float(alignment.score),
        identity=round(alignment.identity, 4),
        reference_aligned=alignment.reference_aligned,
        query_aligned=alignment.query_aligned,
    )


@router.post(
    "/align",
    response_model=AlignResponse,
    responses=ERROR_RESPONSES,
    summary="Align a FASTA sequence to the reference and call variants",
)
def align(record: FastaUpload, reference: Reference) -> AlignResponse:
    """Global alignment of the uploaded sequence against REF_001.

    Returns the score, identity, both full aligned strings and the variants
    (SNVs, insertions, deletions) in 1-based reference coordinates.
    """
    report = align_to_reference(reference.sequence, record.sequence)
    return AlignResponse(
        sequence_id=record.sequence_id,
        reference_id=reference.sequence_id,
        alignment=to_alignment_result(report.alignment),
        variants=[
            VariantOut(position=v.position, ref=v.ref, alt=v.alt, type=v.type)
            for v in report.variants
        ],
    )
