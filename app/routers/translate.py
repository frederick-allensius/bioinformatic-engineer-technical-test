from fastapi import APIRouter

from app.routers.dependencies import FastaUpload, ReadingFrame
from app.schemas.common import ERROR_RESPONSES
from app.schemas.translate import TranslateResponse
from app.services.translation import translate as translate_sequence

router = APIRouter(prefix="/api", tags=["translation"])


@router.post(
    "/translate",
    response_model=TranslateResponse,
    responses=ERROR_RESPONSES,
    summary="Translate a FASTA sequence to protein",
)
def translate(record: FastaUpload, frame: ReadingFrame) -> TranslateResponse:
    """Translate the supplied strand in reading frame 1, 2 or 3 (standard genetic code).

    Stops before the first unambiguous stop codon; N-containing codons become X.
    """
    result = translate_sequence(record.sequence, frame)
    return TranslateResponse(
        sequence_id=record.sequence_id,
        frame=result.frame,
        protein=result.protein,
        length=result.length,
        has_stop_codon=result.has_stop_codon,
    )
