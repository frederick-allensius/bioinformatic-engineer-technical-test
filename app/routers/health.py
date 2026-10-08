from fastapi import APIRouter
from pydantic import BaseModel

from app import __version__
from app.routers.dependencies import Reference

router = APIRouter(prefix="/api", tags=["meta"])


class HealthResponse(BaseModel):
    status: str
    version: str
    reference_id: str
    reference_length: int


@router.get("/health", response_model=HealthResponse, summary="Liveness check")
def health(reference: Reference) -> HealthResponse:
    return HealthResponse(
        status="ok",
        version=__version__,
        reference_id=reference.sequence_id,
        reference_length=len(reference),
    )
