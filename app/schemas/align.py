from typing import Literal

from pydantic import BaseModel, Field


class AlignmentResult(BaseModel):
    score: float = Field(description="Global alignment score under the documented scheme.")
    identity: float = Field(
        ge=0,
        le=1,
        description="Identical A/C/G/T columns / all alignment columns (gaps and N included), "
        "rounded to 4 decimals.",
    )
    reference_aligned: str = Field(description="Full aligned reference, '-' marks gaps.")
    query_aligned: str = Field(description="Full aligned query, '-' marks gaps.")


class VariantOut(BaseModel):
    position: int = Field(
        ge=0,
        description="1-based reference coordinate. Insertions: preceding reference base "
        "(0 = before base 1).",
    )
    ref: str = Field(description="Reference allele, '-' for insertions.")
    alt: str = Field(description="Query allele, '-' for deletions.")
    type: Literal["SNV", "insertion", "deletion"]


class AlignResponse(BaseModel):
    sequence_id: str
    reference_id: str
    alignment: AlignmentResult
    variants: list[VariantOut]

    model_config = {
        "json_schema_extra": {
            "example": {
                "sequence_id": "patient_02",
                "reference_id": "REF_001",
                "alignment": {
                    "score": 294.0,
                    "identity": 0.99,
                    "reference_aligned": "ATGAGTTCG...",
                    "query_aligned": "ATGAGTTCG...",
                },
                "variants": [
                    {"position": 52, "ref": "C", "alt": "A", "type": "SNV"},
                    {"position": 177, "ref": "A", "alt": "T", "type": "SNV"},
                    {"position": 245, "ref": "T", "alt": "G", "type": "SNV"},
                ],
            }
        }
    }
