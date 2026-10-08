from pydantic import BaseModel, Field


class TranslateResponse(BaseModel):
    sequence_id: str
    frame: int = Field(ge=1, le=3, description="Reading frame used (1, 2 or 3).")
    protein: str = Field(description="One-letter protein, X for N-containing codons, no '*'.")
    length: int = Field(ge=0, description="len(protein), X included, stop excluded.")
    has_stop_codon: bool = Field(description="True if an unambiguous stop codon was reached.")

    model_config = {
        "json_schema_extra": {
            "example": {
                "sequence_id": "patient_02",
                "frame": 1,
                "protein": "MSSFVPLLYVSRVGTSNRLFSLLTKFTWVSVTRLRNPAEFAACVRHPYLAHSSQSHEHNLALVITSGQ"
                "STPDSKGEQNWVS",
                "length": 81,
                "has_stop_codon": True,
            }
        }
    }
