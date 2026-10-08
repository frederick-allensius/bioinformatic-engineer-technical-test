from typing import Literal

from pydantic import BaseModel, Field

from app.schemas.align import AlignmentResult, VariantOut

EffectName = Literal[
    "synonymous",
    "missense",
    "nonsense",
    "stop_lost",
    "frameshift",
    "inframe_insertion",
    "inframe_deletion",
    "undetermined",
]


class AnnotatedVariantOut(VariantOut):
    effect: EffectName = Field(description="Predicted consequence in reference frame 1.")
    codon: int | None = Field(
        default=None, description="1-based reference codon where the change starts."
    )
    ref_codon: str | None = Field(default=None, description="Reference codon (SNVs only).")
    alt_codon: str | None = Field(
        default=None, description="Query codon read from the alignment (SNVs only)."
    )
    ref_aa: str | None = Field(default=None, description="Reference residue(s), '-' if none.")
    alt_aa: str | None = Field(default=None, description="Altered residue(s), '-' if none.")
    protein_change: str | None = Field(
        default=None,
        description="Simplified one-letter HGVS-like notation, e.g. p.L82*, p.N36del.",
        examples=["p.L82*"],
    )
    notes: list[str] = Field(
        default_factory=list, description="Caveats, e.g. co-occurring or upstream variants."
    )


class TranslationOut(BaseModel):
    frame: int
    protein: str
    length: int
    has_stop_codon: bool


class AnalyzeResponse(BaseModel):
    sequence_id: str
    reference_id: str
    alignment: AlignmentResult
    variants: list[AnnotatedVariantOut]
    translation: TranslationOut = Field(description="Query translated in frame 1.")
    reference_translation: TranslationOut = Field(
        description="Reference translated in frame 1, for comparison."
    )
    effect_counts: dict[EffectName, int] = Field(description="Number of variants per effect.")
