from pydantic import BaseModel, Field


class ErrorDetail(BaseModel):
    code: str = Field(
        description="Stable machine-readable error code.", examples=["multiple_records"]
    )
    message: str = Field(description="Human-readable explanation of what went wrong.")
    field: str | None = Field(
        default=None, description="Request field the error refers to.", examples=["file"]
    )


class ErrorResponse(BaseModel):
    error: ErrorDetail


ERROR_RESPONSES: dict[int | str, dict] = {
    400: {"model": ErrorResponse, "description": "Missing `file` part."},
    422: {"model": ErrorResponse, "description": "Invalid FASTA content or parameter value."},
}
