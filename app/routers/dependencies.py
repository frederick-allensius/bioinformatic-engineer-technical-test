"""FastAPI dependencies shared by the routers (upload handling, reference access)."""

from __future__ import annotations

from typing import Annotated

from fastapi import Depends, File, Form, Request, UploadFile, status

from app.core.config import Settings
from app.core.errors import APIError
from app.services.fasta import FastaRecord, FastaValidationError, decode_fasta_bytes, parse_fasta


def get_app_settings(request: Request) -> Settings:
    return request.app.state.settings


def get_reference(request: Request) -> FastaRecord:
    return request.app.state.reference


async def fasta_upload(
    settings: Annotated[Settings, Depends(get_app_settings)],
    file: Annotated[
        UploadFile | None,
        File(description="Single-record FASTA file (A/C/G/T/N, 1-1000 nt)."),
    ] = None,
) -> FastaRecord:
    """Read the uploaded ``file`` part and validate it as a single FASTA record."""
    if file is None:
        raise APIError(
            status.HTTP_400_BAD_REQUEST,
            "missing_file",
            "No file uploaded. Send the FASTA file as multipart/form-data in the field 'file'.",
            field="file",
        )
    try:
        # Read one byte past the cap so oversize uploads are detected without
        # reading the whole body into memory.
        data = await file.read(settings.max_upload_bytes + 1)
    finally:
        await file.close()

    if len(data) > settings.max_upload_bytes:
        raise FastaValidationError(
            "file_too_large",
            f"Uploaded file exceeds {settings.max_upload_bytes} bytes; a valid FASTA file "
            f"holds at most {settings.max_sequence_length} nt.",
        )
    return parse_fasta(decode_fasta_bytes(data), max_length=settings.max_sequence_length)


def reading_frame(
    frame: Annotated[
        str | None,
        Form(
            description="Reading frame: 1, 2 or 3 (default 1).",
            json_schema_extra={"enum": ["1", "2", "3"]},
        ),
    ] = None,
) -> int:
    """Parse the optional ``frame`` form field.

    Parsed by hand rather than as ``int`` so that every bad value ("0", "4", "1.5",
    "abc") gets the same clear message instead of a generic type error.
    """
    if frame is None:  # omitted (FastAPI also maps an empty form value to None)
        return 1
    value = frame.strip()
    if value not in {"1", "2", "3"}:
        raise APIError(
            status.HTTP_422_UNPROCESSABLE_CONTENT,
            "invalid_frame",
            f"frame must be 1, 2 or 3; got {frame!r}.",
            field="frame",
        )
    return int(value)


FastaUpload = Annotated[FastaRecord, Depends(fasta_upload)]
Reference = Annotated[FastaRecord, Depends(get_reference)]
ReadingFrame = Annotated[int, Depends(reading_frame)]
