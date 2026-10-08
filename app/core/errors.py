"""Consistent JSON error responses.

Every error returned by the API (validation, 404, 405, unexpected 500) uses the same
envelope, so clients only need one code path:

    {"error": {"code": "invalid_character", "message": "...", "field": "file"}}

Status code policy:

* **400 Bad Request**: the request is missing something required (no ``file`` part).
* **422 Unprocessable Content**: the request is well-formed but its content is invalid
  (malformed FASTA, bad characters, length limit, invalid ``frame``).
"""

from __future__ import annotations

import logging

from fastapi import FastAPI, Request, status
from fastapi.exceptions import RequestValidationError
from fastapi.responses import JSONResponse
from starlette.exceptions import HTTPException as StarletteHTTPException

from app.services.fasta import FastaValidationError

logger = logging.getLogger(__name__)


class APIError(Exception):
    """An error with an explicit HTTP status, raised from the API layer."""

    def __init__(
        self, status_code: int, code: str, message: str, *, field: str | None = None
    ) -> None:
        super().__init__(message)
        self.status_code = status_code
        self.code = code
        self.message = message
        self.field = field


def error_response(
    status_code: int, code: str, message: str, field: str | None = None
) -> JSONResponse:
    body: dict = {"code": code, "message": message}
    if field is not None:
        body["field"] = field
    return JSONResponse(status_code=status_code, content={"error": body})


_HTTP_CODES = {
    status.HTTP_404_NOT_FOUND: "not_found",
    status.HTTP_405_METHOD_NOT_ALLOWED: "method_not_allowed",
}


def register_exception_handlers(app: FastAPI) -> None:
    @app.exception_handler(APIError)
    async def handle_api_error(_: Request, exc: APIError) -> JSONResponse:
        return error_response(exc.status_code, exc.code, exc.message, exc.field)

    @app.exception_handler(FastaValidationError)
    async def handle_fasta_error(_: Request, exc: FastaValidationError) -> JSONResponse:
        return error_response(status.HTTP_422_UNPROCESSABLE_CONTENT, exc.code, exc.message, "file")

    @app.exception_handler(RequestValidationError)
    async def handle_request_validation(_: Request, exc: RequestValidationError) -> JSONResponse:
        first = exc.errors()[0]
        field = str(first["loc"][-1]) if first.get("loc") else None
        if field == "file":
            # Typically `-F file=ACGT` (a text value) instead of `-F file=@seq.fasta`.
            message = "Field 'file' must be an uploaded file, e.g. curl -F \"file=@seq.fasta\"."
        elif field:
            message = f"Invalid value for '{field}': {first['msg']}."
        else:
            message = first["msg"]
        return error_response(
            status.HTTP_422_UNPROCESSABLE_CONTENT, "invalid_request", message, field
        )

    @app.exception_handler(StarletteHTTPException)
    async def handle_http_error(_: Request, exc: StarletteHTTPException) -> JSONResponse:
        code = _HTTP_CODES.get(exc.status_code, "http_error")
        return error_response(exc.status_code, code, str(exc.detail))

    @app.exception_handler(Exception)
    async def handle_unexpected(_: Request, exc: Exception) -> JSONResponse:
        logger.exception("Unhandled error", exc_info=exc)
        return error_response(
            status.HTTP_500_INTERNAL_SERVER_ERROR,
            "internal_error",
            "An unexpected error occurred while processing the request.",
        )
