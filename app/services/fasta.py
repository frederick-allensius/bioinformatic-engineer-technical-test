"""FASTA parsing and validation shared by every endpoint.

The rules follow the API specification:

* exactly one record with a non-empty header;
* 1 to ``max_length`` nucleotides once line breaks and surrounding whitespace are removed;
* only A, C, G, T and N (case-insensitive); the sequence is upper-cased;
* multiline sequences and blank lines are allowed, whitespace *inside* a sequence line is not.

This module knows nothing about HTTP. It raises :class:`FastaValidationError` with a
machine-readable ``code`` and a human-readable ``message``; the API layer decides which
status code to return.
"""

from __future__ import annotations

import re
from dataclasses import dataclass
from pathlib import Path

VALID_NUCLEOTIDES = frozenset("ACGTN")
DEFAULT_MAX_LENGTH = 1000

# Only real line terminators. str.splitlines() would also split on \v, \f, \x1c, ...
# which would let "embedded" whitespace slip through as a line break.
_LINE_BREAK = re.compile(r"\r\n|\r|\n")


class FastaValidationError(ValueError):
    """Raised when input does not satisfy the FASTA rules."""

    def __init__(self, code: str, message: str) -> None:
        super().__init__(message)
        self.code = code
        self.message = message


@dataclass(frozen=True, slots=True)
class FastaRecord:
    sequence_id: str
    description: str
    sequence: str

    def __len__(self) -> int:
        return len(self.sequence)


def decode_fasta_bytes(data: bytes) -> str:
    """Decode an uploaded file as UTF-8 (a leading BOM is tolerated)."""
    try:
        return data.decode("utf-8-sig")
    except UnicodeDecodeError as exc:
        raise FastaValidationError(
            "invalid_encoding",
            "File is not valid UTF-8 text; upload a plain-text FASTA file.",
        ) from exc


def parse_fasta(text: str, *, max_length: int | None = DEFAULT_MAX_LENGTH) -> FastaRecord:
    """Parse and validate a single-record FASTA document.

    Args:
        text: Full file content.
        max_length: Maximum number of nucleotides, or ``None`` for no limit
            (used when loading the trusted reference).

    Returns:
        The validated record with an upper-cased sequence.

    Raises:
        FastaValidationError: if any rule is violated. Errors point at the
            offending line/column where that helps the user fix the file.
    """
    if not text.strip():
        raise FastaValidationError("empty_file", "The uploaded file is empty.")

    # Keep 1-based line numbers so error messages match what a text editor shows.
    lines = [(number, line.strip()) for number, line in enumerate(_LINE_BREAK.split(text), 1)]
    lines = [(number, line) for number, line in lines if line]

    header_line_no, header = lines[0]
    if not header.startswith(">"):
        raise FastaValidationError(
            "malformed_fasta",
            f"Expected a FASTA header starting with '>' on line {header_line_no}, "
            f"found sequence data instead.",
        )

    extra_headers = [number for number, line in lines[1:] if line.startswith(">")]
    if extra_headers:
        raise FastaValidationError(
            "multiple_records",
            f"Expected exactly one FASTA record but found {len(extra_headers) + 1} "
            f"(additional header on line {extra_headers[0]}).",
        )

    description = header[1:].strip()
    if not description:
        raise FastaValidationError(
            "empty_header",
            f"FASTA header on line {header_line_no} is empty; expected an identifier after '>'.",
        )
    sequence_id = description.split()[0]

    chunks: list[str] = []
    for number, line in lines[1:]:
        _check_sequence_line(line, number)
        chunks.append(line)
    sequence = "".join(chunks).upper()

    if not sequence:
        raise FastaValidationError(
            "empty_sequence", f"Record '{sequence_id}' has a header but no sequence."
        )
    if max_length is not None and len(sequence) > max_length:
        raise FastaValidationError(
            "sequence_too_long",
            f"Sequence length is {len(sequence)} nt; the maximum allowed is {max_length} nt.",
        )

    return FastaRecord(sequence_id=sequence_id, description=description, sequence=sequence)


def load_fasta_file(path: str | Path) -> FastaRecord:
    """Load a trusted single-record FASTA file from disk (no length limit)."""
    return parse_fasta(Path(path).read_text(encoding="utf-8-sig"), max_length=None)


def _check_sequence_line(line: str, line_no: int) -> None:
    for column, char in enumerate(line, 1):
        if char.upper() in VALID_NUCLEOTIDES:
            continue
        if char.isspace():
            raise FastaValidationError(
                "embedded_whitespace",
                f"Whitespace inside the sequence at line {line_no}, column {column}.",
            )
        raise FastaValidationError(
            "invalid_character",
            f"Invalid character {char!r} at line {line_no}, column {column}; "
            f"only A, C, G, T and N are allowed.",
        )
