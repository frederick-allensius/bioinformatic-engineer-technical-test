"""Runtime configuration.

Defaults follow the specification. Each value can be overridden with an environment
variable (prefix ``DNA_API_``), which is handy in Docker, e.g.
``DNA_API_REFERENCE_PATH=/data/other_ref.fasta``.
"""

from __future__ import annotations

import os
from dataclasses import dataclass
from functools import lru_cache
from pathlib import Path

PROJECT_ROOT = Path(__file__).resolve().parents[2]
DEFAULT_REFERENCE_PATH = PROJECT_ROOT / "data" / "reference.fasta"


@dataclass(frozen=True, slots=True)
class Settings:
    reference_path: Path = DEFAULT_REFERENCE_PATH
    max_sequence_length: int = 1000
    # Hard cap on raw upload size, checked before parsing so a huge upload is never
    # fully read into memory. 1 MiB is far above any valid 1,000 nt FASTA file.
    max_upload_bytes: int = 1024 * 1024


@lru_cache
def get_settings() -> Settings:
    return Settings(
        reference_path=Path(os.getenv("DNA_API_REFERENCE_PATH", str(DEFAULT_REFERENCE_PATH))),
        max_sequence_length=int(os.getenv("DNA_API_MAX_SEQUENCE_LENGTH", "1000")),
        max_upload_bytes=int(os.getenv("DNA_API_MAX_UPLOAD_BYTES", str(1024 * 1024))),
    )
