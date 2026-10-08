"""Application factory.

Run locally with::

    uvicorn app.main:app --reload
"""

from __future__ import annotations

from fastapi import FastAPI

from app import __version__
from app.core.config import Settings, get_settings
from app.core.errors import register_exception_handlers
from app.routers import align, analyze, health, translate
from app.services.fasta import load_fasta_file

DESCRIPTION = """
DNA sequence analysis against the bundled reference **REF_001** (300 bp synthetic CDS).

* `POST /api/align`: global alignment + variant detection (SNV / insertion / deletion)
* `POST /api/translate`: protein translation in reading frame 1, 2 or 3
* `POST /api/analyze`: alignment + translation + variant effect annotation (bonus)

Upload a single-record FASTA file as `multipart/form-data` in the field `file`.
"""


def create_app(settings: Settings | None = None) -> FastAPI:
    settings = settings or get_settings()

    app = FastAPI(title="DNA Analysis API", version=__version__, description=DESCRIPTION)
    app.state.settings = settings
    # Loaded once at start-up; a missing or invalid reference fails fast here.
    app.state.reference = load_fasta_file(settings.reference_path)

    register_exception_handlers(app)
    app.include_router(health.router)
    app.include_router(align.router)
    app.include_router(translate.router)
    app.include_router(analyze.router)
    return app


app = create_app()
