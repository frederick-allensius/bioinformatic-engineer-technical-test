from pathlib import Path

import pytest
from fastapi.testclient import TestClient

from app.main import create_app
from app.services.fasta import load_fasta_file

ROOT = Path(__file__).resolve().parent.parent
SAMPLES_DIR = ROOT / "samples"
REFERENCE_PATH = ROOT / "data" / "reference.fasta"


@pytest.fixture(scope="session")
def reference_path() -> Path:
    return REFERENCE_PATH


@pytest.fixture(scope="session")
def reference_seq() -> str:
    return load_fasta_file(REFERENCE_PATH).sequence


@pytest.fixture(scope="session")
def samples_dir() -> Path:
    return SAMPLES_DIR


@pytest.fixture(scope="session")
def client() -> TestClient:
    return TestClient(create_app())


@pytest.fixture
def post_fasta(client):
    """POST raw FASTA text/bytes as the multipart ``file`` field.

    ``frame`` is sent as a form field when given; ``filename`` is cosmetic.
    """

    def _post(endpoint: str, content: str | bytes, *, frame=None, filename="query.fasta"):
        data = {} if frame is None else {"frame": str(frame)}
        payload = content.encode() if isinstance(content, str) else content
        return client.post(endpoint, files={"file": (filename, payload)}, data=data)

    return _post


@pytest.fixture
def post_sample(post_fasta, samples_dir):
    def _post(endpoint: str, name: str, **kwargs):
        return post_fasta(endpoint, (samples_dir / f"{name}.fasta").read_bytes(), **kwargs)

    return _post
