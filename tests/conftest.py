from pathlib import Path

import pytest

ROOT = Path(__file__).resolve().parent.parent
SAMPLES_DIR = ROOT / "samples"
REFERENCE_PATH = ROOT / "data" / "reference.fasta"


@pytest.fixture(scope="session")
def reference_path() -> Path:
    return REFERENCE_PATH


@pytest.fixture(scope="session")
def samples_dir() -> Path:
    return SAMPLES_DIR
