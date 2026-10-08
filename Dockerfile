# syntax=docker/dockerfile:1

# ---- base: runtime dependencies + application code ---------------------------
FROM python:3.12-slim AS base

ENV PYTHONDONTWRITEBYTECODE=1 \
    PYTHONUNBUFFERED=1 \
    PIP_NO_CACHE_DIR=1 \
    PIP_DISABLE_PIP_VERSION_CHECK=1

WORKDIR /app
COPY requirements-runtime.txt .
RUN pip install -r requirements-runtime.txt

COPY app ./app
COPY data ./data
COPY samples ./samples

# ---- test: `docker build --target test -t dna-api-test . && docker run --rm dna-api-test`
FROM base AS test
COPY requirements.txt pyproject.toml ./
RUN pip install -r requirements.txt
COPY tests ./tests
CMD ["pytest", "-q"]

# ---- runtime (default target) ------------------------------------------------
FROM base AS runtime
RUN useradd --create-home --uid 10001 appuser
USER appuser
EXPOSE 8000
HEALTHCHECK --interval=30s --timeout=3s --start-period=5s --retries=3 \
    CMD ["python", "-c", "import urllib.request; urllib.request.urlopen('http://127.0.0.1:8000/api/health', timeout=2)"]
CMD ["uvicorn", "app.main:app", "--host", "0.0.0.0", "--port", "8000"]
