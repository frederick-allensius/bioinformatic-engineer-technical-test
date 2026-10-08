"""The JSON error envelope also covers unexpected errors and generic validation errors."""

from fastapi.testclient import TestClient

from app.main import create_app


def make_client_with_extra_routes() -> TestClient:
    app = create_app()

    @app.get("/boom")
    def boom():
        raise RuntimeError("secret internal detail")

    @app.get("/needs-int")
    def needs_int(count: int):
        return {"count": count}

    return TestClient(app, raise_server_exceptions=False)


def test_unexpected_exception_returns_generic_500():
    response = make_client_with_extra_routes().get("/boom")
    assert response.status_code == 500
    error = response.json()["error"]
    assert error["code"] == "internal_error"
    assert "secret" not in error["message"]  # internals are logged, not leaked


def test_generic_validation_error_names_the_field():
    response = make_client_with_extra_routes().get("/needs-int", params={"count": "x"})
    assert response.status_code == 422
    error = response.json()["error"]
    assert error["code"] == "invalid_request"
    assert error["field"] == "count"
    assert "count" in error["message"]
