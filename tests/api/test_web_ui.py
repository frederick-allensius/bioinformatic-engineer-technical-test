"""The bonus web interface is served by the same app."""


def test_index_page(client):
    response = client.get("/")
    assert response.status_code == 200
    assert response.headers["content-type"].startswith("text/html")
    assert "DNA Analysis" in response.text


def test_static_assets(client):
    for path in ["/static/app.js", "/static/app.css", "/static/design-system.css"]:
        assert client.get(path).status_code == 200, path


def test_samples_are_downloadable(client):
    response = client.get("/samples/patient_02.fasta")
    assert response.status_code == 200
    assert response.text.startswith(">patient_02")


def test_unknown_path_uses_json_error_envelope(client):
    response = client.get("/api/does-not-exist")
    assert response.status_code == 404
    assert response.json()["error"]["code"] == "not_found"


def test_health(client):
    assert client.get("/api/health").json() == {
        "status": "ok",
        "version": "1.0.0",
        "reference_id": "REF_001",
        "reference_length": 300,
    }
