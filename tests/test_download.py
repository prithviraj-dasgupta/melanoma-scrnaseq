import pytest
import requests

from melscrna.download import GEO_URL, download, md5sum


def test_existing_file_is_not_downloaded_again(tmp_path, monkeypatch):
    dest = tmp_path / "data.txt.gz"
    dest.write_bytes(b"already here")

    def fail(*args, **kwargs):
        raise AssertionError("should not download")

    monkeypatch.setattr(requests, "get", fail)
    assert download(GEO_URL, dest) == dest
    assert dest.read_bytes() == b"already here"


def test_md5sum(tmp_path):
    path = tmp_path / "x.txt"
    path.write_bytes(b"abc")
    assert md5sum(path) == "900150983cd24fb0d6963f7d28e17f72"


@pytest.mark.network
def test_geo_file_is_reachable():
    response = requests.head(GEO_URL, allow_redirects=True, timeout=30)
    assert response.status_code == 200
