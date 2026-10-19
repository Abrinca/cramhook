"""Shared helpers for the test-data scripts (manifest, sha256, unpacking)."""
import hashlib
import json
import shutil
import tarfile
from pathlib import Path

ROOT = Path(__file__).resolve().parents[2]
TESTDATA = ROOT / "testdata"
MANIFEST = TESTDATA / "manifest.json"
DATA = TESTDATA / "data"          # git-ignored: unpacked tarball + check outputs
EXPECTED = TESTDATA / "expected"  # in git: expected small text outputs


def load_manifest():
    return json.loads(MANIFEST.read_text())


def save_manifest(man):
    MANIFEST.write_text(json.dumps(man, indent=2) + "\n")


def sha256_file(path):
    h = hashlib.sha256()
    with open(path, "rb") as fh:
        for chunk in iter(lambda: fh.read(1 << 20), b""):
            h.update(chunk)
    return h.hexdigest()


def data_dir(man=None):
    """The unpacked test data: ``testdata/data/<name>/``."""
    man = man or load_manifest()
    return DATA / man["name"]


def is_unpacked(man=None):
    man = man or load_manifest()
    marker = data_dir(man) / ".sha256"
    return marker.exists() and marker.read_text().strip() == man["sha256"]


def unpack(tarball, sha256):
    """Verify ``tarball`` and unpack it into ``testdata/data/``."""
    got = sha256_file(tarball)
    if got != sha256:
        raise SystemExit(f"sha256 mismatch for {tarball}:\n  got      {got}\n"
                         f"  expected {sha256} (testdata/manifest.json)")
    man = load_manifest()
    dst = data_dir(man)
    if dst.exists():
        shutil.rmtree(dst)
    DATA.mkdir(parents=True, exist_ok=True)
    with tarfile.open(tarball, "r:gz") as tf:
        tf.extractall(DATA, filter="data")
    if not dst.is_dir():
        raise SystemExit(f"{tarball} did not contain {man['name']}/")
    (dst / ".sha256").write_text(sha256 + "\n")
    print(f"unpacked {tarball} -> {dst.relative_to(ROOT)}")
