#!/usr/bin/env python3
"""Download the public test-data tarball and unpack it into ``testdata/data/``.

The URL and sha256 come from ``testdata/manifest.json`` (the one place they are
defined).  Nothing is downloaded when the unpacked data already match the
manifest sha256.

Usage::

  pixi run fetch-testdata [--force] [--from PATH]   # --from: a local tarball
"""
import argparse
import sys
import urllib.request
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parent))
from common import DATA, ROOT, is_unpacked, load_manifest, sha256_file, unpack  # noqa: E402


def main(argv=None):
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--force", action="store_true", help="download and unpack again")
    ap.add_argument("--from", dest="src", help="use a local tarball instead of the URL")
    args = ap.parse_args(argv)

    man = load_manifest()
    if not man.get("sha256"):
        sys.exit("testdata/manifest.json has no sha256: run `pixi run build-testdata`")
    if is_unpacked(man) and not args.force:
        print(f"{man['name']} already in {DATA.relative_to(ROOT)}/ (sha256 ok)")
        return 0

    DATA.mkdir(parents=True, exist_ok=True)
    if args.src:
        tarball = Path(args.src)
    else:
        tarball = DATA / f"{man['name']}.tar.gz"
        if not (tarball.exists() and sha256_file(tarball) == man["sha256"]):
            if "TODO" in man["url"]:
                sys.exit(f"the download URL is still a placeholder ({man['url']}); "
                         f"use `--from dist/{man['name']}.tar.gz` or `pixi run build-testdata`")
            part = tarball.with_name(tarball.name + ".part")
            print(f"downloading {man['url']}")
            urllib.request.urlretrieve(man["url"], part)
            part.replace(tarball)
    unpack(tarball, man["sha256"])
    return 0


if __name__ == "__main__":
    sys.exit(main())
