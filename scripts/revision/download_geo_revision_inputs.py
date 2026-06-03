from __future__ import annotations

import argparse
import os
import shutil
import time
import urllib.request
from pathlib import Path


FILES = {
    "GSE149614/GSE149614_HCC.scRNAseq.S71915.count.txt.gz": (
        "https://ftp.ncbi.nlm.nih.gov/geo/series/GSE149nnn/GSE149614/suppl/"
        "GSE149614_HCC.scRNAseq.S71915.count.txt.gz"
    ),
    "GSE149614/GSE149614_HCC.metadata.updated.txt.gz": (
        "https://ftp.ncbi.nlm.nih.gov/geo/series/GSE149nnn/GSE149614/suppl/"
        "GSE149614_HCC.metadata.updated.txt.gz"
    ),
    "GSE238264/GSE238264_RAW.tar": (
        "https://ftp.ncbi.nlm.nih.gov/geo/series/GSE238nnn/GSE238264/suppl/"
        "GSE238264_RAW.tar"
    ),
}


def download(url: str, destination: Path) -> None:
    if destination.exists() and destination.stat().st_size > 0:
        print(f"exists={destination} bytes={destination.stat().st_size}")
        return
    destination.parent.mkdir(parents=True, exist_ok=True)
    temporary = destination.with_suffix(destination.suffix + ".part")
    proxy = urllib.request.ProxyHandler(
        {
            "http": os.environ.get("HTTP_PROXY", "http://127.0.0.1:7890"),
            "https": os.environ.get("HTTPS_PROXY", "http://127.0.0.1:7890"),
        }
    )
    opener = urllib.request.build_opener(proxy)
    request = urllib.request.Request(url, headers={"User-Agent": "CSBJ-revision/1.0"})
    last_error = None
    for attempt in range(1, 6):
        try:
            if temporary.exists():
                temporary.unlink()
            with opener.open(request, timeout=120) as response, temporary.open("wb") as output:
                shutil.copyfileobj(response, output, length=1024 * 1024)
            temporary.replace(destination)
            print(f"downloaded={destination} bytes={destination.stat().st_size}")
            return
        except Exception as exc:
            last_error = exc
            print(f"retry={attempt} destination={destination} error={exc}")
            time.sleep(attempt * 3)
    raise RuntimeError(f"Failed to download {url}") from last_error


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--output-dir", type=Path, required=True)
    args = parser.parse_args()
    for relative, url in FILES.items():
        download(url, args.output_dir / relative)


if __name__ == "__main__":
    main()
