#!/usr/bin/env python3
"""
Download the Streptococcus pyogenes MLST scheme and allele FASTA files
from PubMLST.

Outputs:
  outdir/
    profiles.csv
    scheme.json
    loci.txt
    alleles/
      gki.fasta
      gtr.fasta
      murI.fasta
      mutS.fasta
      recP.fasta
      xpt.fasta
      yqiL.fasta
      ...

Usage:
  python download_spyogenes_pubmlst.py
  python download_spyogenes_pubmlst.py --outdir spyogenes_mlst
"""

from __future__ import annotations

import argparse
import json
import sys
from pathlib import Path
from typing import Any
from urllib.error import HTTPError, URLError
from urllib.parse import urlparse
from urllib.request import Request, urlopen
import shutil
import subprocess
import os
from dotenv import load_dotenv

load_dotenv()
API_KEY = os.getenv("PUBMLST_API_KEY")
print(API_KEY)


BASE = "https://rest.pubmlst.org"
DB = "pubmlst_spyogenes_seqdef"
SCHEME_ID = 1
OUTDIR = "spyogenes_mlst"


def get(url: str, accept: str | None = None) -> bytes:
    headers = {"User-Agent": "pubmlst-spyogenes-downloader/1.0"}

    if accept:
        headers["Accept"] = accept

    # Add API key if available
    if API_KEY:
        headers["Authorization"] = f"Bearer {API_KEY}"

    req = Request(url, headers=headers)
    with urlopen(req) as resp:
        return resp.read()


def get_json(url: str) -> Any:
    return json.loads(get(url, accept="application/json").decode("utf-8"))
"""
def get(url: str, accept: str | None = None) -> bytes:
    headers = {"User-Agent": "pubmlst-spyogenes-downloader/1.0"}
    if accept:
        headers["Accept"] = accept
    req = Request(url, headers=headers)
    with urlopen(req) as resp:
        return resp.read()


def get_json(url: str) -> Any:
    data = get(url, accept="application/json")
    return json.loads(data.decode("utf-8"))

"""
def filename_from_locus_url(locus_url: str) -> str:
    path = urlparse(locus_url).path.rstrip("/")
    return path.split("/")[-1]



def make_spyogenes_allele_blastdb(
    spyogenes_dir: str | Path,
    db_name: str = "spyogenes_alleles",
    combined_fasta_name: str = "spyogenes_alleles.fasta",
) -> Path:
    """
    Concatenate all allele FASTA files and build a BLAST database.

    Expected layout:
      spyogenes_dir/
        alleles/
          *.fasta
    """
    spyogenes_dir = Path(spyogenes_dir)
    alleles_dir = spyogenes_dir / "alleles"

    if not alleles_dir.is_dir():
        raise FileNotFoundError(f"Alleles directory not found: {alleles_dir}")

    fasta_files = sorted(alleles_dir.glob("*.fasta"))
    if not fasta_files:
        raise FileNotFoundError(f"No FASTA files found in: {alleles_dir}")

    if shutil.which("makeblastdb") is None:
        raise RuntimeError(
            "makeblastdb not found in PATH. Install BLAST+ and try again."
        )

    combined_fasta = spyogenes_dir / combined_fasta_name

    # Simple concatenation
    with combined_fasta.open("wb") as out_fh:
        for fasta_file in fasta_files:
            with fasta_file.open("rb") as in_fh:
                shutil.copyfileobj(in_fh, out_fh)

            # Ensure newline separation between files
            out_fh.write(b"\n")

    db_prefix = spyogenes_dir / db_name
    cmd = [
        "makeblastdb",
        "-in",
        str(combined_fasta),
        "-dbtype",
        "nucl",
        "-parse_seqids",
        "-out",
        str(db_prefix),
        "-title",
        db_name,
    ]

    result = subprocess.run(cmd, capture_output=True, text=True)

    if result.returncode != 0:
        raise RuntimeError(
            "makeblastdb failed.\n"
            f"Command: {' '.join(cmd)}\n"
            f"stdout:\n{result.stdout}\n"
            f"stderr:\n{result.stderr}"
        )

    return combined_fasta


def update_mlst_db(output_folder: Path = OUTDIR) -> int:

    outdir = Path(output_folder)
    alleles_dir = outdir / "alleles"
    outdir.mkdir(parents=True, exist_ok=True)
    alleles_dir.mkdir(parents=True, exist_ok=True)

    scheme_url = f"{BASE}/db/{DB}/schemes/{SCHEME_ID}"
    profiles_url = f"{scheme_url}/profiles_csv"

    try:
        print(f"[1/4] Fetching scheme metadata: {scheme_url}", file=sys.stderr)
        scheme = get_json(scheme_url)

        print(f"[2/4] Saving scheme metadata to {outdir / 'scheme.json'}", file=sys.stderr)
        with open(outdir / "scheme.json", "w", encoding="utf-8") as fh:
            json.dump(scheme, fh, indent=2, sort_keys=True)

        print(f"[3/4] Downloading profiles: {profiles_url}", file=sys.stderr)
        profiles_csv = get(profiles_url, accept="text/csv")
        with open(outdir / "profiles.csv", "wb") as fh:
            fh.write(profiles_csv)

        loci = scheme.get("loci", [])
        if not loci:
            raise RuntimeError("No loci found in scheme JSON.")

        locus_names = [filename_from_locus_url(url) for url in loci]
        with open(outdir / "loci.txt", "w", encoding="utf-8") as fh:
            for name in locus_names:
                fh.write(name + "\n")

        print(f"[4/4] Downloading {len(loci)} allele FASTA files", file=sys.stderr)
        for i, locus_url in enumerate(loci, start=1):
            locus = filename_from_locus_url(locus_url)
            fasta_url = f"{locus_url}/alleles_fasta"
            outfile = alleles_dir / f"{locus}.fasta"
            print(f"  - ({i}/{len(loci)}) {locus}: {fasta_url}", file=sys.stderr)
            fasta = get(fasta_url, accept="text/plain")
            with open(outfile, "wb") as fh:
                fh.write(fasta)

        print(f"Done. Files written to: {outdir}", file=sys.stderr)
        return 0

    except HTTPError as e:
        print(f"HTTP error: {e.code} {e.reason}", file=sys.stderr)
        print(
            "PubMLST may require authenticated access for some newer records. "
            "If this fails or looks incomplete, check your PubMLST access setup.",
            file=sys.stderr,
        )
        return 1
    except URLError as e:
        print(f"URL error: {e}", file=sys.stderr)
        return 1
    except Exception as e:
        print(f"Error: {e}", file=sys.stderr)
        return 1

def update_mlst_db_cli() -> int:
    parser = argparse.ArgumentParser()
    parser.add_argument("--outdir", default=OUTDIR, help="Output directory")
    args = parser.parse_args()
    update_mlst_db(output_folder = args.outdir)

update_mlst_db("mlst_test_2")