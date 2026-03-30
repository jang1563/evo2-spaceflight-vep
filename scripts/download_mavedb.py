#!/usr/bin/env python3
"""
Download DMS datasets from MaveDB API.

Downloads score sets for our 4 control genes + TERT MPRA:
  - BRCA1: urn:mavedb:00000097-0-2 (Findlay 2018 SGE)
  - TP53: urn:mavedb:00000068-b-1 (Giacomelli 2018)
  - TERT MPRA: urn:mavedb:00000031-b-1 (Kircher 2019 GBM)
  - CHEK2: urn:mavedb:00001203-a-1 (McCarthy-Leo 2024)
  - DNMT3A: Garcia et al. 2025 (bioRxiv 10.1101/2025.09.24.678339)

Usage: python download_mavedb.py
"""

import json
import os
from pathlib import Path

import requests

import sys
sys.path.insert(0, str(Path(__file__).parent))

try:
    from utils.config import DMS_DIR, MPRA_DIR
except ImportError:
    _fallback = Path(os.environ.get("EVO2_ROOT", "/path/to/your/scratch/evo2")) / "project_spaceflight_vep" / "data"
    DMS_DIR = _fallback / "dms"
    MPRA_DIR = _fallback / "mpra"

MAVEDB_API = "https://api.mavedb.org/api/v1"

# Score set URNs — output_dir is now an absolute Path
SCORE_SETS = {
    "brca1_findlay2018": {
        "urn": "urn:mavedb:00000097-0-2",
        "gene": "BRCA1",
        "reference": "Findlay 2018 (SGE)",
        "output_dir": DMS_DIR / "brca1",
        "output_file": "findlay2018_scores.csv",
    },
    "tp53_giacomelli2018": {
        "urn": "urn:mavedb:00000068-b-1",
        "gene": "TP53",
        "reference": "Giacomelli 2018 (nutlin-3 growth)",
        "output_dir": DMS_DIR / "tp53",
        "output_file": "giacomelli2018_scores.csv",
    },
    "tert_kircher2019": {
        "urn": "urn:mavedb:00000031-b-1",
        "gene": "TERT",
        "reference": "Kircher 2019 (TERT GBM promoter MPRA)",
        "output_dir": MPRA_DIR,
        "output_file": "tert_mpra_scores.csv",
    },
}


def download_score_set(name: str, info: dict):
    """Download a single score set from MaveDB."""
    urn = info["urn"]
    output_dir = Path(info["output_dir"])
    output_file = output_dir / info["output_file"]

    output_dir.mkdir(parents=True, exist_ok=True)

    if output_file.exists() and output_file.stat().st_size > 100:
        print(f"  {name}: already downloaded ({output_file.stat().st_size} bytes)")
        return True

    print(f"  {name}: downloading from MaveDB ({urn})...")

    # MaveDB API v1 endpoint for scores
    url = f"{MAVEDB_API}/score-sets/{urn}/scores"

    try:
        resp = requests.get(url, timeout=120)
        resp.raise_for_status()

        # Try JSON first (v1 API returns JSON)
        try:
            data = resp.json()
            # Save as JSON
            json_file = output_dir / info["output_file"].replace(".csv", ".json")
            with open(json_file, 'w') as f:
                json.dump(data, f, indent=2)
            print(f"    Saved JSON: {json_file} ({len(data)} entries)")

            # Also convert to CSV for easier use
            if isinstance(data, list) and len(data) > 0:
                import csv
                keys = data[0].keys()
                with open(output_file, 'w', newline='') as f:
                    writer = csv.DictWriter(f, fieldnames=keys)
                    writer.writeheader()
                    writer.writerows(data)
                print(f"    Saved CSV: {output_file}")
            return True
        except json.JSONDecodeError:
            # Maybe it's plain CSV
            with open(output_file, 'wb') as f:
                f.write(resp.content)
            print(f"    Saved raw: {output_file}")
            return True

    except requests.exceptions.HTTPError as e:
        print(f"    HTTP error: {e}")
        print(f"    Try manual download from https://www.mavedb.org/#/score-sets/{urn}")
        return False
    except requests.exceptions.ConnectionError as e:
        print(f"    Connection error: {e}")
        return False
    except Exception as e:
        print(f"    Error: {e}")
        return False


def search_mavedb(query: str):
    """Search MaveDB for score sets matching a query."""
    url = f"{MAVEDB_API}/score-sets/"
    params = {"search": query, "limit": 10}

    try:
        resp = requests.get(url, params=params, timeout=30)
        resp.raise_for_status()
        data = resp.json()

        if isinstance(data, dict) and "results" in data:
            results = data["results"]
        elif isinstance(data, list):
            results = data
        else:
            print(f"  Unexpected response format for query '{query}'")
            return

        print(f"\nMaveDB search results for '{query}':")
        for item in results[:5]:
            urn = item.get("urn", "?")
            title = item.get("title", "?")
            print(f"  {urn}: {title[:80]}")

    except Exception as e:
        print(f"  Search error: {e}")


def main():
    print("=" * 60)
    print("MaveDB DMS Data Download")
    print("=" * 60)

    # Download known score sets
    for name, info in SCORE_SETS.items():
        print(f"\n{info['gene']} — {info['reference']}")
        download_score_set(name, info)

    # Search for CHEK2 and DNMT3A (not sure of exact URNs)
    print("\n" + "=" * 60)
    print("Searching MaveDB for additional datasets...")
    print("=" * 60)

    search_mavedb("CHEK2")
    search_mavedb("DNMT3A")
    search_mavedb("DNMT3A methylation")

    print("\n" + "=" * 60)
    print("Manual verification needed:")
    print("  1. CHEK2 (McCarthy-Leo 2024): Check paper supplementary data")
    print("  2. DNMT3A (Garcia et al. 2025): Download from bioRxiv 10.1101/2025.09.24.678339 supplementary")
    print("=" * 60)


if __name__ == "__main__":
    main()
