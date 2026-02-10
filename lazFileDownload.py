#!/usr/bin/env python3
import requests
import sys

# ─── USER PARAMETERS ──────────────────────────────────────────────────────────
# Your exact WGS84 bbox:
MIN_LON = -83.35182782546514
MIN_LAT =  33.922604972225635
MAX_LON = -83.35013877033196
MAX_LAT =  33.92044317109077

# How much to expand the bbox (in decimal degrees; ~0.002° ≈ 200 m)
BUFFER = 0.002

# The dataset you want (case-sensitive). 
# You can run with None to list all and pick one.
DATASET_NAME = "3DEP Lidar Point Cloud"

# Output filename
OUTPUT_FILE = "output.laz"
# ──────────────────────────────────────────────────────────────────────────────

def list_datasets():
    """Fetch and print all available TNMAccess datasets."""
    url = "https://tnmaccess.nationalmap.gov/api/v1/datasets"
    resp = requests.get(url, params={"outputFormat": "json"})
    resp.raise_for_status()
    names = [d["datasetName"] for d in resp.json()]
    print("Available datasets:")
    for n in names:
        print(" •", n)
    sys.exit(0)

def download_laz(min_lon, min_lat, max_lon, max_lat, dataset, out_file):
    """Query /products with buffered bbox and download the first LAZ tile."""
    # 1) Buffer the bbox
    min_lon -= BUFFER; max_lon += BUFFER
    min_lat -= BUFFER; max_lat += BUFFER
    print(f"Using buffered bbox: {min_lon:.6f},{min_lat:.6f},{max_lon:.6f},{max_lat:.6f}")

    # 2) Query TNMAccess for products
    url = "https://tnmaccess.nationalmap.gov/api/v1/products"
    params = {
        "datasets": dataset,
        "bbox": f"{min_lon},{min_lat},{max_lon},{max_lat}",
        "outputFormat": "json"
    }
    resp = requests.get(url, params=params)
    resp.raise_for_status()
    items = resp.json().get("items", [])

    if not items:
        raise ValueError(
            "No LAZ files found for this buffered bbox. "
            "Try increasing BUFFER or checking dataset name."
        )

    # 3) Grab the first downloadURL
    download_url = items[0]["downloadURL"]
    print("Found tile:", download_url)

    # 4) Stream-download to disk
    r = requests.get(download_url, stream=True)
    r.raise_for_status()
    with open(out_file, "wb") as f:
        for chunk in r.iter_content(chunk_size=8192):
            f.write(chunk)

    print(f"✅ Downloaded {out_file}")

if __name__ == "__main__":
    if DATASET_NAME is None:
        list_datasets()

    try:
        download_laz(
            MIN_LON, MIN_LAT,
            MAX_LON, MAX_LAT,
            DATASET_NAME,
            OUTPUT_FILE
        )
    except Exception as e:
        print("Error:", e)
        print("If you suspect the dataset name is wrong, re-run with DATASET_NAME = None to list all.") 
        sys.exit(1)
