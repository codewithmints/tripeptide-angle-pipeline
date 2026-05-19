import argparse
import glob
import logging
import os
import sys

import pandas as pd

logging.basicConfig(
    level=logging.INFO,
    format="%(asctime)s [%(levelname)s] %(message)s",
    datefmt="%H:%M:%S",
)
logger = logging.getLogger(__name__)

# Expected schema – used to create an empty DataFrame when no data are found
CSV_COLUMNS = ["pdb_id", "left_res", "central_res", "right_res", "angle", "category"]

# Canonical category ordering (used for a sorted summary table)
CATEGORY_ORDER = ["Tiny", "Small", "Intermediate", "Large", "Bulky"]


def parse_args() -> argparse.Namespace:
    p = argparse.ArgumentParser(
        description="Aggregate per-batch angle CSVs into a single file."
    )
    p.add_argument(
        "--input-dir", required=True,
        help="Directory that contains the per-batch *.csv files.",
    )
    p.add_argument(
        "--output", required=True,
        help="Path for the merged output CSV.",
    )
    return p.parse_args()


def main() -> None:
    args = parse_args()

    # ── Discover batch CSVs ───────────────────────────────────────────────
    pattern   = os.path.join(args.input_dir, "*.csv")
    csv_files = sorted(glob.glob(pattern))

    if not csv_files:
        logger.error("No CSV files found in '%s' (pattern: %s).", args.input_dir, pattern)
        sys.exit(1)

    logger.info("Found %d batch CSV files to aggregate.", len(csv_files))

    # ── Read and concatenate ──────────────────────────────────────────────
    frames: list[pd.DataFrame] = []

    for path in csv_files:
        try:
            df = pd.read_csv(path, dtype={"pdb_id": str, "left_res": str,
                                           "central_res": str, "right_res": str,
                                           "category": str})
            if not df.empty:
                frames.append(df)
        except pd.errors.EmptyDataError:
            logger.debug("Empty file (skipped): %s", path)
        except Exception as exc:
            logger.warning("Could not read '%s': %s – skipping.", path, exc)

    if frames:
        result = pd.concat(frames, ignore_index=True)
    else:
        logger.warning("All batch files were empty – writing empty output.")
        result = pd.DataFrame(columns=CSV_COLUMNS)

    # ── Validate angle column ─────────────────────────────────────────────
    original_len = len(result)
    result = result.dropna(subset=["angle"])
    dropped = original_len - len(result)
    if dropped:
        logger.warning("Dropped %d rows with NaN angles.", dropped)

    # ── Write output ──────────────────────────────────────────────────────
    os.makedirs(os.path.dirname(os.path.abspath(args.output)), exist_ok=True)
    result.to_csv(args.output, index=False)
    logger.info("Aggregated %d records → %s", len(result), args.output)

    # ── Summary table ─────────────────────────────────────────────────────
    if not result.empty and "category" in result.columns:
        counts = (
            result.groupby("category")["angle"]
            .agg(count="count", mean="mean", std="std", median="median")
            .reindex([c for c in CATEGORY_ORDER if c in result["category"].unique()])
        )
        logger.info("\n── Per-category summary ──\n%s\n", counts.to_string())


if __name__ == "__main__":
    main()
