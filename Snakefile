# ─────────────────────────────────────────────────────────────────────────────
# Snakefile  –  Protein Side-Chain Angle Analysis Pipeline
# ─────────────────────────────────────────────────────────────────────────────
# Run with:
#   snakemake --cores <N> --snakefile Snakefile
# ─────────────────────────────────────────────────────────────────────────────

configfile: "config.yaml"

import os
import glob
import math
import sys

# ── Resolve config values ────────────────────────────────────────────────────
PDB_DIR       = config["pdb_dir"]
OUTPUT_DIR    = config["output_dir"]
BATCH_SIZE    = int(config.get("batch_size", 500))
N_WORKERS     = int(config.get("n_workers_per_batch", 4))
DSSP_EXE      = config.get("dssp_executable", "mkdssp")
PLOT_XMIN     = config.get("plot_xmin", -180)
PLOT_XMAX     = config.get("plot_xmax",  180)

BATCHES_DIR   = os.path.join(OUTPUT_DIR, "batches")
BATCH_CSV_DIR = os.path.join(OUTPUT_DIR, "batch_csvs")

# Absolute path to scripts/ — resolves correctly regardless of where
# snakemake is invoked from (workflow.basedir = directory of this Snakefile).
SCRIPTS_DIR = os.path.join(workflow.basedir, "scripts")

# ── Discover PDB files ───────────────────────────────────────────────────────
pdb_files = sorted(set(
    glob.glob(os.path.join(PDB_DIR, "**", "*.pdb"), recursive=True) +
    glob.glob(os.path.join(PDB_DIR, "*.pdb"))
))

if not pdb_files:
    sys.exit(
        f"\n[ERROR] No .pdb files found under '{PDB_DIR}'.\n"
        "Please update 'pdb_dir' in config.yaml and retry.\n"
    )

print(f"[Snakefile] Discovered {len(pdb_files):,} PDB files.", file=sys.stderr)

# ── Create batch manifests at workflow parse time ────────────────────────────
# (lightweight file I/O; idempotent on re-runs)
os.makedirs(BATCHES_DIR,   exist_ok=True)
os.makedirs(BATCH_CSV_DIR, exist_ok=True)

n_batches = max(1, math.ceil(len(pdb_files) / BATCH_SIZE))
BATCHES   = [f"{i:05d}" for i in range(n_batches)]

for idx, batch_id in enumerate(BATCHES):
    start = idx * BATCH_SIZE
    end   = min(start + BATCH_SIZE, len(pdb_files))
    manifest_path = os.path.join(BATCHES_DIR, f"{batch_id}.txt")
    with open(manifest_path, "w") as fh:
        fh.write("\n".join(pdb_files[start:end]) + "\n")

print(f"[Snakefile] Created {n_batches} batch manifests "
      f"(batch_size={BATCH_SIZE}).", file=sys.stderr)

# ── Rules ────────────────────────────────────────────────────────────────────

rule all:
    """Top-level target: final plot + aggregated CSV."""
    input:
        plot = os.path.join(OUTPUT_DIR, "plot.png"),
        csv  = os.path.join(OUTPUT_DIR, "angles.csv"),


rule process_batch:
    """
    Parse one batch of PDB files, extract X–Arg(helix)–X tripeptides,
    compute signed CA→centroid angles, and write results to a CSV.
    """
    input:
        manifest = os.path.join(BATCHES_DIR, "{batch_id}.txt"),
    output:
        csv = os.path.join(BATCH_CSV_DIR, "{batch_id}.csv"),
    threads: N_WORKERS
    params:
        dssp = DSSP_EXE,
    log:
        os.path.join(OUTPUT_DIR, "logs", "parse_{batch_id}.log"),
    shell:
        """
        python {SCRIPTS_DIR}/parse_pdb.py \
            --manifest  {input.manifest} \
            --output    {output.csv} \
            --workers   {threads} \
            --dssp      {params.dssp} \
        2>&1 | tee {log}
        """


rule aggregate:
    """Concatenate all per-batch CSVs into a single angles.csv."""
    input:
        csvs = expand(
            os.path.join(BATCH_CSV_DIR, "{batch_id}.csv"),
            batch_id=BATCHES,
        ),
    output:
        csv = os.path.join(OUTPUT_DIR, "angles.csv"),
    log:
        os.path.join(OUTPUT_DIR, "logs", "aggregate.log"),
    params:
        batch_csv_dir = BATCH_CSV_DIR,
    shell:
        """
        python {SCRIPTS_DIR}/aggregate.py \
            --input-dir {params.batch_csv_dir} \
            --output    {output.csv} \
        2>&1 | tee {log}
        """


rule plot:
    """Generate KDE plot of angle distributions grouped by residue category."""
    input:
        csv = os.path.join(OUTPUT_DIR, "angles.csv"),
    output:
        png = os.path.join(OUTPUT_DIR, "plot.png"),
    log:
        os.path.join(OUTPUT_DIR, "logs", "plot.log"),
    params:
        xmin = PLOT_XMIN,
        xmax = PLOT_XMAX,
    shell:
        """
        python {SCRIPTS_DIR}/plot.py \
            --input  {input.csv} \
            --output {output.png} \
            --xmin   {params.xmin} \
            --xmax   {params.xmax} \
        2>&1 | tee {log}
        """
