# Protein Side-Chain Angle Analysis Pipeline

Snakemake pipeline that analyses **X – Arg(helix) – X** tripeptides in a large
corpus of PDB structures and produces a KDE plot of the signed CA→centroid
angle distribution grouped by left-residue size category.

---

## Directory layout

```
protein_angle_pipeline/
├── Snakefile               # Pipeline orchestration
├── config.yaml             # All user-configurable settings
├── environment.yaml        # Conda environment (Python + DSSP + Snakemake)
├── README.md
└── scripts/
    ├── parse_pdb.py        # Per-batch PDB parsing + angle computation
    ├── aggregate.py        # Merge batch CSVs → single angles.csv
    └── plot.py             # KDE + box-plot saved as plot.png
```

After a run the pipeline writes:

```
results/
├── batches/                # Per-batch PDB path manifests (auto-generated)
├── batch_csvs/             # Per-batch angle CSVs
├── logs/                   # Per-rule log files
├── angles.csv              # Aggregated angle data for all PDB files
└── plot.png                # Final KDE + box-plot figure
```

---

## Prerequisites

### 1 – DSSP binary

DSSP is required for secondary-structure assignment.

| Environment | Install command |
|-------------|-----------------|
| Conda (recommended) | included via `environment.yaml` |
| Ubuntu / Debian | `sudo apt-get install dssp` |
| macOS (Homebrew) | `brew install dssp` |
| Manual | Download from https://swift.cmbi.ru.nl/gv/dssp/ |

Verify with:

```bash
mkdssp --version
# or
dssp --version
```

Update `dssp_executable` in `config.yaml` if the binary is not on your `PATH`
(e.g. `dssp_executable: "/usr/local/bin/mkdssp"`).

### 2 – Python environment

**Option A – Conda (recommended)**

```bash
conda env create -f environment.yaml
conda activate protein_angles
```

**Option B – pip**

```bash
pip install biopython numpy pandas scipy matplotlib seaborn snakemake
# then install DSSP separately (see above)
```

---

## Configuration

Edit **`config.yaml`** before running:

| Key | Description | Default |
|-----|-------------|---------|
| `pdb_dir` | Path to directory containing `.pdb` files (searched recursively) | `/path/to/your/pdb_files` |
| `output_dir` | Where results are written | `results` |
| `batch_size` | PDB files per Snakemake job | `500` |
| `n_workers_per_batch` | Worker processes inside each job | `4` |
| `dssp_executable` | DSSP binary name or path | `mkdssp` |
| `plot_xmin` / `plot_xmax` | KDE x-axis range (degrees) | `-180` / `180` |

---

## ▶ One-line run command

```bash
snakemake --cores all --snakefile Snakefile
```

This is the **single command** that runs the entire pipeline end-to-end.

Snakemake will:
1. Discover all `.pdb` files under `pdb_dir`.
2. Split them into batches and write manifest files.
3. Process all batches in parallel (up to `--cores` limit).
4. Aggregate batch CSVs → `results/angles.csv`.
5. Generate `results/plot.png`.

---

## Performance tuning

For **50 000 PDB files** with the default settings (`batch_size=500`,
`n_workers_per_batch=4`):

| `--cores` | Concurrent jobs | Concurrent workers | ~Wall time† |
|-----------|----------------|--------------------|-------------|
| 8 | 2 | 8 | moderate |
| 16 | 4 | 16 | fast |
| 32 | 8 | 32 | very fast |
| 64 | 16 | 64 | fastest |

† depends heavily on average PDB size and DSSP speed.

**Tune `batch_size` and `n_workers_per_batch` in `config.yaml`:**

- Fewer, larger batches → less Snakemake overhead (recommended for > 50 k files).
- More workers per batch → higher per-job CPU usage.

Total CPU utilisation ≈ `floor(--cores / n_workers_per_batch) × n_workers_per_batch`.

---

## Dry-run (preview without executing)

```bash
snakemake --cores all --dryrun --snakefile Snakefile
```

## Resume after interruption

Snakemake automatically re-uses completed batch CSVs:

```bash
snakemake --cores all --snakefile Snakefile   # just re-run the same command
```

## Re-run from scratch

```bash
snakemake --cores all --snakefile Snakefile --forceall
```

---

## Outputs

| File | Description |
|------|-------------|
| `results/angles.csv` | All extracted records: `pdb_id`, `left_res`, `central_res`, `right_res`, `angle`, `category` |
| `results/plot.png` | Two-panel figure: KDE (top) + box-plot (bottom) |
| `results/logs/` | Per-rule log files for debugging |

### angles.csv schema

| Column | Type | Description |
|--------|------|-------------|
| `pdb_id` | str | PDB file stem (filename without `.pdb`) |
| `left_res` | str | One-letter code of left residue X |
| `central_res` | str | Always `R` (Arginine) |
| `right_res` | str | One-letter code of right residue X |
| `angle` | float | Signed angle in degrees (−180 to 180) |
| `category` | str | Size category of left residue |

### Category mapping

| Category | Residues |
|----------|----------|
| Tiny | G, A |
| Small | V, P, S, T, C |
| Intermediate | D, L, I, N |
| Large | K, E, M, Q, H |
| Bulky | R, F, Y, W |

---

## Troubleshooting

| Symptom | Fix |
|---------|-----|
| `No .pdb files found` | Check `pdb_dir` in `config.yaml`; paths are resolved relative to the working directory. |
| `DSSP failed` / files skipped | Ensure `mkdssp` is on your `PATH`; set `dssp_executable` to the full binary path. |
| `ImportError: No module named 'Bio'` | Activate the conda env or `pip install biopython`. |
| Pipeline hangs on large batches | Reduce `batch_size` or `n_workers_per_batch` to avoid memory pressure. |
| `MemoryError` in aggregate step | The aggregate step reads all CSVs; ensure ≥ 4 GB RAM for very large datasets. |

---

## Citation / dependencies

- **Biopython**: Cock et al. (2009) *Bioinformatics* 25:1422–1423  
- **DSSP**: Kabsch & Sander (1983) *Biopolymers* 22:2577–2637  
- **Snakemake**: Mölder et al. (2021) *F1000Research* 10:33  
