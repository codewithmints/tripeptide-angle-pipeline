# Protein Side Chain Analysis Pipeline

This project analyzes **X–Arg(helix)–X tripeptides** from PDB structures and computes the signed angle between Cα → centroid vectors. The results are grouped by left-residue size category and visualized using a KDE plot.

---

## 📁 Folder Structure


BET-104_tripeptide/
├── Snakefile
├── config.yaml
├── README.md
├── scripts/
│ ├── parse_pdb.py
│ ├── aggregate.py
│ └── plot.py
├── results/
│ ├── angles.csv
│ └── plot.png


---

## ▶ Run Command

snakemake --cores all --configfile config.yaml

This single command runs the full pipeline:

Parses all PDB files
Computes angles
Aggregates results
Generates final plot

## Output Files
results/angles.csv → aggregated angle data
results/plot.png → final KDE plot

## Notes
Input PDB directory is set in config.yaml
DSSP is required for secondary structure assignment
Pipeline is parallelized using Snakemake
