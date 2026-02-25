# the_dataset

Utilities for preparing ligand sets and simple pocket placement.

## Requirements

- Python 3.9+
- Optional: RDKit for robust dedup + rotatable bond counts in `get_unique.py`
- Optional: MUDock converter binary for `convert_to_adtmol2.sh`

## Quick Start

- Deduplicate a multi-mol2 file:
  `python get_unique.py input.mol2 -o unique_mol2_out`
- Generate a sampled dataset (manifest only):
  `python generate_the_dataset.py --ligand_dir unique_mol2_out --summary_csv unique_mol2_out/summary.csv --n 1000 --outdir dataset_out`
- Generate multiple dataset sizes in one shot (aggregate only, no per-ligand copies):
  `./generate_sweep.sh out_sweep data/original_ligands data/original_ligands/summary.csv 50000 200000 1000000`
- Filter a summary CSV by atom/rotor ranges:
  `./filter_summary.sh data/original_ligands/summary.csv tmp_32_64.csv --min-atoms 32 --max-atoms 64`
- Place ligands into a pocket centered on a co-crystal ligand:
  `python place_in_pocket.py --protein protein.pdb --ligand_dir unique_mol2_out --crystal_mol2 crystal.mol2 --outdir placement_out`

## End-to-End Flow From `molecules.mol2`

This example starts from one multi-molecule file (`molecules.mol2`) and ends with a final dataset in either `.mol2` or `.adtmol2`.

Inputs used in the example:
- `data/molecules.mol2` (multi-molecule ligands)
- `data/protein.pdb`
- `data/crystal.mol2`

### 1. Split + deduplicate the source `molecules.mol2`

```bash
python get_unique.py data/molecules.mol2 -o work/unique
```

Outputs:
- `work/unique/*.mol2`
- `work/unique/summary.csv`

### 2. Place ligands into the protein pocket

```bash
python place_in_pocket.py \
  --protein data/protein.pdb \
  --ligand_dir work/unique \
  --outdir work/final_mol2 \
  --radius 10 \
  --crystal_mol2 data/crystal.mol2
```

Outputs:
- `work/final_mol2/placed_ligands/*.mol2`
- `work/final_mol2/placement_summary.csv`
- `work/final_mol2/pocket.pdb`

### 3. Build a sampling summary aligned with placed filenames

```bash
awk -F, 'BEGIN{OFS=","} NR==1{print "ligand,num_atoms,rotatable_bonds"; next} NR>1{print $1 "_placed",$2,0}' \
  work/final_mol2/placement_summary.csv > work/final_mol2/placed_summary.csv
```

Output:
- `work/final_mol2/placed_summary.csv`

### 4. Generate the final sampled MOL2 dataset

```bash
python generate_the_dataset.py \
  --ligand_dir work/final_mol2/placed_ligands \
  --summary_csv work/final_mol2/placed_summary.csv \
  --n 1000 \
  --mode copy \
  --outdir work/dataset_sample
```

Final sampled MOL2 dataset:
- `work/dataset_sample/dataset.csv`
- `work/dataset_sample/files/*.mol2`

### 5. Optional: convert final MOL2 dataset to ADTMOL2

```bash
./convert_to_adtmol2.sh work/dataset_sample/files /path/to/mudock/build
```

Final ADTMOL2 dataset:
- `work/dataset_sample/files/adtmol2/*.adtmol2`

### Minimal directory view after full flow

```text
work/
  unique/
    summary.csv
    *.mol2
  final_mol2/
    pocket.pdb
    placement_summary.csv
    placed_summary.csv
    placed_ligands/
      *.mol2
  dataset_sample/
    dataset.csv
    files/
      *.mol2
      adtmol2/
        *.adtmol2
```

## Scripts

- `get_unique.py`: deduplicate ligands in a multi-mol2 file, write per-ligand mol2 files and `summary.csv`.
- `generate_the_dataset.py`: sample ligands with size bias and write `dataset.csv`; can also copy files or aggregate to a single mol2.
  - Manifest only: `--mode manifest` (default)
  - Copy files: `--mode copy`
  - Aggregate: `--mode aggregate --out_aggregate dataset.mol2`
  - Skip manifest: `--no_manifest`
  - Report bucket distribution: `--report_buckets --bucket_mode LARGE`
  - Custom bucket thresholds: `--bucket_thresholds "32,64,96,128,160,192"`
- `generate_sweep.sh`: generate multiple dataset sizes (aggregate only) using `generate_the_dataset.py`.
- `filter_summary.sh`: filter a summary CSV by atom/rotor ranges (for range-specific sampling).
- `place_in_pocket.py`: extract a protein pocket and translate ligands to the pocket center (no docking).
- `convert_to_adtmol2.sh`: convert `.mol2` to `.adtmol2` using MUDock converter.
  - `./convert_to_adtmol2.sh /path/to/ligands /path/to/mudock/build`

## `generate_the_dataset.py` explained

Current scoring logic is intentionally simple: for the time being, sampling weights are computed only from:
- `num_atoms`
- `rotatable_bonds` (rotamers)

The script reads these values from `--summary_csv`, computes a size score, then gives higher probability to smaller ligands.

Size score used internally:
- `size_score = (num_atoms / atom_scale) + (rotatable_bonds / rotor_scale)`
- `weight = 1 / (size_score ^ bias)` (if `bias=0`, sampling is uniform)

Main input parameters:
- `--ligand_dir`: folder containing input ligand files (`.mol2` or `.adtmol2`)
- `--ext`: ligand extension to read (`mol2` or `adtmol2`)
- `--summary_csv`: CSV with ligand identity + `num_atoms` (+ optional `rotatable_bonds`)
- `--n`: number of samples to draw (with replacement)
- `--outdir`: output directory for `dataset.csv` and optional files
- `--mode`: output mode
  - `manifest`: write only `dataset.csv`
  - `copy`: write `dataset.csv` and copy sampled ligand files
  - `aggregate`: write `dataset.csv` and one combined MOL2/ADTMOL2 file
- `--seed`: random seed for reproducibility
- `--bias`: strength of preference for small ligands (`0` disables size bias)
- `--atom_scale`: how strongly atom count contributes to size score
- `--rotor_scale`: how strongly rotamer count contributes to size score
- `--missing_atoms`: fallback value when `num_atoms` is missing in CSV
- `--missing_rotors`: fallback value when `rotatable_bonds` is missing in CSV
- `--allow_unmatched`: include ligands found in `--ligand_dir` but missing from `--summary_csv` (legacy behavior); by default, unmatched ligands are ignored
- `--no_manifest`: do not write `dataset.csv`
- `--report_buckets`: print bucket distribution summary
- `--bucket_mode`: use predefined bucket thresholds (`DEFAULT`, `MEDIUM`, `LARGE`, `EXTREME`)
- `--bucket_thresholds`: override thresholds with a comma-separated list
- `--ext adtmol2`: the script will match CSV rows against multiple keys; for placed ligands it supports `ligand`, `ligand_placed`, and `out_mol2` basenames so `.adtmol2` inputs map correctly

Minimal example:

```bash
python generate_the_dataset.py \
  --ligand_dir work/final_mol2/placed_ligands \
  --summary_csv work/final_mol2/placed_summary.csv \
  --ext mol2 \
  --n 1000 \
  --mode copy \
  --bias 1.5 \
  --atom_scale 30 \
  --rotor_scale 5 \
  --outdir work/dataset_sample
```

Generate multiple sizes into one output folder (aggregate only, no manifest):

```bash
./generate_sweep.sh out_sweep work/final_mol2/placed_ligands work/final_mol2/placement_summary.csv 50000 200000 1000000
```

## Structure

- `data/`: input data (if used)
- `get_unique.py`, `generate_the_dataset.py`, `place_in_pocket.py`: core scripts
- `convert_to_adtmol2.sh`: mol2 to adtmol2 conversion helper
