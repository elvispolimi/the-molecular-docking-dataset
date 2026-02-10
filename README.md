# the_dataset

Utilities for preparing ligand sets and simple pocket placement.

Requirements

- Python 3.9+
- Optional: RDKit for robust dedup + rotatable bond counts in `get_unique.py`
- Optional: MUDock converter binary for `convert_admol2.sh`

Quick start

- Deduplicate a multi-mol2 file:
  `python get_unique.py input.mol2 -o unique_mol2_out`
- Generate a sampled dataset (manifest only):
  `python generate_the_dataset.py --ligand_dir unique_mol2_out --summary_csv unique_mol2_out/summary.csv --n 1000 --outdir dataset_out`
- Place ligands into a pocket centered on a co-crystal ligand:
  `python place_in_pocket.py --protein protein.pdb --ligand_dir unique_mol2_out --crystal_mol2 crystal.mol2 --outdir placement_out`

Scripts

- `get_unique.py`: deduplicate ligands in a multi-mol2 file, write per-ligand mol2 files and `summary.csv`.
- `generate_the_dataset.py`: sample ligands with size bias and write `dataset.csv`; can also copy files or aggregate to a single mol2.
  - Manifest only: `--mode manifest` (default)
  - Copy files: `--mode copy`
  - Aggregate: `--mode aggregate --out_aggregate dataset.mol2`
- `place_in_pocket.py`: extract a protein pocket and translate ligands to the pocket center (no docking).
- `convert_admol2.sh`: convert `.mol2` to `.adtmol2` using MUDock converter.
  - `./convert_admol2.sh /path/to/ligands /path/to/mudock/build`

Structure

- `data/`: input data (if used)
- `get_unique.py`, `generate_the_dataset.py`, `place_in_pocket.py`: core scripts
- `convert_admol2.sh`: mol2 to adtmol2 conversion helper
