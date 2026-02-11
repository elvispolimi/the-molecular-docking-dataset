# the_dataset

Utilities for preparing ligand sets and simple pocket placement.

Requirements

- Python 3.9+
- Optional: RDKit for robust dedup + rotatable bond counts in `get_unique.py`
- Optional: MUDock converter binary for `convert_to_adtmol2.sh`

Quick start

- Deduplicate a multi-mol2 file:
  `python get_unique.py input.mol2 -o unique_mol2_out`
- Generate a sampled dataset (manifest only):
  `python generate_the_dataset.py --ligand_dir unique_mol2_out --summary_csv unique_mol2_out/summary.csv --n 1000 --outdir dataset_out`
- Place ligands into a pocket centered on a co-crystal ligand:
  `python place_in_pocket.py --protein protein.pdb --ligand_dir unique_mol2_out --crystal_mol2 crystal.mol2 --outdir placement_out`

End-to-end flow from `molecules.mol2`

This example starts from one multi-molecule file (`molecules.mol2`) and ends with a final dataset in either `.mol2` or `.adtmol2`.

Inputs used in the example:
- `data/molecules.mol2` (multi-molecule ligands)
- `data/protein.pdb`
- `data/crystal.mol2`

1. Split + deduplicate the source `molecules.mol2`:

```bash
python get_unique.py data/molecules.mol2 -o work/unique
```

Outputs:
- `work/unique/*.mol2`
- `work/unique/summary.csv`

2. Sample a dataset (copy mode keeps sampled files physically separated):

```bash
python generate_the_dataset.py \
  --ligand_dir work/unique \
  --summary_csv work/unique/summary.csv \
  --n 1000 \
  --mode copy \
  --outdir work/dataset_sample
```

Outputs:
- `work/dataset_sample/dataset.csv`
- `work/dataset_sample/files/*.mol2`

3. Place sampled ligands into the protein pocket:

```bash
python place_in_pocket.py \
  --protein data/protein.pdb \
  --ligand_dir work/dataset_sample/files \
  --outdir work/final_mol2 \
  --radius 10 \
  --crystal_mol2 data/crystal.mol2
```

Final MOL2 dataset:
- `work/final_mol2/placed_ligands/*.mol2`
- `work/final_mol2/placement_summary.csv`
- `work/final_mol2/pocket.pdb`

4. Optional: convert final MOL2 dataset to ADTMOL2:

```bash
./convert_to_adtmol2.sh work/final_mol2/placed_ligands /path/to/mudock/build
```

Final ADTMOL2 dataset:
- `work/final_mol2/placed_ligands/adtmol2/*.adtmol2`

Minimal directory view after full flow:

```text
work/
  unique/
    summary.csv
    *.mol2
  dataset_sample/
    dataset.csv
    files/
      *.mol2
  final_mol2/
    pocket.pdb
    placement_summary.csv
    placed_ligands/
      *.mol2
      adtmol2/
        *.adtmol2
```

Scripts

- `get_unique.py`: deduplicate ligands in a multi-mol2 file, write per-ligand mol2 files and `summary.csv`.
- `generate_the_dataset.py`: sample ligands with size bias and write `dataset.csv`; can also copy files or aggregate to a single mol2.
  - Manifest only: `--mode manifest` (default)
  - Copy files: `--mode copy`
  - Aggregate: `--mode aggregate --out_aggregate dataset.mol2`
- `place_in_pocket.py`: extract a protein pocket and translate ligands to the pocket center (no docking).
- `convert_to_adtmol2.sh`: convert `.mol2` to `.adtmol2` using MUDock converter.
  - `./convert_to_adtmol2.sh /path/to/ligands /path/to/mudock/build`

Structure

- `data/`: input data (if used)
- `get_unique.py`, `generate_the_dataset.py`, `place_in_pocket.py`: core scripts
- `convert_to_adtmol2.sh`: mol2 to adtmol2 conversion helper
