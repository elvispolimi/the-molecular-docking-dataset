#!/usr/bin/env bash
set -euo pipefail

# Usage:
#   ./generate_sweep.sh OUT_BASE_DIR N1 N2 N3 ...
#
# Example:
#   ./generate_sweep.sh out_sweep 50000 200000 1000000

if [[ $# -lt 2 ]]; then
  echo "Usage: $0 OUT_BASE_DIR N1 N2 N3 ..."
  exit 1
fi

OUT_BASE_DIR="$1"
shift

LIGAND_DIR="data/original_ligands"
SUMMARY_CSV="data/original_ligands/summary.csv"

# Mild bias toward smaller ligands (tweak if needed)
BIAS="0.8"
ATOM_SCALE="50"
ROTOR_SCALE="8"

mkdir -p "${OUT_BASE_DIR}"

for N in "$@"; do
  OUTDIR="${OUT_BASE_DIR}/out_${N}"
  mkdir -p "${OUTDIR}"
  python3 generate_the_dataset.py \
    --ligand_dir "${LIGAND_DIR}" \
    --summary_csv "${SUMMARY_CSV}" \
    --n "${N}" \
    --bias "${BIAS}" \
    --atom_scale "${ATOM_SCALE}" \
    --rotor_scale "${ROTOR_SCALE}" \
    --mode aggregate \
    --outdir "${OUTDIR}" \
    --out_aggregate dataset.mol2 \
    --report_buckets \
    --bucket_mode LARGE
done
