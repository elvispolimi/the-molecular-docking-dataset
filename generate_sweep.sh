#!/usr/bin/env bash
set -euo pipefail

# Usage:
#   ./generate_sweep.sh OUT_BASE_DIR LIGAND_DIR SUMMARY_CSV EXT N1 N2 N3 ...
#
# Example:
#   ./generate_sweep.sh out_sweep data/original_ligands data/original_ligands/summary.csv mol2 50000 200000 1000000
#   ./generate_sweep.sh out_sweep data/in_place_ligands/placed_ligands data/in_place_ligands/placement_summary.csv adtmol2 50000 200000 1000000

if [[ $# -lt 5 ]]; then
  echo "Usage: $0 OUT_BASE_DIR LIGAND_DIR SUMMARY_CSV EXT N1 N2 N3 ..."
  exit 1
fi

OUT_BASE_DIR="$1"
LIGAND_DIR="$2"
SUMMARY_CSV="$3"
EXT="$4"
shift 4

# Mild bias toward smaller ligands (tweak if needed)
BIAS="0.8"
ATOM_SCALE="50"
ROTOR_SCALE="8"

mkdir -p "${OUT_BASE_DIR}"

for N in "$@"; do
  OUTDIR="${OUT_BASE_DIR}/out_${N}"
  mkdir -p "${OUTDIR}"
  SUMMARY_FILE="${OUTDIR}/summary.txt"
  REPO_COMMIT="$(git -C . rev-parse HEAD 2>/dev/null || echo unknown)"

  CMD=(python3 generate_the_dataset.py
    --ligand_dir "${LIGAND_DIR}"
    --summary_csv "${SUMMARY_CSV}"
    --ext "${EXT}"
    --n "${N}"
    --bias "${BIAS}"
    --atom_scale "${ATOM_SCALE}"
    --rotor_scale "${ROTOR_SCALE}"
    --mode aggregate
    --no_manifest
    --outdir "${OUTDIR}"
    --out_aggregate dataset.mol2
    --report_buckets
    --bucket_mode LARGE)

  {
    echo "Repo commit: ${REPO_COMMIT}"
    printf "Command:"
    printf " %q" "${CMD[@]}"
    echo
    echo "---"
  } > "${SUMMARY_FILE}"

  "${CMD[@]}" | tee -a "${SUMMARY_FILE}"
done
