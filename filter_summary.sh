#!/usr/bin/env bash
set -euo pipefail

# Filter a summary CSV by atom/rotor ranges.
#
# Usage:
#   ./filter_summary.sh INPUT_CSV OUTPUT_CSV [--min-atoms N] [--max-atoms N] [--min-rotors N] [--max-rotors N]
#
# Notes:
# - Column names are detected from the header. Expected: num_atoms, rotatable_bonds (or rotamers/rotors).
# - If a bound is omitted, it is ignored.

if [[ $# -lt 2 ]]; then
  echo "Usage: $0 INPUT_CSV OUTPUT_CSV [--min-atoms N] [--max-atoms N] [--min-rotors N] [--max-rotors N]"
  exit 1
fi

INPUT_CSV="$1"
OUTPUT_CSV="$2"
shift 2

MIN_ATOMS=""
MAX_ATOMS=""
MIN_ROTORS=""
MAX_ROTORS=""

while [[ $# -gt 0 ]]; do
  case "$1" in
    --min-atoms) MIN_ATOMS="$2"; shift 2 ;;
    --max-atoms) MAX_ATOMS="$2"; shift 2 ;;
    --min-rotors) MIN_ROTORS="$2"; shift 2 ;;
    --max-rotors) MAX_ROTORS="$2"; shift 2 ;;
    *) echo "Unknown option: $1"; exit 1 ;;
  esac
done

awk -F, -v min_a="${MIN_ATOMS}" -v max_a="${MAX_ATOMS}" -v min_r="${MIN_ROTORS}" -v max_r="${MAX_ROTORS}" '
  BEGIN {
    num_atoms_idx = 0;
    rot_idx = 0;
  }
  NR==1 {
    for (i=1; i<=NF; i++) {
      h = $i;
      gsub(/^[ \t]+|[ \t]+$/, "", h);
      if (h == "num_atoms" || h == "atoms" || h == "natoms") num_atoms_idx = i;
      if (h == "rotatable_bonds" || h == "rotamers" || h == "rotors" || h == "n_rotatable_bonds") rot_idx = i;
    }
    if (num_atoms_idx == 0) {
      print "ERROR: num_atoms column not found in header" > "/dev/stderr";
      exit 1;
    }
    print $0;
    next;
  }
  {
    atoms = $num_atoms_idx + 0;
    rot = (rot_idx ? $rot_idx + 0 : 0);

    if (min_a != "" && atoms < min_a) next;
    if (max_a != "" && atoms > max_a) next;
    if (min_r != "" && rot < min_r) next;
    if (max_r != "" && rot > max_r) next;

    print $0;
  }
' "${INPUT_CSV}" > "${OUTPUT_CSV}"
