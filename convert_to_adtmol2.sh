#!/usr/bin/env bash
set -euo pipefail

# Usage:
#   ./convert_to_adtmol2.sh /path/to/ligands /path/to/mudock/build
#
# Outputs go to:
#   /path/to/ligands/adtmol2/<name>.adtmol2

DIR="${1:-.}"
MUDOCK_BUILD="${2:-/opt/mudock/build}"

# Output subfolder name (change if you want)
OUT_SUBDIR="adtmol2"
OUT_DIR="$DIR/$OUT_SUBDIR"
mkdir -p "$OUT_DIR"

shopt -s nullglob
for in_file in "$DIR"/*.mol2; do
  base="$(basename "$in_file" .mol2)"
  out_file="$OUT_DIR/${base}.adtmol2"
  echo "Converting: $in_file -> $out_file"
  "$MUDOCK_BUILD/converter" --input "$in_file" --output "$out_file"
done
