#!/usr/bin/env python3
import os
import csv
import glob
import math
import argparse
import random
import shutil
import re
from typing import Dict, Tuple, Optional, List


def safe_float(x: str) -> Optional[float]:
    if x is None:
        return None
    s = str(x).strip()
    if s == "" or s.upper() == "NA" or s.lower() == "nan":
        return None
    try:
        return float(s)
    except ValueError:
        return None


def load_summary_csv(
    summary_csv: str,
) -> Dict[str, Tuple[Optional[float], Optional[float]]]:
    with open(summary_csv, "r", encoding="utf-8", errors="replace") as f:
        reader = csv.DictReader(f)
        headers = [h.strip() for h in (reader.fieldnames or [])]

        def pick_col(candidates: List[str]) -> Optional[str]:
            for c in candidates:
                if c in headers:
                    return c
            lower = {h.lower(): h for h in headers}
            for c in candidates:
                if c.lower() in lower:
                    return lower[c.lower()]
            return None

        ligand_col = pick_col(["ligand", "name", "mol2_file", "out_mol2", "mol2"])
        atoms_col = pick_col(["num_atoms", "atoms", "natoms"])
        rot_col = pick_col(
            ["rotatable_bonds", "rotamers", "rotors", "n_rotatable_bonds"]
        )

        if ligand_col is None or atoms_col is None:
            raise SystemExit(
                f"CSV must contain ligand/name and num_atoms columns. Found headers: {headers}"
            )

        out: Dict[str, Tuple[Optional[float], Optional[float]]] = {}
        for row in reader:
            lig = (row.get(ligand_col) or "").strip()
            if not lig:
                continue
            lig_base = os.path.splitext(os.path.basename(lig))[0]
            atoms = safe_float(row.get(atoms_col))
            rot = safe_float(row.get(rot_col)) if rot_col else None
            out[lig_base] = (atoms, rot)
        return out


def list_input_files(ligand_dir: str, ext: str) -> List[str]:
    ext = ext.lower().lstrip(".")
    return sorted(glob.glob(os.path.join(ligand_dir, f"*.{ext}")))


def build_weights(
    files, stats, bias, atom_scale, rotor_scale, missing_atoms, missing_rotors, eps=1e-6
):
    weights = []
    meta = []
    for fp in files:
        base = os.path.splitext(os.path.basename(fp))[0]
        atoms, rot = stats.get(base, (None, None))
        if atoms is None:
            atoms = missing_atoms
        if rot is None:
            rot = missing_rotors

        size_score = (atoms / atom_scale) + (rot / rotor_scale)
        size_score = max(size_score, eps)
        w = 1.0 / (size_score**bias) if bias != 0 else 1.0

        weights.append(w)
        meta.append(
            {
                "base": base,
                "path": fp,
                "num_atoms": atoms,
                "rotatable_bonds": rot,
                "size_score": size_score,
                "weight": w,
            }
        )
    return weights, meta


def weighted_choice_index(rng: random.Random, weights: List[float]) -> int:
    total = sum(weights)
    if total <= 0:
        return rng.randrange(len(weights))
    r = rng.random() * total
    acc = 0.0
    for i, w in enumerate(weights):
        acc += w
        if r <= acc:
            return i
    return len(weights) - 1


def split_mol2_blocks(text: str) -> List[str]:
    blocks = re.split(r"(?=@<TRIPOS>MOLECULE\b)", text)
    return [b for b in blocks if b.strip()]


def read_all_blocks(path: str) -> List[str]:
    with open(path, "r", encoding="utf-8", errors="replace") as f:
        txt = f.read()
    blocks = split_mol2_blocks(txt)
    if not blocks:
        # If the file doesn't contain the marker, still keep as a single block
        return [txt.rstrip() + "\n"]
    return [b.rstrip() + "\n" for b in blocks]


def main():
    ap = argparse.ArgumentParser(
        description="Create a random ligand dataset biased toward smaller ligands."
    )
    ap.add_argument("--ligand_dir", required=True, help="Folder containing ligands")
    ap.add_argument(
        "--ext",
        default="mol2",
        choices=["mol2", "adtmol2"],
        help="Input extension to use",
    )
    ap.add_argument(
        "--summary_csv",
        required=True,
        help="CSV with num_atoms and rotatable_bonds/rotamers info",
    )
    ap.add_argument("--outdir", default="dataset_out", help="Output folder")
    ap.add_argument(
        "--n", type=int, required=True, help="Number of samples (with replacement)"
    )
    ap.add_argument("--seed", type=int, default=123, help="Random seed")

    # sampling knobs
    ap.add_argument(
        "--bias",
        type=float,
        default=1.5,
        help="0=uniform; higher => small ligands more likely",
    )
    ap.add_argument(
        "--atom_scale",
        type=float,
        default=30.0,
        help="Atoms scaling (bigger => atoms matter less)",
    )
    ap.add_argument(
        "--rotor_scale",
        type=float,
        default=5.0,
        help="Rotors scaling (bigger => rotors matter less)",
    )
    ap.add_argument(
        "--missing_atoms", type=float, default=50.0, help="Fallback atoms if missing"
    )
    ap.add_argument(
        "--missing_rotors", type=float, default=5.0, help="Fallback rotors if missing"
    )

    # bucket reporting
    ap.add_argument(
        "--report_buckets",
        action="store_true",
        help="Report bucket distribution for sampled ligands",
    )
    ap.add_argument(
        "--bucket_mode",
        choices=["DEFAULT", "MEDIUM", "LARGE", "EXTREME"],
        default="LARGE",
        help="Bucket thresholds for reporting (matches muDock reorder_buffer.hpp)",
    )

    # output mode
    ap.add_argument(
        "--mode",
        choices=["manifest", "copy", "aggregate"],
        default="manifest",
        help="manifest: only dataset.csv; copy: copy sampled files; aggregate: write one combined file",
    )
    ap.add_argument(
        "--out_aggregate",
        default="",
        help="(aggregate mode) output filename, e.g. dataset.mol2",
    )
    args = ap.parse_args()

    os.makedirs(args.outdir, exist_ok=True)

    stats = load_summary_csv(args.summary_csv)
    files = list_input_files(args.ligand_dir, args.ext)
    if not files:
        raise SystemExit(f"No *.{args.ext} files found in {args.ligand_dir}")

    weights, meta = build_weights(
        files,
        stats,
        args.bias,
        args.atom_scale,
        args.rotor_scale,
        args.missing_atoms,
        args.missing_rotors,
    )

    rng = random.Random(args.seed)

    # bucket thresholds for reporting
    bucket_thresholds = {
        "DEFAULT": [256],
        "MEDIUM": [32, 64, 128, 256],
        "LARGE": [32, 64, 128, 160, 192, 256],
        "EXTREME": [16, 32, 58, 64, 96, 128, 160, 192, 256],
    }[args.bucket_mode]
    bucket_counts = [0 for _ in bucket_thresholds]
    bucket_atom_sums = [0 for _ in bucket_thresholds]
    bucket_overflow = 0
    bucket_overflow_atoms = 0

    # Prepare outputs
    manifest_path = os.path.join(args.outdir, "dataset.csv")
    copied_dir = os.path.join(args.outdir, "files")

    agg_path = ""
    agg_handle = None
    if args.mode == "aggregate":
        if not args.out_aggregate.strip():
            raise SystemExit("--out_aggregate is required when --mode aggregate")
        agg_path = os.path.join(args.outdir, args.out_aggregate)
        agg_handle = open(agg_path, "w", encoding="utf-8")

    if args.mode == "copy":
        os.makedirs(copied_dir, exist_ok=True)

    with open(manifest_path, "w", newline="", encoding="utf-8") as f:
        fieldnames = [
            "sample_id",
            "selected_base",
            "source_path",
            "num_atoms",
            "rotatable_bonds",
            "size_score",
            "weight",
            "copied_path",
            "aggregate_file",
        ]
        w = csv.DictWriter(f, fieldnames=fieldnames)
        w.writeheader()

        for k in range(args.n):
            idx = weighted_choice_index(rng, weights)
            m = meta[idx]
            src = m["path"]

            if args.report_buckets:
                atoms = int(m["num_atoms"])
                # count thresholds <= atoms (same logic as reorder_buffer)
                bucket_index = 0
                for t in bucket_thresholds:
                    if t <= atoms:
                        bucket_index += 1
                if bucket_index >= len(bucket_thresholds):
                    bucket_overflow += 1
                    bucket_overflow_atoms += atoms
                else:
                    bucket_counts[bucket_index] += 1
                    bucket_atom_sums[bucket_index] += atoms

            copied = ""
            if args.mode == "copy":
                ext = os.path.splitext(src)[1]
                dst = os.path.join(copied_dir, f"{k:06d}_{m['base']}{ext}")
                shutil.copy2(src, dst)
                copied = os.path.relpath(dst, args.outdir)

            if args.mode == "aggregate":
                # Append all MOL2 blocks from the chosen file
                blocks = read_all_blocks(src)
                # Optional: make molecule names unique by prefixing sample_id
                # (simple text tweak on 2nd line after @<TRIPOS>MOLECULE if present)
                for b in blocks:
                    lines = b.splitlines(True)
                    if len(lines) >= 2 and lines[0].startswith("@<TRIPOS>MOLECULE"):
                        # rename molecule line
                        lines[1] = f"{k:06d}_{m['base']}\n"
                        agg_handle.write("".join(lines))
                    else:
                        agg_handle.write(b.rstrip() + "\n")

            w.writerow(
                {
                    "sample_id": k,
                    "selected_base": m["base"],
                    "source_path": os.path.abspath(src),
                    "num_atoms": m["num_atoms"],
                    "rotatable_bonds": m["rotatable_bonds"],
                    "size_score": f"{m['size_score']:.6f}",
                    "weight": f"{m['weight']:.6e}",
                    "copied_path": copied,
                    "aggregate_file": (
                        args.out_aggregate if args.mode == "aggregate" else ""
                    ),
                }
            )

    if agg_handle:
        agg_handle.close()

    print("DONE")
    print(f"Input ligands: {len(files)} (*.{args.ext})")
    print(f"Samples generated: {args.n} (with replacement)")
    print(f"Manifest: {manifest_path}")
    if args.mode == "copy":
        print(f"Copied files into: {copied_dir}/")
    if args.mode == "aggregate":
        print(f"Aggregate file: {agg_path}")
    if args.report_buckets:
        print(f"Bucket mode: {args.bucket_mode}")
        for i, t in enumerate(bucket_thresholds):
            print(
                f"  bucket <= {t:>3}: ligands={bucket_counts[i]} atoms={bucket_atom_sums[i]}"
            )
        if bucket_overflow:
            print(
                f"  overflow (> {bucket_thresholds[-1]}): ligands={bucket_overflow} atoms={bucket_overflow_atoms}"
            )


if __name__ == "__main__":
    main()
