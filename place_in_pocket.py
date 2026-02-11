#!/usr/bin/env python3
import os
import re
import csv
import glob
import argparse
from typing import List, Tuple, Optional

# ---------- MOL2 helpers ----------


def split_mol2_blocks(text: str):
    blocks = re.split(r"(?=@<TRIPOS>MOLECULE\b)", text)
    # Ignore file preamble/comments before the first MOLECULE marker.
    return [b for b in blocks if b.lstrip().startswith("@<TRIPOS>MOLECULE")]


def mol2_atom_coords(block: str) -> List[Tuple[float, float, float, int]]:
    """
    Return list of (x,y,z, line_index) for atoms in mol2 block.
    We also return the line index so we can rewrite coordinates in-place.
    """
    lines = block.splitlines()
    coords = []
    in_atom = False
    for i, ln in enumerate(lines):
        s = ln.strip()
        if s == "@<TRIPOS>ATOM":
            in_atom = True
            continue
        if s.startswith("@<TRIPOS>") and s != "@<TRIPOS>ATOM":
            in_atom = False
        if in_atom and s and not s.startswith("@<TRIPOS>"):
            parts = ln.split()
            # atom_id atom_name x y z atom_type ...
            if len(parts) >= 6:
                try:
                    x, y, z = float(parts[2]), float(parts[3]), float(parts[4])
                    coords.append((x, y, z, i))
                except ValueError:
                    pass
    if not coords:
        raise ValueError("No ATOM coordinates found in mol2 block.")
    return coords


def centroid_xyz(xyz: List[Tuple[float, float, float]]) -> Tuple[float, float, float]:
    sx = sum(p[0] for p in xyz)
    sy = sum(p[1] for p in xyz)
    sz = sum(p[2] for p in xyz)
    n = len(xyz)
    return (sx / n, sy / n, sz / n)


def translate_mol2_block_to_center(
    block: str, target_center: Tuple[float, float, float]
) -> str:
    lines = block.splitlines()
    coords = mol2_atom_coords(block)
    xyz = [(x, y, z) for (x, y, z, _) in coords]
    cx, cy, cz = centroid_xyz(xyz)
    tx, ty, tz = target_center
    dx, dy, dz = (tx - cx, ty - cy, tz - cz)

    # Rewrite atom lines with shifted coords, preserving most formatting.
    for x, y, z, idx in coords:
        parts = lines[idx].split()
        # parts[2:5] are x y z
        nx, ny, nz = x + dx, y + dy, z + dz
        parts[2] = f"{nx:.4f}"
        parts[3] = f"{ny:.4f}"
        parts[4] = f"{nz:.4f}"
        # Rebuild line with single-space separation (safe mol2)
        lines[idx] = " ".join(parts)

    return "\n".join(lines) + "\n"


def parse_first_mol2_block(path: str) -> str:
    with open(path, "r", encoding="utf-8", errors="replace") as f:
        text = f.read()
    blocks = split_mol2_blocks(text)
    if not blocks:
        raise ValueError(f"No @<TRIPOS>MOLECULE block found in {path}")
    return blocks[0]


# ---------- PDB helpers ----------


def pdb_atom_lines(pdb_path: str) -> List[str]:
    lines = []
    with open(pdb_path, "r", encoding="utf-8", errors="replace") as f:
        for ln in f:
            rec = ln[0:6].strip()
            if rec in ("ATOM", "HETATM"):
                lines.append(ln.rstrip("\n"))
    return lines


def pdb_xyz_from_atom_line(ln: str) -> Optional[Tuple[float, float, float]]:
    # PDB columns: x 31-38, y 39-46, z 47-54 (1-based)
    try:
        x = float(ln[30:38])
        y = float(ln[38:46])
        z = float(ln[46:54])
        return (x, y, z)
    except ValueError:
        return None


def distance2(a, b):
    return (a[0] - b[0]) ** 2 + (a[1] - b[1]) ** 2 + (a[2] - b[2]) ** 2


def pocket_center_from_crystal_mol2(crystal_mol2: str) -> Tuple[float, float, float]:
    blk = parse_first_mol2_block(crystal_mol2)
    coords = mol2_atom_coords(blk)
    xyz = [(x, y, z) for (x, y, z, _) in coords]
    return centroid_xyz(xyz)


def pocket_center_from_pdb_resname(
    pdb_path: str, resname: str
) -> Tuple[float, float, float]:
    atoms = pdb_atom_lines(pdb_path)
    xyz = []
    for ln in atoms:
        # resname is columns 18-20 (0-based 17:20)
        rn = ln[17:20].strip()
        if rn == resname:
            p = pdb_xyz_from_atom_line(ln)
            if p:
                xyz.append(p)
    if not xyz:
        raise ValueError(f"No atoms found for residue name '{resname}' in {pdb_path}")
    return centroid_xyz(xyz)


def write_pocket_pdb(
    protein_pdb: str, out_pdb: str, center: Tuple[float, float, float], radius: float
):
    r2 = radius * radius
    atoms = pdb_atom_lines(protein_pdb)
    kept = []
    for ln in atoms:
        # Keep only protein ATOM records for pocket (typical expectation)
        if ln.startswith("ATOM"):
            p = pdb_xyz_from_atom_line(ln)
            if p and distance2(p, center) <= r2:
                kept.append(ln)
    with open(out_pdb, "w", encoding="utf-8") as f:
        for ln in kept:
            f.write(ln + "\n")
        f.write("END\n")


# ---------- Main ----------


def main():
    ap = argparse.ArgumentParser(
        description="Extract pocket from protein PDB and translate ligands to pocket center (no docking)."
    )
    ap.add_argument("--protein", required=True, help="protein.pdb")
    ap.add_argument(
        "--ligand_dir", required=True, help="Folder of ligand .mol2 files to place"
    )
    ap.add_argument("--outdir", default="placement_out", help="Output folder")
    ap.add_argument(
        "--radius",
        type=float,
        default=8.0,
        help="Pocket radius in Angstroms (for pocket.pdb extraction)",
    )

    # Pocket center definition options:
    g = ap.add_mutually_exclusive_group(required=True)
    g.add_argument(
        "--crystal_mol2", help="Co-crystal ligand mol2 (center = ligand centroid)"
    )
    g.add_argument(
        "--crystal_resname",
        help="Residue name of crystal ligand inside protein.pdb (e.g., LIG)",
    )

    args = ap.parse_args()

    os.makedirs(args.outdir, exist_ok=True)
    pocket_pdb = os.path.join(args.outdir, "pocket.pdb")
    lig_out = os.path.join(args.outdir, "placed_ligands")
    os.makedirs(lig_out, exist_ok=True)

    # 1) Determine pocket center
    if args.crystal_mol2:
        center = pocket_center_from_crystal_mol2(args.crystal_mol2)
        center_source = f"crystal_mol2:{os.path.basename(args.crystal_mol2)}"
    else:
        center = pocket_center_from_pdb_resname(args.protein, args.crystal_resname)
        center_source = f"crystal_resname:{args.crystal_resname}"

    # 2) Write pocket.pdb (protein atoms within radius of center)
    write_pocket_pdb(args.protein, pocket_pdb, center, args.radius)

    # 3) Translate ligands to center
    lig_files = sorted(glob.glob(os.path.join(args.ligand_dir, "*.mol2")))
    if not lig_files:
        raise SystemExit(f"No .mol2 ligands found in {args.ligand_dir}")

    summary_path = os.path.join(args.outdir, "placement_summary.csv")
    rows = []

    for lig_path in lig_files:
        base = os.path.splitext(os.path.basename(lig_path))[0]
        with open(lig_path, "r", encoding="utf-8", errors="replace") as f:
            text = f.read()

        blocks = split_mol2_blocks(text)
        if not blocks:
            print(f"SKIP (no mol2 blocks): {lig_path}")
            continue

        # If a file somehow contains multiple mol2 blocks, place each and write them all back.
        placed_blocks = []
        total_atoms = 0
        for blk in blocks:
            coords = mol2_atom_coords(blk)
            total_atoms += len(coords)
            placed_blocks.append(translate_mol2_block_to_center(blk, center))

        out_mol2 = os.path.join(lig_out, f"{base}_placed.mol2")
        with open(out_mol2, "w", encoding="utf-8") as f:
            f.write("".join(placed_blocks))

        rows.append(
            {
                "ligand": base,
                "num_atoms": total_atoms,
                "out_mol2": os.path.relpath(out_mol2, args.outdir),
            }
        )

    with open(summary_path, "w", newline="", encoding="utf-8") as f:
        w = csv.DictWriter(f, fieldnames=["ligand", "num_atoms", "out_mol2"])
        w.writeheader()
        w.writerows(rows)

    cx, cy, cz = center
    print("DONE")
    print(f"Pocket center ({center_source}): {cx:.3f} {cy:.3f} {cz:.3f}")
    print(f"Pocket PDB written: {pocket_pdb}")
    print(f"Placed ligands written: {lig_out}/")
    print(f"Summary CSV: {summary_path}")


if __name__ == "__main__":
    main()
