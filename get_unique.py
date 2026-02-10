#!/usr/bin/env python3
"""
get_unique.py

Usage:
  python get_unique.py input.mol2 -o unique_mol2_out

Outputs:
  unique_mol2_out/
    ligand_001.mol2
    ligand_002.mol2
    ...
    summary.csv

Notes:
- Uniqueness is detected primarily by RDKit canonical SMILES (if RDKit installed).
- If RDKit is not installed, it falls back to a strict topology-ish hash built from
  elements + bond connectivity/order (works for exact duplicates, but less robust than RDKit).
- Rotamers = number of rotatable bonds (RDKit). If no RDKit, rotamers reported as NA.
"""

import os
import re
import csv
import argparse
import hashlib
from collections import defaultdict


def split_mol2_blocks(text: str):
    # Split into blocks starting at @<TRIPOS>MOLECULE
    # Keep the marker by splitting with a lookahead.
    blocks = re.split(r"(?=@<TRIPOS>MOLECULE\b)", text)
    return [b for b in blocks if b.strip()]


def parse_counts_from_mol2_block(block: str):
    """
    Returns:
      mol_name (str), natoms (int|None), nbonds (int|None)
    Mol2 format:
      @<TRIPOS>MOLECULE
      NAME
      <num_atoms> <num_bonds> ...
    """
    lines = block.splitlines()
    mol_name = None
    natoms = None
    nbonds = None

    # Find MOLECULE header line index
    try:
        i = next(
            idx for idx, ln in enumerate(lines) if ln.strip() == "@<TRIPOS>MOLECULE"
        )
    except StopIteration:
        return mol_name, natoms, nbonds

    if i + 1 < len(lines):
        mol_name = lines[i + 1].strip() or None

    # counts line typically at i+2
    if i + 2 < len(lines):
        parts = lines[i + 2].split()
        if len(parts) >= 2:
            try:
                natoms = int(parts[0])
                nbonds = int(parts[1])
            except ValueError:
                pass

    return mol_name, natoms, nbonds


def parse_atoms_and_bonds_for_fallback(block: str):
    """
    Fallback parser: extract element-ish tokens from ATOM section,
    and bond connectivity/order from BOND section.
    This is used to build a stable hash when RDKit isn't available.

    Returns:
      elements: list[str]
      bonds: list[tuple[int,int,str]]  (a1,a2,order_token)
    """
    lines = block.splitlines()
    elements = []
    bonds = []

    def find_section(start_tag):
        try:
            return next(i for i, ln in enumerate(lines) if ln.strip() == start_tag)
        except StopIteration:
            return None

    atom_i = find_section("@<TRIPOS>ATOM")
    bond_i = find_section("@<TRIPOS>BOND")

    if atom_i is not None:
        end_i = bond_i if bond_i is not None else len(lines)
        for ln in lines[atom_i + 1 : end_i]:
            if not ln.strip() or ln.startswith("@<TRIPOS>"):
                break
            # mol2 atom line format:
            # atom_id atom_name x y z atom_type [subst_id subst_name charge]
            parts = ln.split()
            if len(parts) >= 6:
                atom_type = parts[5]  # e.g. C.3, N.am, O.co2
                # element is before '.', or the first 1-2 letters
                elem = atom_type.split(".")[0]
                elements.append(elem)

    if bond_i is not None:
        end_i = len(lines)
        for ln in lines[bond_i + 1 : end_i]:
            if not ln.strip() or ln.startswith("@<TRIPOS>"):
                break
            # bond line: bond_id origin_atom_id target_atom_id bond_type
            parts = ln.split()
            if len(parts) >= 4:
                try:
                    a1 = int(parts[1])
                    a2 = int(parts[2])
                except ValueError:
                    continue
                btype = parts[3]  # "1", "2", "ar", "am", etc.
                if a1 > a2:
                    a1, a2 = a2, a1
                bonds.append((a1, a2, btype))

    return elements, bonds


def fallback_key(block: str):
    elements, bonds = parse_atoms_and_bonds_for_fallback(block)
    # Normalize:
    # - elements in original order
    # - bonds sorted
    payload = {
        "elements": elements,
        "bonds": sorted(bonds),
    }
    s = repr(payload).encode("utf-8")
    return hashlib.sha256(s).hexdigest()


def try_rdkit_key_and_rotors(block: str):
    """
    If RDKit is installed and can parse mol2, return:
      (key, rotatable_bonds_count, rdkit_num_atoms)
    key is canonical SMILES (isomeric), which is generally robust for duplicates.
    """
    try:
        from rdkit import Chem
        from rdkit import RDLogger
        from rdkit.Chem import rdMolDescriptors
    except Exception:
        return None, None, None

    # Silence RDKit warnings/errors (e.g., kekulize failures) for noisy inputs.
    RDLogger.DisableLog("rdApp.warning")
    RDLogger.DisableLog("rdApp.error")

    mol = Chem.MolFromMol2Block(block, sanitize=True, removeHs=False)
    if mol is None:
        # Try a non-sanitized parse, then sanitize without kekulization.
        mol = Chem.MolFromMol2Block(block, sanitize=False, removeHs=False)
        if mol is None:
            return None, None, None
        try:
            sanitize_ops = Chem.SanitizeFlags.SANITIZE_ALL ^ Chem.SanitizeFlags.SANITIZE_KEKULIZE
            Chem.SanitizeMol(mol, sanitizeOps=sanitize_ops)
        except Exception:
            # If even partial sanitization fails, fall back.
            return None, None, None

    # Canonical isomeric SMILES as identity key
    try:
        key = Chem.MolToSmiles(mol, isomericSmiles=True, canonical=True)
    except Exception:
        return None, None, None

    # Rotamers = rotatable bonds
    try:
        rotors = rdMolDescriptors.CalcNumRotatableBonds(mol)
    except Exception:
        rotors = None

    return key, int(rotors), int(mol.GetNumAtoms())


def safe_filename(s: str, default="ligand"):
    s = s or default
    s = re.sub(r"[^A-Za-z0-9_.-]+", "_", s).strip("_")
    return s[:80] if s else default


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("input_mol2", help="Multi-molecule MOL2 file")
    ap.add_argument(
        "-o", "--outdir", default="unique_mol2_out", help="Output directory"
    )
    args = ap.parse_args()

    os.makedirs(args.outdir, exist_ok=True)

    with open(args.input_mol2, "r", encoding="utf-8", errors="replace") as f:
        text = f.read()

    blocks = split_mol2_blocks(text)
    if not blocks:
        raise SystemExit("No @<TRIPOS>MOLECULE blocks found.")

    # Map key -> data
    rep_block = {}  # key -> representative block text (first occurrence)
    rep_name = {}  # key -> mol name
    rep_atoms = {}  # key -> atom count
    rep_rotors = {}  # key -> rotatable bonds count (or None)
    counts = defaultdict(int)  # key -> repetition count

    used_rdkit = False

    for blk in blocks:
        mol_name, natoms_header, _ = parse_counts_from_mol2_block(blk)

        key, rotors, rdkit_atoms = try_rdkit_key_and_rotors(blk)
        if key is not None:
            used_rdkit = True
            natoms = rdkit_atoms if rdkit_atoms is not None else natoms_header
        else:
            key = fallback_key(blk)
            rotors = None
            natoms = natoms_header

        counts[key] += 1

        if key not in rep_block:
            rep_block[key] = blk
            rep_name[key] = mol_name
            rep_atoms[key] = natoms
            rep_rotors[key] = rotors

    # Write unique mol2 files
    summary_rows = []
    keys_sorted = sorted(counts.keys(), key=lambda k: (-counts[k], k))

    for idx, key in enumerate(keys_sorted, start=1):
        name = rep_name.get(key) or f"ligand_{idx:03d}"
        fname_base = safe_filename(name, default=f"ligand_{idx:03d}")
        out_mol2 = os.path.join(args.outdir, f"{idx:03d}_{fname_base}.mol2")

        with open(out_mol2, "w", encoding="utf-8") as f:
            # keep exact original mol2 formatting
            f.write(rep_block[key].rstrip() + "\n")

        summary_rows.append(
            {
                "unique_index": idx,
                "name": name,
                "identity_key": key if used_rdkit else f"fallback_sha256:{key}",
                "repetitions": counts[key],
                "num_atoms": rep_atoms.get(key),
                "rotatable_bonds": (
                    rep_rotors.get(key) if rep_rotors.get(key) is not None else "NA"
                ),
                "mol2_file": os.path.basename(out_mol2),
            }
        )

    # Write summary.csv
    summary_csv = os.path.join(args.outdir, "summary.csv")
    with open(summary_csv, "w", newline="", encoding="utf-8") as f:
        w = csv.DictWriter(
            f,
            fieldnames=[
                "unique_index",
                "name",
                "identity_key",
                "repetitions",
                "num_atoms",
                "rotatable_bonds",
                "mol2_file",
            ],
        )
        w.writeheader()
        w.writerows(summary_rows)

    print(f"Read molecules: {len(blocks)}")
    print(f"Unique ligands: {len(keys_sorted)}")
    print(f"Wrote: {args.outdir}/ (mol2 files + summary.csv)")
    print(
        f"RDKit used: {'YES' if used_rdkit else 'NO (fallback hashing; rotatable bonds = NA)'}"
    )


if __name__ == "__main__":
    main()
