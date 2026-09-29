#!/usr/bin/env python3
"""
Generate an AlphaFold 3-compatible custom ligand CCD CIF from a SMILES string.

Outputs:
    <COMP_ID>.cif
    <COMP_ID>.png

Example:
    python smiles_to_ccd.py \
        --smiles "CC(=O)Oc1ccccc1C(=O)O" \
        --comp-id ASP \
        --name "Aspirin"
"""

from __future__ import annotations

import argparse
import sys
from pathlib import Path

from rdkit import Chem
from rdkit.Chem import AllChem, Draw
from rdkit.Chem.rdMolDescriptors import CalcExactMolWt, CalcMolFormula


def cif_quote(value: str) -> str:
    """Return a simple quoted CIF string."""
    return "'" + str(value).replace("'", "''") + "'"


def build_3d_molecule(smiles: str, random_seed: int = 0xF00D) -> Chem.Mol:
    """Parse SMILES, add hydrogens, generate a 3D conformer, and optimize it."""
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        raise ValueError(f"RDKit could not parse SMILES: {smiles}")

    mol = Chem.AddHs(mol)

    params = AllChem.ETKDGv3()
    params.randomSeed = random_seed

    status = AllChem.EmbedMolecule(mol, params)
    if status != 0:
        raise RuntimeError("RDKit failed to generate a 3D conformer.")

    # Prefer MMFF when parameters are available; otherwise fall back to UFF.
    if AllChem.MMFFHasAllMoleculeParams(mol):
        AllChem.MMFFOptimizeMolecule(mol)
    elif AllChem.UFFHasAllMoleculeParams(mol):
        AllChem.UFFOptimizeMolecule(mol)

    return mol


def smiles_to_ccd_cif(
    smiles: str,
    comp_id: str,
    comp_name: str = "?",
    formula: str | None = None,
    weight: float | None = None,
    random_seed: int = 0xF00D,
) -> tuple[str, Chem.Mol, dict]:
    """
    Convert a SMILES string to a minimal CCD-style CIF block.

    Returns
    -------
    cif_text : str
        Generated CIF content.
    mol : rdkit.Chem.Mol
        Hydrogen-containing molecule with a 3D conformer.
    metadata : dict
        Formula, exact molecular weight, SMILES, and InChI.
    """
    comp_id = comp_id.strip().upper()
    if not comp_id:
        raise ValueError("comp_id must not be empty.")

    mol = build_3d_molecule(smiles, random_seed=random_seed)

    canonical_smiles = Chem.MolToSmiles(Chem.RemoveHs(mol), canonical=True)
    inchi = Chem.MolToInchi(Chem.RemoveHs(mol))

    if formula is None:
        formula = CalcMolFormula(mol)
    if weight is None:
        weight = CalcExactMolWt(mol)

    # Generate atom names C1, C2, O1, H1, ...
    element_counts: dict[str, int] = {}
    atom_records = []
    conformer = mol.GetConformer()

    for i, atom in enumerate(mol.GetAtoms()):
        symbol = atom.GetSymbol()
        element_counts[symbol] = element_counts.get(symbol, 0) + 1
        atom_name = f"{symbol}{element_counts[symbol]}"
        pos = conformer.GetAtomPosition(i)

        atom_records.append(
            (
                comp_id,
                atom_name,
                symbol,
                atom.GetFormalCharge(),
                "N",
                pos.x,
                pos.y,
                pos.z,
            )
        )

    # Extract bond information.
    bond_records = []
    bond_order_map = {
        Chem.rdchem.BondType.SINGLE: "SING",
        Chem.rdchem.BondType.DOUBLE: "DOUB",
        Chem.rdchem.BondType.TRIPLE: "TRIP",
        Chem.rdchem.BondType.AROMATIC: "AROM",
    }

    for bond in mol.GetBonds():
        atom_1 = atom_records[bond.GetBeginAtomIdx()][1]
        atom_2 = atom_records[bond.GetEndAtomIdx()][1]
        order = bond_order_map.get(bond.GetBondType(), "SING")
        aromatic = "Y" if bond.GetIsAromatic() else "N"
        bond_records.append((comp_id, atom_1, atom_2, order, aromatic))

    lines = [
        f"data_{comp_id}",
        "#",
        f"_chem_comp.id {comp_id}",
        f"_chem_comp.name {cif_quote(comp_name)}",
        "_chem_comp.type non-polymer",
        f"_chem_comp.formula {cif_quote(formula)}",
        "_chem_comp.mon_nstd_parent_comp_id ?",
        "_chem_comp.pdbx_synonyms ?",
        f"_chem_comp.formula_weight {weight:.3f}",
        "#",
        "loop_",
        "_chem_comp_atom.comp_id",
        "_chem_comp_atom.atom_id",
        "_chem_comp_atom.type_symbol",
        "_chem_comp_atom.charge",
        "_chem_comp_atom.pdbx_leaving_atom_flag",
        "_chem_comp_atom.pdbx_model_Cartn_x_ideal",
        "_chem_comp_atom.pdbx_model_Cartn_y_ideal",
        "_chem_comp_atom.pdbx_model_Cartn_z_ideal",
    ]

    lines.extend(
        f"{comp} {atom_id} {symbol} {charge} {leaving} "
        f"{x:.3f} {y:.3f} {z:.3f}"
        for comp, atom_id, symbol, charge, leaving, x, y, z in atom_records
    )

    lines.extend(
        [
            "#",
            "loop_",
            "_chem_comp_bond.comp_id",
            "_chem_comp_bond.atom_id_1",
            "_chem_comp_bond.atom_id_2",
            "_chem_comp_bond.value_order",
            "_chem_comp_bond.pdbx_aromatic_flag",
        ]
    )

    lines.extend(
        f"{comp} {atom_1} {atom_2} {order} {aromatic}"
        for comp, atom_1, atom_2, order, aromatic in bond_records
    )

    lines.append("#")

    metadata = {
        "formula": formula,
        "weight": weight,
        "smiles": canonical_smiles,
        "inchi": inchi,
    }

    return "\n".join(lines) + "\n", mol, metadata


def save_structure_image(mol: Chem.Mol, output_file: Path, comp_id: str) -> None:
    """Save a 2D depiction of the molecule as PNG."""
    mol_2d = Chem.RemoveHs(mol)
    AllChem.Compute2DCoords(mol_2d)

    image = Draw.MolToImage(
        mol_2d,
        size=(800, 600),
        kekulize=True,
        legend=comp_id,
    )
    image.save(output_file)


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(
        description="Convert a SMILES string to a CCD-style CIF for AlphaFold 3."
    )
    parser.add_argument("--smiles", required=True, help="Ligand SMILES string.")
    parser.add_argument("--comp-id", required=True, help="Component ID, e.g. LIG.")
    parser.add_argument("--name", default="?", help="Human-readable compound name.")
    parser.add_argument(
        "--formula",
        default=None,
        help="Optional molecular formula. Calculated automatically if omitted.",
    )
    parser.add_argument(
        "--weight",
        type=float,
        default=None,
        help="Optional molecular weight. Exact mass is calculated if omitted.",
    )
    parser.add_argument(
        "-o",
        "--output-dir",
        type=Path,
        default=Path("."),
        help="Output directory. Default: current directory.",
    )
    parser.add_argument(
        "--seed",
        type=int,
        default=0xF00D,
        help="RDKit ETKDG random seed. Default: 61453.",
    )
    parser.add_argument(
        "--no-image",
        action="store_true",
        help="Do not create a PNG structure image.",
    )
    return parser.parse_args()


def main() -> int:
    args = parse_args()

    comp_id = args.comp_id.strip().upper()
    args.output_dir.mkdir(parents=True, exist_ok=True)

    cif_file = args.output_dir / f"{comp_id}.cif"
    png_file = args.output_dir / f"{comp_id}.png"

    try:
        cif_text, mol, metadata = smiles_to_ccd_cif(
            smiles=args.smiles,
            comp_id=comp_id,
            comp_name=args.name,
            formula=args.formula,
            weight=args.weight,
            random_seed=args.seed,
        )

        cif_file.write_text(cif_text, encoding="utf-8")

        if not args.no_image:
            save_structure_image(mol, png_file, comp_id)

    except Exception as exc:
        print(f"ERROR: {exc}", file=sys.stderr)
        return 1

    print(f"Component : {comp_id}")
    print(f"Formula   : {metadata['formula']}")
    print(f"Exact MW  : {metadata['weight']:.6f}")
    print(f"SMILES    : {metadata['smiles']}")
    print(f"InChI     : {metadata['inchi']}")
    print(f"CIF       : {cif_file.resolve()}")
    if not args.no_image:
        print(f"Image     : {png_file.resolve()}")

    return 0


if __name__ == "__main__":
    raise SystemExit(main())
