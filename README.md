# SMILES to AlphaFold 3 custom ligand CIF

Small standalone Python utility that converts a **SMILES string** into a minimal
CCD-style CIF file that can be used as a **custom ligand in AlphaFold 3**.

The script uses RDKit to:

- parse the SMILES string
- add explicit hydrogens
- generate and optimize a 3D conformer
- calculate the molecular formula and exact molecular mass
- assign atom names such as `C1`, `C2`, `O1`, `H1`, ...
- write atom coordinates and bond information to a CCD-style CIF
- save a 2D PNG depiction for visual verification

## Requirements

The only required Python dependency is RDKit.

Using conda/mamba:

```bash
conda install -c conda-forge rdkit
```

or, in an appropriate Python environment:

```bash
pip install rdkit
```

## Usage

```bash
python smiles_to_ccd.py \
    --smiles "CC(=O)Oc1ccccc1C(=O)O" \
    --comp-id ASPIRIN \
    --name "Aspirin"
```

This creates:

```text
ASPIRIN.cif
ASPIRIN.png
```

A different output directory can be selected with:

```bash
python smiles_to_ccd.py \
    --smiles "CCO" \
    --comp-id ETO \
    --name "Ethanol" \
    --output-dir ligands/
```

## Options

```text
--smiles        SMILES string [required]
--comp-id       Component identifier used in the CIF [required]. Longer IDs are recommended
                for custom ligands to avoid collisions with standard residue/component codes.
--name          Human-readable compound name
--formula       Override automatically calculated molecular formula
--weight        Override automatically calculated exact molecular mass
--output-dir    Output directory
--seed          ETKDG random seed
--no-image      Skip PNG generation
```

## AlphaFold 3

The generated `.cif` file contains a minimal `_chem_comp`, `_chem_comp_atom`,
and `_chem_comp_bond` definition suitable for supplying a custom ligand
definition to AlphaFold 3.

Example workflow:

```text
SMILES
  ↓
smiles_to_ccd.py
  ↓
LIG.cif
  ↓
AlphaFold 3 custom CCD / ligand definition
```

The generated PNG is intended only as a quick visual sanity check of the
interpreted molecular structure.

## Notes

The 3D conformer is generated with RDKit ETKDGv3. Geometry optimization uses
MMFF when parameters are available and otherwise falls back to UFF.

Formal atomic charges are taken directly from the RDKit molecular graph. This
is relevant for charged ligands and is preferable to assigning every atom a
charge of zero.

The script does not attempt to reproduce an existing PDB Chemical Component
Dictionary atom naming scheme. Atom identifiers are generated sequentially by
element (`C1`, `C2`, `N1`, `O1`, ...).

For stereochemically defined ligands, use a SMILES string containing the
appropriate stereochemical annotations.

## Example

```bash
python smiles_to_ccd.py \
    --smiles "C[C@H](N)C(=O)O" \
    --comp-id ALA2 \
    --name "L-alanine"
```

# Reference

If you use this code within your work, please cite us!

[![DOI](https://zenodo.org/badge/1394643775.svg)](https://doi.org/10.5281/zenodo.23034811)

## License

MIT License. The license file is maintained via GitHub.
