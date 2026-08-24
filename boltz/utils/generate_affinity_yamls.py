#!/usr/bin/env python3
# -*- coding: utf-8 -*-

"""
Generate per-ligand YAML files for Boltz-2 affinity prediction.

Features
--------
- Read ligand IDs and SMILES from CSV.
- Protein sequence supplied directly or from FASTA/text file.
- Clean and canonicalize SMILES using RDKit.
- Preserve stereochemistry.
- Generate one YAML per ligand.
- Optional pocket residue constraints.
- Optional precomputed MSA path.
- Protein-only YAML mode for MSA bootstrap.
- Invalid SMILES are logged to invalid_smiles_log.csv.

Example
-------
python utils/generate_affinity_yamls.py \
    --csv input.csv \
    --outdir affinity_yamls \
    --protein-seq "MKV..." \
    --msa-path "/path/to/uniref.a3m"

With pocket constraint:

python utils/generate_affinity_yamls.py \
    --csv input.csv \
    --outdir affinity_yamls \
    --protein-seq "MKV..." \
    --msa-path "/path/to/uniref.a3m" \
    --pocket-residues 383 \
    --max-distance 4.0

Protein-only mode:

python utils/generate_affinity_yamls.py \
    --protein-seq "MKV..." \
    --protein-only-yaml protein_only.yaml
"""

import argparse
from pathlib import Path
from typing import Optional, List

import pandas as pd
import yaml
from rdkit import Chem


# ============================================================
# Protein sequence
# ============================================================

def load_protein_sequence(
    protein_seq: Optional[str],
    protein_seq_file: Optional[str],
) -> str:

    if protein_seq:
        sequence = protein_seq.strip()

        if not sequence:
            raise ValueError("Protein sequence is empty.")

        return sequence

    if protein_seq_file:
        path = Path(protein_seq_file)

        if not path.exists():
            raise FileNotFoundError(
                f"Protein sequence file not found: {path}"
            )

        lines = path.read_text(
            encoding="utf-8"
        ).splitlines()

        # Remove FASTA header lines
        sequence = "".join(
            line.strip()
            for line in lines
            if line.strip() and not line.startswith(">")
        )

        if not sequence:
            raise ValueError(
                f"No protein sequence found in: {path}"
            )

        return sequence

    raise ValueError(
        "Must provide either --protein-seq "
        "or --protein-seq-file."
    )


# ============================================================
# SMILES processing
# ============================================================

def clean_smiles(smiles: str) -> str:
    """
    Replace malformed standalone carbon representations.
    """

    return (
        smiles
        .replace("[C]", "C")
        .replace("[c]", "c")
    )


def canonicalize_smiles(
    smiles: str,
) -> Optional[str]:

    cleaned = clean_smiles(smiles)

    try:

        mol = Chem.MolFromSmiles(cleaned)

        if mol is None:
            return None

        return Chem.MolToSmiles(
            mol,
            canonical=True,
            isomericSmiles=True,
        )

    except Exception:
        return None


# ============================================================
# Build protein + ligand YAML
# ============================================================

def build_yaml_content(
    protein_sequence: str,
    smiles: str,
    msa_path: Optional[str] = None,
    pocket_residues: Optional[List[int]] = None,
    max_distance: float = 4.0,
    force: bool = True,
) -> dict:

    protein_entry = {
        "id": "A",
        "sequence": protein_sequence,
    }

    if msa_path:
        protein_entry["msa"] = msa_path

    content = {
        "version": 1,
        "sequences": [
            {
                "protein": protein_entry
            },
            {
                "ligand": {
                    "id": "B",
                    "smiles": str(smiles),
                }
            },
        ],
    }

    # Optional pocket constraint
    if pocket_residues:

        content["constraints"] = [
            {
                "pocket": {
                    "binder": "B",
                    "contacts": [
                        ["A", residue]
                        for residue in pocket_residues
                    ],
                    "max_distance": max_distance,
                    "force": force,
                }
            }
        ]

    content["properties"] = [
        {
            "affinity": {
                "binder": "B"
            }
        }
    ]

    return content


# ============================================================
# Protein-only YAML
# ============================================================

def write_protein_only_yaml(
    protein_sequence: str,
    outpath: Path,
) -> None:

    content = {
        "version": 1,
        "sequences": [
            {
                "protein": {
                    "id": "A",
                    "sequence": protein_sequence,
                }
            }
        ],
    }

    outpath.parent.mkdir(
        parents=True,
        exist_ok=True,
    )

    with open(
        outpath,
        "w",
        encoding="utf-8",
    ) as f:

        yaml.dump(
            content,
            f,
            sort_keys=False,
            allow_unicode=True,
        )


# ============================================================
# Write ligand YAML
# ============================================================

def write_yaml(
    ligand_name: str,
    smiles: str,
    protein_sequence: str,
    outdir: Path,
    msa_path: Optional[str],
    pocket_residues: Optional[List[int]],
    max_distance: float,
    force: bool,
) -> None:

    yaml_content = build_yaml_content(
        protein_sequence=protein_sequence,
        smiles=smiles,
        msa_path=msa_path,
        pocket_residues=pocket_residues,
        max_distance=max_distance,
        force=force,
    )

    outpath = (
        outdir
        / f"affinity_{ligand_name}.yaml"
    )

    with open(
        outpath,
        "w",
        encoding="utf-8",
    ) as f:

        yaml.dump(
            yaml_content,
            f,
            sort_keys=False,
            allow_unicode=True,
        )


# ============================================================
# Generate ligand YAMLs
# ============================================================

def generate_affinity_yamls(
    csv_file: str,
    output_dir: str,
    ligand_col: str,
    smiles_col: str,
    protein_sequence: str,
    msa_path: Optional[str],
    pocket_residues: Optional[List[int]],
    max_distance: float,
    force: bool,
) -> None:

    csv_path = Path(csv_file)

    if not csv_path.exists():
        raise FileNotFoundError(
            f"Input CSV not found: {csv_path}"
        )

    df = pd.read_csv(csv_path)

    df.columns = [
        column.strip()
        for column in df.columns
    ]

    required = {
        ligand_col,
        smiles_col,
    }

    missing = required - set(df.columns)

    if missing:
        raise ValueError(
            f"Missing required column(s): {missing}"
        )

    outdir = Path(output_dir)

    outdir.mkdir(
        parents=True,
        exist_ok=True,
    )

    print("=" * 70)
    print("BOLTZ AFFINITY YAML GENERATOR")
    print("=" * 70)
    print(f"Input CSV  : {csv_path}")
    print(f"Output     : {outdir}")

    if msa_path:
        print(f"MSA        : {msa_path}")
    else:
        print("MSA        : None")

    if pocket_residues:
        print(
            "Residues   : "
            + ",".join(
                str(x)
                for x in pocket_residues
            )
        )
        print(f"Distance   : {max_distance}")
        print(f"Force      : {force}")
    else:
        print("Constraint : None")

    print("=" * 70)

    generated_count = 0
    error_log = []

    for _, row in df.iterrows():

        ligand = str(
            row[ligand_col]
        ).strip()

        smiles_value = row[
            smiles_col
        ]

        # Empty SMILES
        if (
            pd.isna(smiles_value)
            or str(smiles_value).strip() == ""
        ):

            error_log.append(
                {
                    "Ligand": ligand,
                    "SMILES": smiles_value,
                    "Reason": "NaN/empty",
                }
            )

            continue

        original_smiles = str(
            smiles_value
        ).strip()

        canonical_smiles = (
            canonicalize_smiles(
                original_smiles
            )
        )

        if canonical_smiles is None:

            error_log.append(
                {
                    "Ligand": ligand,
                    "SMILES": original_smiles,
                    "Reason": "Invalid SMILES",
                }
            )

            continue

        write_yaml(
            ligand_name=ligand,
            smiles=canonical_smiles,
            protein_sequence=protein_sequence,
            outdir=outdir,
            msa_path=msa_path,
            pocket_residues=pocket_residues,
            max_distance=max_distance,
            force=force,
        )

        generated_count += 1

    print()
    print("=" * 70)
    print("DONE")
    print("=" * 70)
    print(f"Generated : {generated_count}")
    print(f"Skipped   : {len(error_log)}")
    print(f"Output    : {outdir}")

    if error_log:

        error_path = (
            outdir
            / "invalid_smiles_log.csv"
        )

        pd.DataFrame(
            error_log
        ).to_csv(
            error_path,
            index=False,
        )

        print(
            f"Invalid SMILES log: "
            f"{error_path}"
        )


# ============================================================
# CLI
# ============================================================

if __name__ == "__main__":

    parser = argparse.ArgumentParser(
        description=(
            "Generate Boltz-2 affinity YAML files."
        )
    )

    parser.add_argument(
        "--csv",
        default=None,
        help="Input ligand CSV.",
    )

    parser.add_argument(
        "--outdir",
        default=None,
        help="Output directory for ligand YAML files.",
    )

    parser.add_argument(
        "--ligand-col",
        default="mol_ID",
        help="Ligand ID column. Default: mol_ID",
    )

    parser.add_argument(
        "--smiles-col",
        default="SMILES",
        help="SMILES column. Default: SMILES",
    )

    # Protein sequence
    protein_group = (
        parser.add_mutually_exclusive_group(
            required=True
        )
    )

    protein_group.add_argument(
        "--protein-seq",
        help="Protein sequence.",
    )

    protein_group.add_argument(
        "--protein-seq-file",
        help="Protein FASTA/text file.",
    )

    # MSA
    parser.add_argument(
        "--msa-path",
        default=None,
        help=(
            "Path written into protein.msa "
            "for every generated ligand YAML."
        ),
    )

    # Pocket constraint
    parser.add_argument(
        "--pocket-residues",
        default=None,
        help=(
            "Comma-separated residues, "
            "e.g. 131,194,383."
        ),
    )

    parser.add_argument(
        "--max-distance",
        type=float,
        default=4.0,
        help="Pocket max_distance. Default: 4.0",
    )

    parser.add_argument(
        "--no-force",
        action="store_true",
        help="Set pocket force=false.",
    )

    # Protein-only mode
    parser.add_argument(
        "--protein-only-yaml",
        default=None,
        help=(
            "Create a protein-only YAML "
            "and exit."
        ),
    )

    args = parser.parse_args()

    protein_sequence = (
        load_protein_sequence(
            args.protein_seq,
            args.protein_seq_file,
        )
    )

    # Protein-only mode
    if args.protein_only_yaml:

        outpath = Path(
            args.protein_only_yaml
        )

        write_protein_only_yaml(
            protein_sequence,
            outpath,
        )

        print(
            "[DONE] Protein-only YAML created:"
        )

        print(
            f"  {outpath}"
        )

        raise SystemExit(0)

    # Normal ligand mode
    if not args.csv or not args.outdir:

        raise SystemExit(
            "--csv and --outdir are required "
            "unless using --protein-only-yaml."
        )

    pocket_residues = None

    if args.pocket_residues:

        pocket_residues = [
            int(value.strip())
            for value
            in args.pocket_residues.split(",")
            if value.strip()
        ]

    generate_affinity_yamls(
        csv_file=args.csv,
        output_dir=args.outdir,
        ligand_col=args.ligand_col,
        smiles_col=args.smiles_col,
        protein_sequence=protein_sequence,
        msa_path=args.msa_path,
        pocket_residues=pocket_residues,
        max_distance=args.max_distance,
        force=not args.no_force,
    )