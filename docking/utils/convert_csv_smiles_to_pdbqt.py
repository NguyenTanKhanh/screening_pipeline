#!/usr/bin/env python3

import argparse
import os
import subprocess
from pathlib import Path

import pandas as pd


def convert_csv_smiles_to_pdbqt(csv_path, output_dir):
    pdbqt_dir = Path(output_dir) / "pdbqt_files"
    tmp_dir = Path(output_dir) / "temp_sdf"

    pdbqt_dir.mkdir(parents=True, exist_ok=True)
    tmp_dir.mkdir(parents=True, exist_ok=True)

    df = pd.read_csv(csv_path)

    if "SMILES" not in df.columns or "mol_ID" not in df.columns:
        raise ValueError(
            "CSV file must contain 'SMILES' and 'mol_ID' columns."
        )

    success = 0

    for _, row in df.iterrows():
        smi = str(row["SMILES"])
        mol_id = str(row["mol_ID"])

        sdf_path = tmp_dir / f"{mol_id}.sdf"
        pdbqt_path = pdbqt_dir / f"{mol_id}.pdbqt"

        try:
            print(f"\n[INFO] Processing {mol_id}")

            # Step 1: SMILES -> 3D SDF
            subprocess.run(
                [
                    "obabel",
                    f"-:{smi}",
                    "-O",
                    str(sdf_path),
                    "--gen3D",
                ],
                check=True,
            )

            # Step 2: Minimize structure and convert SDF -> PDBQT
            subprocess.run(
                [
                    "obabel",
                    str(sdf_path),
                    "-O",
                    str(pdbqt_path),
                    "--minimize",
                    "--ff",
                    "MMFF94",
                    "--steps",
                    "500",
                ],
                check=True,
            )

            if pdbqt_path.exists() and pdbqt_path.stat().st_size > 0:
                print(f"[OK] {mol_id}: converted successfully")
                success += 1
            else:
                print(f"[WARNING] {mol_id}: failed to generate PDBQT")

        except subprocess.CalledProcessError as e:
            print(f"[ERROR] {mol_id}: Open Babel failed - {e}")

        except Exception as e:
            print(f"[ERROR] {mol_id}: {e}")

    # Remove temporary SDF files
    for f in tmp_dir.glob("*.sdf"):
        f.unlink()

    try:
        tmp_dir.rmdir()
    except OSError:
        pass

    print("\n========================================")
    print(f"Done: {success}/{len(df)} molecules converted")
    print(f"PDBQT directory: {pdbqt_dir}")
    print("========================================")

    return str(pdbqt_dir)


def main():
    parser = argparse.ArgumentParser(
        description="Convert SMILES from CSV to minimized PDBQT files."
    )

    parser.add_argument(
        "--csv",
        required=True,
        help="Input CSV containing SMILES and mol_ID columns",
    )

    parser.add_argument(
        "--output-dir",
        required=True,
        help="Output directory",
    )

    args = parser.parse_args()

    convert_csv_smiles_to_pdbqt(
        csv_path=args.csv,
        output_dir=args.output_dir,
    )


if __name__ == "__main__":
    main()