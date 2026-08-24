#!/usr/bin/env python3

import argparse
from pathlib import Path

import pandas as pd


def filter_docking_results(
    docking_file,
    smiles_file,
    output_dir,
    cutoff,
    docking_id_col="Ligand",
    docking_score_col="Best Affinity",
    smiles_id_col="mol_ID",
    smiles_smiles_col="SMILES",
):
    docking_file = Path(docking_file)
    smiles_file = Path(smiles_file)
    output_dir = Path(output_dir)

    output_dir.mkdir(parents=True, exist_ok=True)

    out_smi = output_dir / "filtered_molecules.smi"
    out_csv = output_dir / "filtered_molecules_with_scores.csv"

    # ============================================================
    # 1. Load docking results
    # ============================================================

    print("============================================================")
    print("Filtering AutoDock Vina results")
    print("============================================================")
    print(f"Docking results : {docking_file}")
    print(f"SMILES file     : {smiles_file}")
    print(f"Cutoff          : {cutoff}")
    print(f"Condition       : Affinity < {cutoff}")
    print()

    df = pd.read_csv(docking_file)
    df.columns = [c.strip() for c in df.columns]

    required_docking = {docking_id_col, docking_score_col}
    missing = required_docking - set(df.columns)

    if missing:
        raise ValueError(
            f"Missing columns in {docking_file}: {missing}. "
            f"Expected at least {required_docking}"
        )

    # ============================================================
    # 2. Convert docking score to numeric
    # ============================================================

    score_series = (
        df[docking_score_col]
        .astype(str)
        .str.strip()
        .str.replace(",", ".", regex=False)
    )

    df[docking_score_col] = pd.to_numeric(
        score_series,
        errors="coerce",
    )

    before = len(df)

    df = df.dropna(
        subset=[docking_score_col]
    ).copy()

    after = len(df)

    if after < before:
        print(
            f"[WARNING] Dropped {before - after} rows "
            f"with invalid '{docking_score_col}' values."
        )

    print(f"[INFO] Valid docking results: {len(df)}")

    # ============================================================
    # 3. Filter based on affinity cutoff
    # More negative = better
    # ============================================================

    df = df[
        df[docking_score_col] < cutoff
    ].copy()

    print(
        f"[INFO] Molecules passing cutoff "
        f"(< {cutoff} kcal/mol): {len(df)}"
    )

    # Standardize molecule ID
    df.rename(
        columns={docking_id_col: "mol_ID"},
        inplace=True,
    )

    # Convert ID to string
    df["mol_ID"] = df["mol_ID"].astype(str).str.strip()

    # ============================================================
    # 4. Load original SMILES
    # ============================================================

    df_sm = pd.read_csv(smiles_file)
    df_sm.columns = [c.strip() for c in df_sm.columns]

    required_smiles = {
        smiles_id_col,
        smiles_smiles_col,
    }

    missing_sm = required_smiles - set(df_sm.columns)

    if missing_sm:
        raise ValueError(
            f"Missing columns in {smiles_file}: {missing_sm}. "
            f"Expected at least {required_smiles}"
        )

    # Keep only required columns
    df_sm = df_sm[
        [smiles_smiles_col, smiles_id_col]
    ].copy()

    # Standardize names
    df_sm.rename(
        columns={
            smiles_id_col: "mol_ID",
            smiles_smiles_col: "SMILES",
        },
        inplace=True,
    )

    df_sm["mol_ID"] = (
        df_sm["mol_ID"]
        .astype(str)
        .str.strip()
    )

    # ============================================================
    # 5. Remove duplicate molecule IDs
    # ============================================================

    duplicated = df_sm.duplicated(
        subset=["mol_ID"]
    ).sum()

    if duplicated > 0:
        print(
            f"[WARNING] Found {duplicated} duplicated mol_ID(s) "
            "in SMILES file. Keeping first occurrence."
        )

        df_sm = df_sm.drop_duplicates(
            subset=["mol_ID"],
            keep="first",
        ).copy()

    # ============================================================
    # 6. Merge docking results with SMILES
    # ============================================================

    df_merged = pd.merge(
        df,
        df_sm,
        on="mol_ID",
        how="left",
    )

    # Rename score column
    df_merged.rename(
        columns={
            docking_score_col: "Affinity"
        },
        inplace=True,
    )

    # Sort: strongest docking first
    df_merged = df_merged.sort_values(
        by="Affinity",
        ascending=True,
    )

    # ============================================================
    # 7. Check missing SMILES
    # ============================================================

    missing_smiles = (
        df_merged["SMILES"].isna().sum()
    )

    if missing_smiles > 0:
        print(
            f"[WARNING] Missing SMILES for "
            f"{missing_smiles} molecule(s). "
            "Check mol_ID mapping."
        )

    # ============================================================
    # 8. Write .smi
    #
    # Format:
    # SMILES    mol_ID
    # ============================================================

    smi_df = (
        df_merged[
            ["SMILES", "mol_ID"]
        ]
        .dropna(subset=["SMILES"])
        .copy()
    )

    smi_df.to_csv(
        out_smi,
        sep="\t",
        index=False,
        header=False,
    )

    # ============================================================
    # 9. Write CSV
    #
    # Columns:
    # SMILES, mol_ID, Affinity
    # ============================================================

    out_df = df_merged[
        ["SMILES", "mol_ID", "Affinity"]
    ].copy()

    out_df.to_csv(
        out_csv,
        index=False,
    )

    # ============================================================
    # Summary
    # ============================================================

    print()
    print("============================================================")
    print("DONE")
    print("============================================================")
    print(f"Cutoff                 : {cutoff}")
    print(f"Selected molecules     : {len(out_df)}")
    print(f"Molecules with SMILES  : {len(smi_df)}")
    print(f"Missing SMILES         : {missing_smiles}")
    print()
    print("Output files:")
    print(f"  {out_smi}")
    print(f"  {out_csv}")
    print()

    return out_csv, out_smi


def main():
    parser = argparse.ArgumentParser(
        description=(
            "Filter AutoDock Vina docking results "
            "using an affinity cutoff and merge "
            "selected compounds with SMILES."
        )
    )

    parser.add_argument(
        "--docking-file",
        required=True,
        help="CSV containing Vina docking scores.",
    )

    parser.add_argument(
        "--smiles-file",
        required=True,
        help="CSV containing SMILES and mol_ID.",
    )

    parser.add_argument(
        "--output-dir",
        required=True,
        help="Directory for filtered output files.",
    )

    parser.add_argument(
        "--cutoff",
        required=True,
        type=float,
        help=(
            "Affinity cutoff in kcal/mol. "
            "Only compounds with Affinity < cutoff "
            "are retained."
        ),
    )

    args = parser.parse_args()

    filter_docking_results(
        docking_file=args.docking_file,
        smiles_file=args.smiles_file,
        output_dir=args.output_dir,
        cutoff=args.cutoff,
    )


if __name__ == "__main__":
    main()