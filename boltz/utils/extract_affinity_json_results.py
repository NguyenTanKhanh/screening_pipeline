#!/usr/bin/env python3
# -*- coding: utf-8 -*-

"""
Extract Boltz-2 affinity prediction results from all per-ligand result
folders into a single summary CSV.

Usage:
  python utils/extract_affinity_json_results.py
"""

import argparse
import json
import csv
from pathlib import Path


def extract_results(base_dir: str = ".") -> None:

    # Path m?i:
    # boltz_prediction_result/
    # +-- results_affinity/
    #     +-- results_affinity_*/
    #         +-- boltz_results_affinity_*/
    #             +-- predictions/
    #                 +-- affinity_*/
    #                     +-- affinity_affinity_*.json

    pattern = (
        "boltz_prediction_result/"
        "results_affinity/"
        "results_affinity_*/"
        "boltz_results_affinity_*/"
        "predictions/"
        "affinity_*/"
        "affinity_affinity_*.json"
    )

    base_path = Path(base_dir)

    json_files = list(
        base_path.glob(pattern)
    )

    results = []

    for json_file in json_files:

        try:

            with open(
                json_file,
                "r"
            ) as f:

                data = json.load(f)

            mol_id = (
                json_file.stem
                .replace(
                    "affinity_affinity_",
                    ""
                )
            )

            results.append(
                (
                    mol_id,
                    data.get("affinity_pred_value"),
                    data.get("affinity_probability_binary"),
                    data.get("affinity_pred_value1"),
                    data.get("affinity_probability_binary1"),
                    data.get("affinity_pred_value2"),
                    data.get("affinity_probability_binary2"),
                )
            )

        except Exception as e:

            print(
                f"[ERROR] Could not process "
                f"{json_file}: {e}"
            )

    # ========================================================
    # Output m?i n?m trong boltz_prediction_result/
    # ========================================================

    output_file = (
        base_path
        / "boltz_prediction_result"
        / "summary_affinity_results.csv"
    )

    # Ð?m b?o folder t?n t?i
    output_file.parent.mkdir(
        parents=True,
        exist_ok=True
    )

    with open(
        output_file,
        "w",
        newline=""
    ) as f:

        writer = csv.writer(f)

        writer.writerow(
            [
                "mol_ID",
                "affinity_pred_value",
                "affinity_probability_binary",
                "affinity_pred_value1",
                "affinity_probability_binary1",
                "affinity_pred_value2",
                "affinity_probability_binary2",
            ]
        )

        writer.writerows(results)

    print(
        f"[DONE] Extracted "
        f"{len(results)} results "
        f"-> {output_file}"
    )


if __name__ == "__main__":

    parser = argparse.ArgumentParser(
        description=(
            "Extract Boltz-2 affinity results "
            "into a summary CSV."
        )
    )

    parser.add_argument(
        "--base-dir",
        default=".",
        help=(
            "Base directory to search from "
            "(default: current directory)"
        )
    )

    args = parser.parse_args()

    extract_results(
        args.base_dir
    )