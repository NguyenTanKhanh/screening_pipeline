#!/usr/bin/env python3
# -*- coding: utf-8 -*-

"""
Remove NULL bytes from a Boltz-generated A3M MSA file.

The original MSA is backed up before modification.

Usage:

python utils/fix_msa.py \
    --msa-file /path/to/uniref.a3m
"""

import argparse
import shutil
from pathlib import Path


def fix_msa(
    msa_file: str,
) -> None:

    msa_path = Path(msa_file)

    print("=" * 70)
    print("FIX BOLTZ MSA")
    print("=" * 70)
    print(f"MSA: {msa_path}")
    print()

    # ========================================================
    # Check file
    # ========================================================

    if not msa_path.exists():

        raise FileNotFoundError(
            f"MSA file not found: {msa_path}"
        )

    if not msa_path.is_file():

        raise ValueError(
            f"MSA path is not a file: {msa_path}"
        )

    # ========================================================
    # Read
    # ========================================================

    data = msa_path.read_bytes()

    null_count_before = (
        data.count(b"\x00")
    )

    print(
        f"Original size : "
        f"{len(data):,} bytes"
    )

    print(
        f"NULL bytes    : "
        f"{null_count_before}"
    )

    # ========================================================
    # Nothing to fix
    # ========================================================

    if null_count_before == 0:

        print()
        print(
            "[OK] No NULL bytes detected."
        )

        print(
            "[OK] MSA does not require fixing."
        )

        return

    # ========================================================
    # Backup
    # ========================================================

    backup_path = Path(
        str(msa_path)
        + ".backup"
    )

    shutil.copy2(
        msa_path,
        backup_path,
    )

    print()
    print(
        "[INFO] Backup created:"
    )

    print(
        f"  {backup_path}"
    )

    # ========================================================
    # Remove NULL bytes
    # ========================================================

    fixed_data = data.replace(
        b"\x00",
        b"",
    )

    msa_path.write_bytes(
        fixed_data
    )

    # ========================================================
    # Verify
    # ========================================================

    verify_data = (
        msa_path.read_bytes()
    )

    null_count_after = (
        verify_data.count(
            b"\x00"
        )
    )

    print()
    print(
        f"Fixed size    : "
        f"{len(verify_data):,} bytes"
    )

    print(
        f"NULL bytes    : "
        f"{null_count_after}"
    )

    if null_count_after != 0:

        raise RuntimeError(
            "NULL bytes are still present "
            "after fixing MSA."
        )

    # ========================================================
    # Basic text check
    # ========================================================

    try:

        text = verify_data.decode(
            "utf-8"
        )

    except UnicodeDecodeError as error:

        raise RuntimeError(
            "Fixed MSA is not valid UTF-8."
        ) from error

    # First non-empty character should normally be >
    first_nonempty = None

    for line in text.splitlines():

        if line.strip():

            first_nonempty = (
                line.strip()
            )

            break

    if (
        first_nonempty is None
        or not first_nonempty.startswith(">")
    ):

        raise RuntimeError(
            "MSA does not start with a "
            "FASTA/A3M header ('>')."
        )

    print()
    print(
        "[OK] MSA fixed successfully."
    )

    print(
        "[OK] MSA contains no NULL bytes."
    )

    print("=" * 70)


if __name__ == "__main__":

    parser = argparse.ArgumentParser(
        description=(
            "Remove NULL bytes from "
            "a Boltz A3M MSA file."
        )
    )

    parser.add_argument(
        "--msa-file",
        required=True,
        help="Path to the A3M file.",
    )

    args = parser.parse_args()

    fix_msa(
        args.msa_file
    )