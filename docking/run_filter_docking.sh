#!/bin/bash

set -euo pipefail

# ============================================================
# Get folder containing this script
# ============================================================

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"


# ============================================================
# PARAMETERS - EDIT HERE
# ============================================================

CUTOFF="-9.0"


# ============================================================
# Input / Output
# ============================================================

DOCKING_FILE="$SCRIPT_DIR/data_output/docking_results/docking_scores.csv"

SMILES_FILE="$SCRIPT_DIR/data_input/Slc38a9_pKa_filtering.csv"

OUTPUT_DIR="$SCRIPT_DIR/data_output"

PYTHON_SCRIPT="$SCRIPT_DIR/utils/filter_docking_results.py"


# ============================================================
# Check input files
# ============================================================

if [ ! -f "$DOCKING_FILE" ]; then
    echo "[ERROR] Docking results not found:"
    echo "$DOCKING_FILE"
    exit 1
fi

if [ ! -f "$SMILES_FILE" ]; then
    echo "[ERROR] SMILES file not found:"
    echo "$SMILES_FILE"
    exit 1
fi

if [ ! -f "$PYTHON_SCRIPT" ]; then
    echo "[ERROR] Python script not found:"
    echo "$PYTHON_SCRIPT"
    exit 1
fi


# ============================================================
# Create output directory
# ============================================================

mkdir -p "$OUTPUT_DIR"


# ============================================================
# Run
# ============================================================

echo "============================================================"
echo " AutoDock Vina - Filter Docking Results"
echo "============================================================"
echo
echo "Docking results:"
echo "$DOCKING_FILE"
echo
echo "SMILES:"
echo "$SMILES_FILE"
echo
echo "Affinity cutoff: $CUTOFF kcal/mol"
echo "Keeping: Affinity < $CUTOFF"
echo


python "$PYTHON_SCRIPT" \
    --docking-file "$DOCKING_FILE" \
    --smiles-file "$SMILES_FILE" \
    --output-dir "$OUTPUT_DIR" \
    --cutoff "$CUTOFF"


# ============================================================
# Finished
# ============================================================

echo
echo "============================================================"
echo "[DONE] Filtering completed"
echo "============================================================"
echo
echo "Output:"
echo "$OUTPUT_DIR/filtered_molecules.smi"
echo "$OUTPUT_DIR/filtered_molecules_with_scores.csv"
echo