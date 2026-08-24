#!/bin/bash

set -eo pipefail


# ============================================================
# AUTOMATIC AUTODOCK VINA SCREENING PIPELINE
#
# Run:
#   bash run_docking_pipeline.sh
# ============================================================


# ============================================================
# Get directory containing this script
# ============================================================

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"


# ============================================================
# PARAMETERS - EDIT HERE
# ============================================================

# Docking box center
CENTER_X="-50"
CENTER_Y="37"
CENTER_Z="75"

# Docking box size
SIZE_X="30"
SIZE_Y="22.5"
SIZE_Z="30"

# Vina exhaustiveness
EXHAUSTIVENESS="8"

# Filtering cutoff
# Keep molecules with:
# Best Affinity < CUTOFF
CUTOFF="-5.0"


# ============================================================
# CHECK CONDA ENVIRONMENT
# ============================================================

if [ "${CONDA_DEFAULT_ENV:-}" != "drugex" ]; then
    echo
    echo "[ERROR] Conda environment 'drugex' is not active."
    echo
    echo "Please run:"
    echo
    echo "    conda activate drugex"
    echo
    echo "Then run:"
    echo
    echo "    bash run_docking_pipeline.sh"
    echo
    exit 1
fi


# ============================================================
# OPEN BABEL ENVIRONMENT
# ============================================================

export BABEL_DATADIR="/home/khanh/.conda/envs/drugex/share/openbabel/3.1.0"


# ============================================================
# INPUT FILES
# ============================================================

SMILES_FILE="$SCRIPT_DIR/data_input/Slc38a9_pKa_filtering.csv"

RECEPTOR="$SCRIPT_DIR/data_input/slc38a9_human.pdbqt"


# ============================================================
# OUTPUT DIRECTORIES
# ============================================================

OUTPUT_DIR="$SCRIPT_DIR/data_output"

PDBQT_OUTPUT="$OUTPUT_DIR/pdbqt_converted"

LIGAND_DIR="$PDBQT_OUTPUT/pdbqt_files"

DOCKING_OUTPUT="$OUTPUT_DIR/docking_results"

DOCKING_SCORES="$DOCKING_OUTPUT/docking_scores.csv"


# ============================================================
# PYTHON / VINA FILES
# ============================================================

PREPARE_SCRIPT="$SCRIPT_DIR/utils/convert_csv_smiles_to_pdbqt.py"

VINA_SCRIPT="$SCRIPT_DIR/utils/vinascreen.py"

FILTER_SCRIPT="$SCRIPT_DIR/utils/filter_docking_results.py"

VINA_BINARY="$SCRIPT_DIR/utils/vina_1.2.7_linux_x86_64"


# ============================================================
# HEADER
# ============================================================

echo
echo "============================================================"
echo " AUTOMATIC AUTODOCK VINA SCREENING PIPELINE"
echo "============================================================"
echo
echo "Project:"
echo "  $SCRIPT_DIR"
echo
echo "Input SMILES:"
echo "  $SMILES_FILE"
echo
echo "Receptor:"
echo "  $RECEPTOR"
echo
echo "Docking box center:"
echo "  $CENTER_X $CENTER_Y $CENTER_Z"
echo
echo "Docking box size:"
echo "  $SIZE_X $SIZE_Y $SIZE_Z"
echo
echo "Exhaustiveness:"
echo "  $EXHAUSTIVENESS"
echo
echo "Affinity cutoff:"
echo "  $CUTOFF kcal/mol"
echo


# ============================================================
# CHECK INPUT FILES
# ============================================================

if [ ! -f "$SMILES_FILE" ]; then
    echo "[ERROR] SMILES input file not found:"
    echo "$SMILES_FILE"
    exit 1
fi


if [ ! -f "$RECEPTOR" ]; then
    echo "[ERROR] Receptor PDBQT not found:"
    echo "$RECEPTOR"
    exit 1
fi


if [ ! -f "$PREPARE_SCRIPT" ]; then
    echo "[ERROR] Ligand preparation script not found:"
    echo "$PREPARE_SCRIPT"
    exit 1
fi


if [ ! -f "$VINA_SCRIPT" ]; then
    echo "[ERROR] Vina screening script not found:"
    echo "$VINA_SCRIPT"
    exit 1
fi


if [ ! -f "$FILTER_SCRIPT" ]; then
    echo "[ERROR] Filtering script not found:"
    echo "$FILTER_SCRIPT"
    exit 1
fi


if [ ! -f "$VINA_BINARY" ]; then
    echo "[ERROR] Vina binary not found:"
    echo "$VINA_BINARY"
    exit 1
fi


# ============================================================
# CREATE OUTPUT DIRECTORIES
# ============================================================

mkdir -p "$OUTPUT_DIR"
mkdir -p "$PDBQT_OUTPUT"
mkdir -p "$DOCKING_OUTPUT"


# ============================================================
# CHECK OPEN BABEL
# ============================================================

echo "============================================================"
echo " CHECKING SOFTWARE"
echo "============================================================"
echo

if ! command -v obabel >/dev/null 2>&1; then
    echo "[ERROR] Open Babel not found."
    echo
    echo "Make sure you activated:"
    echo
    echo "    conda activate drugex"
    echo
    exit 1
fi

echo "[INFO] Open Babel:"
obabel -V
echo


# ============================================================
# MAKE VINA EXECUTABLE
# ============================================================

chmod +x "$VINA_BINARY"


# ============================================================
# STEP 1 - PREPARE LIGANDS
# ============================================================

echo
echo "============================================================"
echo " STEP 1/3 - PREPARE LIGANDS"
echo "============================================================"
echo

python "$PREPARE_SCRIPT" \
    --csv "$SMILES_FILE" \
    --output-dir "$PDBQT_OUTPUT"


# ============================================================
# CHECK STEP 1 OUTPUT
# ============================================================

if [ ! -d "$LIGAND_DIR" ]; then
    echo
    echo "[ERROR] Ligand directory was not created:"
    echo "$LIGAND_DIR"
    exit 1
fi


PDBQT_COUNT=$(find "$LIGAND_DIR" -maxdepth 1 -type f -name "*.pdbqt" | wc -l)


if [ "$PDBQT_COUNT" -eq 0 ]; then
    echo
    echo "[ERROR] No ligand PDBQT files were generated."
    exit 1
fi


echo
echo "[INFO] Ligands prepared: $PDBQT_COUNT"
echo


# ============================================================
# STEP 2 - AUTODOCK VINA SCREENING
# ============================================================

echo
echo "============================================================"
echo " STEP 2/3 - AUTODOCK VINA SCREENING"
echo "============================================================"
echo

python "$VINA_SCRIPT" \
    --receptor "$RECEPTOR" \
    --ligand_dir "$LIGAND_DIR" \
    --output_dir "$DOCKING_OUTPUT" \
    --center "$CENTER_X" "$CENTER_Y" "$CENTER_Z" \
    --size "$SIZE_X" "$SIZE_Y" "$SIZE_Z" \
    --exhaustiveness "$EXHAUSTIVENESS"


# ============================================================
# CHECK STEP 2 OUTPUT
# ============================================================

if [ ! -f "$DOCKING_SCORES" ]; then
    echo
    echo "[ERROR] Docking score file was not created:"
    echo "$DOCKING_SCORES"
    exit 1
fi


echo
echo "[INFO] Docking completed successfully."
echo
echo "[INFO] Raw docking scores:"
echo "$DOCKING_SCORES"
echo


# ============================================================
# STEP 3 - FILTER DOCKING RESULTS
# ============================================================

echo
echo "============================================================"
echo " STEP 3/3 - FILTER DOCKING RESULTS"
echo "============================================================"
echo

python "$FILTER_SCRIPT" \
    --docking-file "$DOCKING_SCORES" \
    --smiles-file "$SMILES_FILE" \
    --output-dir "$OUTPUT_DIR" \
    --cutoff "$CUTOFF"


# ============================================================
# FINAL OUTPUT
# ============================================================

echo
echo "============================================================"
echo " PIPELINE COMPLETED SUCCESSFULLY"
echo "============================================================"
echo

echo "Raw docking results:"
echo "  $DOCKING_OUTPUT"
echo

echo "Raw docking scores:"
echo "  $DOCKING_SCORES"
echo

echo "Filtered molecules:"
echo "  $OUTPUT_DIR/filtered_molecules.smi"
echo

echo "Filtered molecules with scores:"
echo "  $OUTPUT_DIR/filtered_molecules_with_scores.csv"
echo

echo "============================================================"