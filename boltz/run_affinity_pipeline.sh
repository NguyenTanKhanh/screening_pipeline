#!/bin/bash

set -eo pipefail


# =============================================================================
# BOLTZ-2 AFFINITY SCREENING PIPELINE
#
# STEP 1: Create protein-only YAML
# STEP 2: Generate MSA using Boltz MSA server
# STEP 3: Fix generated MSA
# STEP 4: Generate protein-ligand affinity YAMLs
# STEP 5: Run Boltz-2 for all ligands
# STEP 6: Extract affinity results
#
# Run:
#
#   bash run_affinity_pipeline.sh
#
# =============================================================================


# =============================================================================
# DIRECTORY SETUP
# =============================================================================

# Directory containing this script:
#
# Automatic_screening_pipeline/boltz/

SCRIPT_DIR="$(
    cd "$(dirname "${BASH_SOURCE[0]}")"
    pwd
)"

cd "$SCRIPT_DIR"


# Parent directory:
#
# Automatic_screening_pipeline/

PIPELINE_ROOT="$(
    cd "$SCRIPT_DIR/.."
    pwd
)"


# =============================================================================
# CONFIG
# =============================================================================


# -----------------------------------------------------------------------------
# INPUT CSV
#
# Resolves to:
#
# /media/disk1/khanh/Automatic_screening_pipeline/
# docking/data_output/filtered_molecules_with_scores.csv
#
# Because it is relative to PIPELINE_ROOT, the whole
# Automatic_screening_pipeline folder can be copied elsewhere.
# -----------------------------------------------------------------------------

CSV_FILE="$PIPELINE_ROOT/docking/data_output/filtered_molecules_with_scores.csv"

LIGAND_COL="mol_ID"

SMILES_COL="SMILES"


# -----------------------------------------------------------------------------
# PROTEIN SEQUENCE
# -----------------------------------------------------------------------------

PROTEIN_SEQUENCE="DDKMAHHTLLLGSGHVGLRNLGNTCFLNAVLQCLSSTRPLRDFCLRRDFRQEVPGGGRAQELTEAFADVIGALWHPDSCEAVNPTRFRAVFQKYVPSFSGYSQQDAQEFLKLLMERLHLEINRRGRRAPPILANGPVPSPPRRGGALLEEPELSDDDRANLMWKRYLEREDSKIVDLFVGQLKSCLKCQACGYRSTTFEVFCDLSLPIPKKGFAGGKVSLRDCFNLFTKEEELESENAPVCDRCRQKTRSTKKLTVQRFPRILVLHLNRFSASRGSIKKSSVGVDFPLQRLSLGDFASDKAGSPVYQLYALCNHSGSVHYGHYTALCRCQTGWHVYNDSRVSPVSENQVASSEGYVLFYQLMQEPPRCL"


# -----------------------------------------------------------------------------
# OPTIONAL POCKET RESIDUE CONSTRAINT
#
# No constraint:
#   RESIDUES=""
#
# One residue:
#   RESIDUES="383"
#
# Multiple residues:
#   RESIDUES="131,194,383"
# -----------------------------------------------------------------------------

RESIDUES=""

MAX_DISTANCE="4.0"


# -----------------------------------------------------------------------------
# GPU
# -----------------------------------------------------------------------------

GPU_ID="1"


# =============================================================================
# OUTPUT PATHS
# =============================================================================

RESULT_ROOT="$SCRIPT_DIR/boltz_prediction_result"

PROTEIN_ONLY_YAML="$RESULT_ROOT/protein_only.yaml"

MSA_OUTPUT_DIR="$RESULT_ROOT/msa"

YAML_DIR="$RESULT_ROOT/affinity_yamls"

RESULTS_DIR="$RESULT_ROOT/results_affinity"

SUMMARY_FILE="$RESULT_ROOT/summary_affinity_results.csv"


# =============================================================================
# MSA FILE
# =============================================================================

MSA_FILE="$MSA_OUTPUT_DIR/boltz_results_protein_only/msa/protein_only_unpaired_tmp_env/uniref.a3m"


# =============================================================================
# UTILITY SCRIPTS
# =============================================================================

UTILS_DIR="$SCRIPT_DIR/utils"

GENERATE_SCRIPT="$UTILS_DIR/generate_affinity_yamls.py"

FIX_MSA_SCRIPT="$UTILS_DIR/fix_msa.py"

EXTRACT_SCRIPT="$UTILS_DIR/extract_affinity_json_results.py"


# =============================================================================
# CREATE OUTPUT ROOT
# =============================================================================

mkdir -p "$RESULT_ROOT"


# =============================================================================
# FUNCTION: EXTRACT CURRENT RESULTS
# =============================================================================

extract_current_results() {

    echo
    echo "======================================================================"
    echo " EXTRACTING CURRENTLY AVAILABLE RESULTS"
    echo "======================================================================"
    echo

    if [[ ! -f "$EXTRACT_SCRIPT" ]]; then

        echo "[WARNING] Extraction script not found:"
        echo "  $EXTRACT_SCRIPT"

        return 0

    fi

    python "$EXTRACT_SCRIPT" \
        --base-dir "$SCRIPT_DIR" \
        || true

    if [[ -f "$SUMMARY_FILE" ]]; then

        echo
        echo "[INFO] Summary:"
        echo "  $SUMMARY_FILE"

    fi
}


# =============================================================================
# FUNCTION: HANDLE ERROR / CTRL+C
# =============================================================================

pipeline_interrupted() {

    EXIT_CODE=$?

    trap - ERR INT TERM

    echo
    echo "======================================================================"
    echo " PIPELINE STOPPED"
    echo "======================================================================"
    echo
    echo "[INFO] Extracting completed predictions..."

    extract_current_results

    echo
    echo "[INFO] Existing results were preserved."
    echo
    echo "[INFO] Partial summary:"
    echo "  $SUMMARY_FILE"
    echo

    exit "$EXIT_CODE"
}


trap pipeline_interrupted ERR INT TERM


# =============================================================================
# HEADER
# =============================================================================

echo
echo "======================================================================"
echo " BOLTZ-2 AFFINITY SCREENING PIPELINE"
echo "======================================================================"

echo
echo "Pipeline root:"
echo "  $PIPELINE_ROOT"

echo
echo "Boltz directory:"
echo "  $SCRIPT_DIR"

echo
echo "Ligand CSV:"
echo "  $CSV_FILE"

echo
echo "Boltz output:"
echo "  $RESULT_ROOT"

echo
echo "GPU:"
echo "  $GPU_ID"

echo

if [[ -n "$RESIDUES" ]]; then

    echo "Pocket constraint:"
    echo "  Residues     : $RESIDUES"
    echo "  Max distance : $MAX_DISTANCE"

else

    echo "Pocket constraint:"
    echo "  None"

fi

echo
echo "======================================================================"


# =============================================================================
# CHECK INPUT CSV
# =============================================================================

if [[ ! -f "$CSV_FILE" ]]; then

    echo
    echo "[ERROR] Ligand CSV not found:"
    echo "  $CSV_FILE"

    false

fi


echo
echo "[OK] Ligand CSV found:"
echo "  $CSV_FILE"


# =============================================================================
# SHOW NUMBER OF INPUT LIGANDS
# =============================================================================

echo

python - <<PY
import pandas as pd

csv_file = r"$CSV_FILE"

df = pd.read_csv(csv_file)

print("[INFO] Input CSV rows :", len(df))
print("[INFO] Columns        :", df.columns.tolist())

if "$LIGAND_COL" in df.columns:
    print("[INFO] Ligand IDs:")
    for x in df["$LIGAND_COL"]:
        print("       ", x)
PY


# =============================================================================
# CHECK PROTEIN
# =============================================================================

if [[ -z "$PROTEIN_SEQUENCE" ]]; then

    echo
    echo "[ERROR] PROTEIN_SEQUENCE is empty."

    false

fi


# =============================================================================
# CHECK UTILITIES
# =============================================================================

if [[ ! -f "$GENERATE_SCRIPT" ]]; then

    echo
    echo "[ERROR] Missing:"
    echo "  $GENERATE_SCRIPT"

    false

fi


if [[ ! -f "$FIX_MSA_SCRIPT" ]]; then

    echo
    echo "[ERROR] Missing:"
    echo "  $FIX_MSA_SCRIPT"

    false

fi


if [[ ! -f "$EXTRACT_SCRIPT" ]]; then

    echo
    echo "[ERROR] Missing:"
    echo "  $EXTRACT_SCRIPT"

    false

fi


# =============================================================================
# CHECK BOLTZ
# =============================================================================

if ! command -v boltz >/dev/null 2>&1; then

    echo
    echo "[ERROR] 'boltz' command not found."
    echo
    echo "Activate boltz_env first."

    false

fi


# =============================================================================
# GPU
# =============================================================================

export CUDA_VISIBLE_DEVICES="$GPU_ID"


# =============================================================================
# STEP 1/6
# CREATE PROTEIN-ONLY YAML
# =============================================================================

echo
echo "======================================================================"
echo " STEP 1/6 - CREATE PROTEIN-ONLY YAML"
echo "======================================================================"
echo


python "$GENERATE_SCRIPT" \
    --protein-seq "$PROTEIN_SEQUENCE" \
    --protein-only-yaml "$PROTEIN_ONLY_YAML"


if [[ ! -f "$PROTEIN_ONLY_YAML" ]]; then

    echo
    echo "[ERROR] protein_only.yaml was not created."

    false

fi


echo
echo "[OK] Created:"
echo "  $PROTEIN_ONLY_YAML"


# =============================================================================
# STEP 2/6
# GENERATE MSA
# =============================================================================

echo
echo "======================================================================"
echo " STEP 2/6 - GENERATE PROTEIN MSA"
echo "======================================================================"
echo


# Delete previous MSA because protein sequence may have changed.

if [[ -d "$MSA_OUTPUT_DIR" ]]; then

    echo "[INFO] Removing old MSA:"
    echo "  $MSA_OUTPUT_DIR"

    rm -rf "$MSA_OUTPUT_DIR"

fi


mkdir -p "$MSA_OUTPUT_DIR"


CUDA_VISIBLE_DEVICES="$GPU_ID" \
boltz predict "$PROTEIN_ONLY_YAML" \
    --use_msa_server \
    --use_potentials \
    --out_dir "$MSA_OUTPUT_DIR"


# =============================================================================
# CHECK MSA
# =============================================================================

if [[ ! -f "$MSA_FILE" ]]; then

    echo
    echo "[ERROR] MSA was not generated."
    echo
    echo "Expected:"
    echo "  $MSA_FILE"

    false

fi


if [[ ! -s "$MSA_FILE" ]]; then

    echo
    echo "[ERROR] MSA exists but is empty:"
    echo "  $MSA_FILE"

    false

fi


echo
echo "[OK] MSA generated:"
echo "  $MSA_FILE"


# =============================================================================
# STEP 3/6
# FIX MSA
# =============================================================================

echo
echo "======================================================================"
echo " STEP 3/6 - FIX MSA"
echo "======================================================================"
echo


python "$FIX_MSA_SCRIPT" \
    --msa-file "$MSA_FILE"


# =============================================================================
# VERIFY MSA
# =============================================================================

NULL_COUNT=$(
python - <<PY
from pathlib import Path

p = Path(r"$MSA_FILE")

data = p.read_bytes()

print(data.count(b"\x00"))
PY
)


if [[ "$NULL_COUNT" != "0" ]]; then

    echo
    echo "[ERROR] MSA still contains NULL byte(s):"
    echo "  $NULL_COUNT"

    false

fi


echo
echo "[OK] MSA verification passed."
echo "[OK] NULL bytes: 0"


# =============================================================================
# STEP 4/6
# GENERATE AFFINITY YAMLs
# =============================================================================

echo
echo "======================================================================"
echo " STEP 4/6 - GENERATE AFFINITY YAML FILES"
echo "======================================================================"
echo


# IMPORTANT:
#
# Completely remove old YAML directory.
#
# Therefore old ligand YAMLs cannot remain when the input CSV changes.

if [[ -d "$YAML_DIR" ]]; then

    echo "[INFO] Removing old affinity YAMLs:"
    echo "  $YAML_DIR"

    rm -rf "$YAML_DIR"

fi


mkdir -p "$YAML_DIR"


# =============================================================================
# OPTIONAL POCKET CONSTRAINT
# =============================================================================

RESIDUE_ARGS=()


if [[ -n "$RESIDUES" ]]; then

    RESIDUE_ARGS=(
        --pocket-residues "$RESIDUES"
        --max-distance "$MAX_DISTANCE"
    )

fi


# =============================================================================
# GENERATE YAMLs
# =============================================================================

python "$GENERATE_SCRIPT" \
    --csv "$CSV_FILE" \
    --outdir "$YAML_DIR" \
    --ligand-col "$LIGAND_COL" \
    --smiles-col "$SMILES_COL" \
    --protein-seq "$PROTEIN_SEQUENCE" \
    --msa-path "$MSA_FILE" \
    "${RESIDUE_ARGS[@]}"


# =============================================================================
# COUNT GENERATED YAMLs
# =============================================================================

YAML_COUNT=$(
    find "$YAML_DIR" \
        -maxdepth 1 \
        -type f \
        -name "affinity_*.yaml" \
        | wc -l
)


echo
echo "[INFO] Number of generated YAMLs:"
echo "  $YAML_COUNT"


if [[ "$YAML_COUNT" -eq 0 ]]; then

    echo
    echo "[ERROR] No affinity YAML files were generated."

    false

fi


echo
echo "[INFO] Generated YAML files:"

find "$YAML_DIR" \
    -maxdepth 1 \
    -type f \
    -name "affinity_*.yaml" \
    -printf "  %f\n"


# =============================================================================
# STEP 5/6
# RUN BOLTZ AFFINITY PREDICTIONS
# =============================================================================

echo
echo "======================================================================"
echo " STEP 5/6 - RUN BOLTZ-2 AFFINITY PREDICTIONS"
echo "======================================================================"
echo


mkdir -p "$RESULTS_DIR"


CURRENT=0


for YAML_FILE in "$YAML_DIR"/affinity_*.yaml
do

    CURRENT=$((CURRENT + 1))

    FILENAME="$(
        basename "$YAML_FILE"
    )"

    LIGAND_ID="${FILENAME%.yaml}"

    OUTDIR="$RESULTS_DIR/results_${LIGAND_ID}"


    echo
    echo "----------------------------------------------------------------------"
    echo " Ligand $CURRENT / $YAML_COUNT"
    echo "----------------------------------------------------------------------"

    echo
    echo "YAML:"
    echo "  $YAML_FILE"

    echo
    echo "Output:"
    echo "  $OUTDIR"

    echo


    # =========================================================================
    # RESUME
    # =========================================================================

    if [[ -d "$OUTDIR" ]]; then

        echo "[INFO] Existing result directory found."
        echo "[INFO] Skipping:"
        echo "  $FILENAME"

        continue

    fi


    # =========================================================================
    # BOLTZ PREDICTION
    # =========================================================================

    CUDA_VISIBLE_DEVICES="$GPU_ID" \
    boltz predict "$YAML_FILE" \
        --use_potentials \
        --out_dir "$OUTDIR"


    echo
    echo "[OK] Finished:"
    echo "  $FILENAME"


    # =========================================================================
    # UPDATE SUMMARY AFTER EACH LIGAND
    # =========================================================================

    python "$EXTRACT_SCRIPT" \
        --base-dir "$SCRIPT_DIR"

done


# =============================================================================
# STEP 6/6
# FINAL EXTRACTION
# =============================================================================

echo
echo "======================================================================"
echo " STEP 6/6 - FINAL RESULT EXTRACTION"
echo "======================================================================"
echo


extract_current_results


# =============================================================================
# NORMAL COMPLETION
# =============================================================================

trap - ERR INT TERM


echo
echo "======================================================================"
echo " PIPELINE COMPLETED SUCCESSFULLY"
echo "======================================================================"

echo
echo "Input CSV:"
echo "  $CSV_FILE"

echo
echo "Protein-only YAML:"
echo "  $PROTEIN_ONLY_YAML"

echo
echo "MSA:"
echo "  $MSA_FILE"

echo
echo "Affinity YAMLs:"
echo "  $YAML_DIR"

echo
echo "Affinity predictions:"
echo "  $RESULTS_DIR"

echo
echo "Summary:"
echo "  $SUMMARY_FILE"

echo
echo "======================================================================"