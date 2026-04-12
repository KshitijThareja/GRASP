#!/bin/bash

# This script automates the construction of a DROSOPHILA pangenome graph using Minigraph.
# It will automatically create the initial graph if it's missing or invalid.

# --- Configuration ---
MINIGRAPH_EXEC="/home/kshitij/FY_Project/GRASP/minigraph/minigraph"
GENOMES_DIR="/home/kshitij/FY_Project/GRASP/drosophila_dataset"

# The reference genome is in a subdirectory
REF_GENOME_REL_PATH="GCA_000001215.4/GCA_000001215.4_Release_6_plus_ISO1_MT_genomic.fna"
REF_GENOME_FILE="$GENOMES_DIR/$REF_GENOME_REL_PATH"
INITIAL_GRAPH="drosophila_graph_step1.gfa"

# --- Script Start ---
echo "Starting DROSOPHILA pangenome graph construction..."

# Navigate into the directory with all the files
cd "$GENOMES_DIR" || { echo "Error: Could not find directory $GENOMES_DIR"; exit 1; }

# --- FIX: Automatically create/validate the initial graph ---
# Check if the reference genome exists first.
if [ ! -f "$REF_GENOME_FILE" ]; then
    echo "Error: Reference genome file '$REF_GENOME_FILE' not found!"
    exit 1
fi

# Check if the initial graph is missing or too small (e.g., less than 1KB).
if [ ! -f "$INITIAL_GRAPH" ] || [ $(stat -c%s "$INITIAL_GRAPH") -lt 1024 ]; then
    echo "----------------------------------------------------"
    echo "Initial graph '$INITIAL_GRAPH' is missing or invalid."
    echo "Recreating it now..."
    echo "----------------------------------------------------"
    rm -f "$INITIAL_GRAPH" # Remove old invalid file if it exists
    
    # Create the base graph using the -x ggs preset.
    "$MINIGRAPH_EXEC" -cxggs "$REF_GENOME_FILE" "$REF_GENOME_FILE" > "$INITIAL_GRAPH"

    # Check that the new graph was created successfully.
    if [ $? -ne 0 ] || [ ! -s "$INITIAL_GRAPH" ]; then
        echo "!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!"
        echo "Error: Failed to create the initial graph."
        echo "Please check your minigraph installation and reference file."
        echo "!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!"
        exit 1
    fi
    echo "Initial graph created successfully."
fi


# The 'current_graph.gfa' will be the input for each step.
cp "$INITIAL_GRAPH" current_graph.gfa

# Collect all genome files excluding the reference.
# Since files are in respective folders we use `find` instead of `*.fna`.
mapfile -t STRAIN_FILES < <(find . -type f -name "*.fna" | grep -v "$REF_GENOME_REL_PATH")

TOTAL_GENOMES=${#STRAIN_FILES[@]}
COUNT=0

for strain_file in "${STRAIN_FILES[@]}"; do
    COUNT=$((COUNT+1))
    echo "----------------------------------------------------"
    echo "Processing genome $COUNT of $TOTAL_GENOMES: $strain_file"
    echo "----------------------------------------------------"

    # Run minigraph
    "$MINIGRAPH_EXEC" -cxggs current_graph.gfa "$strain_file" > next_graph.gfa

    # Added critical error checking.
    # This checks if the last command failed OR if the output file is empty (0 bytes).
    if [ $? -ne 0 ] || [ ! -s "next_graph.gfa" ]; then
        echo "!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!"
        echo "Error: Minigraph failed or produced an empty file for '$strain_file'."
        echo "Stopping script."
        echo "!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!"
        rm -f next_graph.gfa # Clean up the empty file
        exit 1
    fi

    # The new graph becomes the input for the next iteration
    mv next_graph.gfa current_graph.gfa
done

# MODIFIED: Renamed the final graph
mv current_graph.gfa drosophila_pangenome_final.gfa

echo "===================================================="
echo "✅ Pangenome graph construction complete!"
echo "Final graph saved as: drosophila_pangenome_final.gfa"
echo "===================================================="
