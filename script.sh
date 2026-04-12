#!/bin/bash

# This script automates the construction of a YEAST pangenome graph using Minigraph.
# It will automatically create the initial graph if it's missing or invalid.

# --- Configuration ---
# MODIFIED: Updated paths for the Yeast Genome project
MINIGRAPH_EXEC="/home/kshitij/FY_Project/GRASP/minigraph/minigraph"
GENOMES_DIR="/home/kshitij/FY_Project/GRASP/YeastGenomes"

# MODIFIED: Set the R64 genome as the reference file
REF_GENOME_FILE="GCF_000146045.2_R64_genomic.fna"
INITIAL_GRAPH="graph_step1.gfa"

# --- Script Start ---
echo "Starting YEAST pangenome graph construction..."

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

# Use a safer way to find files to get an accurate count for the progress message.
# Note: The glob '*.fna' will match all the '...ic.fna' files in your screenshot.
TOTAL_GENOMES=$(find . -maxdepth 1 -name "*.fna" | grep -v "$REF_GENOME_FILE" | wc -l)
COUNT=0

# Use a safer for-loop to iterate through genome files instead of parsing 'ls'.
for strain_file in *.fna; do
    # Skip the reference genome file itself
    if [[ "$strain_file" == "$REF_GENOME_FILE" ]]; then
        continue
    fi
    
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
mv current_graph.gfa yeast_pangenome_final.gfa

echo "===================================================="
echo "✅ Pangenome graph construction complete!"
echo "Final graph saved as: yeast_pangenome_final.gfa"
echo "===================================================="
