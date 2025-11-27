#!/bin/bash

# Usage: ./run_parallel.sh <GFA_FILE> <READS_FILE> [THREADS]

GFA=$1
READS=$2
THREADS=${3:-4}

if [ -z "$GFA" ] || [ -z "$READS" ]; then
    echo "Usage: $0 <GFA_FILE> <READS_FILE> [THREADS]"
    exit 1
fi

echo "[Pipeline] Step 1: Partitioning Graph and Dispatching Reads..."
# Run PanAligner in Partition Mode
./PanAligner --partition -t $THREADS $GFA $READS > /dev/null 2>&1

echo "[Pipeline] Step 1 Complete. Checking for partitions..."

# Find all partition GFA files
PARTITIONS=$(ls partition_*.gfa 2>/dev/null)

if [ -z "$PARTITIONS" ]; then
    echo "[Error] No partition files found. Partitioning failed?"
    exit 1
fi

echo "[Pipeline] Step 2: Running Parallel Alignment..."

# Function to run alignment for a single partition
run_alignment() {
    export LD_LIBRARY_PATH=/home/kshitij/local/lib:$LD_LIBRARY_PATH
    part_gfa=$1
    # Extract partition ID from filename (partition_X.gfa)
    pid=$(echo $part_gfa | sed 's/partition_\([0-9]*\).gfa/\1/')
    reads_fq="reads_part_${pid}.fq"
    
    if [ -f "$reads_fq" ]; then
        echo "[Aligner] Processing Partition $pid with $reads_fq..."
        # Run in Normal Mode (no --partition)
        ./PanAligner -c -t 2 $part_gfa $reads_fq > "alignment_part_${pid}.gaf" 2> "log_part_${pid}.txt"
    else
        echo "[Aligner] No reads for Partition $pid. Skipping."
    fi
}

export -f run_alignment

# Use GNU parallel if available, otherwise simple loop (backgrounding)
if command -v parallel >/dev/null 2>&1; then
    echo "$PARTITIONS" | parallel -j $THREADS run_alignment {}
else
    echo "[Warning] GNU parallel not found. Running in simple background loop."
    for gfa in $PARTITIONS; do
        run_alignment $gfa &
    done
    wait
fi

echo "[Pipeline] Step 3: Merging Results..."
cat alignment_part_*.gaf > final_alignment.gaf
echo "[Pipeline] All alignments complete. Final output: final_alignment.gaf"
