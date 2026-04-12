#!/bin/bash
# =============================================================================
# benchmark.sh  -  GRASP data-collection script for paper figures
#
# Figures produced:
#   Fig 1 - Accuracy vs K (# partitions) + cut-edge fraction vs K
#   Fig 2 - Correct-mapping rate vs halo depth d, compared with PA baseline
#   Fig 3 - Resource usage (RAM, wall time): PA vs GRASP Align-Only vs GRASP Full
#   Fig 4 - Scalability: speedup / RAM / time vs K on Drosophila
#
# Usage:
#   cd /home/kshitij/FY_Project/GRASP
#   bash benchmark.sh [yeast|drosophila|all]  [fig1|fig2|fig3|fig4|all]
#
# Output: results/ directory with TSV files per figure.
# =============================================================================

set -euo pipefail
export LD_LIBRARY_PATH=/home/kshitij/local/lib:${LD_LIBRARY_PATH:-}

# ── Paths ──────────────────────────────────────────────────────────────────
GRASP_DIR="$(cd "$(dirname "$0")" && pwd)"
PA_ORIGINAL="$GRASP_DIR/PanAlignerOriginal/PanAligner"   # baseline (unmodified)
PA_GRASP="$GRASP_DIR/PanAligner/PanAligner"              # GRASP-modified

YEAST_GFA="$GRASP_DIR/PanAligner/yeast_pangenome_final.gfa"
DROSO_GFA="$GRASP_DIR/PanAligner/drosophila_pangenome_final.gfa"

# Real reads (for Fig 3 / Fig 4 resource benchmarks)
YEAST_READS="$GRASP_DIR/YeastGenomes/SRR36120465.fasta"
DROSO_READS="$GRASP_DIR/PanAligner/drosophila_test.fasta"

# Reference genomes (for simulated reads used in Fig 1 / Fig 2 accuracy)
YEAST_REF="$GRASP_DIR/YeastGenomes/GCF_000146045.2_R64_genomic.fna"
DROSO_REF="$GRASP_DIR/drosophila_dataset/GCA_000001215.4/GCA_000001215.4_Release_6_plus_ISO1_MT_genomic.fna"

RESULTS="$GRASP_DIR/results"
SIM_READS_YEAST="$RESULTS/sim_reads_yeast.fasta"
SIM_READS_DROSO="$RESULTS/sim_reads_droso.fasta"

N_THREADS=12    # match available CPU cores
SIM_N=3000      # number of simulated reads per dataset
SIM_LEN=10000   # simulated read length (long reads)

# K used for the Fig 3 resource-usage run.  Lower = fewer parallel workers =
# less peak RAM.  Override via env: FIG3_K=4 bash benchmark.sh ...
FIG3_K="${FIG3_K:-4}"

# K values swept in Fig 4 scalability experiment.
# Override via env: FIG4_K_VALUES="1 2 4 6" bash benchmark.sh ...
FIG4_K_VALUES="${FIG4_K_VALUES:-1 2 4 6}"

# ── Arguments ─────────────────────────────────────────────────────────────
DATASET="${1:-all}"   # yeast | drosophila | all
FIG="${2:-all}"       # fig1 | fig2 | fig3 | fig4 | all

# ── Helpers ────────────────────────────────────────────────────────────────
mkdir -p "$RESULTS"

log()  { echo "[bench] $*"; }
die()  { echo "[ERROR] $*" >&2; exit 1; }

# Run a command, capture peak RSS (KB->GB) and wall time (s).
# Writes <prefix>.time with two fields: wall_s  peak_rss_gb
run_timed() {
    local prefix="$1"; shift
    local tmp; tmp=$(mktemp)
    /usr/bin/time -v "$@" 2>"$tmp" || true
    local wall
    wall=$(grep "Elapsed (wall clock)" "$tmp" | awk '{print $NF}' | awk -F: '{
        if (NF==3) print $1*3600+$2*60+$3
        else        print $1*60+$2 }')
    local rss
    rss=$(grep "Maximum resident set" "$tmp" | awk '{print $NF}')
    local rss_gb
    rss_gb=$(awk "BEGIN{printf \"%.3f\", ${rss:-0}/1048576}")
    echo "${wall:-0} $rss_gb" > "${prefix}.time"
    # /usr/bin/time -v output lines are all tab-indented; the command's own
    # stderr lines are not – so strip the tab-prefixed timing lines to recover
    # the program's stderr and save it as <prefix>.log
    grep -v "^	" "$tmp" > "${prefix}.log" 2>/dev/null || true
    rm -f "$tmp"
}

# Extract cut-edge percentage from PanAligner stderr log
cut_frac_from_log() {
    local logfile="$1"
    grep "Cut edges:" "$logfile" 2>/dev/null \
        | tail -1 \
        | grep -oP '[\d.]+(?=%)' \
        | awk '{printf "%.4f", $1/100}' \
        || echo "NA"
}

# Count alignments in a GAF file
mapped_reads() { grep -c '.' "$1" 2>/dev/null || echo 0; }

# ── Simulate reads ──────────────────────────────────────────────────────────
simulate_reads() {
    log "Simulating reads for accuracy experiments..."
    if [ ! -f "$SIM_READS_YEAST" ]; then
        log "  Yeast: ${SIM_N} reads x ${SIM_LEN}bp"
        python3 "$GRASP_DIR/simulate_reads.py" \
            --fasta "$YEAST_REF" --output "$SIM_READS_YEAST" \
            --length "$SIM_LEN" --num "$SIM_N" --seed 42
    else
        log "  Yeast simulated reads already exist, skipping."
    fi
    if [ ! -f "$SIM_READS_DROSO" ]; then
        log "  Drosophila: ${SIM_N} reads x ${SIM_LEN}bp"
        python3 "$GRASP_DIR/simulate_reads.py" \
            --fasta "$DROSO_REF" --output "$SIM_READS_DROSO" \
            --length "$SIM_LEN" --num "$SIM_N" --seed 42
    else
        log "  Drosophila simulated reads already exist, skipping."
    fi
}

# ── Run one GRASP partition+align cycle ────────────────────────────────────
# run_grasp <workdir> <gfa> <reads> <K> <depth>
run_grasp() {
    local wdir="$1" gfa="$2" reads="$3" K="$4" D="$5"
    mkdir -p "$wdir"
    # Remove stale partition files so a fresh run is clean
    rm -f "$wdir"/partition_*.gfa "$wdir"/reads_part_*.fq \
          "$wdir"/alignment_part_*.gaf "$wdir"/aligned.gaf
    pushd "$wdir" > /dev/null

    log "    GRASP K=$K d=$D in $(basename $wdir)"

    # Step 1 - Partition the graph and dispatch reads
    # stdout discarded; run_timed saves stderr to partition.log automatically
    run_timed "partition" \
        "$PA_GRASP" --partition --nparts "$K" --halo-depth "$D" \
            -t "$N_THREADS" "$gfa" "$reads" \
        > /dev/null
    local partition_wall partition_rss
    read -r partition_wall partition_rss < partition.time 2>/dev/null || { partition_wall=0; partition_rss=0; }

    local nparts_found
    nparts_found=$(ls partition_*.gfa 2>/dev/null | wc -l)
    if [ "$nparts_found" -eq 0 ]; then
        log "    WARNING: no partition files found"
        popd > /dev/null
        return
    fi
    log "    $nparts_found partition files found"

    # Step 2 - Parallel alignment (one job per partition, wall-clock timed)
    local align_wall_start=$SECONDS
    local max_align_rss="0"
    for pgfa in partition_*.gfa; do
        pid="${pgfa#partition_}"; pid="${pid%.gfa}"
        rfq="reads_part_${pid}.fq"
        if [ -f "$rfq" ]; then
            if [ "${GRASP_SEQUENTIAL:-0}" = "1" ]; then
                run_timed "align_${pid}" \
                    "$PA_GRASP" -cx lr -t "$N_THREADS" "$pgfa" "$rfq" \
                    > "alignment_part_${pid}.gaf"
            else
                run_timed "align_${pid}" \
                    "$PA_GRASP" -cx lr -t 2 "$pgfa" "$rfq" \
                    > "alignment_part_${pid}.gaf" &
            fi
        fi
    done
    [ "${GRASP_SEQUENTIAL:-0}" = "1" ] || wait
    local align_wall=$(( SECONDS - align_wall_start ))

    # Collect max RAM across align workers
    for f in align_*.time; do
        [ -f "$f" ] || continue
        r=$(awk '{print $2}' "$f" 2>/dev/null || echo 0)
        if awk "BEGIN{exit !($r > $max_align_rss)}"; then max_align_rss=$r; fi
    done

    # Step 3 - Merge
    cat alignment_part_*.gaf > aligned.gaf 2>/dev/null || touch aligned.gaf

    # Write align-only timing summary
    echo "$align_wall $max_align_rss" > align_only.time

    popd > /dev/null
}

# ── Run alignment-only on pre-existing partition files ─────────────────────
# run_grasp_align_only <workdir>
# Keeps partition_*.gfa and reads_part_*.fq intact.
# Deletes old alignment outputs, re-runs workers, writes align_only.time.
run_grasp_align_only() {
    local wdir="$1"
    pushd "$wdir" > /dev/null

    local nparts_found
    nparts_found=$(ls partition_*.gfa 2>/dev/null | wc -l)
    if [ "$nparts_found" -eq 0 ]; then
        log "    WARNING: no partition files in $(basename $wdir) – skipping align-only"
        popd > /dev/null; return
    fi

    # Wipe only the alignment outputs, leave partition + read files untouched
    rm -f alignment_part_*.gaf log_part_*.txt aligned.gaf align_*.time align_only.time

    log "    GRASP Align-Only: $nparts_found partitions in $(basename $wdir)"

    local align_start=$SECONDS
    for pgfa in partition_*.gfa; do
        pid="${pgfa#partition_}"; pid="${pid%.gfa}"
        rfq="reads_part_${pid}.fq"
        if [ -f "$rfq" ]; then
            if [ "${GRASP_SEQUENTIAL:-0}" = "1" ]; then
                run_timed "align_${pid}" \
                    "$PA_GRASP" -cx lr -t "$N_THREADS" "$pgfa" "$rfq" \
                    > "alignment_part_${pid}.gaf"
            else
                run_timed "align_${pid}" \
                    "$PA_GRASP" -cx lr -t 2 "$pgfa" "$rfq" \
                    > "alignment_part_${pid}.gaf" &
            fi
        fi
    done
    [ "${GRASP_SEQUENTIAL:-0}" = "1" ] || wait
    local align_wall=$(( SECONDS - align_start ))

    local max_align_rss="0"
    for f in align_*.time; do
        [ -f "$f" ] || continue
        r=$(awk '{print $2}' "$f" 2>/dev/null || echo 0)
        if awk "BEGIN{exit !($r > $max_align_rss)}"; then max_align_rss=$r; fi
    done

    cat alignment_part_*.gaf > aligned.gaf 2>/dev/null || touch aligned.gaf
    echo "$align_wall $max_align_rss" > align_only.time
    log "    GRASP Align-Only done: wall=${align_wall}s ram=${max_align_rss}GB"

    popd > /dev/null
}

# ── PanAligner baseline ─────────────────────────────────────────────────────
run_baseline() {
    local wdir="$1" gfa="$2" reads="$3"
    mkdir -p "$wdir"
    if [ -f "$wdir/baseline.gaf" ]; then
        log "  PA baseline already exists in $(basename $wdir), skipping."
        return
    fi
    log "  PA baseline in $(basename $wdir)"
    # stdout = alignments; run_timed saves stderr to $wdir/baseline.log
    run_timed "$wdir/baseline" \
        "$PA_ORIGINAL" -cx lr -t "$N_THREADS" "$gfa" "$reads" \
        > "$wdir/baseline.gaf"
}

# ── Eval accuracy helper ────────────────────────────────────────────────────
eval_agreement() {
    local baseline_gaf="$1" query_gaf="$2"
    python3 "$GRASP_DIR/eval_accuracy.py" \
        --mode agreement \
        --baseline "$baseline_gaf" \
        --query    "$query_gaf" \
        2>/dev/null | grep '^AGREEMENT_RATE=' | cut -d= -f2 || echo "NA"
}

# ===========================================================================
# Fig 1 - Alignment accuracy vs K + cut-edge fraction vs K
# ===========================================================================
run_fig1() {
    local dataset="$1"
    local gfa reads
    if [ "$dataset" = "yeast" ]; then
        gfa="$YEAST_GFA"; reads="$SIM_READS_YEAST"
    else
        gfa="$DROSO_GFA"; reads="$SIM_READS_DROSO"
    fi

    local outtsv="$RESULTS/fig1_${dataset}.tsv"
    printf "K\tcut_edge_frac\tagreement_rate\tmapping_rate_grasp\tmapping_rate_baseline\n" \
        > "$outtsv"

    local bdir="$RESULTS/baseline_${dataset}_fig1"
    run_baseline "$bdir" "$gfa" "$reads"
    local b_mapped; b_mapped=$(mapped_reads "$bdir/baseline.gaf")

    for K in 2 4 6 8 10 12; do
        local wdir="$RESULTS/fig1_${dataset}_K${K}"
        run_grasp "$wdir" "$gfa" "$reads" "$K" 1

        local cef; cef=$(cut_frac_from_log "$wdir/partition.log")
        local q_mapped=0
        [ -f "$wdir/aligned.gaf" ] && q_mapped=$(mapped_reads "$wdir/aligned.gaf")
        local agree="NA"
        if [ -f "$wdir/aligned.gaf" ] && [ -f "$bdir/baseline.gaf" ]; then
            agree=$(eval_agreement "$bdir/baseline.gaf" "$wdir/aligned.gaf")
        fi
        printf "%s\t%s\t%s\t%s\t%s\n" "$K" "$cef" "$agree" "$q_mapped" "$b_mapped" >> "$outtsv"
        log "  Fig1 $dataset K=$K cut=$cef agree=$agree"
    done
    log "Fig 1 ($dataset) written to $outtsv"
}

# ===========================================================================
# Fig 2 - Correct-mapping rate vs halo depth d
# ===========================================================================
run_fig2() {
    local dataset="$1"
    local gfa reads
    if [ "$dataset" = "yeast" ]; then
        gfa="$YEAST_GFA"; reads="$SIM_READS_YEAST"
    else
        gfa="$DROSO_GFA"; reads="$SIM_READS_DROSO"
    fi

    local outtsv="$RESULTS/fig2_${dataset}.tsv"
    printf "depth\tagreement_rate\tmapping_rate\n" > "$outtsv"

    local bdir="$RESULTS/baseline_${dataset}_fig2"
    run_baseline "$bdir" "$gfa" "$reads"
    local b_mapped; b_mapped=$(mapped_reads "$bdir/baseline.gaf")
    # Record baseline as depth=baseline row
    printf "baseline\t1.0000\t%s\n" "$b_mapped" >> "$outtsv"

    for D in 0 1 2 3 5 8 10 15 20; do
        local wdir="$RESULTS/fig2_${dataset}_d${D}"
        run_grasp "$wdir" "$gfa" "$reads" 8 "$D"   # fixed K=8

        local q_mapped=0
        [ -f "$wdir/aligned.gaf" ] && q_mapped=$(mapped_reads "$wdir/aligned.gaf")
        local agree="NA"
        if [ -f "$wdir/aligned.gaf" ] && [ -f "$bdir/baseline.gaf" ]; then
            agree=$(eval_agreement "$bdir/baseline.gaf" "$wdir/aligned.gaf")
        fi
        printf "%s\t%s\t%s\n" "$D" "$agree" "$q_mapped" >> "$outtsv"
        log "  Fig2 $dataset d=$D agree=$agree mapped=$q_mapped"
    done
    log "Fig 2 ($dataset) written to $outtsv"
}

# ===========================================================================
# Fig 3 - Resource usage comparison
# ===========================================================================
run_fig3() {
    local dataset="$1"
    local gfa reads
    if [ "$dataset" = "yeast" ]; then
        gfa="$YEAST_GFA"; reads="$YEAST_READS"
    else
        gfa="$DROSO_GFA"; reads="$DROSO_READS"
    fi

    local outtsv="$RESULTS/fig3_${dataset}.tsv"
    printf "config\twall_s\tpeak_ram_gb\n" > "$outtsv"

    # 1. PanAligner baseline
    local bdir="$RESULTS/fig3_${dataset}_baseline"
    mkdir -p "$bdir"
    log "  Fig3 $dataset: running PanAligner baseline"
    run_timed "$bdir/run" \
        "$PA_ORIGINAL" -cx lr -t "$N_THREADS" "$gfa" "$reads" \
        > "$bdir/out.gaf"
    local b_wall b_rss
    read -r b_wall b_rss < "$bdir/run.time"
    printf "PanAligner\t%s\t%s\n" "$b_wall" "$b_rss" >> "$outtsv"
    log "  Fig3 $dataset PanAligner: wall=${b_wall}s ram=${b_rss}GB"

    # 2. GRASP Full Run (partition + dispatch + align)
    log "  Fig3 $dataset: running GRASP Full Run (K=${FIG3_K} d=1)"
    local fdir="$RESULTS/fig3_${dataset}_grasp_full"
    local full_start=$SECONDS
    run_grasp "$fdir" "$gfa" "$reads" "$FIG3_K" 1
    local full_wall=$(( SECONDS - full_start ))

    # peak RAM = max of partition step and per-worker align step
    local part_rss="0"
    [ -f "$fdir/partition.time" ] && read -r _ part_rss < "$fdir/partition.time"
    local align_rss="0"
    for f in "$fdir"/align_*.time; do
        [ -f "$f" ] || continue
        r=$(awk '{print $2}' "$f")
        if awk "BEGIN{exit !($r > $align_rss)}"; then align_rss=$r; fi
    done
    # True peak = max(partition step RAM, per-worker align RAM)
    local full_rss
    full_rss=$(awk "BEGIN{print ($part_rss > $align_rss) ? $part_rss : $align_rss}")
    printf "GRASP_Full\t%s\t%s\n" "$full_wall" "$full_rss" >> "$outtsv"
    log "  Fig3 $dataset GRASP Full: wall=${full_wall}s ram=${full_rss}GB"

    # 3. GRASP Align-Only – separate timed run on the partition files from step 2
    log "  Fig3 $dataset: running GRASP Align-Only (reusing partition files)"
    run_grasp_align_only "$fdir"
    local ao_wall="0" ao_rss="0"
    [ -f "$fdir/align_only.time" ] && read -r ao_wall ao_rss < "$fdir/align_only.time"
    printf "GRASP_AlignOnly\t%s\t%s\n" "$ao_wall" "$ao_rss" >> "$outtsv"
    log "  Fig3 $dataset GRASP Align-Only: wall=${ao_wall}s ram=${ao_rss}GB"

    log "Fig 3 ($dataset) written to $outtsv"
}

# ===========================================================================
# Fig 4 - Scalability vs K (run on Drosophila by default)
# ===========================================================================
run_fig4() {
    local dataset="${1:-drosophila}"
    local gfa reads
    if [ "$dataset" = "yeast" ]; then
        gfa="$YEAST_GFA"; reads="$YEAST_READS"
    else
        gfa="$DROSO_GFA"; reads="$DROSO_READS"
    fi

    local outtsv="$RESULTS/fig4_${dataset}.tsv"
    printf "K\tspeedup\tpeak_ram_gb\talign_only_s\ttotal_s\n" > "$outtsv"

    # Baseline timing
    local bdir="$RESULTS/fig4_${dataset}_baseline"
    mkdir -p "$bdir"
    log "  Fig4 $dataset: PA baseline"
    run_timed "$bdir/run" \
        "$PA_ORIGINAL" -cx lr -t "$N_THREADS" "$gfa" "$reads" \
        > "$bdir/out.gaf"
    local b_wall b_rss
    read -r b_wall b_rss < "$bdir/run.time"
    log "  Fig4 $dataset baseline: wall=${b_wall}s"

    for K in $FIG4_K_VALUES; do
        local wdir="$RESULTS/fig4_${dataset}_K${K}"
        local start_total=$SECONDS
        run_grasp "$wdir" "$gfa" "$reads" "$K" 1
        local total_wall=$(( SECONDS - start_total ))

        local ao_wall="0" ao_rss="0"
        [ -f "$wdir/align_only.time" ] && read -r ao_wall ao_rss < "$wdir/align_only.time"

        local speedup
        speedup=$(awk "BEGIN{if ($total_wall>0) printf \"%.3f\", $b_wall/$total_wall; else print \"NA\"}")
        printf "%s\t%s\t%s\t%s\t%s\n" "$K" "$speedup" "$ao_rss" "$ao_wall" "$total_wall" >> "$outtsv"
        log "  Fig4 $dataset K=$K speedup=$speedup align=${ao_wall}s total=${total_wall}s"
    done
    log "Fig 4 ($dataset) written to $outtsv"
}

# ===========================================================================
# Main dispatcher
# ===========================================================================
simulate_reads

run_for_dataset() {
    local ds="$1"
    case "$FIG" in
        fig1) run_fig1 "$ds" ;;
        fig2) run_fig2 "$ds" ;;
        fig3) run_fig3 "$ds" ;;
        fig4) run_fig4 "$ds" ;;
        all)
            run_fig1 "$ds"
            run_fig2 "$ds"
            run_fig3 "$ds"
            ;;
    esac
}

case "$DATASET" in
    yeast)
        run_for_dataset yeast
        if [ "$FIG" = "fig4" ] || [ "$FIG" = "all" ]; then run_fig4 drosophila; fi
        ;;
    drosophila)
        run_for_dataset drosophila
        ;;
    all)
        run_for_dataset yeast
        run_for_dataset drosophila
        if [ "$FIG" = "fig4" ] || [ "$FIG" = "all" ]; then run_fig4 drosophila; fi
        ;;
    *)
        die "Unknown dataset '$DATASET'. Use: yeast | drosophila | all"
        ;;
esac

log "================================================================"
log "All done. Results in: $RESULTS/"
log "TSV files:"
ls "$RESULTS"/*.tsv 2>/dev/null | sed 's/^/  /'
log "================================================================"
