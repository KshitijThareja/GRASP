#!/bin/bash
# =============================================================================
# run_sequential.sh  -  Run all benchmark cases one at a time to avoid OOM.
#
# Each case (dataset x figure x K/d value) runs fully before the next starts.
# Alignment workers within each GRASP run are also serialized (no parallel
# subprocesses), so only one PanAligner process holds memory at a time.
#
# Usage:
#   bash run_sequential.sh              # all datasets, all figures
#   bash run_sequential.sh yeast        # yeast only, all figures
#   bash run_sequential.sh yeast fig1   # single figure
# =============================================================================

set -euo pipefail

SCRIPT_DIR="$(cd "$(dirname "$0")" && pwd)"
BENCH="$SCRIPT_DIR/benchmark.sh"

DATASET="${1:-all}"
FIG="${2:-all}"

# Fig 1 and Fig 2 measure accuracy only (no wall-time in TSV).
# Sequential workers are safe there and keep peak RAM low.
# Fig 3 and Fig 4 measure wall time and peak RAM – must stay parallel.
run_case() {
    local ds="$1" fig="$2"
    echo ""
    echo "========================================================"
    echo "  Starting: $ds $fig"
    echo "========================================================"
    case "$fig" in
        fig1|fig2) GRASP_SEQUENTIAL=1 bash "$BENCH" "$ds" "$fig" ;;
        fig3|fig4) GRASP_SEQUENTIAL=0 bash "$BENCH" "$ds" "$fig" ;;
        *)         GRASP_SEQUENTIAL=0 bash "$BENCH" "$ds" "$fig" ;;
    esac
    echo "  Done: $ds $fig"
}

case "$DATASET" in
    yeast)
        case "$FIG" in
            all)
                run_case yeast fig1
                run_case yeast fig2
                run_case yeast fig3
                ;;
            *)  run_case yeast "$FIG" ;;
        esac
        ;;
    drosophila)
        case "$FIG" in
            all)
                run_case drosophila fig1
                run_case drosophila fig2
                run_case drosophila fig3
                run_case drosophila fig4
                ;;
            *)  run_case drosophila "$FIG" ;;
        esac
        ;;
    all)
        case "$FIG" in
            all)
                run_case yeast      fig1
                run_case yeast      fig2
                run_case yeast      fig3
                run_case drosophila fig1
                run_case drosophila fig2
                run_case drosophila fig3
                run_case drosophila fig4
                ;;
            *)
                run_case yeast      "$FIG"
                run_case drosophila "$FIG"
                if [ "$FIG" = "fig4" ]; then : ; else
                    run_case drosophila fig4
                fi
                ;;
        esac
        ;;
    *)
        echo "Usage: $0 [yeast|drosophila|all] [fig1|fig2|fig3|fig4|all]" >&2
        exit 1
        ;;
esac

echo ""
echo "========================================================"
echo "  All cases complete."
echo "========================================================"
