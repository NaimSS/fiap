#!/usr/bin/env bash
# Main trunk DP against tree DP on every instance, for m indicators.
# Produces data/results/comparison.csv.
#
#   bash scripts/bench/comparison_run_all.sh [--ms=1,2,3,...] [--reps=3]
#       [--instances='glob'] [--only=linear|general] [--force] [out.csv]
#
# Output (default data/results/comparison.csv, kept unless --force):
#   approach,instance,n_nodes,n_edges,diam_arcs,slots,m,value,time_ms,status
#
# Uses the bin/nocap/ builds; CAPPED=1 switches to bin/ with MEM_CAP_MB (1024).
set -uo pipefail
HERE="$(cd "$(dirname "$0")" && pwd)"
REPO="$(cd "$HERE/../.." && pwd)"
M_LIST="1,2,3,5,8,10,15,20,25,30,40,50,75,100"
REPS=3
GLOB="$REPO/data/instances/*.fi"
ONLY=""
FORCE=0
OUT="$REPO/data/results/comparison.csv"
for arg in "$@"; do
    case "$arg" in
        --ms=*) M_LIST="${arg#*=}" ;;
        --reps=*) REPS="${arg#*=}" ;;
        --instances=*) GLOB="${arg#*=}" ;;
        --only=*) ONLY="${arg#*=}" ;;
        --force) FORCE=1 ;;
        -*) echo "unknown option: $arg" >&2; exit 1 ;;
        *) OUT="$arg" ;;
    esac
done
if [[ -e "$OUT" && "$FORCE" != 1 ]]; then
    echo "error: $OUT exists; pass another output path or --force" >&2
    exit 1
fi
if [[ "${CAPPED:-0}" == 1 ]]; then
    bash "$REPO/scripts/build.sh" >/dev/null || exit 1
else
    bash "$REPO/scripts/build.sh" --no-cap >/dev/null || exit 1
fi
mkdir -p "$(dirname "$OUT")"
echo "approach,instance,n_nodes,n_edges,diam_arcs,slots,m,value,time_ms,status" > "$OUT"
for FI in $GLOB; do
    for M in ${M_LIST//,/ }; do
        [[ "$ONLY" == "general" ]] || bash "$HERE/comparison_run_linear.sh"  "$FI" "$M" "$REPS" | tee -a "$OUT"
        [[ "$ONLY" == "linear"  ]] || bash "$HERE/comparison_run_general.sh" "$FI" "$M" "$REPS" | tee -a "$OUT"
    done
done
echo "wrote $OUT"
