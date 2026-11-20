#!/bin/bash
# Draws the main trunk and the tree placement, m = 10, for alt03, JAB04 and
# ATH06 into out/placement_figs/.
#
#   bash scripts/viz/make_placement_figs.sh

set -euo pipefail

SCRIPT_DIR="$(cd "$(dirname "$0")" && pwd)"
REPO="$(cd "$SCRIPT_DIR/../.." && pwd)"
OUT_DIR="$REPO/out/placement_figs"
INSTANCES="$REPO/data/instances"

for arg in "$@"; do
    case "$arg" in
        -h|--help) sed -n '2,/^set -euo/p' "$0" | sed '$d' | sed 's/^# \{0,1\}//'; exit 0 ;;
        *) echo "unknown argument: $arg" >&2; exit 1 ;;
    esac
done

# instance:m:extra flags for the linear picture (alt03 uses the tree's colours)
FIGS=(
    "rede_mt_feeder_alt03:10:--general-colors"
    "JaboticabalJAB04:10:"
    "JaboticabalATH06:10:"
)

mkdir -p "$OUT_DIR"
for entry in "${FIGS[@]}"; do
    IFS=: read -r inst m lin_flags <<< "$entry"
    fi="$INSTANCES/$inst.fi"
    echo "== $inst m=$m"
    lin="$OUT_DIR/${inst}_m${m}_linear.png"
    gen="$OUT_DIR/${inst}_m${m}_general.png"
    # shellcheck disable=SC2086
    bash "$SCRIPT_DIR/linear/run.sh" "$fi" "$m" --no-title $lin_flags --out "$lin"
    bash "$SCRIPT_DIR/general/run.sh" "$fi" "$m" --no-title --out "$gen"
done

echo "figures in $OUT_DIR"
