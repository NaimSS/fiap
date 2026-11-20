#!/usr/bin/env bash
# One instance with lagrange_solver, m indicators. Prints one CSV row in the
# comparison_run_all.sh schema, with diam_arcs and slots empty and the
# fastest of [reps] (3) runs as time.
#
#   bash scripts/bench/comparison_run_general.sh <instance.fi> <m> [reps]
#
# status is the solver's own: exact, exact-collinear, exact-squeeze, or GAP
# (the regex below stops at the underscore of GAP_OR_COLLINEAR).
set -uo pipefail
REPO="$(cd "$(dirname "$0")/../.." && pwd)"
if [[ "${CAPPED:-0}" == 1 ]]; then BIN="$REPO/bin"; export MEM_CAP_MB="${MEM_CAP_MB:-1024}"
else BIN="$REPO/bin/nocap"; fi
[[ -x "$BIN/lagrange_solver" ]] || { echo "error: $BIN/lagrange_solver missing; run scripts/build.sh" >&2; exit 1; }
FI="$1"; M="$2"; REPS="${3:-3}"
INST="$(basename "$FI" .fi)"
read -r N_NODES N_EDGES < <(awk '$1=="p"{print $3, $4; exit}' "$FI")

TMP="$(mktemp "${TMPDIR:-/tmp}/conv_XXXXXX")"
trap 'rm -f "$TMP"' EXIT
# an empty tree would give value 0 with status exact
if ! ( "$BIN/converter" < "$FI" > "$TMP"; exit $? ) 2>/dev/null || [[ ! -s "$TMP" ]]; then
    echo "general,$INST,$N_NODES,$N_EDGES,,,$M,,,ERROR"; exit 0
fi

BEST=""; VAL=""; STATUS=""; CUTS=""
for ((r = 0; r < REPS; r++)); do
    OUT=$("$BIN/lagrange_solver" "$M" 1 1 < "$TMP" 2>/dev/null) || { echo "general,$INST,$N_NODES,$N_EDGES,,,$M,,,ERROR"; exit 0; }
    LINE=$(echo "$OUT" | head -1)
    V=${LINE%% *}
    T=$(echo "$LINE" | sed -E 's/.* time: ([0-9.e-]+)s.*/\1/')
    STATUS=$(echo "$LINE" | sed -E 's/.*status: ([A-Za-z-]+).*/\1/')
    CUTS=$(echo "$LINE" | sed -E 's/.*cuts: ([0-9]+).*/\1/')
    T_MS=$(awk "BEGIN{printf \"%.4f\", $T * 1000}")
    VAL="$V"
    if [[ -z "$BEST" ]] || awk "BEGIN{exit !($T_MS < $BEST)}"; then BEST="$T_MS"; fi
done
echo "general,$INST,$N_NODES,$N_EDGES,,,$M,$VAL,$BEST,$STATUS"
