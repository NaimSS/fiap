#!/usr/bin/env bash
# One instance with dp_nlogW on the main trunk, m indicators (so K = m + 1
# sections). Prints one CSV row in the comparison_run_all.sh schema, with
# the fastest of [reps] (3) runs as time.
#
#   bash scripts/bench/comparison_run_linear.sh <instance.fi> <m> [reps]
#
# slots is the number of non-zero entries of the linear array. The trunk
# can't take more than slots - 1 indicators, so larger m get status SATURATED.
set -uo pipefail
REPO="$(cd "$(dirname "$0")/../.." && pwd)"
if [[ "${CAPPED:-0}" == 1 ]]; then BIN="$REPO/bin"; export MEM_CAP_MB="${MEM_CAP_MB:-1024}"
else BIN="$REPO/bin/nocap"; fi
[[ -x "$BIN/dp_nlogW" ]] || { echo "error: $BIN/dp_nlogW missing; run scripts/build.sh" >&2; exit 1; }
FI="$1"; M="$2"; REPS="${3:-3}"
INST="$(basename "$FI" .fi)"
DIAM="$FI.diameter"
[[ -f "$DIAM" ]] || "$BIN/fi_to_diameter" "$FI" >/dev/null

read -r N_NODES N_EDGES < <(awk '$1=="p"{print $3, $4; exit}' "$FI")
DIAM_ARCS=$(awk '$1=="diameter"{print NF-1; exit}' "$DIAM")
SLOTS=$("$BIN/convert_to_line" < "$DIAM" | sed -n 2p | tr ' ' '\n' | grep -vc '^0$')

K=$((M + 1))
BEST=""; VAL=""
for ((r = 0; r < REPS; r++)); do
    OUT=$("$BIN/dp_nlogW" "$DIAM" "$K" 2>/dev/null) || { echo "linear,$INST,$N_NODES,$N_EDGES,$DIAM_ARCS,$SLOTS,$M,,,ERROR"; exit 0; }
    V=$(echo "$OUT" | awk -F': ' '/^k = /{print $2}')
    T=$(echo "$OUT" | awk '/^time_ms:/{print $2}')
    VAL="$V"
    if [[ -z "$BEST" ]] || awk "BEGIN{exit !($T < $BEST)}"; then BEST="$T"; fi
done
STATUS=OK
(( K > SLOTS )) && STATUS=SATURATED
echo "linear,$INST,$N_NODES,$N_EDGES,$DIAM_ARCS,$SLOTS,$M,$VAL,$BEST,$STATUS"
