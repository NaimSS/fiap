#!/usr/bin/env bash
# Runs lagrange_solver ("lagrange") and default_solver ("default") on Campinas
# over a K sweep. data/results/bench_campinas.csv was made with the defaults.
#
#   bash scripts/bench/bench_campinas.sh [OUT.csv]
#
# Output (default out/bench_campinas.csv):
#   solver,instance,K,time_s,peak_heap_mb,search_peak_mb,value,tag,status,timeout_s
# search_peak_mb is the heap peak of the lambda search alone, before cut
# recovery. tag is what the solver printed as status; status is OK only when
# tag is exact*.
#
# Env: MEM_CAP_MB (3072), TIMEOUT_S (600), KS (1,2,5,...,800),
# SOLVERS (lagrange,default), INSTANCE (the Campinas .fi).
set -uo pipefail

ROOT="$(cd "$(dirname "$0")/../.." && pwd)"
WORK="$ROOT/out/bench_campinas"
OUT_CSV="${1:-$ROOT/out/bench_campinas.csv}"
MEM_CAP_MB="${MEM_CAP_MB:-3072}"
TIMEOUT_S="${TIMEOUT_S:-600}"
IFS=',' read -ra KS <<< "${KS:-1,2,5,10,25,50,100,200,400,800}"
IFS=',' read -ra SOLVERS <<< "${SOLVERS:-lagrange,default}"
INSTANCE="${INSTANCE:-$ROOT/data/instances/rede_mt_municipality_campinas.fi}"
NAME="$(basename "$INSTANCE" .fi)"; NAME="${NAME#rede_mt_municipality_}"

log(){ printf '%s\n' "$*" >&2; }
[[ -f "$INSTANCE" ]] || { log "ERROR: no such instance: $INSTANCE"; exit 1; }
mkdir -p "$WORK" "$(dirname "$OUT_CSV")"

bash "$ROOT/scripts/build.sh" >&2 || { log "ERROR: build failed"; exit 1; }
bin_of(){ case "$1" in lagrange) echo "$ROOT/bin/lagrange_solver" ;; default) echo "$ROOT/bin/default_solver" ;; *) return 1 ;; esac; }
for b in "${SOLVERS[@]}"; do bin_of "$b" >/dev/null || { log "ERROR: unknown solver '$b'"; exit 1; }; done

GEN="$WORK/$NAME.gen"
if [[ ! -s "$GEN" ]]; then
    log "[prep] converter $NAME"
    "$ROOT/bin/converter" < "$INSTANCE" > "$GEN" || { log "ERROR: converter failed"; exit 1; }
fi

TMP_OUT="$WORK/.run.out"; TMP_ERR="$WORK/.run.err"
trap 'rm -f "$TMP_OUT" "$TMP_ERR"' EXIT

run_one(){ # stdin_file timeout binary args...
    local stdin_file="$1" tmo="$2"; shift 2
    R_STATUS=""; R_VALUE=""; R_TIME=""; R_MEM=""; R_SMEM=""; R_TAG=""
    MEM_CAP_MB="$MEM_CAP_MB" MEM_CAP_REPORT=1 perl -e 'alarm shift @ARGV; exec @ARGV' \
        "$tmo" "$@" < "$stdin_file" > "$TMP_OUT" 2> "$TMP_ERR"
    local rc=$?
    if   (( rc == 0 ));   then R_STATUS=OK
    elif (( rc == 142 )); then R_STATUS=TIMEOUT
    elif grep -q "memory cap exceeded" "$TMP_ERR"; then R_STATUS=OOM
    elif (( rc == 134 )); then R_STATUS=ASSERT
    else R_STATUS="ERROR($rc)"
    fi
    local kb; kb=$(grep -o 'mem_cap_peak_kb: [0-9]*' "$TMP_ERR" | tail -1 | awk '{print $2}')
    [[ -n "$kb" ]] && R_MEM=$(awk -v k="$kb" 'BEGIN{printf "%.1f", k/1024}')
    if [[ "$R_STATUS" == OK ]]; then
        R_VALUE=$(grep -- '- k=' "$TMP_OUT" | head -1 | awk '{print $1}')
        R_TIME=$(grep -o 'time: [0-9.eE+-]*' "$TMP_OUT" | head -1 | awk '{print $2}')
        R_TAG=$(grep -o 'status: [a-zA-Z_-]*' "$TMP_OUT" | head -1 | awk '{print $2}')
        [[ "$R_TAG" == exact* ]] || R_STATUS="${R_TAG:-ERROR(parse)}"
        local skb; skb=$(grep -o 'search_peak_kb: [0-9]*' "$TMP_OUT" | head -1 | awk '{print $2}')
        [[ -n "$skb" ]] && R_SMEM=$(awk -v k="$skb" 'BEGIN{printf "%.1f", k/1024}')
    fi
}

echo "solver,instance,K,time_s,peak_heap_mb,search_peak_mb,value,tag,status,timeout_s" > "$OUT_CSV"
for b in "${SOLVERS[@]}"; do
    for k in "${KS[@]}"; do
        log "[run] $b $NAME K=$k"
        run_one "$GEN" "$TIMEOUT_S" "$(bin_of "$b")" "$k"
        log "    -> $R_STATUS $R_VALUE t=${R_TIME}s mem=${R_MEM}MB (search ${R_SMEM:--}MB) tag=$R_TAG"
        echo "$b,$NAME,$k,$R_TIME,$R_MEM,$R_SMEM,$R_VALUE,$R_TAG,$R_STATUS,$TIMEOUT_S" >> "$OUT_CSV"
    done
done
log "results written to $OUT_CSV"
