#!/usr/bin/env bash
# Runs lagrange_solver and the six linear DPs on Campinas and Paulinia.
#
#   bash scripts/bench/bench_campinas_paulinia.sh [OUT.csv]
#
# Output (default out/campinas_paulinia_bench.csv):
#   solver,instance,K,time_s,peak_heap_mb,value,status
# For lagrange_solver K is the number of indicators; for the linear DPs it is
# the number of sections, taken as 10..100% of the trunk length.
#
# Env: MEM_CAP_MB (3072), TIMEOUT_S (120), GEN_KS (1,2,5,...,800),
# INSTANCES (paulinia,campinas).
set -uo pipefail

ROOT="$(cd "$(dirname "$0")/../.." && pwd)"
INST_DIR="$ROOT/data/instances"
WORK="$ROOT/out/bench_cp"
BIN="$ROOT/bin"
OUT_CSV="${1:-$ROOT/out/campinas_paulinia_bench.csv}"

MEM_CAP_MB="${MEM_CAP_MB:-3072}"
TIMEOUT_S="${TIMEOUT_S:-120}"
PCTS=(10 25 50 75 100)
IFS=',' read -ra GEN_KS <<< "${GEN_KS:-1,2,5,10,25,50,100,200,400,800}"
IFS=',' read -ra INSTANCES <<< "${INSTANCES:-paulinia,campinas}"

GEN_BINS=(lagrange_solver)
LIN_BINS=(dp_nlogW dp_nk dp_n2k dp_nk_smawk dp_nk_larsch dp_nlogW_larsch)
FI_paulinia="$INST_DIR/rede_mt_municipality_paulinia.fi"
FI_campinas="$INST_DIR/rede_mt_municipality_campinas.fi"

mkdir -p "$WORK" "$(dirname "$OUT_CSV")"
log(){ printf '%s\n' "$*" >&2; }
fi_of(){ eval "printf '%s' \"\$FI_$1\""; }

bash "$ROOT/scripts/build.sh" >&2 || { log "ERROR: build failed"; exit 1; }

for name in "${INSTANCES[@]}"; do
    fi="$(fi_of "$name")"
    [[ -f "$fi" ]] || { log "ERROR: unknown instance '$name'"; exit 1; }
    dia="$fi.diameter"
    gen="$WORK/$name.gen"
    if [[ ! -f "$dia" ]]; then
        log "[prep] fi_to_diameter $name ..."
        "$BIN/fi_to_diameter" "$fi" || { log "ERROR: fi_to_diameter failed ($name)"; exit 1; }
    fi
    if [[ ! -s "$gen" ]]; then
        log "[prep] converter $name -> $name.gen ..."
        "$BIN/converter" < "$fi" > "$gen" || { log "ERROR: converter failed ($name)"; exit 1; }
    fi
done

TMP_OUT="$WORK/.run.out"; TMP_ERR="$WORK/.run.err"
trap 'rm -f "$TMP_OUT" "$TMP_ERR"' EXIT

# run_one <stdin file or ""> <binary> [args...]
# sets R_STATUS R_VALUE R_TIME R_MEM. perl's alarm survives exec, so a
# timeout kills the solver itself with exit code 142.
run_one(){
    local stdin_file="$1"; shift
    R_STATUS=""; R_VALUE=""; R_TIME=""; R_MEM=""
    local rc
    if [[ -n "$stdin_file" ]]; then
        MEM_CAP_MB="$MEM_CAP_MB" MEM_CAP_REPORT=1 perl -e 'alarm shift @ARGV; exec @ARGV' \
            "$TIMEOUT_S" "$@" < "$stdin_file" > "$TMP_OUT" 2> "$TMP_ERR"
    else
        MEM_CAP_MB="$MEM_CAP_MB" MEM_CAP_REPORT=1 perl -e 'alarm shift @ARGV; exec @ARGV' \
            "$TIMEOUT_S" "$@" > "$TMP_OUT" 2> "$TMP_ERR"
    fi
    rc=$?

    if   (( rc == 0 ));   then R_STATUS=OK
    elif (( rc == 142 )); then R_STATUS=TIMEOUT
    elif grep -q "memory cap exceeded" "$TMP_ERR"; then R_STATUS=OOM
    elif (( rc == 134 )); then R_STATUS=ASSERT
    else R_STATUS="ERROR($rc)"
    fi

    local kb; kb=$(grep -o 'mem_cap_peak_kb: [0-9]*' "$TMP_ERR" | tail -1 | awk '{print $2}')
    [[ -n "$kb" ]] && R_MEM=$(awk -v k="$kb" 'BEGIN{printf "%.1f", k/1024}')

    if [[ "$R_STATUS" == OK ]]; then
        if grep -q -- '- k=' "$TMP_OUT"; then                 # tree DP
            local line; line=$(grep -m1 -- '- k=' "$TMP_OUT")
            R_VALUE=$(awk '{print $1}' <<< "$line")
            R_TIME=$(sed -nE 's/.*, time: ([0-9.eE+-]+)s.*/\1/p' <<< "$line")
            # anything but exact* (i.e. GAP_OR_COLLINEAR) is only an upper bound
            local st; st=$(sed -nE 's/.*, status: ([A-Za-z_-]+).*/\1/p' <<< "$line")
            [[ -n "$st" && "$st" != exact* ]] && R_STATUS="$st"
        elif grep -q '^k = ' "$TMP_OUT"; then                 # linear DP
            R_VALUE=$(grep '^k = ' "$TMP_OUT" | head -1 | sed -E 's/.*: //')
            local ms; ms=$(grep -o 'time_ms: [0-9.]*' "$TMP_OUT" | head -1 | awk '{print $2}')
            R_TIME=$(awk -v m="$ms" 'BEGIN{ if(m=="") print ""; else printf "%.6f", m/1000 }')
        fi
        [[ -z "$R_VALUE" || -z "$R_TIME" ]] && R_STATUS="ERROR(parse)"
    fi
}

maxk_of(){ awk '/^diameter/{print NF-1; exit}' "$1"; }   # number of trunk arcs

echo "solver,instance,K,time_s,peak_heap_mb,value,status" > "$OUT_CSV"

for name in "${INSTANCES[@]}"; do
    fi="$(fi_of "$name")"
    gen="$WORK/$name.gen"
    dia="$fi.diameter"

    log ""
    log "================ $name ================"

    for b in "${GEN_BINS[@]}"; do
        for k in "${GEN_KS[@]}"; do
            log "[run] $b $name K=$k ..."
            run_one "$gen" "$BIN/$b" "$k"
            log "      -> $R_STATUS  value=$R_VALUE  time=${R_TIME}s  mem=${R_MEM}MB"
            echo "$b,$name,$k,$R_TIME,$R_MEM,$R_VALUE,$R_STATUS" >> "$OUT_CSV"
        done
    done

    maxk="$(maxk_of "$dia")"
    [[ -z "$maxk" || "$maxk" -lt 1 ]] && { log "WARN: no diameter for $name, skipping linear"; continue; }
    KS=()
    for p in "${PCTS[@]}"; do
        k=$(( (p * maxk + 50) / 100 ))
        (( k >= 1 )) && KS+=("$k")
    done
    KS=($(printf '%s\n' "${KS[@]}" | sort -n -u))
    log "[run] $name linear: diameter maxK=$maxk, K set = ${KS[*]}"
    for b in "${LIN_BINS[@]}"; do
        for k in "${KS[@]}"; do
            log "[run] $b $name K=$k ..."
            run_one "" "$BIN/$b" "$dia" "$k"
            log "      -> $R_STATUS  value=$R_VALUE  time=${R_TIME}s  mem=${R_MEM}MB"
            echo "$b,$name,$k,$R_TIME,$R_MEM,$R_VALUE,$R_STATUS" >> "$OUT_CSV"
        done
    done
done

log ""
log "=== results written to $OUT_CSV ==="
if command -v column >/dev/null 2>&1; then column -t -s, "$OUT_CSV" >&2; else cat "$OUT_CSV" >&2; fi
