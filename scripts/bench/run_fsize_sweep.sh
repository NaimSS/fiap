#!/usr/bin/env bash
# Envelope sizes (F_max, sum |F|) and timings of the stats builds of
# lagrange_solver ("lagrange") and default_solver ("nolag").
#
#   ./run_fsize_sweep.sh [-o OUT.csv] [-k K1,K2,...] [-s lagrange,nolag] [FILE.fi ...]
#
# Defaults: all data/instances/*.fi, K = 1,2,3,4,5,10,...,50, both solvers,
# OUT = data/results/fsize_sweep.csv (it won't overwrite, so pass -o).
# Binaries, converted trees and stderr of failed runs go under
# out/fsize_sweep/<OUT name>/.
#
# Once nolag times out or runs out of memory on an instance, its larger K are
# marked SKIPPED, since its cost only grows with K.
#
# Env: TIMEOUT_S (600), MEM_CAP_MB (1024).
#
# fsize_sweep_bigK.csv was made with
#   ./run_fsize_sweep.sh -o data/results/fsize_sweep_bigK.csv \
#       -s lagrange -k 1,2,5,10,25,50,75,100,150,200,300,400 \
#       data/instances/rede_mt_feeder_{spd31,alt03,gtu22,aux12,gav23}.fi \
#       data/instances/rede_mt_municipality_{paulinia,campinas}.fi

set -uo pipefail

SCRIPT_DIR="$(cd "$(dirname "$0")" && pwd)"
REPO="$(cd "$SCRIPT_DIR/../.." && pwd)"

TIMEOUT_S="${TIMEOUT_S:-600}"
CAP_MB="${MEM_CAP_MB:-1024}"
OUT="$REPO/data/results/fsize_sweep.csv"
KS=(1 2 3 4 5 10 15 20 25 30 35 40 50)
SOLVERS=(lagrange nolag)

while (($# > 0)); do
    case "$1" in
        -o) OUT="$2"; shift 2 ;;
        -k) IFS=',' read -ra KS <<< "$2"; shift 2 ;;
        -s) IFS=',' read -ra SOLVERS <<< "$2"; shift 2 ;;
        -h) sed -n '2,/^set -uo/p' "$0" | sed '$d' | sed 's/^# \{0,1\}//'; exit 0 ;;
        --) shift; break ;;
        -*) echo "unknown option: $1" >&2; exit 1 ;;
        *) break ;;
    esac
done

FILES=("$@")

for s in "${SOLVERS[@]}"; do
    [[ "$s" == lagrange || "$s" == nolag ]] || { echo "error: unknown solver '$s'" >&2; exit 1; }
done
for k in "${KS[@]}"; do
    [[ "$k" =~ ^[1-9][0-9]*$ ]] || { echo "error: bad K '$k' (need an integer >= 1)" >&2; exit 1; }
done
((${#KS[@]} > 0)) || { echo "error: empty K list" >&2; exit 1; }
KS=($(printf '%s\n' "${KS[@]}" | sort -n -u))
for f in ${FILES[@]+"${FILES[@]}"}; do
    [[ "$f" == *.fi && -f "$f" ]] || { echo "error: not a .fi file: $f" >&2; exit 1; }
done
if ((${#FILES[@]} == 0)); then
    shopt -s nullglob
    FILES=("$REPO"/data/instances/*.fi)
    shopt -u nullglob
fi
if ((${#FILES[@]} == 0)); then echo "error: no instances found" >&2; exit 1; fi

if [[ -e "$OUT" ]]; then
    echo "error: $OUT exists; pass -o NEW.csv or move it away" >&2
    exit 1
fi
if [[ ! -d "$(dirname "$OUT")" || ! -w "$(dirname "$OUT")" ]]; then
    echo "error: cannot write $OUT (directory missing or read-only)" >&2
    exit 1
fi

WORK="$REPO/out/fsize_sweep/$(basename "$OUT" .csv)"
BIN_DIR="$WORK/bin"
TREE_DIR="$WORK/trees"
LOG_DIR="$WORK/logs"
if [[ -e "$WORK" ]]; then
    echo "error: work directory $WORK exists (from an earlier run with this OUT name); move it away" >&2
    exit 1
fi

source "$REPO/src/helpers/run_tracked.sh"
trap rt_cleanup EXIT

# copy the binaries so a rebuild during the sweep can't change them
mkdir -p "$BIN_DIR" "$TREE_DIR" "$LOG_DIR"
bash "$REPO/scripts/build.sh" >&2 || { echo "error: build failed" >&2; exit 1; }
cp "$REPO/bin/converter" "$BIN_DIR/converter"
cp "$REPO/bin/lagrange_solver_stats" "$BIN_DIR/lagrange"
cp "$REPO/bin/default_solver_stats" "$BIN_DIR/nolag"

# convert, then sort by n
TREES=()
for f in "${FILES[@]}"; do
    tree="$TREE_DIR/$(basename "$f" .fi)"
    [[ -e "$tree" ]] && { echo "error: two instances named $(basename "$f" .fi)" >&2; exit 1; }
    "$BIN_DIR/converter" < "$f" > "$tree.tmp" 2>"$LOG_DIR/converter_$(basename "$f" .fi).err" \
        || { echo "error: converter failed on $f" >&2; rm -f "$tree.tmp"; exit 1; }
    mv "$tree.tmp" "$tree"
    TREES+=("$tree")
done
sorted=()
while IFS=$'\t' read -r _ path; do sorted+=("$path"); done < <(
    for t in "${TREES[@]}"; do printf '%s\t%s\n' "$(head -1 "$t" | awk '{print $1}')" "$t"; done | sort -n)
TREES=("${sorted[@]}")

COLS="solver,instance,n,K,optimal_value,time_s,value_time_s,recover_time_s,cuts,bs_steps,max_segs,merge_calls,merge_seg_work,const_calls,const_walk,seg_samples,seg_sum,rec_merge_calls,rec_merge_seg_work,rec_const_calls,rec_const_walk,rec_seg_samples,rec_seg_sum,rec_builds,rec_positions,peak_heap_mb,solver_status,status"
echo "$COLS" > "$OUT"
NCOLS=$(awk -F, '{print NF}' <<< "$COLS")

# fail_row <solver status> <status>: a row with the numbers left empty
fail_row(){
    local empties; empties=$(printf '%*s' $((NCOLS - 6)) '' | tr ' ' ',')
    printf '%s,%s,%s,%s%s,%s,%s\n' "$solver" "$name" "$n" "$k" "$empties" "${1%%,*}" "$2" >> "$OUT"
}

total=$(( ${#SOLVERS[@]} * ${#TREES[@]} * ${#KS[@]} )); count=0
for solver in "${SOLVERS[@]}"; do
    BIN="$BIN_DIR/$solver"
    for t in "${TREES[@]}"; do
        name="$(basename "$t")"
        n="$(head -1 "$t" | awk '{print $1}')"
        dead=""
        for k in "${KS[@]}"; do
            count=$((count + 1))
            printf '\r[%d/%d] %-8s %-32s n=%-7s K=%-4s' "$count" "$total" "$solver" "$name" "$n" "$k" >&2
            if [[ -n "$dead" ]]; then fail_row "" SKIPPED; continue; fi

            if [[ "$solver" == lagrange ]]; then
                run_tracked "$TIMEOUT_S" "$CAP_MB" "$BIN" "$k" 1 < "$t"
            else
                run_tracked "$TIMEOUT_S" "$CAP_MB" "$BIN" "$k" < "$t"
            fi
            rt_parse_tree_dp
            err="$LOG_DIR/${solver}_${name}_K${k}.err"
            cp "$RT_ERR" "$err"

            if [[ "$RT_STATUS" != OK && -z "$TD_VALUE" ]]; then
                fail_row "$TD_SOLVER_STATUS" "$RT_STATUS"
                [[ "$solver" == nolag && ( "$RT_STATUS" == TIMEOUT || "$RT_STATUS" == OOM ) ]] \
                    && dead="$RT_STATUS"
                echo >&2; echo "  $solver $name K=$k: $RT_STATUS  [stderr: $err]" >&2
                continue
            fi

            row=$(python3 - "$RT_OUT" "$solver" "$name" "$n" "$k" "$RT_MEM_MB" "$RT_STATUS" 2>>"$err" <<'PY'
import re, sys
out, solver, name, n, k, mem, status = sys.argv[1:8]
text = open(out).read()

def grab(pat):
    m = re.search(pat, text, re.M)
    return m.group(1) if m else ""

def field(tag, key):
    m = re.search(r"^" + tag + r":(.*)$", text, re.M)
    if not m:
        return ""
    m = re.search(r"\b" + key + r"=(\d+)", m.group(1))
    return m.group(1) if m else ""

num = r"([0-9.]+(?:e[-+]?\d+)?)"
value    = grab(r"^(-?\d+)\s+- k=")
time_s   = grab(r", time:\s*" + num + "s")
sstat    = grab(r"status:\s*([\w-]+)")
bs_steps = grab(r"bs_steps:\s*(\d+)") or field("hlstats", "probes")
stats = [field("hlstats", c) for c in ("maxSegs", "mergeCalls", "mergeSegWork",
                                        "constCalls", "constWalk", "segSamples", "segSum")]
rec = [field("hlrecover", c) for c in ("mergeCalls", "mergeSegWork", "constCalls",
                                        "constWalk", "segSamples", "segSum",
                                        "recBuilds", "recPositions")]
missing = [c for c, v in (("value", value), ("time", time_s), ("status", sstat),
                          ("hlstats", stats[0]), ("hlrecover", rec[0])) if v == ""]
if missing:
    sys.exit("parse: missing " + ", ".join(missing))
print(",".join([solver, name, n, k, value, time_s,
                grab(r"value_time:\s*" + num + "s"), grab(r"recover_time:\s*" + num + "s"),
                grab(r"cuts:\s*(\d+)"), bs_steps] + stats + rec + [mem, sstat, status]))
PY
            )
            if (($? != 0)) || [[ -z "$row" ]]; then
                fail_row "$TD_SOLVER_STATUS" "ERROR(parse)"
                echo >&2; echo "  $solver $name K=$k: could not parse  [stderr: $err]" >&2
            else
                echo "$row" >> "$OUT"
                [[ "$RT_STATUS" == OK ]] && rm -f "$err"
            fi
        done
    done
done
echo >&2
echo "Done. Results in $OUT" >&2
