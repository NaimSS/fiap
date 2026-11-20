#!/bin/bash
# Runs the six linear DPs over a set of .fi.diameter files and K values.
#
#   bash scripts/bench/linear_bench_all.sh <instances> <ks> [out.csv] [--all-k]
#
# <instances> is a glob ('data/instances/*.fi.diameter') or a .txt with one
# path per line. <ks> is a list (10,20,30), a range (1-1080) or a .txt with
# one K per line. --all-k runs every K from 1 to the length of the linear
# array instead. K is the number of sections, so K - 1 indicators.
#
# Output (default out/linear_bench_all.csv): dp,instance,K,time_ms,value,status
# The dp column uses the names dpnk and dpnlogW for dp_nk and dp_nlogW, as in
# the result CSVs.
#
# Env: TIMEOUT_S (600), MEM_CAP_MB (3072).

set -uo pipefail

ROOT="$(cd "$(dirname "$0")/../.." && pwd)"

TIMEOUT_S="${TIMEOUT_S:-600}"
MEM_CAP_MB="${MEM_CAP_MB:-3072}"
source "$ROOT/src/helpers/run_tracked.sh"
trap rt_cleanup EXIT

if [[ $# -lt 2 ]]; then
    echo "Usage: $0 <instances_glob_or_txt> <k_values> [output.csv] [--all-k]"
    echo "  e.g. $0 'data/instances/*.fi.diameter' '10,20,30' out/results.csv"
    echo "  e.g. $0 'data/instances/1080nodes.fi.diameter' '1-207' out/results.csv"
    exit 1
fi

INSTANCES_ARG="$1"
K_ARG="$2"
OUTPUT="$ROOT/out/linear_bench_all.csv"
ALL_K=false
for arg in "${@:3}"; do
    case "$arg" in
        --all-k) ALL_K=true ;;
        *)       OUTPUT="$arg" ;;
    esac
done
mkdir -p "$(dirname "$OUTPUT")"

# binary:name in the dp column
DP_TABLE=(
    "dp_n2k:dp_n2k"
    "dp_nk:dpnk"
    "dp_nlogW:dpnlogW"
    "dp_nk_smawk:dp_nk_smawk"
    "dp_nk_larsch:dp_nk_larsch"
    "dp_nlogW_larsch:dp_nlogW_larsch"
)

bash "$ROOT/scripts/build.sh" > /dev/null || { echo "Error: scripts/build.sh failed"; exit 1; }

declare -a DP_BINS=() DP_NAMES=()
for entry in "${DP_TABLE[@]}"; do
    DP_BINS+=("$ROOT/bin/${entry%%:*}")
    DP_NAMES+=("${entry##*:}")
done

declare -a INSTANCE_FILES=()
if [[ "$INSTANCES_ARG" == *.txt && -f "$INSTANCES_ARG" ]]; then
    while IFS= read -r fpath || [[ -n "$fpath" ]]; do
        [[ -z "$fpath" || "$fpath" == \#* ]] && continue
        INSTANCE_FILES+=("$fpath")
    done < "$INSTANCES_ARG"
else
    for f in $INSTANCES_ARG; do
        [[ -f "$f" ]] && INSTANCE_FILES+=("$f")
    done
fi

if [[ ${#INSTANCE_FILES[@]} -eq 0 ]]; then
    echo "Error: no .fi.diameter files found for: $INSTANCES_ARG"
    exit 1
fi

declare -a K_VALUES=()
if [[ "$K_ARG" == *.txt && -f "$K_ARG" ]]; then
    while IFS= read -r k || [[ -n "$k" ]]; do
        [[ -z "$k" || "$k" == \#* ]] && continue
        K_VALUES+=("$k")
    done < "$K_ARG"
elif [[ "$K_ARG" =~ ^([0-9]+)(-|\.\.)([0-9]+)$ ]]; then
    start="${BASH_REMATCH[1]}"
    end="${BASH_REMATCH[3]}"
    if (( 10#$start > 10#$end )); then
        echo "Error: empty K range: $K_ARG"
        exit 1
    fi
    for k in $( seq "$start" "$end" ); do
        K_VALUES+=("$k")
    done
else
    IFS=',' read -ra K_VALUES <<< "$K_ARG"
fi

if [[ "$ALL_K" != true ]]; then
    if [[ ${#K_VALUES[@]} -eq 0 ]]; then
        echo "Error: no K values found for: $K_ARG"
        exit 1
    fi
    for k in "${K_VALUES[@]}"; do
        if ! [[ "$k" =~ ^[0-9]+$ ]] || (( 10#$k < 1 )); then
            echo "Error: K must be a positive integer, got '$k' (from: $K_ARG)"
            exit 1
        fi
    done
fi

# --all-k: K goes up to 2 * (trunk arcs) + 1, the length of the linear array
declare -a INSTANCE_MAX_K=()
if [[ "$ALL_K" == true ]]; then
    for fpath in "${INSTANCE_FILES[@]}"; do
        mk=$(grep "^diameter" "$fpath" | awk '{print 2*(NF-1)+1}')
        INSTANCE_MAX_K+=("${mk:-0}")
    done
fi

echo "dp,instance,K,time_ms,value,status" > "$OUTPUT"

if [[ "$ALL_K" == true ]]; then
    total=0
    for i in "${!INSTANCE_FILES[@]}"; do
        total=$(( total + INSTANCE_MAX_K[$i] * ${#DP_BINS[@]} ))
    done
else
    total=$(( ${#INSTANCE_FILES[@]} * ${#K_VALUES[@]} * ${#DP_BINS[@]} ))
fi
done_count=0

for i in "${!INSTANCE_FILES[@]}"; do
    fpath="${INSTANCE_FILES[$i]}"
    iname="$(basename "$fpath")"

    if [[ "$ALL_K" == true ]]; then
        maxK="${INSTANCE_MAX_K[$i]}"
        K_ITER=( $(seq 1 "$maxK") )
    else
        K_ITER=("${K_VALUES[@]}")
    fi

    for K in "${K_ITER[@]}"; do
        for j in "${!DP_BINS[@]}"; do
            dp_bin="${DP_BINS[$j]}" dp_name="${DP_NAMES[$j]}"
            done_count=$(( done_count + 1 ))
            printf "[%d/%d] %-16s  %-40s  K=%s\n" "$done_count" "$total" "$dp_name" "$iname" "$K"

            run_tracked "$TIMEOUT_S" "$MEM_CAP_MB" "$dp_bin" "$fpath" "$K" < /dev/null
            value="" time_ms="" status="$RT_STATUS"
            if [[ "$status" == OK ]]; then
                value=$( grep -m1 "^k = " "$RT_OUT" | sed 's/.*: //' )
                time_ms=$( grep -m1 "^time_ms:" "$RT_OUT" | sed 's/time_ms: //' )
                [[ -z "$value" || -z "$time_ms" ]] && status="ERROR(parse)"
            fi
            if [[ "$status" != OK ]]; then
                value="" time_ms=""
                printf "    -> %s: %s\n" "$status" "$( grep -v -e '^mem_cap_peak_kb:' -e 'resident set size' -e 'terminated abnormally' "$RT_ERR" | head -1 )"
            fi
            echo "${dp_name},${iname},${K},${time_ms},${value},${status}" >> "$OUTPUT"
        done
    done
done

echo ""
echo "Results written to: $OUTPUT"
echo ""
column -t -s ',' "$OUTPUT" 2>/dev/null || cat "$OUTPUT"
