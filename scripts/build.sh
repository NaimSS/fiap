#!/usr/bin/env bash
# Compiles the helpers and all solvers into bin/.
#
#   bash scripts/build.sh [--no-cap] [--force]
#
# By default solvers are built with src/helpers/mem_cap.h, so MEM_CAP_MB caps
# their heap and MEM_CAP_REPORT=1 prints the peak. --no-cap builds them
# without it into bin/nocap/ (used by the comparison scripts and
# bench_reps.py). --force rebuilds everything. CXX picks the compiler.
set -euo pipefail

ROOT="$(cd "$(dirname "$0")/.." && pwd)"
CXX="${CXX:-g++}"
CAP=1
FORCE=0
for arg in "$@"; do
    case "$arg" in
        --no-cap) CAP=0 ;;
        --force)  FORCE=1 ;;
        -h|--help) sed -n '2,/^set -euo/p' "$0" | sed '$d' | sed 's/^# \{0,1\}//'; exit 0 ;;
        *) echo "build.sh: unknown option '$arg'" >&2; exit 1 ;;
    esac
done

MEM_CAP_H="$ROOT/src/helpers/mem_cap.h"
if ((CAP)); then
    BIN="$ROOT/bin"
    SOLVER_FLAGS=(-std=c++17 -O2 -include "$MEM_CAP_H")
else
    BIN="$ROOT/bin/nocap"
    SOLVER_FLAGS=(-std=c++17 -O2)
fi
mkdir -p "$BIN"

# build <name> <src> <deps...> -- <flags...>
build(){
    local name="$1" src="$2"; shift 2
    local deps=()
    while (($# > 0)) && [[ "$1" != -- ]]; do deps+=("$1"); shift; done
    shift
    local out="$BIN/$name" stale=$FORCE d
    [[ -x "$out" ]] || stale=1
    for d in "$src" ${deps[@]+"${deps[@]}"}; do [[ "$d" -nt "$out" ]] && stale=1; done
    if ((stale == 0)); then echo "up to date: ${out#$ROOT/}"; return 0; fi
    echo "building ${out#$ROOT/}"
    "$CXX" "$@" -I "$(dirname "$src")" "$src" -o "$out"
}

H="$ROOT/src/helpers"
L="$ROOT/src/linear_dp"
G="$ROOT/src/general_dp"

build converter       "$H/converter.cpp"       -- -std=c++14 -O2
build fi_to_diameter  "$H/fi_to_diameter.cpp"  -- -std=c++17 -O2
build convert_to_line "$H/convert_to_line.cpp" -- -std=c++17 -O2

LIN_DEPS=("$L/convert_utils.h" "$L/swamk.h" "$L/larsch.h" "$MEM_CAP_H")
for d in dp_n2k dp_nk dp_nk_smawk dp_nk_larsch dp_nlogW dp_nlogW_larsch; do
    build "$d" "$L/$d.cpp" "${LIN_DEPS[@]}" -- "${SOLVER_FLAGS[@]}"
done

TREE_DEPS=("$G/slope_min_helper.h" "$MEM_CAP_H")
build lagrange_solver       "$G/lagrange_solver.cpp"       "${TREE_DEPS[@]}" -- "${SOLVER_FLAGS[@]}"
build default_solver        "$G/default_solver.cpp"        "${TREE_DEPS[@]}" -- "${SOLVER_FLAGS[@]}"
build lagrange_solver_stats "$G/lagrange_solver_stats.cpp" "${TREE_DEPS[@]}" -- "${SOLVER_FLAGS[@]}" -DTREE_STATS
build default_solver_stats  "$G/default_solver_stats.cpp"  "${TREE_DEPS[@]}" -- "${SOLVER_FLAGS[@]}"
build lagrange_solver_without_hld "$G/lagrange_solver_without_hld.cpp" \
    "$G/slope_min_helper_without_hld.h" "$G/debug.h" "$MEM_CAP_H" -- "${SOLVER_FLAGS[@]}"

echo "all binaries in ${BIN#$ROOT/}/"
