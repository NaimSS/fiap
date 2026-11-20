#!/usr/bin/env python3
"""Re-runs the linear DP benchmark several times over the same grid.

    python3 scripts/bench/bench_reps.py GRID.csv out/linear_reps [--reps 5]

The (dp, instance, K) grid is read from an existing linear_bench_all.sh CSV.
Each pass goes over the grid in a different random order, so a burst of
system load doesn't land on every repetition of the same run. Writes every
sample to <prefix>_raw.csv and median/min/max to <prefix>.csv. Exits 1 if a
run returns different values across repetitions.

Uses the bin/nocap/ builds; --capped uses bin/ instead.
"""

import argparse
import csv
import os
import random
import re
import statistics
import subprocess
import sys
import time

HERE = os.path.dirname(os.path.abspath(__file__))
ROOT = os.path.dirname(os.path.dirname(HERE))
INSTANCE_DIR = os.path.join(ROOT, "data", "instances")
BIN_DIR = os.path.join(ROOT, "bin", "nocap")

# dp column name -> binary in bin/nocap/
SOLVERS = {
    "dp_n2k": "dp_n2k",
    "dpnk": "dp_nk",
    "dp_nk_smawk": "dp_nk_smawk",
    "dp_nk_larsch": "dp_nk_larsch",
    "dpnlogW": "dp_nlogW",
    "dp_nlogW_larsch": "dp_nlogW_larsch",
}

TIME_RE = re.compile(r"time_ms:\s*([\d.]+)")
VAL_RE = re.compile(r"k = \d+:\s*(-?\d+)")


def compile_all(capped):
    subprocess.run(
        ["bash", os.path.join(ROOT, "scripts", "build.sh")] + ([] if capped else ["--no-cap"]),
        check=True, stdout=subprocess.DEVNULL,
    )


def run_once(binary, instance, k):
    """(time_ms, value), or (None, None) if the run failed."""
    path = os.path.join(INSTANCE_DIR, instance)
    p = subprocess.run(
        [os.path.join(BIN_DIR, SOLVERS[binary]), path, str(k)],
        capture_output=True,
        text=True,
    )
    if p.returncode != 0:
        return None, None
    mt, mv = TIME_RE.search(p.stdout), VAL_RE.search(p.stdout)
    if not mt or not mv:
        return None, None
    return float(mt.group(1)), mv.group(1)


def main():
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("grid_csv", help="existing bench CSV defining the pair grid")
    ap.add_argument("out_prefix", help="writes <prefix>_raw.csv and <prefix>.csv")
    ap.add_argument("--reps", type=int, default=5)
    ap.add_argument("--capped", action="store_true",
                    help="use the mem_cap.h builds in bin/ (honour MEM_CAP_MB) instead of bin/nocap/")
    args = ap.parse_args()
    global BIN_DIR
    if args.capped:
        BIN_DIR = os.path.join(ROOT, "bin")

    with open(args.grid_csv) as fh:
        grid = sorted({(r["dp"], r["instance"], int(r["K"])) for r in csv.DictReader(fh)})
    unknown = sorted({d for d, _, _ in grid if d not in SOLVERS})
    if unknown:
        sys.exit(f"unknown solver names in {args.grid_csv}: {unknown} (expected {sorted(SOLVERS)})")
    missing = {i for _, i, _ in grid if not os.path.exists(os.path.join(INSTANCE_DIR, i))}
    if missing:
        sys.exit(f"missing instances in {INSTANCE_DIR}: {sorted(missing)}")

    compile_all(args.capped)
    print(f"{len(grid)} pairs x {args.reps} reps = {len(grid) * args.reps} runs", flush=True)

    raw_path = f"{args.out_prefix}_raw.csv"
    samples = {}
    values = {}
    t_start = time.time()

    with open(raw_path, "w", newline="") as fh:
        w = csv.writer(fh)
        w.writerow(["rep", "dp", "instance", "K", "time_ms", "value"])
        for rep in range(args.reps):
            order = list(grid)
            random.Random(1234 + rep).shuffle(order)
            for n, key in enumerate(order):
                dp, inst, k = key
                t, v = run_once(dp, inst, k)
                w.writerow([rep, dp, inst, k, "ERROR" if t is None else f"{t:.3f}", v or ""])
                if t is not None:
                    samples.setdefault(key, []).append(t)
                    values.setdefault(key, set()).add(v)
                if n % 2000 == 0:
                    fh.flush()
                    el = time.time() - t_start
                    done = rep * len(grid) + n + 1
                    tot = len(grid) * args.reps
                    print(
                        f"rep {rep + 1}/{args.reps}  {n}/{len(grid)}  "
                        f"elapsed {el / 60:.1f} min  eta {(el / done * (tot - done)) / 60:.1f} min",
                        flush=True,
                    )

    mismatched = sorted(k for k, vs in values.items() if len(vs) > 1)

    out_path = f"{args.out_prefix}.csv"
    with open(out_path, "w", newline="") as fh:
        w = csv.writer(fh)
        w.writerow(["dp", "instance", "K", "time_ms", "time_min", "time_max", "n_ok", "value"])
        for key in grid:
            dp, inst, k = key
            s = samples.get(key, [])
            if not s:
                w.writerow([dp, inst, k, "ERROR", "ERROR", "ERROR", 0, ""])
                continue
            w.writerow([
                dp, inst, k,
                f"{statistics.median(s):.3f}",
                f"{min(s):.3f}", f"{max(s):.3f}", len(s),
                "MISMATCH" if key in mismatched else next(iter(values[key])),
            ])

    print(f"\nwrote {raw_path} and {out_path} in {(time.time() - t_start) / 60:.1f} min", flush=True)
    if mismatched:
        for key in mismatched:
            print(f"value mismatch across reps: {key}: {sorted(values[key])}", file=sys.stderr)
        sys.exit(f"{len(mismatched)} pair(s) returned different values across repetitions")


if __name__ == "__main__":
    main()
