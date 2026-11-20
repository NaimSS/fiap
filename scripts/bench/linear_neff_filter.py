#!/usr/bin/env python3
"""Drops the rows of a linear_bench_all.sh CSV where K is past n_eff.

    python3 scripts/bench/linear_neff_filter.py IN.csv OUT_PREFIX

The linear array of an instance has only n_eff non-zero entries, so the
optimum stops improving at K = n_eff and every larger K repeats the same
problem. n_eff is taken as the first K with f(K) = f(K + 1).

Writes <prefix>_neff.csv and <prefix>_labeled_neff.csv, the second with the
solver names used in the paper's figure. This is how
linear_bench_min5_neff.csv and linear_bench_labeled_neff.csv were made, from
the minimum over 5 runs of linear_bench_all.sh.
"""

import collections
import sys

LABELS = {
    "dp_n2k": "Naive $O(n^2m)$",
    "dpnk": "CHT",
    "dpnlogW": "CHT + Lagrange",
    "dp_nk_smawk": "SMAWK",
    "dp_nk_larsch": "LARSCH",
    "dp_nlogW_larsch": "LARSCH + Lagrange",
}
HEADER = "dp,instance,K,time_ms,value,status"


def main():
    if len(sys.argv) != 3:
        sys.exit(__doc__)
    src, prefix = sys.argv[1], sys.argv[2]
    with open(src, newline="") as fh:
        lines = fh.read().splitlines(True)
    if not lines or lines[0].strip() != HEADER:
        sys.exit(f"{src}: expected header {HEADER}")
    rows = [ln.rstrip("\r\n").split(",") for ln in lines[1:]]
    if any(len(r) != 6 for r in rows):
        sys.exit(f"{src}: every row must have 6 fields")

    values = collections.defaultdict(dict)  # instance -> K -> set of values
    for dp, inst, k, _, val, status in rows:
        if status == "OK":
            values[inst].setdefault(int(k), set()).add(val)
    neff = {}
    for inst, by_k in values.items():
        for k in sorted(by_k):
            if len(by_k[k]) > 1:
                sys.exit(f"{inst} K={k}: solvers disagree on the value {sorted(by_k[k])}")
            if k + 1 in by_k and by_k[k] == by_k[k + 1]:
                neff[inst] = k
                break
        else:
            print(f"warning: {inst}: f(K) never repeats, n_eff not reached; keeping all rows",
                  file=sys.stderr)
            neff[inst] = max(by_k)

    kept = [(ln, r) for ln, r in zip(lines[1:], rows) if int(r[2]) <= neff.get(r[1], -1)]
    with open(prefix + "_neff.csv", "w", newline="") as fh:
        fh.write(lines[0])
        fh.writelines(ln for ln, _ in kept)
    with open(prefix + "_labeled_neff.csv", "w", newline="") as fh:
        fh.write(lines[0])
        for ln, r in kept:
            if r[0] not in LABELS:
                sys.exit(f"no label for solver {r[0]}")
            fh.write(LABELS[r[0]] + ln[len(r[0]):])
    print(f"kept {len(kept)} of {len(rows)} rows; n_eff per instance: "
          + ", ".join(f"{i}={neff[i]}" for i in sorted(neff)), file=sys.stderr)


if __name__ == "__main__":
    main()
