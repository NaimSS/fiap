#!/usr/bin/env python3
"""Random tree generator in converter format: 'n k' / weights / parents(2..n).

usage: gen_tree.py n [seed] [mode] [maxw] [zero_frac]
modes: random      - uniform random attachment
       netlike     - long paths with occasional branches (distribution-like)
       caterpillar - spine with hanging leaves
       path        - single path
       binary      - balanced-ish binary tree
       star        - all nodes hang off the root
Weights are random in [0, maxw] with a fraction zero_frac (default 0.35) of
forced zeros (mirrors real .fi conversions). A small maxw (1-3) makes many
cut counts tie at the critical lambda, which is what exposes false
GAP_OR_COLLINEAR results; maxw=1 zero_frac=0 is the
all-unit-weight case.
"""
import random
import sys

n = int(sys.argv[1])
seed = int(sys.argv[2]) if len(sys.argv) > 2 else 1
mode = sys.argv[3] if len(sys.argv) > 3 else "random"
maxw = int(sys.argv[4]) if len(sys.argv) > 4 else 10000
zero_frac = float(sys.argv[5]) if len(sys.argv) > 5 else 0.35
rng = random.Random(seed)

par = [0] * (n + 1)
if mode == "random":
    for i in range(2, n + 1):
        par[i] = rng.randint(1, i - 1)
elif mode == "netlike":
    # grow paths; occasionally start a new branch from a random old node
    cur = 1
    for i in range(2, n + 1):
        if rng.random() < 0.05:
            cur = rng.randint(1, i - 1)
        par[i] = cur
        cur = i
elif mode == "caterpillar":
    spine = 1
    for i in range(2, n + 1):
        if i % 2 == 0:
            par[i] = spine
            spine = i
        else:
            par[i] = spine
elif mode == "path":
    for i in range(2, n + 1):
        par[i] = i - 1
elif mode == "binary":
    for i in range(2, n + 1):
        par[i] = i // 2
elif mode == "star":
    for i in range(2, n + 1):
        par[i] = 1
else:
    raise SystemExit("unknown mode " + mode)

w = [0] * (n + 1)
for i in range(1, n + 1):
    w[i] = 0 if rng.random() < zero_frac else rng.randint(0 if zero_frac > 0 else 1, maxw)

out = []
out.append(f"{n} {n}")
out.append(" ".join(str(w[i]) for i in range(1, n + 1)))
out.append(" ".join(str(par[i]) for i in range(2, n + 1)))
print("\n".join(out))
