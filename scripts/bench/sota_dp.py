#!/usr/bin/env python3
"""lagrange_solver on the 19 networks and 9 budgets of Usberti et al. (2024).
Produces data/results/sota_dp.csv.

    python3 scripts/bench/sota_dp.py [--reps=5] [--force] [out.csv]

Needs bin/nocap/lagrange_solver (scripts/build.sh --no-cap).

The solver works on lengths rounded down to 10^-5 km, but the Tabares et al.
networks have lengths with up to 15 decimals. To get the cost of our cut set
on the original lengths we need the map from tree vertices back to network
edges, which bin/converter doesn't output, so this script builds the same
tree itself.

Columns: instance,n_nodes,n_edges,length_km,m,value,cost_exact_km2,time_ms,status
value is the optimum on the rounded lengths, in (10^-5 km)^2, so value / 1e10
is a lower bound for the original instance. cost_exact_km2 is our cut set on
the original lengths, an upper bound. time_ms is the fastest of the reps.
"""

import csv
import os
import re
import subprocess
import sys
from collections import deque
from decimal import ROUND_FLOOR, Decimal

REPO = os.path.dirname(os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
BIN = os.path.join(REPO, "bin", "nocap", "lagrange_solver")
INSTANCES = [
    "bus_13_3", "bus_29_1", "bus_32_1", "bus_83_11", "bus_135_8", "bus_201_3", "bus_873_7",
    "85nodes", "137nodes", "417nodes", "1080nodes",
    "JaboticabalJAB05", "JaboticabalATH05", "JaboticabalATH04", "JaboticabalATH03",
    "JaboticabalJAB03", "JaboticabalJAB04", "JaboticabalATH06", "JaboticabalJAB02",
]
MS = [1, 2, 3, 4, 5, 10, 15, 20, 25]
SCALE = Decimal(100000)


def read_fi(path):
    """(n_nodes, [feeders], [(u, v, length as Decimal km)])"""
    feeders, edges, n = [], [], 0
    with open(path) as fh:
        for line in fh:
            t = line.split()
            if not t:
                continue
            if t[0] == "p":
                n = int(t[2])
            elif t[0] == "f":
                feeders += [int(x) for x in t[1:]]
            elif t[0] == "e":
                edges.append((int(t[1]), int(t[2]), Decimal(t[3])))
    return n, feeders, edges


def expand(feeders, edges):
    """Same tree as bin/converter: one vertex per node and per edge, in BFS
    order from the feeder (or from a weight-0 root above several feeders).
    Returns (parent, exact), 0-based, exact[v] = original length on v."""
    adj = {}
    for j, (u, v, _) in enumerate(edges):
        adj.setdefault(u, []).append((v, j))
        adj.setdefault(v, []).append((u, j))
    parent, exact = [-1], [Decimal(0)]
    if len(feeders) == 1:
        index = {feeders[0]: 0}
    else:
        index = {}
        for f in feeders:
            index[f] = len(parent)
            parent.append(0)
            exact.append(Decimal(0))
    queue = deque(feeders)
    while queue:
        u = queue.popleft()
        for w, j in adj.get(u, []):
            if w in index:
                continue
            e = len(parent)              # vertex of edge j
            parent.append(index[u])
            exact.append(edges[j][2])
            index[w] = len(parent)
            parent.append(e)
            exact.append(Decimal(0))
            queue.append(w)
    if len(parent) != 2 * len(edges) + 1 + (len(feeders) if len(feeders) > 1 else 0):
        sys.exit("network is not a forest rooted at its feeders")
    return parent, exact


def cost_of_cuts(parent, weight, cuts):
    """Sum of squared component weights, cutting the edge above each v in cuts."""
    acc = list(weight)
    total = 0
    for v in range(len(parent) - 1, 0, -1):
        if v in cuts:
            total += acc[v] * acc[v]
        else:
            acc[parent[v]] += acc[v]
    return total + acc[0] * acc[0]


def solve(tree_text, m):
    out = subprocess.run([BIN, str(m), "1", "1"], input=tree_text, text=True,
                         capture_output=True, check=True).stdout.splitlines()
    head = out[0]
    value = int(head.split()[0])
    t = float(re.search(r" time: ([0-9.e+-]+)s", head).group(1))
    status = re.search(r"status: ([A-Za-z_-]+)", head).group(1)
    cuts = {int(x) - 1 for x in out[1].split(":")[1].split()}
    return value, t * 1000, status, cuts


def main():
    reps, force, out = 5, False, os.path.join(REPO, "data", "results", "sota_dp.csv")
    for a in sys.argv[1:]:
        if a.startswith("--reps="):
            reps = int(a.split("=")[1])
        elif a == "--force":
            force = True
        elif a in ("-h", "--help"):
            print(__doc__)
            return
        else:
            out = a
    if not os.path.exists(BIN):
        sys.exit(f"{BIN} missing; run scripts/build.sh --no-cap")
    if os.path.exists(out) and not force:
        sys.exit(f"{out} exists; pass --force or another path")

    rows = []
    for inst in INSTANCES:
        n, feeders, edges = read_fi(os.path.join(REPO, "data", "instances", inst + ".fi"))
        parent, exact = expand(feeders, edges)
        rounded = [int((w * SCALE).to_integral_value(rounding=ROUND_FLOOR)) for w in exact]
        tree = "%d %d\n%s\n%s\n" % (len(parent), len(parent),
                                    " ".join(map(str, rounded)),
                                    " ".join(str(p + 1) for p in parent[1:]))
        length = sum(exact)
        for m in MS:
            best = None
            for _ in range(reps):
                value, t, status, cuts = solve(tree, m)
                best = t if best is None else min(best, t)
            assert len(cuts) <= m and cost_of_cuts(parent, rounded, cuts) == value, (inst, m)
            cost = cost_of_cuts(parent, exact, cuts)
            rows.append([inst, n, len(edges), f"{length:.6f}", m, value,
                         f"{cost:.9f}", f"{best:.4f}", status])
            print(*rows[-1], sep=",")

    with open(out, "w", newline="") as fh:
        w = csv.writer(fh)
        w.writerow(["instance", "n_nodes", "n_edges", "length_km", "m", "value",
                    "cost_exact_km2", "time_ms", "status"])
        w.writerows(rows)
    print("wrote", os.path.relpath(out, REPO), file=sys.stderr)


if __name__ == "__main__":
    main()
