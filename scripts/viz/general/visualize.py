#!/usr/bin/env python3
"""Draws the network with lagrange_solver's placement of at most K indicators,
one colour per section.

    python visualize.py <instance.fi> <K> [--style ...] [--out ...] [--show]
                        [--no-title] [--allow-uncertified]

Runs the converter and lagrange_solver (compiled into out/bin/viz/), then maps
the cut vertices back to the edges of the .fi file. Refuses to draw a
GAP_OR_COLLINEAR result unless --allow-uncertified.
"""
import argparse
import colorsys
import os
import subprocess
import sys
from collections import deque
from pathlib import Path

SCRIPT_DIR = Path(__file__).resolve().parent
REPO = SCRIPT_DIR.parent.parent.parent  # scripts/viz/general -> repo root
HELPERS = REPO / "src" / "helpers"
GENERAL_DP = REPO / "src" / "general_dp"
BIN = REPO / "out" / "bin" / "viz"


# grey for a segment that dominates the drawing
_MUTED_DOMINANT = (0.74, 0.74, 0.77)


def _distinct_palette(n):
    """n distinct colours: golden-angle hue steps, cycling through a few
    saturation/value pairs so hues that come back close still differ."""
    if n <= 0:
        return []
    golden = 0.6180339887498949
    shades = [
        (0.90, 0.93), (0.62, 0.82), (0.96, 0.60), (0.55, 0.98),
        (0.82, 0.72), (0.72, 0.91), (1.00, 0.80), (0.60, 0.67),
        (0.88, 0.52), (0.68, 0.78), (0.95, 0.70), (0.58, 0.88),
    ]
    out, h = [], 0.045
    for i in range(n):
        s, v = shades[i % len(shades)]
        out.append(colorsys.hsv_to_rgb(h % 1.0, s, v))
        h += golden
    return out


def _color_dist2(a, b):
    """Squared redmean RGB distance."""
    rmean = (a[0] + b[0]) / 2.0
    dr, dg, db = a[0] - b[0], a[1] - b[1], a[2] - b[2]
    return (2 + rmean) * dr * dr + 4 * dg * dg + (3 - rmean) * db * db


def _segment_adjacency(n_comp, comp, parent, n, nodes, comp_of_orig, used_orig,
                       seg_size):
    """Pairs of components that should get different colours: the two sides
    of each cut, components drawn close to each other (same 3x3 grid block),
    and any two of the 22 largest."""
    adj = [set() for _ in range(n_comp)]
    for v in range(2, n + 1):
        a, b = comp[v], comp[parent[v]]
        if a != b:
            adj[a].add(b)
            adj[b].add(a)
    if seg_size is not None and n_comp > 1:
        big = sorted(range(n_comp), key=lambda c: -seg_size[c])[:min(n_comp, 22)]
        for i in range(len(big)):
            for j in range(i + 1, len(big)):
                adj[big[i]].add(big[j])
                adj[big[j]].add(big[i])
    if used_orig:
        xs = [nodes[o][0] for o in used_orig]
        ys = [nodes[o][1] for o in used_orig]
        minx, miny = min(xs), min(ys)
        ncell = max(1, int(len(used_orig) ** 0.5))
        cw = (max(xs) - minx) / ncell or 1.0
        ch = (max(ys) - miny) / ncell or 1.0
        grid = {}
        for o in used_orig:
            x, y = nodes[o]
            cell = (int((x - minx) / cw), int((y - miny) / ch))
            grid.setdefault(cell, set()).add(comp_of_orig[o])
        for (cx, cy), _here in grid.items():
            near = set()
            for dx in (-1, 0, 1):
                for dy in (-1, 0, 1):
                    near |= grid.get((cx + dx, cy + dy), set())
            near = list(near)
            if len(near) > 24:
                continue
            for i in range(len(near)):
                for j in range(i + 1, len(near)):
                    adj[near[i]].add(near[j])
                    adj[near[j]].add(near[i])
    return adj


def _assign_segment_colors(n_comp, adj, palette, seg_size):
    """Greedy colouring, largest components first: each one takes the unused
    palette colour farthest from its coloured neighbours."""
    assigned, used = {}, [False] * len(palette)
    order = sorted(range(n_comp),
                   key=lambda c: (-seg_size[c], -len(adj[c]), c))
    for c in order:
        ncols = [palette[assigned[nb]] for nb in adj[c] if nb in assigned]
        best_idx, best_score = None, None
        for idx in range(len(palette)):
            if used[idx]:
                continue
            # with no coloured neighbour, take the palette in order
            score = min((_color_dist2(palette[idx], nc) for nc in ncols),
                        default=-idx)
            if best_score is None or score > best_score:
                best_idx, best_score = idx, score
        used[best_idx] = True
        assigned[c] = best_idx
    return {c: palette[i] for c, i in assigned.items()}


def compute_segment_colors(n_comp, comp, parent, n, nodes, comp_of_orig,
                           used_orig, seg_size, mute_dominant=True):
    """Component -> RGB. seg_size[c] is the number of drawn edges of c. With
    mute_dominant, a segment much larger than all others is drawn grey."""
    palette = _distinct_palette(n_comp)
    if n_comp > 600:
        # too many for the quadratic colouring
        cmap = {c: palette[c] for c in range(n_comp)}
    else:
        adj = _segment_adjacency(n_comp, comp, parent, n, nodes, comp_of_orig,
                                 used_orig, seg_size)
        cmap = _assign_segment_colors(n_comp, adj, palette, seg_size)
    if mute_dominant and n_comp >= 6:
        dom = max(range(n_comp), key=lambda c: seg_size[c])
        total = sum(seg_size)
        second = sorted(seg_size, reverse=True)[1]
        if total and seg_size[dom] >= 0.12 * total and seg_size[dom] >= 1.5 * second:
            cmap[dom] = _MUTED_DOMINANT
    return cmap


def build_if_needed(src: Path, binary: Path, deps=()):
    """Compile src unless binary is newer than it and its deps."""
    if binary.exists() and all(binary.stat().st_mtime >= p.stat().st_mtime
                               for p in (src, *deps)):
        return
    binary.parent.mkdir(parents=True, exist_ok=True)
    print(f"building {binary.name} ...", file=sys.stderr)
    subprocess.run(["g++", "-O2", "-std=c++17", str(src), "-o", str(binary)],
                   check=True)


def parse_fi(path: Path):
    """(node id -> (x, y), feeder ids)"""
    nodes = {}
    feeders = []
    with open(path) as f:
        for line in f:
            s = line.strip()
            if not s or s[0] == "#":
                continue
            tok = s.split()
            key = tok[0]
            if key == "f":
                feeders.extend(int(t) for t in tok[1:])
            elif key == "v" and len(tok) >= 4:
                nodes[int(tok[1])] = (float(tok[2]), float(tok[3]))
    return nodes, feeders


def run_converter(fi_path: Path):
    """Converter output with EMIT_META=1: (text, n, weights, parent, meta),
    1-based, where meta[v] is the .fi node of tree vertex v or -1 for an edge
    vertex or the super-root."""
    build_if_needed(HELPERS / "converter.cpp", BIN / "converter")
    text = fi_path.read_text()
    out = subprocess.run(
        [str(BIN / "converter")], input=text, capture_output=True,
        text=True, check=True, env={**os.environ, "EMIT_META": "1"},
    ).stdout
    lines = out.splitlines()
    if len(lines) < 4:
        sys.exit("error: converter did not emit a metadata line "
                 "(is EMIT_META supported by helpers/converter.cpp?)")
    n = int(lines[0].split()[0])
    weights = [0] + [int(x) for x in lines[1].split()]
    parents_raw = [int(x) for x in lines[2].split()]
    meta = [0] + [int(x) for x in lines[3].split()]
    if not (len(weights) == n + 1 == len(meta) and len(parents_raw) == n - 1):
        sys.exit("error: converter output is inconsistent")
    parent = [0, 0] + parents_raw
    return out, n, weights, parent, meta


def run_dp(conv_text: str, K: int):
    """(value, cut vertices, status)"""
    binary = BIN / "lagrange_solver"
    build_if_needed(GENERAL_DP / "lagrange_solver.cpp", binary,
                    deps=(GENERAL_DP / "slope_min_helper.h",))
    # K, skip brute, print the cuts, root at node 1 (the converter's root)
    out = subprocess.run(
        [str(binary), str(K), "1", "1", "0"],
        input=conv_text, capture_output=True, text=True, check=True,
    ).stdout
    value, cuts, status = None, None, None
    for line in out.splitlines():
        s = line.strip()
        if s.startswith("cut edge child-vertices:"):
            cuts = [int(x) for x in s[len("cut edge child-vertices:"):].split()]
        elif value is None and " - k=" in s:
            value = int(s.split()[0])
            if ", status: " in s:
                status = s.split(", status: ", 1)[1].split(",", 1)[0].strip()
    if value is None or cuts is None or status is None:
        sys.exit(f"error: could not parse DP output:\n{out}")
    return value, cuts, status


def main():
    ap = argparse.ArgumentParser(
        description="Visualize the general tree DP's optimal FI placement.")
    ap.add_argument("instance", help="path to a .fi instance file")
    ap.add_argument("K", type=int,
                    help="number of fault indicators (yields K+1 segments)")
    ap.add_argument("--style", choices=["auto", "small", "big"], default="auto",
                    help="rendering style (auto: small if <=60 nodes, else big)")
    ap.add_argument("--out", default=None, help="output PNG path")
    ap.add_argument("--show", action="store_true", help="also display the figure")
    ap.add_argument("--no-title", action="store_true",
                    help="omit the title (for figures that carry a LaTeX caption)")
    ap.add_argument("--allow-uncertified", action="store_true",
                    help="draw even when lagrange_solver reports GAP_OR_COLLINEAR "
                         "(the value is then an upper bound, not the optimum)")
    ap.add_argument("--no-mute-dominant", action="store_true",
                    help="keep the largest segment a full color instead of "
                         "fading it to grey")
    args = ap.parse_args()

    instance = Path(args.instance).resolve()
    if not instance.exists():
        sys.exit(f"error: file not found: {instance}")
    if not instance.name.endswith(".fi"):
        sys.exit(f"error: expected a .fi file, got '{instance.name}'")
    if args.K < 1:
        sys.exit("error: K must be >= 1")

    conv_text, n, weights, parent, meta = run_converter(instance)
    if args.K > n - 1:
        sys.exit(f"error: K={args.K} exceeds the {n - 1} cuts available on "
                 f"this instance")

    value, cuts, status = run_dp(conv_text, args.K)
    certified = status.startswith("exact")
    if not certified and not args.allow_uncertified:
        sys.exit(f"error: lagrange_solver status is {status} for K={args.K}: the value "
                 f"{value} is an upper bound, not a certified optimum "
                 f"(--allow-uncertified to draw it anyway)")
    if len(cuts) > args.K:
        sys.exit(f"error: the DP returned {len(cuts)} cuts for a budget of {args.K}")
    if len(cuts) < args.K:
        print(f"note: the optimum uses {len(cuts)} of the {args.K} allowed cuts",
              file=sys.stderr)

    nodes, feeders = parse_fi(instance)

    children = {v: [] for v in range(1, n + 1)}
    for v in range(2, n + 1):
        children[parent[v]].append(v)

    cut_set = set(cuts)
    # the converter numbers parents before children
    comp = [0] * (n + 1)
    next_comp = 1
    for v in range(2, n + 1):
        if v in cut_set:
            comp[v] = next_comp
            next_comp += 1
        else:
            comp[v] = comp[parent[v]]
    n_comp = next_comp

    # the drawn partition must have the solver's cost
    totals = [0] * n_comp
    for v in range(1, n + 1):
        totals[comp[v]] += weights[v]
    recomputed = sum(t * t for t in totals)
    if recomputed != value:
        sys.exit(f"error: self-check failed: partition cost {recomputed} "
                 f"does not match the DP objective {value}")

    root_is_super = meta[1] == -1  # several feeders

    comp_of_orig = {}
    for v in range(1, n + 1):
        if meta[v] != -1:
            comp_of_orig[meta[v]] = comp[v]

    # each edge vertex of the tree is one .fi edge
    draw_edges = []  # (origA, origB, component)
    for v in range(1, n + 1):
        if v == 1 or meta[v] != -1:
            continue
        p = parent[v]
        kids = children[v]
        if meta[p] == -1 or len(kids) != 1 or meta[kids[0]] == -1:
            continue
        draw_edges.append((meta[p], meta[kids[0]], comp[v]))

    fi_specs = []  # (node the indicator sits next to, other end of the edge)
    structural_cuts = 0
    for c in cuts:
        if meta[c] != -1:
            vp = parent[c]
            if vp == 1 and root_is_super:
                structural_cuts += 1
                continue
            if meta[vp] != -1:
                continue
            gp = parent[vp]
            fi_specs.append((meta[c], meta[gp]))
        else:
            kids = children[c]
            if meta[parent[c]] == -1 or len(kids) != 1:
                structural_cuts += 1
                continue
            fi_specs.append((meta[parent[c]], meta[kids[0]]))

    used_orig = set(comp_of_orig)
    missing = sorted(o for o in used_orig if o not in nodes)
    if missing:
        sys.exit(f"error: {len(missing)} node(s) referenced by the network have "
                 f"no 'v' coordinate line; cannot draw")
    if len({nodes[o] for o in used_orig}) <= 1:
        sys.exit("error: instance has no usable coordinates (all nodes "
                 "coincide); only instances with real coordinates are supported")

    style = args.style
    if style == "auto":
        style = "small" if len(used_orig) <= 60 else "big"

    import matplotlib
    if not args.show:
        matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    from matplotlib.collections import LineCollection

    seg_size = [0] * n_comp
    for _a, _b, cidx in draw_edges:
        seg_size[cidx] += 1

    seg_color_map = compute_segment_colors(
        n_comp, comp, parent, n, nodes, comp_of_orig, used_orig, seg_size,
        mute_dominant=not args.no_mute_dominant)
    seg_color = lambda s: seg_color_map[s]

    # small segments on top
    draw_edges.sort(key=lambda e: -seg_size[e[2]])

    # shrink dots and lines on big networks
    if style == "small":
        node_s, edge_lw = 230.0, 2.6
    else:
        scale = (1500.0 / max(len(used_orig), 1)) ** 0.5
        node_s = max(1.0, min(14.0, 14.0 * scale))
        edge_lw = max(0.3, min(0.9, 0.9 * scale))

    fig, ax = plt.subplots(figsize=(14, 10))

    segs, colors, lws = [], [], []
    for a, b, cidx in draw_edges:
        segs.append([nodes[a], nodes[b]])
        colors.append(seg_color(cidx))
        lws.append(edge_lw)
    ax.add_collection(LineCollection(segs, colors=colors, linewidths=lws,
                                     capstyle="round", zorder=1))

    if style == "small":
        for o in used_orig:
            x, y = nodes[o]
            ax.scatter([x], [y], s=230, color=seg_color(comp_of_orig[o]),
                       edgecolors="black", linewidths=0.8, zorder=3)
            ax.text(x, y, str(o), ha="center", va="center", fontsize=6,
                    color="white", zorder=4)
    else:
        ordered = sorted(used_orig, key=lambda o: -seg_size[comp_of_orig[o]])
        xs = [nodes[o][0] for o in ordered]
        ys = [nodes[o][1] for o in ordered]
        cs = [seg_color(comp_of_orig[o]) for o in ordered]
        ax.scatter(xs, ys, s=node_s, color=cs, edgecolors="none", zorder=3)

    for fid in feeders:
        if fid in nodes:
            fx, fy = nodes[fid]
            ax.scatter([fx], [fy], s=620 if style == "small" else 90,
                       marker="s", facecolors="none", edgecolors="black",
                       linewidths=1.8 if style == "small" else 1.0, zorder=5)

    fi_size = 150 if style == "small" else 55
    t = 0.2
    for boundary, other in fi_specs:
        bx, by = nodes[boundary]
        ox, oy = nodes[other]
        mx, my = bx + t * (ox - bx), by + t * (oy - by)
        ax.scatter([mx], [my], s=fi_size, color="#39FF14",
                   edgecolors="black", linewidths=1.0, zorder=6)

    ax.set_aspect("equal")
    ax.autoscale()
    ax.margins(0.05)
    ax.axis("off")
    extra = (f"   ({structural_cuts} structural)" if structural_cuts else "")
    if not certified:
        extra += f"   [{status}: upper bound]"
    if not args.no_title:
        ax.set_title(
            f"{instance.name}   |   K={args.K} fault indicator"
            f"{'' if args.K == 1 else 's'},  {n_comp} segments{extra}\n"
            f"objective = {value}   ({value / 1e10:.6g} km²)",
            fontsize=11,
        )

    if args.out:
        out_path = Path(args.out).resolve()
    else:
        base = instance.name[:-len(".fi")]
        out_path = REPO / "out" / "viz" / f"{base}_m{args.K}_general.png"
        out_path.parent.mkdir(parents=True, exist_ok=True)
    fig.savefig(out_path, dpi=(220 if style == "big" else 150), bbox_inches="tight")
    print(f"saved {out_path}")
    print(f"objective = {value} ({value / 1e10:.6g} km^2), m = {args.K}, " + f"{n_comp} segments, status {status}")
    if args.show:
        plt.show()


if __name__ == "__main__":
    main()
