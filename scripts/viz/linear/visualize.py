#!/usr/bin/env python3
"""Draws the network with the optimal placement of m indicators on the main
trunk, one colour per section.

    python visualize.py <instance.fi|.fi.diameter> <m> [--style ...] [--out ...]
                        [--no-title] [--draw-all-edges] [--show]

Runs fi_to_diameter, convert_to_line and dp_n2k (compiled into out/bin/viz/),
then maps the cuts back to trunk edges.

Edges whose coordinates make them over 1 km and 20x longer than their stated
length are not drawn (--draw-all-edges keeps them). In our instances that is
only the two edges at the feeder of rede_mt_feeder_spd31, which would
otherwise cross the whole map.
"""
import argparse
import colorsys
import math
import os
import subprocess
import sys
from collections import deque
from pathlib import Path

SCRIPT_DIR = Path(__file__).resolve().parent
REPO = SCRIPT_DIR.parent.parent.parent  # scripts/viz/linear -> repo root
HELPERS = REPO / "src" / "helpers"
LINEAR_DP = REPO / "src" / "linear_dp"
BIN = REPO / "out" / "bin" / "viz"


def build_if_needed(src: Path, binary: Path, deps=()):
    """Compile src unless binary is newer than it and its deps."""
    if binary.exists() and all(binary.stat().st_mtime >= p.stat().st_mtime
                               for p in (src, *deps)):
        return
    binary.parent.mkdir(parents=True, exist_ok=True)
    print(f"building {binary.name} ...", file=sys.stderr)
    subprocess.run(["g++", "-O2", "-std=c++17", str(src), "-o", str(binary)],
                   check=True)


def ensure_diameter(instance: Path) -> Path:
    """The .fi.diameter of an instance, made with fi_to_diameter if missing."""
    name = instance.name
    if name.endswith(".fi.diameter"):
        return instance
    if name.endswith(".fi"):
        diam = instance.parent / (name + ".diameter")
        if not diam.exists():
            build_if_needed(HELPERS / "fi_to_diameter.cpp", BIN / "fi_to_diameter")
            print(f"generating {diam.name} ...", file=sys.stderr)
            subprocess.run([str(BIN / "fi_to_diameter"), str(instance)], check=True)
        return diam
    sys.exit(f"error: expected a .fi or .fi.diameter file, got '{name}'")


def parse_diameter_file(path: Path):
    """(nodes, edges, lengths, feeder, trunk arc indices)"""
    nodes, edges, lengths, feeder, diam_arcs = {}, [], [], None, []
    reading_diam = False
    with open(path) as f:
        for line in f:
            s = line.strip()
            if not s or s[0] in "#%":
                continue
            tok = s.split()
            key = tok[0]
            if key in ("diameter", "DIAMETER"):
                reading_diam = True
                diam_arcs.extend(int(t) for t in tok[1:])
                continue
            if reading_diam and all(t.lstrip("-").isdigit() for t in tok):
                diam_arcs.extend(int(t) for t in tok)
                continue
            reading_diam = False
            if key == "f" and len(tok) >= 2:
                feeder = int(tok[1])
            elif key == "v" and len(tok) >= 4:
                nodes[int(tok[1])] = (float(tok[2]), float(tok[3]))
            elif key == "e" and len(tok) >= 4:
                edges.append((int(tok[1]), int(tok[2])))
                lengths.append(float(tok[3]))
    if not diam_arcs:
        sys.exit(f"error: no diameter section found in {path.name}")
    return nodes, edges, lengths, feeder, diam_arcs


def section_color(s):
    """Golden-ratio hue walk, alternating lightness. Skips the greens so the
    indicator markers stand out."""
    u = (0.2 + s * 0.6180339887) % 1.0
    h = (0.46 + 0.74 * u) % 1.0
    lum = 0.42 if s % 2 == 0 else 0.58
    r, g, b = colorsys.hls_to_rgb(h, lum, 0.75)
    return (r, g, b, 1.0)


def general_section_colors(fi_path: Path, m, node_seg, n_sections):
    """Gives each section the colour that ../general/visualize.py gives the
    tree segment sharing most of its nodes, so the two figures can be put
    side by side."""
    import importlib.util
    spec = importlib.util.spec_from_file_location(
        "general_viz", SCRIPT_DIR.parent / "general" / "visualize.py")
    gv = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(gv)

    conv_text, n, _w, parent, meta = gv.run_converter(fi_path)
    _value, cuts, _status = gv.run_dp(conv_text, m)
    gnodes, _feeders = gv.parse_fi(fi_path)
    children = {v: [] for v in range(1, n + 1)}
    for v in range(2, n + 1):
        children[parent[v]].append(v)
    cut_set, comp, n_comp = set(cuts), [0] * (n + 1), 1
    for v in range(2, n + 1):
        if v in cut_set:
            comp[v], n_comp = n_comp, n_comp + 1
        else:
            comp[v] = comp[parent[v]]
    comp_of_orig = {meta[v]: comp[v] for v in range(1, n + 1) if meta[v] != -1}
    seg_size = [0] * n_comp
    for v in range(2, n + 1):
        p, kids = parent[v], children[v]
        if (meta[v] == -1 and meta[p] != -1 and len(kids) == 1
                and meta[kids[0]] != -1):
            seg_size[comp[v]] += 1
    cmap = gv.compute_segment_colors(n_comp, comp, parent, n, gnodes,
                                     comp_of_orig, set(comp_of_orig), seg_size)

    votes = [dict() for _ in range(n_sections)]
    for nid, s in node_seg.items():
        c = comp_of_orig.get(nid)
        if c is not None:
            votes[s][c] = votes[s].get(c, 0) + 1
    colors, taken = {}, set()
    for s in sorted(range(n_sections), key=lambda s: -sum(votes[s].values())):
        for c, _cnt in sorted(votes[s].items(), key=lambda kv: -kv[1]):
            if c not in taken:
                taken.add(c)
                colors[s] = cmap[c]
                break
    spare = [col for col in gv._distinct_palette(n_comp + n_sections)
             if col not in {cmap[c] for c in taken}]
    for s in range(n_sections):
        if s not in colors:
            colors[s] = spare.pop(0)
    return {s: (*colors[s][:3], 1.0) for s in colors}


def is_geographic(nodes):
    """Coordinates look like lon/lat of a small region."""
    xs = [p[0] for p in nodes.values()]
    ys = [p[1] for p in nodes.values()]
    if not xs:
        return False
    in_range = all(-180 <= x <= 180 for x in xs) and all(-90 <= y <= 90 for y in ys)
    small = (max(xs) - min(xs)) < 5 and (max(ys) - min(ys)) < 5
    fractional = any(x != int(x) for x in xs)
    return in_range and small and fractional


def inconsistent_edges(nodes, edges, lengths):
    """Edges drawn over 1 km and 20x their stated length."""
    if not is_geographic(nodes):
        return set()
    bad = set()
    for eidx, ((u, v), length_km) in enumerate(zip(edges, lengths)):
        (x1, y1), (x2, y2) = nodes[u], nodes[v]
        kx = 111.32 * math.cos(math.radians((y1 + y2) / 2))
        drawn_km = math.hypot((x2 - x1) * kx, (y2 - y1) * 110.57)
        if drawn_km > 1.0 and drawn_km > 20 * length_km:
            bad.add(eidx)
    return bad


def linearize(diam_file: Path):
    """convert_to_line with EMIT_META=1: (weights, meta)"""
    build_if_needed(HELPERS / "convert_to_line.cpp", BIN / "convert_to_line")
    text = diam_file.read_text()
    out = subprocess.run(
        [str(BIN / "convert_to_line")], input=text, capture_output=True,
        text=True, check=True, env={**os.environ, "EMIT_META": "1"},
    ).stdout
    lines = out.splitlines()
    n = int(lines[0])
    weights = [int(x) for x in lines[1].split()]
    meta = [int(x) for x in lines[2].split()]
    if not (len(weights) == n == len(meta)):
        sys.exit("error: convert_to_line output is inconsistent")
    return weights, meta


def run_dp_stdin(cleaned, K):
    """(value, cut positions) from dp_n2k with PRINT_CUTS=1"""
    build_if_needed(LINEAR_DP / "dp_n2k.cpp", BIN / "dp_n2k",
                    deps=(LINEAR_DP / "convert_utils.h",))
    stdin = f"{len(cleaned)}\n{' '.join(map(str, cleaned))}\n"
    out = subprocess.run(
        [str(BIN / "dp_n2k"), str(K)], input=stdin, capture_output=True,
        text=True, check=True, env={**os.environ, "PRINT_CUTS": "1"},
    ).stdout
    value, cuts = None, []
    for line in out.splitlines():
        if line.startswith("k = "):
            value = int(line.split(":", 1)[1])
        elif line.startswith("cuts:"):
            cuts = [int(x) for x in line[len("cuts:"):].split()]
    if value is None:
        sys.exit(f"error: could not parse dp_n2k output:\n{out}")
    return value, cuts


def dp_benchmark_value(diam_file: Path, K):
    """dp_n2k on the file itself, to cross-check run_dp_stdin."""
    try:
        out = subprocess.run(
            [str(BIN / "dp_n2k"), str(diam_file), str(K)],
            capture_output=True, text=True, check=True,
        ).stdout
    except subprocess.CalledProcessError:
        return None
    for line in out.splitlines():
        if line.startswith("k = "):
            return int(line.split(":", 1)[1])
    return None


def main():
    ap = argparse.ArgumentParser(description="Visualize FIAP optimal FI placement.")
    ap.add_argument("instance", help="path to a .fi or .fi.diameter file")
    ap.add_argument("m", type=int, help="number of fault indicators (yields m+1 sections)")
    ap.add_argument("--style", choices=["auto", "small", "big"], default="auto",
                    help="rendering style (auto: small if <=60 nodes, else big)")
    ap.add_argument("--out", default=None, help="output PNG path")
    ap.add_argument("--show", action="store_true", help="also display the figure")
    ap.add_argument("--no-title", action="store_true",
                    help="omit the in-image title (paper figures carry a caption)")
    ap.add_argument("--draw-all-edges", action="store_true",
                    help="also draw edges whose drawn length contradicts their stated length")
    ap.add_argument("--general-colors", action="store_true",
                    help="draw in the general visualizer's style and palette, giving each "
                         "section the colour of the matching general-DP segment "
                         "(for side-by-side comparison figures; needs a .fi)")
    args = ap.parse_args()

    instance = Path(args.instance).resolve()
    if not instance.exists():
        sys.exit(f"error: file not found: {instance}")
    if args.m < 0:
        sys.exit("error: m must be >= 0")
    K = args.m + 1  # sections

    diam_file = ensure_diameter(instance)
    nodes, edges, lengths, feeder, diam_arcs = parse_diameter_file(diam_file)

    edge_node_ids = {x for e in edges for x in e}
    missing = edge_node_ids - nodes.keys()
    if missing:
        sys.exit(f"error: {len(missing)} node(s) referenced by edges have no "
                 f"'v' coordinate line; cannot draw")
    if len({nodes[i] for i in nodes}) <= 1:
        sys.exit("error: instance has no usable coordinates (all nodes coincide); "
                 "only instances with real coordinates are supported")

    weights, meta = linearize(diam_file)

    # drop the zeros, remembering where each entry came from
    cleaned, orig_pos = [], []
    for i, w in enumerate(weights):
        if w != 0:
            cleaned.append(w)
            orig_pos.append(i)
    if len(cleaned) < K:
        sys.exit(f"error: m={args.m} needs {K} sections, more than the "
                 f"{len(cleaned)} non-zero segment slots available on the diameter")

    value, cuts = run_dp_stdin(cleaned, K)
    if len(cuts) != args.m:
        sys.exit(f"error: expected {args.m} cuts, dp_n2k returned {len(cuts)}")

    bench = dp_benchmark_value(diam_file, K)
    if bench is not None and bench != value:
        sys.exit(f"error: optimal value mismatch: stdin run {value} vs "
                 f"dp_n2k benchmark {bench}; the converter and convert_utils.h disagree")

    diam_arc_set = set(diam_arcs)
    cut_orig = sorted(orig_pos[c] for c in cuts)

    def seg_of(pos):
        return sum(1 for cp in cut_orig if cp < pos)

    diam_node_seg, arc_seg = {}, {}
    for q, mv in enumerate(meta):
        if q % 2 == 0:
            diam_node_seg[mv] = seg_of(q)   # trunk node
        else:
            arc_seg[mv] = seg_of(q)         # trunk arc

    # the marker goes near the end of the cut arc that starts the new section
    fi_markers = []
    for c in cuts:
        p = orig_pos[c]
        if p % 2 == 0:
            boundary_node, arc_eidx = meta[p], meta[p + 1]
        else:
            boundary_node, arc_eidx = meta[p + 1], meta[p]
        u, v = edges[arc_eidx]
        other = v if u == boundary_node else u
        bx, by = nodes[boundary_node]
        ox, oy = nodes[other]
        t = 0.2
        fi_markers.append((bx + t * (ox - bx), by + t * (oy - by)))

    # off-trunk nodes take the section of their trunk node
    adj = {nid: [] for nid in nodes}
    for eidx, (u, v) in enumerate(edges):
        adj[u].append((v, eidx))
        adj[v].append((u, eidx))
    node_seg = dict(diam_node_seg)
    dq = deque(diam_node_seg.keys())
    while dq:
        u = dq.popleft()
        for v, eidx in adj[u]:
            if eidx in diam_arc_set or v in node_seg:
                continue
            node_seg[v] = node_seg[u]
            dq.append(v)

    def edge_segment(eidx, u, v):
        if eidx in diam_arc_set:
            return arc_seg[eidx]
        return node_seg.get(u, node_seg.get(v, 0))

    style = args.style
    if style == "auto":
        style = "small" if len(nodes) <= 60 else "big"

    import matplotlib
    if not args.show:
        matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    from matplotlib.collections import LineCollection

    seg_color = section_color
    gstyle = args.general_colors
    if gstyle:
        fi_path = diam_file.parent / diam_file.name[:-len(".diameter")]
        if not fi_path.exists():
            sys.exit(f"error: --general-colors needs {fi_path.name} next to the .diameter file")
        gmap = general_section_colors(fi_path, args.m, node_seg, K)
        seg_color = lambda s: gmap[s]

    skip = set() if args.draw_all_edges else inconsistent_edges(nodes, edges, lengths)
    for eidx in sorted(skip):
        u, v = edges[eidx]
        print(f"not drawing edge {eidx} ({u}-{v}, stated {lengths[eidx]} km): its "
              f"coordinates put it far longer than stated (--draw-all-edges to keep)",
              file=sys.stderr)

    fig, ax = plt.subplots(figsize=(14, 10))

    if gstyle:
        # same look as the general figures, off-trunk parts half transparent
        if style == "small":
            node_s, edge_lw = 230.0, 2.6
        else:
            scale = (1500.0 / max(len(nodes), 1)) ** 0.5
            node_s = max(1.0, min(14.0, 14.0 * scale))
            edge_lw = max(0.3, min(0.9, 0.9 * scale))
        segs, colors = [], []
        for eidx, (u, v) in enumerate(edges):
            if eidx in skip:
                continue
            r, g, b, _ = seg_color(edge_segment(eidx, u, v))
            segs.append([nodes[u], nodes[v]])
            colors.append((r, g, b, 1.0 if eidx in diam_arc_set else 0.5))
        ax.add_collection(LineCollection(segs, colors=colors, linewidths=edge_lw,
                                         capstyle="round", zorder=1))
        drawn = [nid for nid in nodes if nid in node_seg]
        ax.scatter([nodes[i][0] for i in drawn], [nodes[i][1] for i in drawn],
                   s=node_s, edgecolors="none", zorder=3,
                   color=[(*seg_color(node_seg[i])[:3],
                           1.0 if i in diam_node_seg else 0.5) for i in drawn])
        if feeder is not None and feeder in nodes:
            fx, fy = nodes[feeder]
            ax.scatter([fx], [fy], s=620 if style == "small" else 90, marker="s",
                       facecolors="none", edgecolors="black",
                       linewidths=1.8 if style == "small" else 1.0, zorder=5)
        fi_size = 150 if style == "small" else 55
        for mx, my in fi_markers:
            ax.scatter([mx], [my], s=fi_size, color="#39FF14",
                       edgecolors="black", linewidths=1.0, zorder=6)
    else:
        draw_default(ax, style, nodes, edges, skip, diam_arc_set, diam_node_seg,
                     node_seg, edge_segment, seg_color, feeder, fi_markers,
                     LineCollection)

    ax.set_aspect("equal")
    ax.autoscale()
    ax.margins(0.05)
    ax.axis("off")
    if not args.no_title:
        ax.set_title(
            f"{instance.name}   |   m={args.m} fault indicator"
            f"{'' if args.m == 1 else 's'},  {K} sections\n"
            f"objective = {value}   ({value / 1e10:.6g} km²)",
            fontsize=11,
        )

    if args.out:
        out_path = Path(args.out).resolve()
    else:
        base = instance.name
        for suf in (".fi.diameter", ".fi"):
            if base.endswith(suf):
                base = base[:-len(suf)]
                break
        out_path = REPO / "out" / "viz" / f"{base}_m{args.m}_linear.png"
        out_path.parent.mkdir(parents=True, exist_ok=True)
    fig.savefig(out_path, dpi=220 if (gstyle and style == "big") else 150,
                bbox_inches="tight")
    print(f"saved {out_path}")
    print(f"objective = {value} ({value / 1e10:.6g} km^2), m = {args.m}, " + f"{K} sections")
    if args.show:
        plt.show()


def draw_default(ax, style, nodes, edges, skip, diam_arc_set, diam_node_seg,
                 node_seg, edge_segment, seg_color, feeder, fi_markers,
                 LineCollection):
    """Thick trunk, thin translucent branches, gold star on the feeder."""
    segs, colors, lws = [], [], []
    for eidx, (u, v) in enumerate(edges):
        if eidx in skip:
            continue
        on_diam = eidx in diam_arc_set
        r, g, b, _ = seg_color(edge_segment(eidx, u, v))
        segs.append([nodes[u], nodes[v]])
        colors.append((r, g, b, 1.0 if on_diam else 0.5))
        if style == "small":
            lws.append(5.5 if on_diam else 1.6)
        else:
            lws.append(3.0 if on_diam else 0.9)
    ax.add_collection(LineCollection(segs, colors=colors, linewidths=lws,
                                     capstyle="round", zorder=1))

    if style == "small":
        for nid, (x, y) in nodes.items():
            on_diam = nid in diam_node_seg
            ax.scatter([x], [y], s=430 if on_diam else 200,
                       color=seg_color(node_seg.get(nid, 0)),
                       edgecolors="black", linewidths=2.4 if on_diam else 0.6,
                       zorder=3)
            ax.text(x, y, str(nid), ha="center", va="center",
                    fontsize=8 if on_diam else 6,
                    fontweight="bold" if on_diam else "normal",
                    color="white", zorder=4)
    else:
        # no labels in the big style, so mark the trunk nodes
        dx = [nodes[nid][0] for nid in diam_node_seg]
        dy = [nodes[nid][1] for nid in diam_node_seg]
        dc = [seg_color(diam_node_seg[nid]) for nid in diam_node_seg]
        ax.scatter(dx, dy, s=20, color=dc, edgecolors="black",
                   linewidths=0.5, zorder=3)

    if feeder is not None and feeder in nodes:
        fx, fy = nodes[feeder]
        ax.scatter([fx], [fy], s=900 if style == "small" else 360,
                   marker="*", color="#FFD700",
                   edgecolors="black", linewidths=1.4, zorder=7)

    fi_size = 150 if style == "small" else 60
    for mx, my in fi_markers:
        ax.scatter([mx], [my], s=fi_size, color="#39FF14",
                   edgecolors="black", linewidths=1.0, zorder=6)


if __name__ == "__main__":
    main()
