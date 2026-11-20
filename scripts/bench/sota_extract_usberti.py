#!/usr/bin/env python3
"""Pulls the published results of Usberti et al. (2024) out of the authors'
spreadsheet (Comparacao_Flows_vs_EdgePartition_vs_Literature.xlsx, not
included). Produces data/results/sota_usberti.csv.

    python3 -I scripts/bench/sota_extract_usberti.py <results.xlsx> [out.csv]

Their runs: Gurobi 9.5.1 and BRKGA on an i7-3930K with 32 GB, 1 h per run.
Costs are in hours of fault location time, sum of squared section lengths
over (10 km/h x network length).

Columns: instance,m,lm_ub_h,lm_lb_h,lm_time_s,cepm_ub_h,cepm_lb_h,cepm_time_s,
brkga_ub_h,brkga_time_s,article_gap_pct
LM is empty for 1080nodes (out of memory). bus_13_3 has no m = 20, 25.
article_gap_pct is the gap of the article's Table 1.
"""

import csv
import os
import re
import sys
import xml.etree.ElementTree as ET
import zipfile

REPO = os.path.dirname(os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
NS = {"m": "http://schemas.openxmlformats.org/spreadsheetml/2006/main"}
REL = "{http://schemas.openxmlformats.org/officeDocument/2006/relationships}id"
SHEET = "Resultados_Artigo_Editado"
GAP_SHEET, GAP_COL = "Resultados_Artigo_Editado_Horar", "O"
NAMES = {n: n for n in ("bus_13_3", "bus_29_1", "bus_32_1", "bus_83_11", "bus_135_8",
                        "bus_201_3", "bus_873_7", "85nodes", "137nodes", "417nodes",
                        "1080nodes")}
# the spreadsheet names the Jaboticabal feeders by size
NAMES.update({"jab_222nodes": "JaboticabalJAB05", "jab_417nodes": "JaboticabalATH05",
              "jab_446nodes": "JaboticabalATH04", "jab_553nodes": "JaboticabalATH03",
              "jab_598nodes": "JaboticabalJAB03", "jab_620nodes": "JaboticabalJAB04",
              "jab_647nodes": "JaboticabalATH06", "jab_678nodes": "JaboticabalJAB02"})
EXPECTED_ROWS = 18 * 9 + 7
COLS = {"C": "lm_ub_h", "H": "lm_lb_h", "O": "lm_time_s",
        "P": "cepm_ub_h", "U": "cepm_lb_h", "AB": "cepm_time_s",
        "AC": "brkga_ub_h", "AJ": "brkga_time_s"}


def sheet_rows(path, sheet):
    z = zipfile.ZipFile(path)
    strings = ["".join(t.text or "" for t in si.iter("{%s}t" % NS["m"]))
               for si in ET.fromstring(z.read("xl/sharedStrings.xml")).findall("m:si", NS)]
    rels = {r.get("Id"): r.get("Target")
            for r in ET.fromstring(z.read("xl/_rels/workbook.xml.rels"))}
    target = None
    for s in ET.fromstring(z.read("xl/workbook.xml")).find("m:sheets", NS):
        if s.get("name") == sheet:
            target = rels[s.get(REL)]
    if target is None:
        sys.exit(f"sheet {sheet} not found")
    root = ET.fromstring(z.read("xl/" + target))
    for row in root.find("m:sheetData", NS).findall("m:row", NS):
        cells = {}
        for c in row.findall("m:c", NS):
            v = c.find("m:v", NS)
            if v is None:
                continue
            val = strings[int(v.text)] if c.get("t") == "s" else v.text
            cells[re.sub(r"\d+", "", c.get("r"))] = val
        yield cells


def number(s):
    try:
        return repr(float(s))
    except (TypeError, ValueError):
        return ""


def main():
    if len(sys.argv) < 2 or sys.argv[1] in ("-h", "--help"):
        print(__doc__)
        return
    out = sys.argv[2] if len(sys.argv) > 2 else os.path.join(
        REPO, "data", "results", "sota_usberti.csv")
    gaps = {}
    for cells in sheet_rows(sys.argv[1], GAP_SHEET):
        if cells.get("A") in NAMES:
            gaps[(NAMES[cells["A"]], int(float(cells["B"])))] = number(cells.get(GAP_COL))
    rows = []
    for cells in sheet_rows(sys.argv[1], SHEET):
        inst = NAMES.get(cells.get("A"))
        if inst is None:
            continue
        m = int(float(cells["B"]))
        rows.append([inst, m] + [number(cells.get(c)) for c in COLS] + [gaps[(inst, m)]])
    order = list(NAMES.values())
    rows.sort(key=lambda r: (order.index(r[0]), r[1]))
    if len(rows) != EXPECTED_ROWS:
        sys.exit(f"expected {EXPECTED_ROWS} rows, found {len(rows)}")
    with open(out, "w", newline="") as fh:
        w = csv.writer(fh)
        w.writerow(["instance", "m"] + list(COLS.values()) + ["article_gap_pct"])
        w.writerows(rows)
    print("wrote", out, file=sys.stderr)


if __name__ == "__main__":
    main()
