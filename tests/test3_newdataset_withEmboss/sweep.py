#!/usr/bin/env python3
"""
Sweep all combinations of structure_type and hairpin_selection_type by
importing the binomial modules directly (no subprocess / no config.yaml rewrites).

Run from any directory:
    python /path/to/tests/test3/sweep.py
"""
import sys
import os
import io
import math
import contextlib
import re
from itertools import product
from pathlib import Path

BINOMIAL_DIR = Path(__file__).parent.parent.parent / "binomial"
OUT_PATH     = Path(__file__).parent / "results.txt"

# ── set up path and working directory ─────────────────────────────────────────
# load.py uses a relative path for the genome file, so we must be in binomial/
os.chdir(BINOMIAL_DIR)
sys.path.insert(0, str(BINOMIAL_DIR))

# ── import modules (genome, mutations, hairpins are loaded once here) ─────────
import load
import functions
import all_hairpins as ah
import hairpin_groups as hg

print(f"Genome length : {load.genome.length:,} bp")
print(f"Mutations     : {load.mut_cnt:,}")
print(f"Hairpins      : {len(load.all_hairpins_list):,}")
print()

# load.py exports the hairpin list as 'all_hairpins_list', but the updated
# calculate_pval_* functions reference it internally as 'hairpins_list'.
# Inject the name into both module namespaces so the functions can find it.
ah.hairpins_list = load.all_hairpins_list
hg.hairpins_list = load.all_hairpins_list

# ── helper: parse captured stdout from calculate_pval_* functions ──────────────
def parse_output(text):
    """Extract hits, frac_pct, pvalue from printed function output."""
    hits    = float("nan")
    frac    = float("nan")
    pvalue  = float("nan")
    m = re.search(r"Mutations in structures:\s*(\d+)", text)
    if m:
        hits = int(m.group(1))
    m = re.search(r"Fraction of targets in structures:\s*([\d.]+)%", text)
    if m:
        frac = float(m.group(1))
    # printed as pvalue*100 with %
    m = re.search(r"p-value:\s*([\d.eE+\-]+)%", text)
    if m:
        pvalue = float(m.group(1)) / 100.0
    return hits, frac, pvalue

# ── sweep parameters ──────────────────────────────────────────────────────────
STRUCTURE_TYPES = ["hairpin", "spacer", "ct_end", "c_end"]
SELECTION_TYPES = ["most_stable", "greedy", "max_cov", "min_cov", "all"]

# ── run all combinations ───────────────────────────────────────────────────────
results = []
total = len(STRUCTURE_TYPES) * len(SELECTION_TYPES)
done  = 0

for structure_type, selection_type in product(STRUCTURE_TYPES, SELECTION_TYPES):
    done += 1
    print(f"[{done:>2}/{total}] structure={structure_type:<8}  selection={selection_type}", flush=True)

    # patch the module-level hit_type that hit_or_not() reads at call time
    functions.hit_type = structure_type

    buf = io.StringIO()
    with contextlib.redirect_stdout(buf):
        if selection_type == "all":
            log10p = ah.calculate_pval_all(load.all_hairpins_list, structure_type)
        else:
            log10p = hg.calculate_pval_groups(load.all_hairpins_list, selection_type, structure_type)
    captured = buf.getvalue()
    print(captured, end="")   # echo to terminal so user can follow progress

    hits, frac, pvalue = parse_output(captured)
    results.append(dict(
        structure_type=structure_type,
        selection_type=selection_type,
        mut_cnt=load.mut_cnt,
        hits=hits,
        frac_pct=frac,
        pvalue=pvalue,
        log10p=log10p,
    ))
    log10p_str = "inf" if math.isinf(log10p) else f"{log10p:.4f}"
    print(f"  → -log10(p) = {log10p_str}\n")

# ── format table ──────────────────────────────────────────────────────────────
col_w = dict(structure_type=14, selection_type=12, mut_cnt=8, hits=6,
             frac_pct=14, pvalue=14, log10p=14)

header = (f"{'structure_type':<{col_w['structure_type']}}  "
          f"{'selection_type':<{col_w['selection_type']}}  "
          f"{'mut_cnt':>{col_w['mut_cnt']}}  "
          f"{'hits':>{col_w['hits']}}  "
          f"{'frac_targets_%':>{col_w['frac_pct']}}  "
          f"{'p_value':>{col_w['pvalue']}}  "
          f"{'-log10(p)':>{col_w['log10p']}}")
sep = "-" * len(header)

lines = [header, sep]
for r in results:
    log10p_str = ("inf" if math.isinf(r["log10p"]) else
                  "NA"  if math.isnan(r["log10p"]) else
                  f"{r['log10p']:.4f}")
    pval_str   = ("NA" if math.isnan(r["pvalue"]) else
                  f"{r['pvalue']:.6f}" if r["pvalue"] >= 1e-3 else
                  f"{r['pvalue']:.6e}")
    hits_str   = "ERR" if (r["hits"] is None or (isinstance(r["hits"], float) and math.isnan(r["hits"]))) else str(int(r["hits"]))
    frac_str   = "NA"  if math.isnan(r["frac_pct"]) else f"{r['frac_pct']:.4f}"

    lines.append(
        f"{r['structure_type']:<{col_w['structure_type']}}  "
        f"{r['selection_type']:<{col_w['selection_type']}}  "
        f"{r['mut_cnt']:>{col_w['mut_cnt']}}  "
        f"{hits_str:>{col_w['hits']}}  "
        f"{frac_str:>{col_w['frac_pct']}}  "
        f"{pval_str:>{col_w['pvalue']}}  "
        f"{log10p_str:>{col_w['log10p']}}"
    )
    # blank line between structure type blocks
    if r is results[-1] or results[results.index(r)+1]["structure_type"] != r["structure_type"]:
        lines.append("")

table_str = "\n".join(lines)
print("\n" + table_str)

with open(OUT_PATH, "w") as f:
    f.write(table_str + "\n")
print(f"Results saved → {OUT_PATH}")
