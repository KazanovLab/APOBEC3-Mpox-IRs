"""
Test 10: Fixed-stem, fixed-loop heatmaps (mismatches = 0).

Layout: 5-row × 4-column grid of heatmaps.
  rows    → selection_type: most_stable, greedy, max_cov, min_cov, all
  columns → structure_type: hairpin, spacer, ct_end, c_end

Each heatmap:
  X axis: loop length in [2, 3, 4, 5, 6, 7, 8, 9, 10]  (exact spacer length)
  Y axis: stem length in [4, 5, 6, 7, 8, 9, 10, 11, 12]  (stem_min == stem_max)
  Cell color: −log10(p-value)
  Cell text:  "hits/targets" on line 1, p-value on line 2

EMBOSS is run with gaplimit=loop_len (max spacer), then hairpins are filtered
to retain only those with:
  1. spacer_length == loop_len  (exact loop length)
  2. stem is non-extendable: flanking nucleotides outside the stem are not
     Watson-Crick complementary (i.e. the stem cannot be extended further).

Saves: one combined PDF/PNG + one PDF/PNG per (selection_type, structure_type) pair.
"""

import sys
import os
import io
import contextlib
import math
import re
import pickle
import argparse
from itertools import product
from pathlib import Path

import numpy as np
import matplotlib.pyplot as plt
from matplotlib.colors import LinearSegmentedColormap, TwoSlopeNorm

# ── paths ─────────────────────────────────────────────────────────────────────
BINOMIAL_DIR = Path(__file__).parent.parent.parent / "binomial"
sys.path.insert(0, str(BINOMIAL_DIR))
os.chdir(BINOMIAL_DIR)

import load
from emboss import get_palindrome
from hairpin_groups import calculate_pval_groups
import all_hairpins as ah

# ── command-line arguments ─────────────────────────────────────────────────────
parser = argparse.ArgumentParser(description="Test 10: fixed-stem/loop heatmaps.")
parser.add_argument("--mm", type=int, default=0,
                    help="Number of mismatches allowed in stem (default: 0)")
args = parser.parse_args()
MM = args.mm
print(f"Running with mismatches = {MM}")

# ── grid parameters ────────────────────────────────────────────────────────────
LOOP_VALUES     = [2, 3, 4, 5, 6, 7, 8, 9, 10]
STEM_VALUES     = [4, 5, 6, 7, 8, 9, 10, 11, 12]

STRUCTURE_TYPES = ["hairpin", "spacer", "ct_end", "c_end"]
SELECTION_TYPES = ["most_stable", "greedy", "max_cov", "min_cov", "all"]

# ── parse captured stdout ──────────────────────────────────────────────────────
def parse_output(text):
    hits              = 0
    targets_in_struct = 0
    targets_p         = 0.0
    pvalue            = 1.0
    m = re.search(r"Mutations in structures:\s*(\d+)", text)
    if m:
        hits = int(m.group(1))
    m = re.search(r"Targets in structures:\s*(\d+)", text)
    if m:
        targets_in_struct = int(m.group(1))
    m = re.search(r"Fraction of targets in structures:\s*([\d.eE+\-]+)%", text)
    if m:
        targets_p = float(m.group(1)) / 100.0
    m = re.search(r"p-value:\s*([\d.eE+\-]+)%", text)
    if m:
        pvalue = float(m.group(1)) / 100.0
    return hits, targets_in_struct, targets_p, pvalue

# ── maximal-stem filter ────────────────────────────────────────────────────────
_COMPLEMENT = {'A': 'T', 'T': 'A', 'G': 'C', 'C': 'G'}
_genome_seq  = load.genome.sequence

# genome-wide denominator sizes (used for panel-level combined stats)
_GENOME_DENOM = {
    "hairpin": len(re.findall("TC", _genome_seq)) + len(re.findall("GA", _genome_seq)),
    "spacer":  len(re.findall("TC", _genome_seq)) + len(re.findall("GA", _genome_seq)),
    "ct_end":  len(re.findall("TC", _genome_seq)) + len(re.findall("GA", _genome_seq)),
    "c_end":   len(re.findall("C",  _genome_seq)) + len(re.findall("G",  _genome_seq)),
}
_MUT_CNT = load.mut_cnt

def is_extendable(h):
    """Return True if the nucleotides flanking the stem are Watson-Crick
    complementary, meaning the stem could be extended by one more base pair."""
    if h.start == 0 or h.end >= len(_genome_seq) - 1:
        return False
    left  = _genome_seq[h.start - 1].upper()
    right = _genome_seq[h.end + 1].upper()
    return _COMPLEMENT.get(left) == right

# ── pre-compute hairpin lists (one per stem/loop combo) ───────────────────────
print("Pre-computing hairpin lists...")
hairpin_cache = {}
for stem_len in STEM_VALUES:
    for loop_len in LOOP_VALUES:
        pa_file      = get_palindrome(stem_len, stem_len, loop_len, MM)
        all_hairpins = load.load_hairpins(str(pa_file))
        # keep only hairpins with exact spacer length AND non-extendable stem
        hairpins = [h for h in all_hairpins
                    if h.spacer_length == loop_len and not is_extendable(h)]
        hairpin_cache[(stem_len, loop_len)] = hairpins
        print(f"  stem={stem_len}, loop={loop_len}: {len(hairpins)}/{len(all_hairpins)} hairpins (exact loop, non-extendable)")

# ── compute all results ────────────────────────────────────────────────────────
# key: (sel_type, struct_type, stem_len, loop_len)
# val: (log10p, hits, targets_in_struct, targets_p, pvalue)
results = {}

for sel_type in SELECTION_TYPES:
    for struct_type in STRUCTURE_TYPES:
        print(f"\n[{sel_type} / {struct_type}]")
        for stem_len in STEM_VALUES:
            for loop_len in LOOP_VALUES:
                hairpins = hairpin_cache[(stem_len, loop_len)]
                if not hairpins:
                    results[(sel_type, struct_type, stem_len, loop_len)] = (0.0, 0, 0, 0.0, 1.0)
                    continue
                buf = io.StringIO()
                with contextlib.redirect_stdout(buf):
                    if sel_type == "all":
                        log10p = ah.calculate_pval_all(hairpins, struct_type)
                    else:
                        log10p = calculate_pval_groups(hairpins, sel_type, struct_type)
                captured = buf.getvalue()
                hits, targets_in_struct, targets_p, pvalue = parse_output(captured)
                results[(sel_type, struct_type, stem_len, loop_len)] = (
                    log10p, hits, targets_in_struct, targets_p, pvalue
                )
                print(f"  stem={stem_len} loop={loop_len}: "
                      f"hits={hits}, tgt={targets_in_struct}, "
                      f"p={pvalue:.5f}, -log10p={log10p:.3f}")

# ── panel-level combined stats ────────────────────────────────────────────────
from scipy.stats import binomtest as _binomtest
panel_stats = {}
for _sel in SELECTION_TYPES:
    for _str in STRUCTURE_TYPES:
        _non_empty = [
            results[(_sel, _str, sl, ll)]
            for sl in STEM_VALUES for ll in LOOP_VALUES
            if results[(_sel, _str, sl, ll)][2] > 0
        ]
        _hits       = sum(v[1] for v in _non_empty)
        _targets    = sum(v[2] for v in _non_empty)
        _gdenom     = _GENOME_DENOM[_str]
        _p0         = min(_targets / _gdenom, 1.0) if _gdenom > 0 else 0.0
        _k          = min(_hits, _MUT_CNT)
        if _p0 > 0:
            _pv = _binomtest(_k, _MUT_CNT, _p0, alternative='greater').pvalue
        else:
            _pv = 1.0
        _lp = -math.log10(_pv) if _pv > 0 else float('inf')
        panel_stats[(_sel, _str)] = (_hits, _targets, _p0, _lp)

# ── save results for standalone plotting ──────────────────────────────────────
out_dir = Path(__file__).parent
with open(out_dir / f"results_mm{MM}.pkl", "wb") as f:
    pickle.dump({
        "results":         results,
        "panel_stats":     panel_stats,
        "LOOP_VALUES":     LOOP_VALUES,
        "STEM_VALUES":     STEM_VALUES,
        "STRUCTURE_TYPES": STRUCTURE_TYPES,
        "SELECTION_TYPES": SELECTION_TYPES,
    }, f)
print(f"Results saved to results_mm{MM}.pkl")

# ── colormap (shared across all panels) ───────────────────────────────────────
SIG_THRESH = 1.3  # −log10(0.05)

finite_vals = [v[0] for v in results.values() if not math.isinf(v[0])]
vmin = 0.0
vmax = math.ceil(max(finite_vals)) if finite_vals else 5

n_nonsig      = 128
n_sig         = 128
colors_nonsig = plt.cm.Greys(np.linspace(0.0, 0.55, n_nonsig))
colors_sig    = plt.cm.YlOrRd(np.linspace(0.2,  1.0, n_sig))
cmap_custom   = LinearSegmentedColormap.from_list(
    "sig_cmap", np.vstack([colors_nonsig, colors_sig])
)
norm = TwoSlopeNorm(vmin=vmin, vcenter=SIG_THRESH, vmax=vmax)

x_labels = [str(v) for v in LOOP_VALUES]
y_labels  = [str(v) for v in STEM_VALUES]

# ── helper: draw one heatmap into an axes ──────────────────────────────────────
def draw_panel(ax, fig, sel_type, struct_type, fontsize_tick=7, fontsize_ann=5):
    n_rows = len(STEM_VALUES)
    n_cols = len(LOOP_VALUES)

    mat = np.zeros((n_rows, n_cols))
    for r, stem_len in enumerate(STEM_VALUES):
        for c, loop_len in enumerate(LOOP_VALUES):
            log10p = results[(sel_type, struct_type, stem_len, loop_len)][0]
            mat[r, c] = min(log10p, vmax)

    im = ax.imshow(mat, aspect="auto", origin="lower", cmap=cmap_custom, norm=norm)

    ax.set_xticks(range(n_cols))
    ax.set_xticklabels(x_labels, fontsize=fontsize_tick)
    ax.set_yticks(range(n_rows))
    ax.set_yticklabels(y_labels, fontsize=fontsize_tick)
    ax.set_xlabel("Loop length (nt)", fontsize=fontsize_tick)
    ax.set_ylabel("Stem length (bp)", fontsize=fontsize_tick)
    ph, pt, pp0, plp = panel_stats[(sel_type, struct_type)]
    plp_str = f"{plp:.2f}" if not math.isinf(plp) else "inf"
    ax.set_title(
        f"{struct_type} / {sel_type}\n"
        f"All: {ph}/{pt}  p\u2080={pp0:.3g}  \u2212log\u2081\u2080p={plp_str}",
        fontsize=fontsize_tick + 1, fontweight="bold"
    )

    for r, stem_len in enumerate(STEM_VALUES):
        for c, loop_len in enumerate(LOOP_VALUES):
            log10p, hits, targets_in_struct, targets_p, pvalue = results[
                (sel_type, struct_type, stem_len, loop_len)
            ]
            log10p_str = f"{log10p:.2f}" if not math.isinf(log10p) else "inf"
            ann = f"{hits}/{targets_in_struct}\n{log10p_str}"
            text_color = (
                "white"
                if log10p > SIG_THRESH + (vmax - SIG_THRESH) * 0.5
                else "black"
            )
            ax.text(
                c, r, ann,
                ha="center", va="center",
                fontsize=fontsize_ann, fontweight="bold",
                color=text_color, linespacing=1.2,
            )

    cbar = fig.colorbar(im, ax=ax, fraction=0.046, pad=0.04)
    cbar.set_label("−log₁₀(p)", fontsize=fontsize_tick - 1)
    cbar.ax.tick_params(labelsize=fontsize_tick - 1)
    cbar.ax.axhline(y=SIG_THRESH, color="black", linewidth=1.0, linestyle="--")
    cbar.ax.text(2.6, SIG_THRESH, "p=0.05", va="center",
                 fontsize=fontsize_tick - 2, color="black")

# ── combined 5×4 figure ────────────────────────────────────────────────────────
fig, axes = plt.subplots(
    len(SELECTION_TYPES), len(STRUCTURE_TYPES),
    figsize=(18, 22),
)
for row_i, sel_type in enumerate(SELECTION_TYPES):
    for col_i, struct_type in enumerate(STRUCTURE_TYPES):
        draw_panel(axes[row_i, col_i], fig, sel_type, struct_type)

fig.suptitle(
    f"Binomial test significance  (mm={MM}, fixed stem & loop length)\n"
    "Cell: hits/targets  |  −log₁₀(p)      Color: −log₁₀(p-value)",
    fontsize=12, y=1.002,
)
plt.tight_layout()
plt.savefig(out_dir / f"heatmap_all_mm{MM}.pdf", dpi=150, bbox_inches="tight")
plt.savefig(out_dir / f"heatmap_all_mm{MM}.png", dpi=150, bbox_inches="tight")
plt.close()
print(f"\nSaved: heatmap_all_mm{MM}.pdf / .png")

# ── individual figures (one per panel) ────────────────────────────────────────
for sel_type in SELECTION_TYPES:
    for struct_type in STRUCTURE_TYPES:
        fig, ax = plt.subplots(figsize=(9, 7))
        draw_panel(ax, fig, sel_type, struct_type, fontsize_tick=9, fontsize_ann=7)
        fig.suptitle(
            f"Binomial test: {struct_type} / {sel_type}  (mm={MM}, fixed stem & loop)\n"
            "Cell: hits/targets  |  −log₁₀(p)      Color: −log₁₀(p-value)",
            fontsize=10, y=1.01,
        )
        plt.tight_layout()
        fname = f"heatmap_mm{MM}_{sel_type}_{struct_type}"
        plt.savefig(out_dir / f"{fname}.pdf", dpi=150, bbox_inches="tight")
        plt.savefig(out_dir / f"{fname}.png", dpi=150, bbox_inches="tight")
        plt.close()
        print(f"Saved: {fname}.pdf / .png")
