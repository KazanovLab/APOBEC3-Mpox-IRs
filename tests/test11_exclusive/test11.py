"""
Test 11: Exclusive counting — fixed stem & loop, four mutually non-overlapping categories.

Structure types (exclusive, no overlap):
  tc_end          — TC/GA at spacer boundary (C of TC at right boundary r,
                    G of GA at left boundary l)
  c_end_excl      — C/G at spacer boundary NOT in TC/GA context
  spacer_interior — TC/GA in interior spacer positions [l+1, r-1]
  stem            — TC/GA in stem arm positions only

Denominators:
  tc_end, spacer_interior, stem  → total TC+GA positions in genome
  c_end_excl                     → total (C not in TC) + (G not in GA) in genome

Layout: 5-row × 4-column grid of heatmaps.
  rows    → selection_type: most_stable, greedy, max_cov, min_cov, all
  columns → structure_type: tc_end, c_end_excl, spacer_interior, stem

Each heatmap:
  X axis: loop length in [2, 3, 4, 5, 6, 7, 8, 9, 10]  (exact spacer length)
  Y axis: stem length in [4, 5, 6, 7, 8, 9, 10, 11, 12] (stem_min == stem_max)
  Cell color: −log10(p-value)
  Cell text:  "hits/targets" on line 1, −log10(p) on line 2
"""

import sys
import os
import re
import math
import pickle
import argparse
from pathlib import Path
from copy import deepcopy

import numpy as np
import matplotlib.pyplot as plt
from matplotlib.colors import LinearSegmentedColormap, TwoSlopeNorm
from scipy.stats import binomtest

# ── paths ─────────────────────────────────────────────────────────────────────
BINOMIAL_DIR = Path(__file__).parent.parent.parent / "binomial"
sys.path.insert(0, str(BINOMIAL_DIR))
os.chdir(BINOMIAL_DIR)

import load
from emboss import get_palindrome
from functions import find_groups, max_coverage, max_coverage_spacer

# ── command-line arguments ─────────────────────────────────────────────────────
parser = argparse.ArgumentParser(description="Test 11: exclusive structure-type heatmaps.")
parser.add_argument("--mm", type=int, default=0,
                    help="Number of mismatches allowed in stem (default: 0)")
args = parser.parse_args()
MM = args.mm
print(f"Running with mismatches = {MM}")

# ── grid parameters ────────────────────────────────────────────────────────────
LOOP_VALUES     = [2, 3, 4, 5, 6, 7, 8, 9, 10]
STEM_VALUES     = [4, 5, 6, 7, 8, 9, 10, 11, 12]

STRUCTURE_TYPES = ["tc_end", "c_end_excl", "spacer_interior", "stem"]
SELECTION_TYPES = ["most_stable", "greedy", "max_cov", "min_cov", "all"]

# ── genome / mutation data (loaded via load module) ───────────────────────────
_genome_seq   = load.genome.sequence
_mut_list     = load.mutations_list   # 0-based, with multiplicity
_mut_cnt      = load.mut_cnt

# ── pre-compute genome-wide target sets (done once) ───────────────────────────
print("Pre-computing genome-wide target sets...")
_tc_pos    = set(m.end() - 1 for m in re.finditer("TC", _genome_seq))   # C of TC
_ga_pos    = set(m.start()   for m in re.finditer("GA", _genome_seq))   # G of GA
_all_c     = set(m.start()   for m in re.finditer("C",  _genome_seq))
_all_g     = set(m.start()   for m in re.finditer("G",  _genome_seq))
_c_not_tc  = _all_c - _tc_pos   # C not preceded by T
_g_not_ga  = _all_g - _ga_pos   # G not followed by A

# denominator sizes (constant across all cells)
_denom_tc_ga      = len(_tc_pos | _ga_pos)
_denom_c_not_tc   = len(_c_not_tc | _g_not_ga)
print(f"  TC+GA positions in genome    : {_denom_tc_ga}")
print(f"  C_not_TC+G_not_GA in genome  : {_denom_c_not_tc}")

# ── maximal-stem filter ────────────────────────────────────────────────────────
_COMPLEMENT = {'A': 'T', 'T': 'A', 'G': 'C', 'C': 'G'}

def is_extendable(h):
    """Return True if flanking nucleotides outside the stem are WC complementary."""
    if h.start == 0 or h.end >= len(_genome_seq) - 1:
        return False
    left  = _genome_seq[h.start - 1].upper()
    right = _genome_seq[h.end + 1].upper()
    return _COMPLEMENT.get(left) == right

# ── selection helpers ──────────────────────────────────────────────────────────
_priority_is_max = False  # matches load.priority_is_max

def _most_stable(group):
    k = -1 if _priority_is_max else 1
    return [min(group, key=lambda x: (x.score * k, k * x.length, x.start))]

def _min_coverage(group):
    return [min(group, key=lambda x: (x.length, x.score, x.start))]

def _greedy_choose(hairpins_group, taken=None, skip_sort=False):
    if taken is None:
        taken = []
    hairpins_group_copy = deepcopy(hairpins_group)
    taken_copy = taken.copy()

    if not skip_sort:
        k = -1 if _priority_is_max else 1
        hairpins_group_copy.sort(key=lambda x: x.score * k)
        n = len(hairpins_group_copy)
        i = 0
        equal_groups = []
        while i < n:
            current_group = [hairpins_group_copy[i]]
            j = i + 1
            while j < n and hairpins_group_copy[i].score == hairpins_group_copy[j].score:
                current_group.append(hairpins_group_copy[j])
                j += 1
            current_group.sort(key=lambda x: x.length * -1)
            equal_groups.append(current_group)
            i = j
        hairpins_group_copy = equal_groups

    for p in range(len(hairpins_group_copy)):
        equal_group = hairpins_group_copy[p]
        if equal_group == []:
            continue
        usable = [h for h in equal_group
                  if all(h.can_exist(j) for j in taken_copy)]
        n = len(usable)
        if n == 0:
            continue
        elif n == 1:
            taken_copy.append(usable[0])
        elif n <= 11:
            variants = []
            for hairpin in usable:
                new_hg = deepcopy(hairpins_group_copy)
                new_hg[p] = [h for h in new_hg[p] if h is not hairpin]
                temp = _greedy_choose(new_hg[p:], taken_copy + [hairpin], True)
                variants.append(temp)
            return max(variants, key=lambda x: len(x))
        else:
            hairpin = usable[0]
            new_hg = deepcopy(hairpins_group_copy)
            new_hg[p] = [h for h in new_hg[p] if h is not hairpin]
            return _greedy_choose(new_hg[p:], taken_copy + [hairpin], True)

    return taken_copy


def select_hairpins(hairpins, sel_type):
    """Return a list of (non-overlapping) hairpins according to sel_type."""
    if sel_type == "all":
        return list(hairpins)

    groups = find_groups(list(hairpins))
    selected = []
    for group in groups:
        if sel_type == "most_stable":
            choice = _most_stable(group)
        elif sel_type == "greedy":
            choice = _greedy_choose(group)
        elif sel_type == "max_cov":
            choice = max_coverage(group)
        elif sel_type == "min_cov":
            choice = _min_coverage(group)
        else:
            choice = _most_stable(group)
        selected.extend(choice)
    return selected

# ── exclusive target computation ───────────────────────────────────────────────
def compute_exclusive(hairpins, struct_type):
    """
    Return (hits, n_struct_targets, targets_p) for the given exclusive struct_type.

    hits        — number of mutation observations (with multiplicity) in struct_targets
    n_struct    — number of distinct genomic positions that are structure targets
    targets_p   — n_struct / genome_denominator
    """
    if struct_type == "tc_end":
        struct_targets = set()
        for h in hairpins:
            if h.spacer_length < 2:
                continue
            l, r = h.spacer_index
            if r in _tc_pos:
                struct_targets.add(r)
            if l in _ga_pos:
                struct_targets.add(l)
        targets_p = len(struct_targets) / _denom_tc_ga if _denom_tc_ga else 0.0

    elif struct_type == "c_end_excl":
        struct_targets = set()
        for h in hairpins:
            if h.spacer_length < 1:
                continue
            l, r = h.spacer_index
            if r in _c_not_tc:
                struct_targets.add(r)
            if l in _g_not_ga:
                struct_targets.add(l)
        targets_p = len(struct_targets) / _denom_c_not_tc if _denom_c_not_tc else 0.0

    elif struct_type == "spacer_interior":
        struct_targets = set()
        for h in hairpins:
            l, r = h.spacer_index
            for pos in range(l + 1, r):   # interior: [l+1, r-1]
                if pos in _tc_pos or pos in _ga_pos:
                    struct_targets.add(pos)
        targets_p = len(struct_targets) / _denom_tc_ga if _denom_tc_ga else 0.0

    elif struct_type == "stem":
        struct_targets = set()
        for h in hairpins:
            ls, le = h.stem_indexes[0]      # left arm: [start, start+stem_len-1]
            rs, re_ = h.stem_indexes[1]     # right arm: [end-stem_len+1, end]
            for pos in range(ls, le + 1):
                if pos in _tc_pos or pos in _ga_pos:
                    struct_targets.add(pos)
            for pos in range(rs, re_ + 1):
                if pos in _tc_pos or pos in _ga_pos:
                    struct_targets.add(pos)
        targets_p = len(struct_targets) / _denom_tc_ga if _denom_tc_ga else 0.0

    else:
        raise ValueError(f"Unknown struct_type: {struct_type}")

    hits = sum(1 for m in _mut_list if m in struct_targets)
    return hits, len(struct_targets), targets_p

# ── pre-compute hairpin lists (one per stem/loop combo) ───────────────────────
print("Pre-computing hairpin lists...")
hairpin_cache = {}
for stem_len in STEM_VALUES:
    for loop_len in LOOP_VALUES:
        pa_file      = get_palindrome(stem_len, stem_len, loop_len, MM)
        all_hairpins = load.load_hairpins(str(pa_file))
        hairpins = [h for h in all_hairpins
                    if h.spacer_length == loop_len and not is_extendable(h)]
        hairpin_cache[(stem_len, loop_len)] = hairpins
        print(f"  stem={stem_len}, loop={loop_len}: "
              f"{len(hairpins)}/{len(all_hairpins)} hairpins "
              f"(exact loop, non-extendable)")

# ── compute all results ────────────────────────────────────────────────────────
# key: (sel_type, struct_type, stem_len, loop_len)
# val: (log10p, hits, n_struct, targets_p, pvalue)
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

                selected = select_hairpins(hairpins, sel_type)

                hits, n_struct, targets_p = compute_exclusive(selected, struct_type)

                if targets_p > 0:
                    result = binomtest(hits, _mut_cnt, targets_p, alternative='greater')
                    pvalue = result.pvalue
                else:
                    pvalue = 1.0

                log10p = -math.log10(pvalue) if pvalue > 0 else float("inf")
                results[(sel_type, struct_type, stem_len, loop_len)] = (
                    log10p, hits, n_struct, targets_p, pvalue
                )
                print(f"  stem={stem_len} loop={loop_len}: "
                      f"hits={hits}, tgt={n_struct}, "
                      f"p={pvalue:.5f}, -log10p={log10p:.3f}")

# ── panel-level combined stats ────────────────────────────────────────────────
_GENOME_DENOM = {
    "tc_end":          _denom_tc_ga,
    "c_end_excl":      _denom_c_not_tc,
    "spacer_interior": _denom_tc_ga,
    "stem":            _denom_tc_ga,
}
panel_stats = {}
for _sel in SELECTION_TYPES:
    for _str in STRUCTURE_TYPES:
        _non_empty = [
            results[(_sel, _str, sl, ll)]
            for sl in STEM_VALUES for ll in LOOP_VALUES
            if results[(_sel, _str, sl, ll)][2] > 0
        ]
        _hits    = sum(v[1] for v in _non_empty)
        _targets = sum(v[2] for v in _non_empty)
        _gdenom  = _GENOME_DENOM[_str]
        _p0      = min(_targets / _gdenom, 1.0) if _gdenom > 0 else 0.0
        _k       = min(_hits, _mut_cnt)
        if _p0 > 0:
            _pv = binomtest(_k, _mut_cnt, _p0, alternative='greater').pvalue
        else:
            _pv = 1.0
        _lp = -math.log10(_pv) if _pv > 0 else float('inf')
        panel_stats[(_sel, _str)] = (_hits, _targets, _p0, _lp)

# ── save results ───────────────────────────────────────────────────────────────
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
print(f"\nResults saved to results_mm{MM}.pkl")

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
            log10p, hits, n_struct, targets_p, pvalue = results[
                (sel_type, struct_type, stem_len, loop_len)
            ]
            log10p_str = f"{log10p:.2f}" if not math.isinf(log10p) else "inf"
            ann = f"{hits}/{n_struct}\n{log10p_str}"
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
    f"Binomial test significance  (mm={MM}, exclusive structure categories)\n"
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
            f"Binomial test: {struct_type} / {sel_type}  (mm={MM}, exclusive)\n"
            "Cell: hits/targets  |  −log₁₀(p)      Color: −log₁₀(p-value)",
            fontsize=10, y=1.01,
        )
        plt.tight_layout()
        fname = f"heatmap_mm{MM}_{sel_type}_{struct_type}"
        plt.savefig(out_dir / f"{fname}.pdf", dpi=150, bbox_inches="tight")
        plt.savefig(out_dir / f"{fname}.png", dpi=150, bbox_inches="tight")
        plt.close()
        print(f"Saved: {fname}.pdf / .png")
