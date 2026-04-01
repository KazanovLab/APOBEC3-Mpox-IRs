"""
Standalone script: load results_mm{MM}.pkl produced by test11.py and regenerate
the combined 5×4 heatmap figure (heatmap_all_mm{MM}.pdf / heatmap_all_mm{MM}.png).

Run from any directory:
    python plot_combined.py [--mm 0]
"""

import math
import pickle
import argparse
from pathlib import Path

import numpy as np
import matplotlib.pyplot as plt
from matplotlib.colors import LinearSegmentedColormap, TwoSlopeNorm

parser = argparse.ArgumentParser(description="Plot combined 5×4 heatmap from test11 results.")
parser.add_argument("--mm", type=int, default=0,
                    help="Number of mismatches used in test11 run (default: 0)")
args = parser.parse_args()
MM = args.mm

out_dir = Path(__file__).parent

# ── load results ──────────────────────────────────────────────────────────────
with open(out_dir / f"results_mm{MM}.pkl", "rb") as f:
    data = pickle.load(f)

results         = data["results"]
panel_stats     = data["panel_stats"]
LOOP_VALUES     = data["LOOP_VALUES"]
STEM_VALUES     = data["STEM_VALUES"]
STRUCTURE_TYPES = data["STRUCTURE_TYPES"]
SELECTION_TYPES = data["SELECTION_TYPES"]

# ── colormap ──────────────────────────────────────────────────────────────────
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

# ── draw one panel ─────────────────────────────────────────────────────────────
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
print(f"Saved: heatmap_all_mm{MM}.pdf / heatmap_all_mm{MM}.png")
