"""
Test 6: Number of hairpins and genome coverage for stem_min_length = 4, 5, 6.

All other parameters taken from binomial/config.yaml (stem_max=30, loop=10,
mismatches=1). Results presented as a bar chart with a secondary y-axis for
genome coverage percentage.
"""

import sys
import os
import subprocess
from pathlib import Path

# ── paths ─────────────────────────────────────────────────────────────────────
BINOMIAL_DIR = Path(__file__).parent.parent.parent / "binomial"
sys.path.insert(0, str(BINOMIAL_DIR))

import yaml
import matplotlib.pyplot as plt
import matplotlib.ticker as ticker
import numpy as np

os.chdir(BINOMIAL_DIR)

# ── load config ───────────────────────────────────────────────────────────────
with open("config.yaml") as f:
    params = yaml.safe_load(f)

GENOME_PATH    = params["genome_path"]
STEM_MAX       = params["stem_max_length"]
LOOP_LENGTH    = params["loop_length"]
NUM_MISMATCHES = params["number_mismatches"]

# ── genome length ──────────────────────────────────────────────────────────────
with open(GENOME_PATH) as f:
    f.readline()
    GENOME_LENGTH = len(f.read().replace("\n", ""))

# ── parse EMBOSS output ────────────────────────────────────────────────────────
def parse_emboss(filepath):
    """Return (count, intervals) from an EMBOSS palindrome output file.
    intervals: list of (start, end) tuples, 0-based with end exclusive."""
    with open(filepath) as f:
        lines = f.readlines()
    intervals = []
    for i in range(12, len(lines), 4):
        if lines[i].strip():
            pos = int(lines[i].split()[0])
            end = int(lines[i + 2].split()[0])
            intervals.append((pos - 1, end))   # 0-based, end exclusive
    return len(intervals), intervals

def genome_coverage_pct(intervals, genome_length):
    """Fraction of genome covered by merged intervals, as percent."""
    if not intervals:
        return 0.0
    sorted_ivs = sorted(intervals)
    covered = 0
    cur_s, cur_e = sorted_ivs[0]
    for s, e in sorted_ivs[1:]:
        if s < cur_e:
            cur_e = max(cur_e, e)
        else:
            covered += cur_e - cur_s
            cur_s, cur_e = s, e
    covered += cur_e - cur_s
    return covered * 100 / genome_length

# ── main loop ─────────────────────────────────────────────────────────────────
stem_min_values = [4, 5, 6]
emboss_dir      = BINOMIAL_DIR / "emboss"
emboss_dir.mkdir(exist_ok=True)

print(f"Fixed: stem_max={STEM_MAX}, loop={LOOP_LENGTH}, mismatches={NUM_MISMATCHES}")
print(f"{'stem_min':>8}  {'count':>9}  {'cov%':>7}")
print("-" * 32)

counts    = []
coverages = []
for stem_min in stem_min_values:
    tag         = f"mm{NUM_MISMATCHES}_stemin{stem_min}_stemax{STEM_MAX}_spacer{LOOP_LENGTH}"
    emboss_file = emboss_dir / f"emboss_{tag}.txt"

    if not emboss_file.exists():
        subprocess.run(
            [
                "palindrome",
                "-sequence",      str(GENOME_PATH),
                "-minpallen",     str(stem_min),
                "-maxpallen",     str(STEM_MAX),
                "-gaplimit",      str(LOOP_LENGTH),
                "-nummismatches", str(NUM_MISMATCHES),
                "-outfile",       str(emboss_file),
                "-overlap",
            ],
            check=True, capture_output=True,
        )

    n, intervals = parse_emboss(str(emboss_file))
    cov = genome_coverage_pct(intervals, GENOME_LENGTH)
    print(f"{stem_min:>8}  {n:>9,}  {cov:>6.2f}%")
    counts.append(n)
    coverages.append(cov)

# ── plot ──────────────────────────────────────────────────────────────────────
x      = np.arange(len(stem_min_values))
labels = [f"stem_min={v}" for v in stem_min_values]
bar_w  = 0.45
color_count = "#2166AC"
color_cov   = "#D6604D"

fig, ax = plt.subplots(figsize=(6, 4))
ax2 = ax.twinx()

bars = ax.bar(x, counts, width=bar_w, color=color_count, alpha=0.85,
              label="Hairpin count")
ax2.bar(x + bar_w, coverages, width=bar_w, color=color_cov, alpha=0.85,
        label="Genome coverage (%)")

# value labels on bars
for bar, n in zip(bars, counts):
    ax.text(bar.get_x() + bar.get_width() / 2, bar.get_height() + 150,
            f"{n:,}", ha="center", va="bottom", fontsize=9, color=color_count)
for xi, cov in zip(x + bar_w, coverages):
    ax2.text(xi, cov + 0.8, f"{cov:.1f}%",
             ha="center", va="bottom", fontsize=9, color=color_cov)

ax.set_xticks(x + bar_w / 2)
ax.set_xticklabels(labels, fontsize=11)
ax.set_ylabel("Number of hairpins", fontsize=12)
ax.yaxis.set_major_formatter(ticker.FuncFormatter(lambda v, _: f"{int(v):,}"))
ax.tick_params(axis="y", labelsize=10)
ax.set_ylim(0, max(counts) * 1.15)

ax2.set_ylim(0, 100)
ax2.set_ylabel("Genome coverage (%)", fontsize=12, color=color_cov)
ax2.tick_params(axis="y", labelcolor=color_cov, labelsize=10)
ax2.yaxis.set_major_formatter(ticker.FuncFormatter(lambda v, _: f"{v:.0f}%"))

ax.set_title(f"Hairpin count and genome coverage vs minimum stem length\n"
             f"(stem_max={STEM_MAX}, loop={LOOP_LENGTH}, mismatches={NUM_MISMATCHES})",
             fontsize=11)

# combined legend
h1, l1 = ax.get_legend_handles_labels()
h2, l2 = ax2.get_legend_handles_labels()
ax.legend(h1 + h2, l1 + l2, fontsize=10, loc="upper right")

ax.grid(axis="y", linestyle="--", alpha=0.4)
plt.tight_layout()

out_dir = Path(__file__).parent
plt.savefig(out_dir / "hairpins_vs_stem_min.pdf", dpi=150)
plt.savefig(out_dir / "hairpins_vs_stem_min.png", dpi=150)
print("\nPlot saved: hairpins_vs_stem_min.pdf / .png")
