"""
Test 5: Number of hairpins produced by EMBOSS palindrome vs loop_length.

Loops loop_length from 2 to 50, all other parameters taken from
binomial/config.yaml. For each value, hairpins are counted by parsing
the EMBOSS output file directly (no energy calculation), which avoids
the nn_energy() limitation of loops <= 30 nt.
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

os.chdir(BINOMIAL_DIR)

# ── load config ───────────────────────────────────────────────────────────────
with open("config.yaml") as f:
    params = yaml.safe_load(f)

GENOME_PATH    = params["genome_path"]
STEM_MIN       = params["stem_min_lenght"]       # note: typo preserved from config
STEM_MAX       = params["stem_max_length"]
NUM_MISMATCHES = params["number_mismatches"]

# ── count hairpins directly from EMBOSS output (no energy calculation) ────────
def count_hairpins_emboss(filepath):
    """Count hairpins by parsing EMBOSS palindrome output file directly."""
    with open(filepath) as f:
        lines = f.readlines()
    count = 0
    for i in range(12, len(lines), 4):
        if lines[i].strip():
            count += 1
    return count

# ── run EMBOSS and count ───────────────────────────────────────────────────────
loop_values = list(range(2, 51))
emboss_dir  = BINOMIAL_DIR / "emboss"
emboss_dir.mkdir(exist_ok=True)

print(f"Parameters: stem_min={STEM_MIN}, stem_max={STEM_MAX}, mismatches={NUM_MISMATCHES}")
print(f"{'loop_len':>8}  {'count':>9}")
print("-" * 22)

counts = []
for loop_len in loop_values:
    tag         = f"mm{NUM_MISMATCHES}_stemin{STEM_MIN}_stemax{STEM_MAX}_spacer{loop_len}"
    emboss_file = emboss_dir / f"emboss_{tag}.txt"

    if not emboss_file.exists():
        subprocess.run(
            [
                "palindrome",
                "-sequence",      str(GENOME_PATH),
                "-minpallen",     str(STEM_MIN),
                "-maxpallen",     str(STEM_MAX),
                "-gaplimit",      str(loop_len),
                "-nummismatches", str(NUM_MISMATCHES),
                "-outfile",       str(emboss_file),
                "-overlap",
            ],
            check=True, capture_output=True,
        )

    n = count_hairpins_emboss(str(emboss_file))
    print(f"{loop_len:>8}  {n:>9}")
    counts.append(n)

# ── plot ──────────────────────────────────────────────────────────────────────
fig, ax = plt.subplots(figsize=(7, 4))

ax.plot(loop_values, counts,
        marker="o", markersize=4, linewidth=1.5, color="#2166AC")

ax.set_xlabel("Maximum loop length (nt)", fontsize=12)
ax.set_ylabel("Number of hairpins", fontsize=12)
ax.set_title("EMBOSS palindrome hairpin count vs maximum loop length\n"
             f"(stem_min={STEM_MIN}, stem_max={STEM_MAX}, mismatches={NUM_MISMATCHES})",
             fontsize=11)

ax.xaxis.set_major_locator(ticker.MultipleLocator(10))
ax.xaxis.set_minor_locator(ticker.MultipleLocator(2))
ax.yaxis.set_major_formatter(ticker.FuncFormatter(lambda x, _: f"{int(x):,}"))
ax.grid(axis="y", linestyle="--", alpha=0.4)
ax.tick_params(axis="both", labelsize=10)

plt.tight_layout()

out_dir = Path(__file__).parent
plt.savefig(out_dir / "hairpins_vs_loop_len.pdf", dpi=150)
plt.savefig(out_dir / "hairpins_vs_loop_len.png", dpi=150)
print("\nPlot saved: hairpins_vs_loop_len.pdf / .png")
