"""
Extract all ct_end hairpins (mutated and unmutated) for specified (stem, loop) pairs.
Compare energy and trinucleotide context (third nucleotide after TC) between groups.

Selection: most_stable (one per overlapping group). mm = 0.
Only hairpins with a TC or GA spacer boundary are included.

Saves to out_dir:
  hairpins_all.tsv      — per-hairpin table across all pairs
  summary.txt           — energy + trinucleotide summary per pair
  energy_comparison.pdf — boxplots + stripplots of energy by group
"""

import sys, os, re
from pathlib import Path
from collections import defaultdict, Counter

import numpy as np
import matplotlib.pyplot as plt
from scipy.stats import mannwhitneyu, ttest_ind, permutation_test, fisher_exact

BINOMIAL_DIR = Path(__file__).parent.parent.parent / "binomial"
sys.path.insert(0, str(BINOMIAL_DIR))
os.chdir(BINOMIAL_DIR)

import load
from emboss import get_palindrome
from functions import find_groups

MM      = 0
PAIRS   = [(5, 3), (6, 3), (7, 3), (5, 4), (6, 4)]
OUT_DIR = Path(__file__).parent
OUT_DIR.mkdir(exist_ok=True)

_seq        = load.genome.sequence
_mut_counts = Counter(load.mutations_list)   # pos → number of mutation observations
_COMP    = {'A': 'T', 'T': 'A', 'G': 'C', 'C': 'G'}

def rc(base):
    return _COMP.get(base.upper(), 'N')

def is_extendable(h):
    if h.start == 0 or h.end >= len(_seq) - 1:
        return False
    return _COMP.get(_seq[h.start - 1].upper()) == _seq[h.end + 1].upper()

def most_stable(group):
    return min(group, key=lambda x: (x.score, x.length, x.start))

_tc_pos = set(m.end() - 1 for m in re.finditer("TC", _seq))
_ga_pos = set(m.start()   for m in re.finditer("GA", _seq))

# ── stat helpers ───────────────────────────────────────────────────────────────
def mean_diff_stat(x, y):
    """Test statistic: mean(x) - mean(y)."""
    return np.mean(x) - np.mean(y)

def run_stats(mut_scores, unm_scores, label=""):
    """Return dict of statistical results; print summary."""
    n_m, n_u = len(mut_scores), len(unm_scores)
    results = dict(n_mut=n_m, n_unm=n_u)

    if n_m == 0 or n_u == 0:
        results.update(mw_p=None, mw_stat=None, tt_p=None, tt_stat=None,
                       perm_p=None, mean_diff=None)
        return results

    m_arr = np.array(mut_scores)
    u_arr = np.array(unm_scores)
    results['mean_diff'] = float(np.mean(m_arr) - np.mean(u_arr))

    # Mann-Whitney U (two-sided — small samples make direction unclear)
    if n_m >= 2 and n_u >= 2:
        mw = mannwhitneyu(m_arr, u_arr, alternative='two-sided')
        results['mw_stat'] = float(mw.statistic)
        results['mw_p']    = float(mw.pvalue)
        # rank-biserial correlation as effect size
        results['rbc'] = float(1 - 2 * mw.statistic / (n_m * n_u))
    else:
        results['mw_stat'] = None
        results['mw_p']    = None
        results['rbc']     = None

    # Welch t-test (two-sided)
    if n_m >= 2 and n_u >= 2:
        tt = ttest_ind(m_arr, u_arr, equal_var=False, alternative='two-sided')
        results['tt_stat'] = float(tt.statistic)
        results['tt_p']    = float(tt.pvalue)
    else:
        results['tt_stat'] = None
        results['tt_p']    = None

    # Permutation test on mean difference (two-sided, 9999 permutations)
    if n_m >= 2 and n_u >= 2:
        perm = permutation_test(
            (m_arr, u_arr), mean_diff_stat,
            permutation_type='independent',
            alternative='two-sided',
            n_resamples=9999,
            random_state=42,
        )
        results['perm_p'] = float(perm.pvalue)
    else:
        results['perm_p'] = None

    return results

def fmt_p(p):
    if p is None: return "   n/a  "
    if p < 0.001: return f"{p:.2e}"
    return f"{p:.4f}"

# ── collect rows ───────────────────────────────────────────────────────────────
all_rows    = []
pair_data   = {}   # (stem, loop) -> {'mut': [scores], 'unm': [scores], 'rows': [...]}
summary_lines = []

for stem_len, loop_len in PAIRS:
    pa_file  = get_palindrome(stem_len, stem_len, loop_len, MM)
    all_hp   = load.load_hairpins(str(pa_file))
    filtered = [h for h in all_hp
                if h.spacer_length == loop_len and not is_extendable(h)]
    groups   = find_groups(filtered)
    selected = sorted([most_stable(g) for g in groups], key=lambda h: h.start)

    rows = []
    for h in selected:
        l, r = h.spacer_index
        if (r not in _tc_pos) and (l not in _ga_pos):
            continue

        tc_target = r in _tc_pos
        ga_target = l in _ga_pos
        tc_third  = _seq[r + 1].upper() if (tc_target and r + 1 < len(_seq)) else ""
        ga_third  = rc(_seq[l - 1])     if (ga_target and l > 0)             else ""
        n_tc = _mut_counts.get(r, 0) if tc_target else 0
        n_ga = _mut_counts.get(l, 0) if ga_target else 0
        n_mut = n_tc + n_ga   # total mutation observations at ct_end positions

        mut_label = ""
        if n_tc: mut_label += f"TC({n_tc})"
        if n_ga: mut_label += ("+" if n_tc else "") + f"GA({n_ga})"
        if not mut_label: mut_label = "-"

        rows.append(dict(
            stem=stem_len, loop=loop_len,
            pos=h.start + 1, end=h.end + 1,
            score=h.score,
            tc_target="yes" if tc_target else "no", tc_third=tc_third,
            ga_target="yes" if ga_target else "no", ga_third=ga_third,
            n_mut=n_mut, mutated=mut_label, sequence=h.palindrome,
        ))
        all_rows.append(rows[-1])

    mut_rows = [r for r in rows if r['mutated'] != "-"]
    unm_rows = [r for r in rows if r['mutated'] == "-"]
    # weight each mutated hairpin by its number of mutation observations
    mut_scores = [r['score'] for r in mut_rows for _ in range(r['n_mut'])]
    unm_scores = [r['score'] for r in unm_rows]
    def _trinucs(group_rows, weighted=True):
        d = defaultdict(int)
        for row in group_rows:
            w = row['n_mut'] if (weighted and row['n_mut'] > 0) else 1
            if row['tc_third']: d[row['tc_third']] += w
            if row['ga_third']: d[row['ga_third']] += w
        return d

    pair_data[(stem_len, loop_len)] = dict(
        mut=mut_scores, unm=unm_scores, rows=rows,
        trinuc_mut=_trinucs(mut_rows, weighted=True),
        trinuc_unm=_trinucs(unm_rows, weighted=False),
    )

    stats = run_stats(mut_scores, unm_scores)

    lines = []
    n_mut_obs = sum(r['n_mut'] for r in mut_rows)
    lines.append("\n" + "=" * 80)
    lines.append(f"  stem={stem_len}  loop={loop_len}  "
                 f"({len(rows)} ct_end hairpins: "
                 f"{len(mut_rows)} mutated [{n_mut_obs} observations], "
                 f"{len(unm_rows)} unmutated)")
    lines.append("=" * 80)

    lines.append(f"\n{'pos':>8}  {'score':>7}  {'mut':>10}  "
                 f"{'TC?':>4}  {'3rd':>4}  {'GA?':>4}  {'3rd':>4}  sequence")
    lines.append("-" * 80)
    for row in rows:
        lines.append(
            f"{row['pos']:>8}  {row['score']:>7.2f}  {row['mutated']:>10}  "
            f"{row['tc_target']:>4}  {row['tc_third'] or '-':>4}  "
            f"{row['ga_target']:>4}  {row['ga_third'] or '-':>4}  "
            f"{row['sequence']}"
        )

    # ── energy detail ──────────────────────────────────────────────────────────
    lines.append(f"\n── Energy (ΔG) {'─'*60}")
    for label, scores, hairpins in [
        ("Mutated",   mut_scores, mut_rows),
        ("Unmutated", unm_scores, unm_rows),
    ]:
        if not scores:
            lines.append(f"  {label:10}: none")
            continue
        n_obs = len(scores)
        n_hp  = len(hairpins)
        lines.append(
            f"  {label:10}: {n_hp} hairpins / {n_obs} observations  "
            f"mean={np.mean(scores):7.3f}  median={np.median(scores):7.3f}  "
            f"sd={np.std(scores, ddof=1) if n_obs>1 else 0.0:6.3f}  "
            f"min={min(scores):6.2f}  max={max(scores):6.2f}"
        )
        lines.append(f"  {'':10}  values: {[round(s,2) for s in sorted(scores)]}")

    lines.append(f"\n── Statistical tests (mutated vs unmutated energy) {'─'*28}")
    lines.append(f"  mean difference (mut − unm): "
                 f"{stats['mean_diff']:+.3f}" if stats['mean_diff'] is not None
                 else "  mean difference: n/a")
    lines.append(f"  Mann-Whitney U:  stat={stats['mw_stat'] if stats['mw_stat'] is not None else 'n/a':>8}  "
                 f"p={fmt_p(stats['mw_p'])}  "
                 f"rank-biserial r={stats['rbc']:+.3f}" if stats['rbc'] is not None
                 else "  Mann-Whitney U:  n/a (n<2)")
    lines.append(f"  Welch t-test:    stat={stats['tt_stat']:.3f}  p={fmt_p(stats['tt_p'])}"
                 if stats['tt_p'] is not None
                 else "  Welch t-test:    n/a (n<2)")
    lines.append(f"  Permutation:     p={fmt_p(stats['perm_p'])}"
                 if stats['perm_p'] is not None
                 else "  Permutation:     n/a (n<2)")

    # trinucleotide (mutated weighted by observation count)
    lines.append(f"\n── Third nucleotide after TC (T[C]x) {'─'*40}")
    for label, group in [("Mutated", mut_rows), ("Unmutated", unm_rows)]:
        thirds = defaultdict(int)
        for row in group:
            w = row['n_mut'] if row['n_mut'] > 0 else 1
            if row['tc_third']: thirds[row['tc_third']] += w
            if row['ga_third']: thirds[row['ga_third']] += w
        lines.append(f"  {label:10}: {dict(sorted(thirds.items()))}")

    summary_lines.extend(lines)
    for line in lines:
        print(line)

# ── combined analysis across all pairs ────────────────────────────────────────
all_mut = [r['score'] for r in all_rows if r['n_mut'] > 0 for _ in range(r['n_mut'])]
all_unm = [r['score'] for r in all_rows if r['n_mut'] == 0]
stats_all = run_stats(all_mut, all_unm)

comb_lines = []
comb_lines.append("\n" + "=" * 80)
comb_lines.append("  COMBINED (all pairs pooled)")
comb_lines.append("=" * 80)
n_mut_hp = sum(1 for r in all_rows if r['n_mut'] > 0)
comb_lines.append(f"  Mutated  : {n_mut_hp} hairpins / {len(all_mut)} observations  "
                  f"mean={np.mean(all_mut):.3f}  median={np.median(all_mut):.3f}  "
                  f"sd={np.std(all_mut, ddof=1):.3f}")
comb_lines.append(f"  Unmutated: {len(all_unm)} hairpins  "
                  f"mean={np.mean(all_unm):.3f}  median={np.median(all_unm):.3f}  "
                  f"sd={np.std(all_unm, ddof=1):.3f}")
comb_lines.append(f"  mean diff (mut − unm): {stats_all['mean_diff']:+.3f}")
comb_lines.append(f"  Mann-Whitney U: stat={stats_all['mw_stat']:.1f}  "
                  f"p={fmt_p(stats_all['mw_p'])}  "
                  f"rank-biserial r={stats_all['rbc']:+.3f}")
comb_lines.append(f"  Welch t-test:   stat={stats_all['tt_stat']:.3f}  p={fmt_p(stats_all['tt_p'])}")
comb_lines.append(f"  Permutation:    p={fmt_p(stats_all['perm_p'])}")

summary_lines.extend(comb_lines)
for line in comb_lines:
    print(line)

# ── save summary.txt ───────────────────────────────────────────────────────────
with open(OUT_DIR / "summary.txt", "w") as f:
    f.write("\n".join(summary_lines) + "\n")

# ── save hairpins_all.tsv ──────────────────────────────────────────────────────
tsv_fields = ["stem", "loop", "pos", "end", "score",
              "tc_target", "tc_third", "ga_target", "ga_third",
              "n_mut", "mutated", "sequence"]
with open(OUT_DIR / "hairpins_all.tsv", "w") as f:
    f.write("\t".join(tsv_fields) + "\n")
    for row in all_rows:
        f.write("\t".join(str(row[k]) for k in tsv_fields) + "\n")

# ── energy boxplot figure ──────────────────────────────────────────────────────
n_pairs  = len(PAIRS)
fig, axes = plt.subplots(1, n_pairs + 1, figsize=(4 * (n_pairs + 1), 5))

def draw_energy_ax(ax, mut_rows, unm_rows, title, stats):
    mut_sc = [r['score'] for r in mut_rows for _ in range(r['n_mut'])]
    unm_sc = [r['score'] for r in unm_rows]
    data   = [mut_sc, unm_sc]
    n_mut_hp  = len(mut_rows)
    n_mut_obs = sum(r['n_mut'] for r in mut_rows)
    labels = [f"Mutated\n({n_mut_hp} hp / {n_mut_obs} obs)",
              f"Unmutated\n({len(unm_rows)} hp)"]
    colors = ["#d62728", "#1f77b4"]

    bp = ax.boxplot(data, patch_artist=True, widths=0.4,
                    medianprops=dict(color="black", linewidth=2))
    for patch, color in zip(bp['boxes'], colors):
        patch.set_facecolor(color)
        patch.set_alpha(0.6)

    # strip: mutated hairpins — dot size proportional to mutation count
    rng = np.random.default_rng(0)
    for row in mut_rows:
        jitter = rng.uniform(-0.12, 0.12)
        ax.scatter(1 + jitter, row['score'],
                   color=colors[0], alpha=0.85,
                   s=20 + 25 * row['n_mut'], zorder=3)
    for sc in unm_sc:
        jitter = rng.uniform(-0.12, 0.12)
        ax.scatter(2 + jitter, sc, color=colors[1], alpha=0.7, s=20, zorder=3)

    ax.set_xticks([1, 2])
    ax.set_xticklabels(labels, fontsize=7)
    ax.set_ylabel("Energy (ΔG)", fontsize=9)
    ax.set_title(title, fontsize=9, fontweight="bold")

    # annotate p-values
    p_lines = []
    if stats['mw_p']   is not None: p_lines.append(f"MWU  p={fmt_p(stats['mw_p']).strip()}")
    if stats['tt_p']   is not None: p_lines.append(f"t    p={fmt_p(stats['tt_p']).strip()}")
    if stats['perm_p'] is not None: p_lines.append(f"perm p={fmt_p(stats['perm_p']).strip()}")
    if stats['mean_diff'] is not None:
        p_lines.insert(0, f"Δmean={stats['mean_diff']:+.2f}")
    ax.text(0.97, 0.97, "\n".join(p_lines),
            transform=ax.transAxes, fontsize=7,
            va="top", ha="right",
            bbox=dict(boxstyle="round,pad=0.3", facecolor="white", alpha=0.8))

for ax, (stem_len, loop_len) in zip(axes[:-1], PAIRS):
    d = pair_data[(stem_len, loop_len)]
    st = run_stats(d['mut'], d['unm'])
    draw_energy_ax(ax, [r for r in d['rows'] if r['n_mut'] > 0],
                        [r for r in d['rows'] if r['n_mut'] == 0],
                        f"stem={stem_len}, loop={loop_len}", st)

# combined panel
all_mut_rows = [r for r in all_rows if r['n_mut'] > 0]
all_unm_rows = [r for r in all_rows if r['n_mut'] == 0]
draw_energy_ax(axes[-1], all_mut_rows, all_unm_rows, "All pairs combined", stats_all)

fig.suptitle("Energy (ΔG) — mutated vs unmutated ct_end hairpins (mm=0, most_stable)",
             fontsize=11, y=1.01)
plt.tight_layout()
plt.savefig(OUT_DIR / "energy_comparison.pdf", dpi=150, bbox_inches="tight")
plt.savefig(OUT_DIR / "energy_comparison.png", dpi=150, bbox_inches="tight")
plt.close()

# ── trinucleotide figure ───────────────────────────────────────────────────────
NUCS   = ['A', 'C', 'G', 'T']
COLORS = {"Mutated": "#d62728", "Unmutated": "#1f77b4"}

def trinuc_counts(rows, weighted):
    d = defaultdict(int)
    for row in rows:
        w = row['n_mut'] if (weighted and row['n_mut'] > 0) else 1
        if row['tc_third']: d[row['tc_third']] += w
        if row['ga_third']: d[row['ga_third']] += w
    return d

def draw_trinuc_ax(ax, tc_mut, tc_unm, title):
    """Grouped bar chart of trinucleotide proportions; annotate with counts."""
    x     = np.arange(len(NUCS))
    width = 0.35

    for offset, (label, tc) in enumerate([("Mutated", tc_mut), ("Unmutated", tc_unm)]):
        counts = np.array([tc.get(n, 0) for n in NUCS])
        total  = counts.sum()
        props  = counts / total if total > 0 else counts
        bars   = ax.bar(x + (offset - 0.5) * width, props,
                        width, label=label, color=COLORS[label], alpha=0.75)
        for bar, cnt in zip(bars, counts):
            if cnt > 0:
                ax.text(bar.get_x() + bar.get_width() / 2,
                        bar.get_height() + 0.01,
                        str(int(cnt)), ha='center', va='bottom', fontsize=7)

    # Fisher's exact test: G vs non-G
    g_mut  = tc_mut.get('G', 0);  ng_mut  = sum(tc_mut.get(n, 0)  for n in NUCS if n != 'G')
    g_unm  = tc_unm.get('G', 0);  ng_unm  = sum(tc_unm.get(n, 0) for n in NUCS if n != 'G')
    table  = [[g_mut, ng_mut], [g_unm, ng_unm]]
    if g_mut + ng_mut > 0 and g_unm + ng_unm > 0:
        _, p_fish = fisher_exact(table, alternative='two-sided')
        fish_str = f"Fisher G vs rest: p={fmt_p(p_fish).strip()}"
    else:
        fish_str = "Fisher: n/a"

    ax.set_xticks(x)
    ax.set_xticklabels([f"TC{n}" for n in NUCS], fontsize=8)
    ax.set_ylabel("Proportion", fontsize=9)
    ax.set_title(title, fontsize=9, fontweight="bold")
    ax.set_ylim(0, 1.15)
    ax.legend(fontsize=7)
    ax.text(0.97, 0.97, fish_str, transform=ax.transAxes,
            fontsize=7, va='top', ha='right',
            bbox=dict(boxstyle="round,pad=0.3", facecolor="white", alpha=0.8))

fig2, axes2 = plt.subplots(1, len(PAIRS) + 1, figsize=(4 * (len(PAIRS) + 1), 5))
axes2_flat  = axes2.flatten()

for ax, (stem_len, loop_len) in zip(axes2_flat, PAIRS):
    d = pair_data[(stem_len, loop_len)]
    draw_trinuc_ax(ax, d['trinuc_mut'], d['trinuc_unm'],
                   f"stem={stem_len}, loop={loop_len}")

# combined panel
all_mut_rows_t = [r for r in all_rows if r['n_mut'] > 0]
all_unm_rows_t = [r for r in all_rows if r['n_mut'] == 0]
draw_trinuc_ax(axes2_flat[5],
               trinuc_counts(all_mut_rows_t, weighted=True),
               trinuc_counts(all_unm_rows_t, weighted=False),
               "All pairs combined")

fig2.suptitle(
    "Third nucleotide after TC  (T[C]x) — mutated (obs-weighted) vs unmutated\n"
    "Bar height = proportion within group; numbers = counts",
    fontsize=11, y=1.01,
)
plt.tight_layout()
plt.savefig(OUT_DIR / "trinuc_comparison.pdf", dpi=150, bbox_inches="tight")
plt.savefig(OUT_DIR / "trinuc_comparison.png", dpi=150, bbox_inches="tight")
plt.close()

print(f"\nSaved: summary.txt")
print(f"Saved: hairpins_all.tsv  ({len(all_rows)} hairpins)")
print(f"Saved: energy_comparison.pdf / .png")
print(f"Saved: trinuc_comparison.pdf / .png")
