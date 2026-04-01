"""
Extract all ct_end hairpins (mutated and unmutated) for specified (stem, loop) pairs.
Compare energy and trinucleotide context (third nucleotide after TC) between groups.

Selection: most_stable (one per overlapping group).
mm = 0.
"""

import sys, os, re
from pathlib import Path
from collections import defaultdict

BINOMIAL_DIR = Path(__file__).parent.parent.parent / "binomial"
sys.path.insert(0, str(BINOMIAL_DIR))
os.chdir(BINOMIAL_DIR)

import load
from emboss import get_palindrome
from functions import find_groups

MM    = 0
PAIRS = [(5, 3), (6, 3), (7, 3), (5, 4), (6, 4)]

_seq        = load.genome.sequence
_mut_set    = set(load.mutations_list)   # 0-based
_COMP       = {'A': 'T', 'T': 'A', 'G': 'C', 'C': 'G'}

def rc(base):
    return _COMP.get(base.upper(), 'N')

def is_extendable(h):
    if h.start == 0 or h.end >= len(_seq) - 1:
        return False
    return _COMP.get(_seq[h.start - 1].upper()) == _seq[h.end + 1].upper()

def most_stable(group):
    return min(group, key=lambda x: (x.score, x.length, x.start))

_tc_pos = set(m.end() - 1 for m in re.finditer("TC", _seq))   # C of TC (0-based)
_ga_pos = set(m.start()   for m in re.finditer("GA", _seq))   # G of GA (0-based)

# ── per-pair analysis ──────────────────────────────────────────────────────────
for stem_len, loop_len in PAIRS:
    pa_file  = get_palindrome(stem_len, stem_len, loop_len, MM)
    all_hp   = load.load_hairpins(str(pa_file))
    filtered = [h for h in all_hp
                if h.spacer_length == loop_len and not is_extendable(h)]
    groups   = find_groups(filtered)
    selected = sorted([most_stable(g) for g in groups], key=lambda h: h.start)

    rows = []
    for h in selected:
        l, r = h.spacer_index   # 0-based: l = left spacer end, r = right spacer end

        tc_target = r in _tc_pos   # C at right boundary, preceded by T
        ga_target = l in _ga_pos   # G at left boundary, followed by A

        # trinucleotide T[C]x on the relevant strand
        # TC_right: seq[r-1]=T, seq[r]=C, seq[r+1]=x  → third nuc = seq[r+1]
        # GA_left:  minus strand T[C]x where x = RC(seq[l-1])
        tc_third = _seq[r + 1].upper()  if (tc_target and r + 1 < len(_seq)) else None
        ga_third = rc(_seq[l - 1])      if (ga_target and l > 0)             else None

        mut_tc = tc_target and (r in _mut_set)
        mut_ga = ga_target and (l in _mut_set)

        rows.append(dict(
            pos      = h.start + 1,           # 1-based
            end      = h.end   + 1,
            score    = h.score,
            seq      = h.palindrome,          # stem SPACER stem with spaces
            tc_target= tc_target,
            ga_target= ga_target,
            tc_third = tc_third,
            ga_third = ga_third,
            mut_tc   = mut_tc,
            mut_ga   = mut_ga,
            mutated  = mut_tc or mut_ga,
        ))

    mutated   = [r for r in rows if     r['mutated']]
    unmutated = [r for r in rows if not r['mutated']]

    hdr = f"  stem={stem_len}  loop={loop_len}  ({len(selected)} hairpins: {len(mutated)} mutated, {len(unmutated)} unmutated)"
    print("\n" + "=" * 80)
    print(hdr)
    print("=" * 80)

    # ── full table ──────────────────────────────────────────────────────────────
    print(f"\n{'pos':>8}  {'score':>7}  {'mut':>5}  {'TC?':>4}  {'3rd':>4}  {'GA?':>4}  {'3rd':>4}  sequence")
    print("-" * 80)
    for row in rows:
        mut_str  = ("TC" if row['mut_tc'] else "  ") + ("+" if row['mut_tc'] and row['mut_ga'] else "") + ("GA" if row['mut_ga'] else "  ")
        mut_flag = "*" if row['mutated'] else " "
        tc_str   = "yes" if row['tc_target'] else " no"
        ga_str   = "yes" if row['ga_target'] else " no"
        print(f"{row['pos']:>8}  {row['score']:>7.2f}  {mut_flag:>1}{mut_str:>5}  "
              f"{tc_str:>4}  {row['tc_third'] or '-':>4}  "
              f"{ga_str:>4}  {row['ga_third'] or '-':>4}  "
              f"{row['seq']}")

    # ── energy comparison ───────────────────────────────────────────────────────
    print(f"\n── Energy (ΔG) ─────────────────────────────────────────────────────")
    for label, group in [("Mutated", mutated), ("Unmutated", unmutated)]:
        if not group:
            print(f"  {label:10}: none")
            continue
        scores = [r['score'] for r in group]
        print(f"  {label:10}: n={len(group):2}  mean={sum(scores)/len(scores):7.3f}"
              f"  min={min(scores):7.3f}  max={max(scores):7.3f}"
              f"  values: {[round(s,2) for s in sorted(scores)]}")

    # ── trinucleotide comparison ────────────────────────────────────────────────
    print(f"\n── Third nucleotide after TC (T[C]x) ──────────────────────────────")
    for label, group in [("Mutated", mutated), ("Unmutated", unmutated)]:
        thirds = defaultdict(int)
        for row in group:
            if row['tc_third']: thirds[row['tc_third']] += 1
            if row['ga_third']: thirds[row['ga_third']] += 1
        print(f"  {label:10}: {dict(sorted(thirds.items()))}")
