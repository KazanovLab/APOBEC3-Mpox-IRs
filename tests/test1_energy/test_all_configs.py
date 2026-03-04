"""
Run binomial/main.py for every combination of structure_type × hairpin_selection_type
and collect p-values into a tab-separated table.

Output: tests/pvalue_results.txt
"""

import subprocess
import sys
import os
import re
from itertools import product

BINOMIAL_DIR = os.path.join(os.path.dirname(os.path.dirname(os.path.dirname(os.path.abspath(__file__)))), "binomial")
CONFIG_PATH = os.path.join(BINOMIAL_DIR, "config.yaml")
OUTPUT_PATH = os.path.join(os.path.dirname(os.path.abspath(__file__)), "pvalue_results.txt")

STRUCTURE_TYPES = ["hairpin", "spacer", "ct_end", "c_end"]
SELECTION_TYPES = ["most_stable", "greedy", "max_cov", "min_cov", "all"]


def write_config(structure_type, selection_type):
    with open(CONFIG_PATH, "w") as f:
        f.write(f'structure_type: "{structure_type}" # options: hairpin, spacer, ct_end, c_end\n')
        f.write(f'hairpin_selection_type: "{selection_type}" # options: most_stable, greedy, max_cov, min_cov, all\n')


def run_main():
    env = os.environ.copy()
    env["MPLBACKEND"] = "Agg"  # suppress interactive plot windows
    result = subprocess.run(
        [sys.executable, "main.py"],
        cwd=BINOMIAL_DIR,
        capture_output=True,
        text=True,
        env=env,
    )
    return result.stdout, result.stderr, result.returncode


def parse_pvalue(stdout):
    match = re.search(r"p-value:\s*([\d.]+)%", stdout)
    if match:
        return match.group(1)
    return "NA"


def main():
    with open(CONFIG_PATH, "r") as f:
        original_config = f.read()

    # rows: structure_type, columns: selection_type
    results = {}

    try:
        total = len(STRUCTURE_TYPES) * len(SELECTION_TYPES)
        done = 0
        for structure_type, selection_type in product(STRUCTURE_TYPES, SELECTION_TYPES):
            done += 1
            label = f"{structure_type} / {selection_type}"
            print(f"[{done}/{total}] {label} ...", end=" ", flush=True)

            write_config(structure_type, selection_type)
            stdout, stderr, returncode = run_main()
            pvalue = parse_pvalue(stdout)

            if returncode != 0:
                print(f"ERROR (rc={returncode})")
                if stderr.strip():
                    # show only first meaningful error line
                    first_err = next((l for l in stderr.splitlines() if l.strip()), "")
                    print(f"   {first_err[:120]}")
                pvalue = "ERROR"
            else:
                print(f"p = {pvalue}%")

            results.setdefault(structure_type, {})[selection_type] = pvalue

    finally:
        with open(CONFIG_PATH, "w") as f:
            f.write(original_config)

    # ── write tab-separated table ──────────────────────────────────────────────
    col_w = 12  # fixed column width for aligned plain-text view

    header_cols = ["structure_type \\ selection_type"] + SELECTION_TYPES
    sep = "\t"

    lines = []
    lines.append(sep.join(header_cols))
    for st in STRUCTURE_TYPES:
        row = [st] + [results.get(st, {}).get(sel, "NA") for sel in SELECTION_TYPES]
        lines.append(sep.join(row))

    table_text = "\n".join(lines) + "\n"

    with open(OUTPUT_PATH, "w") as f:
        f.write(table_text)

    # also print a readable version to stdout
    print()
    print("=" * 70)
    print("p-value table (%):")
    print("=" * 70)
    # header
    print(f"{'':20}", end="")
    for sel in SELECTION_TYPES:
        print(f"{sel:>12}", end="")
    print()
    print("-" * 70)
    for st in STRUCTURE_TYPES:
        print(f"{st:<20}", end="")
        for sel in SELECTION_TYPES:
            val = results.get(st, {}).get(sel, "NA")
            print(f"{val:>12}", end="")
        print()

    print()
    print(f"Results saved to: {OUTPUT_PATH}")


if __name__ == "__main__":
    main()
