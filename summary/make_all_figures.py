#!/usr/bin/env python3
"""Build every figure and summary table for the report in one command.

Replaces the hand-maintained list of per-script invocations, so adding or
reordering a figure is a one-line change here instead of a change in whatever
shell history produced the last run.

Usage
-----
    conda activate MASH
    python summary/make_all_figures.py              # everything
    python summary/make_all_figures.py manhattan    # one stage
    python summary/make_all_figures.py --list       # show stages

Why one driver
--------------
Figure layout used to be decided independently in each script (figsize, dpi,
font sizes, colour palette), which is how the Manhattan panels ended up with
mismatched x-axes while every other figure looked fine. The shared geometry now
lives in summary/plotstyle.py; this driver just fixes the execution order.

Environment
-----------
Python stages need the MASH env (nilearn, nibabel, matplotlib, pandas). The R
stages need the gdc env, so they are NOT run from here -- run them with

    conda activate gdc
    Rscript summary/04d_percent_SNP.R
    Rscript summary/03_plot_sa_horizontal.R
    Rscript summary/04_compare_r2.R
    Rscript summary/07_method_diagnostics.R

then render with quarto. Run `python summary/make_all_figures.py --list` to see
which stages still need R.
"""
from __future__ import annotations

import argparse
import subprocess
import sys
import time
from pathlib import Path

SUMMARY = Path(__file__).resolve().parent
ROOT = SUMMARY.parent

# Ordered: raw estimates -> derived tables -> per-panel figures. Later stages
# read files written by earlier ones, so do not reorder without checking.
STAGES = [
    ("manhattan", "05_manhattan.py", "Sys-Sys Manhattan panels, ordered by AdjHE-RE"),
    ("sa-correlation", "04b_h2_SA_correlation.py", "SA-FC correlations + per-network detail"),
    ("heatmap", "04c_FC_heatmap.py", "FC heatmaps (Twin + AdjHE-RE)"),
    ("within-outside", "04e_within_outside.py", "Within vs outside-network percentiles"),
    ("fc-distance", "04f_FC_distance.py", "FC vs Euclidean/geodesic distance"),
    ("gordon-brain", "02a_brain_gordon.py", "Gordon surfaces + chords"),
    ("proba-brain", "02b_brain_probaconns.py", "ProbaConns surfaces + chords"),
    ("sa-brain-wo", "02c_brain_sa.py", "SA wo_total surfaces"),
    ("sa-brain-wtot", "02d_brain_sa_total.py", "SA w_total surfaces"),
    ("trim", "trim_png_margins.py", "crop white margins from brain figures"),
]

# Produced by the R scripts, run separately in the gdc env.
R_STAGES = [
    ("r-percent-snp", "04d_percent_SNP.R", "Percent of twin h2 explained by SNP"),
    ("r-sa-barplot", "03_plot_sa_horizontal.R", "SA horizontal barplots (w_total + wo_total)"),
    ("r-compare-r2", "04_compare_r2.R", "Twin-vs-MASH R2, SA w_total vs wo_total"),
    ("r-diagnostics", "07_method_diagnostics.R", "Method diagnostics incl. GCTA/FE/HEreg"),
]


def _list_stages() -> None:
    print("Python stages (run in the MASH env):\n")
    for name, script, desc in STAGES:
        print(f"  {name:18} {script:32} {desc}")
    print("\nR stages (run separately in the gdc env):\n")
    for name, script, desc in R_STAGES:
        print(f"  {name:18} {script:32} {desc}")
    print("\nThen:  quarto render summary/06_results_catalog.qmd")


def _run(name: str, script: str, desc: str) -> int:
    path = SUMMARY / script
    if not path.exists():
        print(f"  SKIP {name}: {script} not found")
        return 1
    print(f"\n=== {name}: {desc}")
    # The trim stage edits the files the brain stages just wrote, so it runs
    # in place and keeps a .orig copy.
    argv = [sys.executable, str(path)]
    if name == "trim":
        argv += ["--dir", str(ROOT / "results" / "summary" / "brain"),
                 "--inplace", "--backup"]
    print(f"    {' '.join(argv[1:])}")
    t0 = time.time()
    proc = subprocess.run(argv, cwd=str(ROOT))
    dt = time.time() - t0
    if proc.returncode != 0:
        print(f"    FAILED ({proc.returncode}) after {dt:.1f}s")
    else:
        print(f"    ok ({dt:.1f}s)")
    return proc.returncode


def main() -> int:
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("stages", nargs="*", help="stage names to run (default: all)")
    ap.add_argument("--list", action="store_true", help="list stages and exit")
    ap.add_argument("--keep-going", action="store_true",
                    help="continue after a failing stage (default: stop)")
    args = ap.parse_args()

    if args.list:
        _list_stages()
        return 0

    all_stages = dict((n, (n, s, d)) for n, s, d in STAGES)
    wanted = [all_stages[n] for n in args.stages if n in all_stages]
    unknown = [n for n in args.stages if n not in all_stages]
    if unknown:
        print(f"Unknown stage(s): {', '.join(unknown)}", file=sys.stderr)
        _list_stages()
        return 2
    if not wanted:
        wanted = [(n, s, d) for n, s, d in STAGES]

    if any(n.startswith("r-") for n, _, _ in wanted):
        print("Note: R stages are not run by this driver; use Rscript in the gdc env.\n")

    failures = []
    for name, script, desc in wanted:
        rc = _run(name, script, desc)
        if rc != 0:
            failures.append(name)
            if not args.keep_going:
                print(f"\nStopped at {name}. Re-run with --keep-going to continue anyway.")
                break

    print("\n" + "=" * 60)
    if failures:
        print(f"FAILED stages: {', '.join(failures)}")
        return 1
    print("All requested Python stages completed.")
    print("Next: R stages in the gdc env, then quarto render.")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())