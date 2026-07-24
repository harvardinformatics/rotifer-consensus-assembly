#!/usr/bin/env python
"""Plot the purge_dups read-depth histogram (PB.stat) with the calcuts cutoffs
overlaid, so a human can judge whether calcuts placed the haploid/diploid boundary
correctly -- the key check for a possibly-bimodal / degenerate-tetraploid depth
distribution before committing to scaffolding.

Best-effort: on any parse failure it still writes a placeholder PNG carrying the
error text, so the workflow never dies on a plotting hiccup.

Usage: purge_hist.py --stat PB.stat --cuts cutoffs --title STR --out hist.png
  PB.stat : whitespace-separated, first two numeric columns = <depth> <count>
  cutoffs : calcuts output line(s); all integers found are drawn as vertical cuts
"""
import sys, argparse
import matplotlib; matplotlib.use("Agg")
import matplotlib.pyplot as plt

p = argparse.ArgumentParser()
p.add_argument("--stat", required=True)
p.add_argument("--cuts", required=True)
p.add_argument("--title", default="purge_dups depth")
p.add_argument("--out", required=True)
a = p.parse_args()

try:
    depth, count = [], []
    for l in open(a.stat):
        f = l.split()
        if len(f) < 2:
            continue
        try:
            d = float(f[0]); c = float(f[1])
        except ValueError:
            continue
        depth.append(d); count.append(c)

    cuts = []
    for l in open(a.cuts):
        if l.startswith("#") or not l.strip():
            continue
        for tok in l.replace(",", " ").split():
            try:
                cuts.append(int(float(tok)))
            except ValueError:
                pass
        break
    cuts = sorted({c for c in cuts if c > 0})

    fig, ax = plt.subplots(figsize=(9, 5))
    if depth:
        ax.fill_between(depth, count, step="mid", color="#4c72b0", alpha=0.6)
        ax.plot(depth, count, color="#274b7a", lw=0.8)
    ymax = max(count) if count else 1
    for x in cuts:
        ax.axvline(x, color="#d62728", ls="--", lw=1)
        ax.text(x, ymax * 0.97, str(x), color="#d62728", rotation=90,
                va="top", ha="right", fontsize=8)
    # focus the x-range: out to ~2.5x the largest cutoff (drops the long high-depth tail)
    if depth:
        hi = (max(cuts) * 2.5) if cuts else max(depth)
        ax.set_xlim(0, max(hi, 10))
    ax.set_xlabel("read depth"); ax.set_ylabel("bases at this depth")
    ax.set_title(a.title + "   (red dashed = calcuts cutoffs)")
    plt.tight_layout(); plt.savefig(a.out, dpi=150); plt.close()
    sys.stderr.write(f"purge_hist: {len(depth)} depth bins, cutoffs={cuts} -> {a.out}\n")
except Exception as e:
    fig, ax = plt.subplots(figsize=(9, 5))
    ax.text(0.5, 0.5, f"purge_hist failed:\n{e}\n\ninspect PB.stat + cutoffs by hand",
            ha="center", va="center", wrap=True)
    ax.axis("off"); plt.savefig(a.out, dpi=150); plt.close()
    sys.stderr.write(f"purge_hist ERROR: {e}\n")
