#!/usr/bin/env python
"""Blobtools-style contamination report for one primary assembly.

Combines per-contig GC, read coverage, windowed-BLAST taxonomy, and (optionally)
a fresh FCS-GX report into:
  --report      TSV: name, length, gc, coverage, cov_ratio, fcsgx, blast_kingdom,
                     candidate, reason        (the "data file with reasons")
  --candidates  contig names auto-flagged for removal (SUPERSET, for review)
  --plotprefix  writes {prefix}.blob.png (GC vs coverage, colored by call),
                {prefix}.gc_hist.png, {prefix}.cov_hist.png

A contig is a removal CANDIDATE if ANY of: FCS-GX flagged it; GC > gcflag;
dominant BLAST kingdom is non-metazoan (Bacteria/Fungi/Viruses/Archaea). The
reason column records which. This is a SUPERSET for human review, not an
auto-removal — the user curates {iso}.remove.txt from it.
"""
import sys, argparse
from collections import Counter, defaultdict
import matplotlib; matplotlib.use("Agg")
import matplotlib.pyplot as plt

p = argparse.ArgumentParser()
p.add_argument("--gc", required=True)         # asm_stats --per-contig: name length gc
p.add_argument("--cov", required=True)        # samtools coverage
p.add_argument("--fcsgx", required=True)      # FCS-GX report (or placeholder)
p.add_argument("--blastn", default=None)      # windowed blastn tsv (may be absent)
p.add_argument("--gcflag", type=float, default=0.45)
p.add_argument("--report", required=True)
p.add_argument("--candidates", required=True)
p.add_argument("--plotprefix", required=True)
args = p.parse_args()

# GC + length  (seqkit fx2tab -n -l -g: "name<TAB>length<TAB>GC%", no header)
gc, length = {}, {}
for l in open(args.gc):
    if l.startswith("#") or not l.strip(): continue
    f = l.split("\t"); g = float(f[2])
    if g > 1.5: g /= 100.0            # seqkit reports GC as a percentage -> fraction
    length[f[0]] = int(f[1]); gc[f[0]] = g

# coverage (samtools coverage: col1 rname, col7 meandepth)
cov = {}
for l in open(args.cov):
    if l.startswith("#") or not l.strip(): continue
    f = l.split("\t")
    if len(f) >= 7: cov[f[0]] = float(f[6])

# FCS-GX calls (skip placeholder / comment lines)
fcs = {}
for l in open(args.fcsgx):
    if l.startswith("#") or not l.strip(): continue
    f = l.split("\t")
    if len(f) >= 6: fcs[f[0]] = f"{f[4]}:{f[5]}"   # action:div

# BLAST dominant kingdom per contig (best hit per window)
btax = {}
if args.blastn:
    best = {}
    try:
        for l in open(args.blastn):
            f = l.rstrip("\n").split("\t")
            if len(f) >= 10 and f[0] not in best: best[f[0]] = f
    except FileNotFoundError:
        best = {}
    perc = defaultdict(Counter)
    for q, f in best.items():
        perc[q.split(":")[0]][f[8] or "?"] += 1
    for c, cnt in perc.items():
        k, v = cnt.most_common(1)[0]
        btax[c] = (k, v, sum(cnt.values()))

# host coverage baseline = median meandepth over big low-GC contigs
import statistics
host = [cov.get(c, 0) for c in gc if gc[c] < 0.40 and length[c] > 100_000 and c in cov]
host_med = statistics.median(host) if host else 1.0

NONMETA = ("Bacteria", "Fungi", "Viruses", "Archaea")
rows, cands = [], []
for c in sorted(gc, key=lambda x: -length[x]):
    g = gc[c]; L = length[c]; cv = cov.get(c, 0.0); r = cv / host_med if host_med else 0
    fk = fcs.get(c, "-")
    bk = btax.get(c, ("-", 0, 0))
    reasons = []
    if c in fcs: reasons.append(f"FCS-GX[{fk}]")
    if g > args.gcflag: reasons.append(f"GC={g:.2f}")
    if bk[0] in NONMETA: reasons.append(f"BLAST={bk[0]}({bk[1]}/{bk[2]})")
    if r > 2.5 or (r < 0.3 and cv > 0): reasons.append(f"cov_ratio={r:.1f}")
    cand = bool(c in fcs or g > args.gcflag or bk[0] in NONMETA)
    rows.append((c, L, g, cv, r, fk, f"{bk[0]}({bk[1]}/{bk[2]})", cand, ";".join(reasons) or "-"))
    if cand: cands.append(c)

with open(args.report, "w") as o:
    o.write("name\tlength\tgc\tcoverage\tcov_ratio\tfcsgx\tblast_kingdom\tcandidate\treason\n")
    for c, L, g, cv, r, fk, bk, cand, why in rows:
        o.write(f"{c}\t{L}\t{g:.3f}\t{cv:.1f}\t{r:.2f}\t{fk}\t{bk}\t{cand}\t{why}\n")
with open(args.candidates, "w") as o:
    for c in cands: o.write(c + "\n")

# ---- plots ----
def color(row):
    c, L, g, cv, r, fk, bk, cand, why = row
    if c in fcs or bk.startswith(NONMETA): return "#d62728"   # contaminant (red)
    if cand: return "#ff7f0e"                                  # candidate/unknown (orange)
    return "#7f7f7f"                                           # host (grey)

xs = [g for (_, _, g, *_ ) in rows]
ys = [max(cv, 0.1) for (_, _, _, cv, *_ ) in rows]
ss = [max(4, (L ** 0.5) / 30) for (_, L, *_ ) in rows]
cs = [color(row) for row in rows]

fig, ax = plt.subplots(figsize=(9, 7))
ax.scatter(xs, ys, s=ss, c=cs, alpha=0.6, edgecolors="none")
ax.set_yscale("log"); ax.set_xlabel("GC fraction"); ax.set_ylabel("ONT mean depth (log)")
ax.axvline(args.gcflag, color="k", ls=":", lw=0.8, label=f"GC flag={args.gcflag}")
ax.axhline(host_med, color="b", ls=":", lw=0.8, label=f"host median={host_med:.0f}x")
ax.set_title(f"blob: {args.plotprefix.split('/')[-1]}  (grey=host, orange=candidate, red=contaminant)")
ax.legend(loc="best", fontsize=8)
plt.tight_layout(); plt.savefig(args.plotprefix + ".blob.png", dpi=150); plt.close()

for vals, lab, fn in [([g for (_,_,g,*_) in rows], "GC fraction", ".gc_hist.png"),
                      (ys, "ONT mean depth", ".cov_hist.png")]:
    fig, ax = plt.subplots(figsize=(8, 4))
    ax.hist(vals, bins=60, color="#4c72b0")
    if lab.startswith("ONT"): ax.set_xscale("log")
    ax.set_xlabel(lab); ax.set_ylabel("contigs")
    plt.tight_layout(); plt.savefig(args.plotprefix + fn, dpi=150); plt.close()

sys.stderr.write(f"blob_report: {len(rows)} contigs, {len(cands)} candidates flagged "
                 f"(host median depth {host_med:.0f}x) -> {args.report}\n")
