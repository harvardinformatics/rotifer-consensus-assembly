#!/usr/bin/env python
"""Whole-genome dotplot from a PAF (query on X, target on Y).

One line segment per alignment; scaffolds concatenated along each axis in size
order with faint boundary gridlines. This is a plain line-based dotplot (the kind
you actually want) -- no midpoint-scatter tricks.

Usage: paf_dotplot.py MAxMM.paf out.png [MINLEN_bp=20000]
"""
import sys
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt

paf, out = sys.argv[1], sys.argv[2]
MIN = int(sys.argv[3]) if len(sys.argv) > 3 else 20000

rows, qlen, tlen = [], {}, {}
for l in open(paf):
    f = l.split('\t')
    q, ql, qs, qe, strand = f[0], int(f[1]), int(f[2]), int(f[3]), f[4]
    t, tl, ts, te = f[5], int(f[6]), int(f[7]), int(f[8])
    qlen[q] = ql; tlen[t] = tl
    if qe - qs >= MIN:
        rows.append((q, qs, qe, strand, t, ts, te))

qo = sorted(qlen, key=lambda x: -qlen[x])
to = sorted(tlen, key=lambda x: -tlen[x])
qoff, o = {}, 0
for s in qo: qoff[s] = o; o += qlen[s]
QT = o
toff, o = {}, 0
for s in to: toff[s] = o; o += tlen[s]
TT = o

fig, ax = plt.subplots(figsize=(12, 12))
for q, qs, qe, strand, t, ts, te in rows:
    x0, x1 = qoff[q] + qs, qoff[q] + qe
    if strand == '+':
        y0, y1 = toff[t] + ts, toff[t] + te
    else:
        y0, y1 = toff[t] + te, toff[t] + ts
    ax.plot([x0, x1], [y0, y1], '-', lw=0.4, color='navy')
for s in qo: ax.axvline(qoff[s], color='0.85', lw=0.3)
for s in to: ax.axhline(toff[s], color='0.85', lw=0.3)
ax.set_xlim(0, QT); ax.set_ylim(0, TT)
ax.set_xlabel('MA (bp; scaffolds by size)')
ax.set_ylabel('MM (bp; scaffolds by size)')
ax.set_title('MA x MM synteny')
plt.tight_layout()
plt.savefig(out, dpi=150)
