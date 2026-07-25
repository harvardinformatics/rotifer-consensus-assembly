#!/usr/bin/env python
"""Apply de-chimerization cuts to a scaffold FASTA.

Usage: apply_cuts.py in.fa cuts.tsv out.fa [snap_bp=20000]
  cuts.tsv : <scaffold>\\t<pos_mb>   (one per line; a scaffold may appear >once)

Each listed scaffold is split at the given position(s); a cut is snapped to the
nearest N-gap within snap_bp (scaffolder joins are N-gaps, so a legitimate cut
lands in a gap). Pieces are renamed <scaffold>_p1, _p2, ... in coordinate order.
Scaffolds with no cut pass through unchanged.
"""
import sys, re
from collections import defaultdict

infa, cutf, outfa = sys.argv[1:4]
SNAP = int(sys.argv[4]) if len(sys.argv) > 4 else 20000

cuts = defaultdict(list)
for line in open(cutf):
    if not line.strip():
        continue
    s, mb = line.split()[:2]
    cuts[s].append(int(round(float(mb) * 1e6)))

def read_fa(fn):
    name = None; seq = []
    for l in open(fn):
        if l.startswith('>'):
            if name: yield name, ''.join(seq)
            name = l[1:].split()[0]; seq = []
        else:
            seq.append(l.strip())
    if name: yield name, ''.join(seq)

def snap_to_gap(seq, pos, snap):
    lo, hi = max(0, pos - snap), min(len(seq), pos + snap)
    best, at = snap + 1, None
    for g in re.finditer('N{10,}', seq[lo:hi]):
        c = lo + (g.start() + g.end()) // 2
        if abs(c - pos) < best:
            best, at = abs(c - pos), c
    return at if at is not None else pos

with open(outfa, 'w') as o:
    for name, seq in read_fa(infa):
        short = name.split(",")[0]      # cuts.tsv / cut_sites use the short id (before the comma)
        if short not in cuts:
            o.write(f">{short}\n{seq}\n"); continue
        pts = sorted(set(min(max(1, snap_to_gap(seq, p, SNAP)), len(seq) - 1) for p in cuts[short]))
        bounds = [0] + pts + [len(seq)]
        for i in range(len(bounds) - 1):
            o.write(f">{short}_p{i+1}\n{seq[bounds[i]:bounds[i+1]]}\n")
        sys.stderr.write(f"{short}: cut into {len(bounds)-1} at {[round(p/1e6,2) for p in pts]} Mb\n")
