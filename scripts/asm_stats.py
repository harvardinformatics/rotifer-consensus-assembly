#!/usr/bin/env python
"""Assembly statistics.

Usage:
  asm_stats.py asm.fa                 # summary: #seqs, total, N50, largest, GC%
  asm_stats.py asm.fa --per-contig    # TSV: name  length  gc_fraction
"""
import sys

fa = sys.argv[1]
per_contig = "--per-contig" in sys.argv[2:]

def read_fa(fn):
    name = None; buf = []
    for l in open(fn):
        if l.startswith('>'):
            if name: yield name, ''.join(buf)
            name = l[1:].split()[0]; buf = []
        else:
            buf.append(l.strip())
    if name: yield name, ''.join(buf)

def gc(s):
    s = s.upper()
    at_gc = sum(s.count(b) for b in "ACGT")
    g = s.count("G") + s.count("C")
    return g / at_gc if at_gc else 0.0

recs = [(n, s) for n, s in read_fa(fa)]

if per_contig:
    print("name\tlength\tgc")
    for n, s in recs:
        print(f"{n}\t{len(s)}\t{gc(s):.4f}")
    sys.exit(0)

lens = sorted((len(s) for _, s in recs), reverse=True)
tot = sum(lens)
half, run, n50 = tot / 2, 0, 0
for L in lens:
    run += L
    if run >= half:
        n50 = L; break
allseq = "".join(s for _, s in recs)
print(f"sequences   : {len(recs)}")
print(f"total_bp    : {tot} ({tot/1e6:.2f} Mb)")
print(f"largest_bp  : {lens[0] if lens else 0}")
print(f"N50_bp      : {n50} ({n50/1e6:.2f} Mb)")
print(f"num_>=1Mb   : {sum(1 for L in lens if L>=1_000_000)}")
print(f"GC_percent  : {100*gc(allseq):.2f}")
