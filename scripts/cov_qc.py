#!/usr/bin/env python
"""Coverage / collapse / telomere QC on the consensus reference.

Maps BOTH isolates' ONT to the same reference (done in the rule); here we read the
two `samtools coverage` outputs and report per-scaffold depth, a collapse flag
(depth/median > collapse_ratio), a low-coverage flag (< lowcov_frac x median), and
an assumption-free telomere scan (dominant terminal 6-mer fraction + canonical
TTAGGG-family density). Median is taken over chromosome-named scaffolds (chr*).

Usage: cov_qc.py ref.fa cov_MA.txt cov_MM.txt \\
                 [collapse_ratio=1.7] [lowcov_frac=0.3] [telo_topk=0.25] [telo_canon=0.30]
"""
import sys, statistics
from collections import Counter

ref, covMA, covMM = sys.argv[1], sys.argv[2], sys.argv[3]
COLLAPSE = float(sys.argv[4]) if len(sys.argv) > 4 else 1.7
LOWCOV   = float(sys.argv[5]) if len(sys.argv) > 5 else 0.3
TK       = float(sys.argv[6]) if len(sys.argv) > 6 else 0.25
TC       = float(sys.argv[7]) if len(sys.argv) > 7 else 0.30

def readcov(fn):
    d = {}
    for l in open(fn):
        if l.startswith('#') or not l.strip(): continue
        f = l.split('\t'); d[f[0]] = float(f[6])   # col7 = meandepth
    return d

cMA, cMM = readcov(covMA), readcov(covMM)
seqs = {}; n = None; b = []
for line in open(ref):
    if line.startswith('>'):
        if n: seqs[n] = ''.join(b)
        n = line[1:].split()[0]; b = []
    else:
        b.append(line.strip().upper())
if n: seqs[n] = ''.join(b)

chrs = [s for s in seqs if s.startswith('chr')]
medMA = statistics.median([cMA.get(s, 0) for s in chrs]) or 1
medMM = statistics.median([cMM.get(s, 0) for s in chrs]) or 1

def revc(s): return s.translate(str.maketrans('ACGT', 'TGCA'))[::-1]
CANON = ['TTAGGG', 'TTAGGGG', 'TTAGG']

def top_kmer(sub, k=6):
    if len(sub) < k * 3: return ('', 0.0)
    c = Counter(sub[i:i+k] for i in range(len(sub) - k + 1))
    kmer, cnt = c.most_common(1)[0]
    return (kmer, cnt / (len(sub) - k + 1))

def canon_density(sub):
    return max((sub.count(m) + sub.count(revc(m))) * len(m) / max(1, len(sub)) for m in CANON)

W = 6000
print(f"median chromosome ONT meandepth: MA={medMA:.0f}x  MM={medMM:.0f}x   "
      f"(COLLAPSE if ratio>{COLLAPSE} ; LOWCOV if <{LOWCOV}x median)")
print(f"{'scaffold':9} {'Mb':>5} {'MAx':>5} {'MMx':>5} {'MA/m':>5} {'MM/m':>5}  "
      f"{'5p topk:frac/canon':>22} {'3p topk:frac/canon':>22}  flag")
for s in sorted(seqs, key=lambda x: -len(seqs[x])):
    seq = seqs[s]; L = len(seq); e5, e3 = seq[:W], seq[-W:]
    k5, f5 = top_kmer(e5); k3, f3 = top_kmer(e3)
    c5, c3 = canon_density(e5), canon_density(e3)
    rMA = cMA.get(s, 0) / medMA; rMM = cMM.get(s, 0) / medMM
    flag = ''
    if max(rMA, rMM) > COLLAPSE: flag += 'COLLAPSE? '
    if max(cMA.get(s, 0), cMM.get(s, 0)) < LOWCOV * max(medMA, medMM): flag += 'LOWCOV? '
    t5 = (f5 > TK or c5 > TC); t3 = (f3 > TK or c3 > TC)
    if t5 or t3: flag += 'TELO[' + ('5' if t5 else '') + ('3' if t3 else '') + ']'
    d5 = f"{k5}:{f5:.2f}/c{c5:.2f}"; d3 = f"{k3}:{f3:.2f}/c{c3:.2f}"
    print(f"{s:9} {L/1e6:5.1f} {cMA.get(s,0):5.0f} {cMM.get(s,0):5.0f} "
          f"{rMA:5.2f} {rMM:5.2f}  {d5:>22} {d3:>22}  {flag}")
