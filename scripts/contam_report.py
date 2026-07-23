#!/usr/bin/env python
"""Contamination report for the consensus reference.

Per-contig GC + (optional) windowed-BLAST taxonomy summary for GC-flagged contigs
-> a TSV suggesting which contigs to drop. The final removal list is a REVIEWED
config value (contam.drop_contigs); this only produces the evidence.

Usage: contam_report.py consensus.fa gc_flag_high <contam_dir>
  <contam_dir> may contain flagged.txt, blastn.tsv, blastx.tsv (from the rule).
BLAST tsv outfmt: qseqid pident length evalue bitscore staxids sscinames
                  sblastnames sskingdoms stitle
"""
import sys, os
from collections import Counter, defaultdict

fa, gc_hi, d = sys.argv[1], float(sys.argv[2]), sys.argv[3]

def read_fa(fn):
    name = None; buf = []
    for l in open(fn):
        if l.startswith('>'):
            if name: yield name, ''.join(buf)
            name = l[1:].split()[0]; buf = []
        else: buf.append(l.strip())
    if name: yield name, ''.join(buf)

def gc(s):
    s = s.upper(); acgt = sum(s.count(b) for b in "ACGT")
    return (s.count("G") + s.count("C")) / acgt if acgt else 0.0

def load_blast(fn):
    best = {}  # qseqid -> row (first = top bitscore)
    if os.path.exists(fn):
        for l in open(fn):
            f = l.rstrip("\n").split("\t")
            if len(f) >= 10 and f[0] not in best:
                best[f[0]] = f
    per_contig = defaultdict(Counter)  # contig -> Counter(sskingdom)
    for q, f in best.items():
        per_contig[q.split(":")[0]][f[8] or "?"] += 1
    return per_contig

bn = load_blast(os.path.join(d, "blastn.tsv"))
bx = load_blast(os.path.join(d, "blastx.tsv"))

def dom(counter):
    if not counter: return "-"
    k, v = counter.most_common(1)[0]
    return f"{k}({v}/{sum(counter.values())})"

print("contig\tlength\tgc\tgc_flag\tblastn_kingdom\tblastx_kingdom\tcandidate_drop")
for name, seq in sorted(read_fa(fa), key=lambda x: -len(x[1])):
    g = gc(seq); flag = "HIGH_GC" if g > gc_hi else ""
    kn, kx = dom(bn.get(name, Counter())), dom(bx.get(name, Counter()))
    non_metazoan = any(x in (kn + kx) for x in ("Bacteria", "Fungi", "Viruses", "Archaea"))
    cand = "DROP?" if (flag and (non_metazoan or (not bn and not bx))) else ""
    print(f"{name}\t{len(seq)}\t{g:.3f}\t{flag}\t{kn}\t{kx}\t{cand}")
