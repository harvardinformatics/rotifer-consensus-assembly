#!/usr/bin/env python
"""Classify cross-isolate scaffold correspondence from a MA x MM PAF.

Edges = scaffold pairs sharing >= MINALN aligned bp. Nodes are ISOLATE-PREFIXED
(MA:<scaf> / MM:<scaf>) so that scaffolds sharing a name across isolates do NOT
get merged into one component (this bug produced false "complex" clusters before).
Connected components are classified:
  1:1           one MA + one MM  (clean)
  simple-split  1:N              (benign -- underassembly in one isolate)
  complex       N:M (>1 each)    (the real problem to fix)

Usage: synteny_classify.py MAxMM.paf out.tsv [MINALN_bp=1000000]
PAF query = MA, target = MM.
"""
import sys
from collections import defaultdict

paf, out = sys.argv[1], sys.argv[2]
MIN = int(sys.argv[3]) if len(sys.argv) > 3 else 1_000_000

w = defaultdict(int)
for l in open(paf):
    f = l.split('\t')
    q, t, nmatch = "MA:" + f[0], "MM:" + f[5], int(f[9])
    w[(q, t)] += nmatch
edges = {k: v for k, v in w.items() if v >= MIN}

adj = defaultdict(set)
for (q, t) in edges:
    adj[q].add(t); adj[t].add(q)

seen, comps = set(), []
for n in list(adj):
    if n in seen:
        continue
    stack, comp = [n], set()
    while stack:
        x = stack.pop()
        if x in seen:
            continue
        seen.add(x); comp.add(x)
        stack += [y for y in adj[x] if y not in seen]
    comps.append(comp)

def comp_aln(c):
    return sum(v for (a, b), v in edges.items() if a in c and b in c)

with open(out, 'w') as o:
    o.write("class\tn_MA\tn_MM\tMA_scaffolds\tMM_scaffolds\ttotal_aln_bp\n")
    n11 = nsplit = ncomplex = 0
    for c in sorted(comps, key=lambda c: -comp_aln(c)):
        ma = sorted(x[3:] for x in c if x.startswith("MA:"))
        mm = sorted(x[3:] for x in c if x.startswith("MM:"))
        if len(ma) == 1 and len(mm) == 1:
            cls = "1:1"; n11 += 1
        elif (len(ma) == 1) != (len(mm) == 1):
            cls = "simple-split"; nsplit += 1
        else:
            cls = "complex"; ncomplex += 1
        o.write(f"{cls}\t{len(ma)}\t{len(mm)}\t{','.join(ma)}\t{','.join(mm)}\t{comp_aln(c)}\n")
sys.stderr.write(f"components: {n11} 1:1, {nsplit} simple-split, {ncomplex} complex\n")
