#!/usr/bin/env python
"""Cross-isolate scaffold correspondence + per-homolog %id/length from a BASE-LEVEL
MA x MM PAF (minimap2 -c: carries the de:f: gap-compressed divergence tag).

PAF query = MA (col1), target = MM (col6)  -- the rule runs `minimap2 MM MA`.

For each MA<->MM scaffold pair we accumulate aligned length (sum of alignment-block
lengths) and a length-weighted gap-compressed %id = 100*(1 - sum(de*alen)/sum(alen)).
Pairs are written to edges.tsv. Only pairs with aln >= MINALN AND %id >= MINID
(true orthologs) are used to build correspondence components -- diverged
paralogs/ohnologs are kept in edges.tsv but must NOT merge chromosomes into false
"complex" clusters. Nodes are isolate-prefixed so a shared scaffold name across
isolates is never merged.

Components:  1:1 (one MA + one MM) | simple-split (1:N, under-assembly) | complex (N:M).

Usage: synteny_classify.py MAxMM.paf correspondence.tsv edges.tsv [MINALN=1e6] [MINID=0.90]
"""
import sys
from collections import defaultdict

paf, out, edgeout = sys.argv[1:4]
MIN   = int(float(sys.argv[4])) if len(sys.argv) > 4 else 1_000_000
MINID = float(sys.argv[5]) if len(sys.argv) > 5 else 0.90

alen = defaultdict(int); dew = defaultdict(float)
for l in open(paf):
    f = l.rstrip("\n").split('\t')
    if len(f) < 12:
        continue
    q, t = "MA:" + f[0], "MM:" + f[5]          # query=MA, target=MM
    a = int(f[10])
    de = 0.0
    for x in f[12:]:
        if x.startswith("de:f:"):
            de = float(x[5:]); break
    alen[(q, t)] += a
    dew[(q, t)] += de * a

def pid(k):                                     # gap-compressed identity fraction
    return 1 - dew[k] / alen[k] if alen[k] else 0.0

edges = {k: v for k, v in alen.items() if v >= MIN}

with open(edgeout, 'w') as e:
    e.write("MA\tMM\taln_bp\tpct_id\ttype\n")
    for (q, t), v in sorted(edges.items(), key=lambda kv: -kv[1]):
        typ = "ortholog" if pid((q, t)) >= MINID else "paralog"
        e.write(f"{q[3:]}\t{t[3:]}\t{v}\t{100*pid((q,t)):.2f}\t{typ}\n")

# correspondence components built on ORTHOLOG edges only
orth = {k: v for k, v in edges.items() if pid(k) >= MINID}
adj = defaultdict(set)
for (q, t) in orth:
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

def c_edges(c):
    return [k for k in orth if k[0] in c and k[1] in c]
def c_aln(c):
    return sum(orth[k] for k in c_edges(c))
def c_id(c):
    tot = c_aln(c); d = sum(dew[k] for k in c_edges(c))
    return 1 - d / tot if tot else 0.0

with open(out, 'w') as o:
    o.write("class\tn_MA\tn_MM\tMA_scaffolds\tMM_scaffolds\taln_bp\tpct_id\n")
    n11 = nsplit = ncomplex = 0
    for c in sorted(comps, key=lambda c: -c_aln(c)):
        ma = sorted(x[3:] for x in c if x.startswith("MA:"))
        mm = sorted(x[3:] for x in c if x.startswith("MM:"))
        if len(ma) == 1 and len(mm) == 1:
            cls = "1:1"; n11 += 1
        elif (len(ma) == 1) != (len(mm) == 1):
            cls = "simple-split"; nsplit += 1
        else:
            cls = "complex"; ncomplex += 1
        o.write(f"{cls}\t{len(ma)}\t{len(mm)}\t{','.join(ma)}\t{','.join(mm)}\t{c_aln(c)}\t{100*c_id(c):.2f}\n")
sys.stderr.write(f"components: {n11} 1:1, {nsplit} simple-split, {ncomplex} complex "
                 f"(MINID={MINID}, MINALN={MIN}); {len(edges)-len(orth)} paralog edges excluded\n")
