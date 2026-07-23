#!/usr/bin/env python
"""Build the single cross-strain consensus reference from the per-isolate cut
assemblies, following the human-chosen spec (JSON emitted from config.consensus).

Usage: build_consensus.py consensus_spec.json out.fa

spec keys:
  cut:          {ISO: path_to_{ISO}.cut.fa}      # source assemblies
  gap_len:      int                              # N-run inserted at joins
  chromosomes:  [{name, take: "ISO:contig"}]     # 1 representative per chromosome
  joins:        [{name, pieces: ["ISO:c1","ISO:c2"]}]  # optional, concatenated with a gap
  unplaced:     ["ISO:contig", ...]              # carried through, renamed unpl01.. by size

Contig references are ISOLATE-PREFIXED ("MA:scaffold3") because the two isolates
share scaffold names.
"""
import sys, json

spec = json.load(open(sys.argv[1]))
outfa = sys.argv[2]
GAP = "N" * int(spec.get("gap_len", 100))

def read_fa(fn):
    d = {}; name = None; buf = []
    for l in open(fn):
        if l.startswith('>'):
            if name: d[name] = ''.join(buf)
            name = l[1:].split()[0]; buf = []
        else:
            buf.append(l.strip())
    if name: d[name] = ''.join(buf)
    return d

seqs = {iso: read_fa(path) for iso, path in spec["cut"].items()}

def get(ref):
    iso, ctg = ref.split(":", 1)
    if iso not in seqs or ctg not in seqs[iso]:
        sys.exit(f"ERROR: consensus spec references missing contig {ref}")
    return seqs[iso][ctg]

records = []   # (name, seq, source_label)
for c in spec.get("chromosomes", []):
    records.append((c["name"], get(c["take"]), c["take"]))
for j in spec.get("joins", []):
    seq = GAP.join(get(p) for p in j["pieces"])
    records.append((j["name"], seq, "+".join(j["pieces"])))

unpl = sorted(spec.get("unplaced", []), key=lambda r: -len(get(r)))
for i, ref in enumerate(unpl, 1):
    records.append((f"unpl{i:02d}", get(ref), ref))

with open(outfa, "w") as o, open(outfa + ".source.tsv", "w") as s:
    s.write("name\tsource\tlength\n")
    for name, seq, src in records:
        o.write(f">{name}\n")
        for k in range(0, len(seq), 80):
            o.write(seq[k:k+80] + "\n")
        s.write(f"{name}\t{src}\t{len(seq)}\n")

tot = sum(len(seq) for _, seq, _ in records)
sys.stderr.write(f"consensus: {len(records)} sequences, {tot/1e6:.2f} Mb "
                 f"({len(spec.get('chromosomes',[]))+len(spec.get('joins',[]))} chromosomes + "
                 f"{len(unpl)} unplaced) -> {outfa}\n")
