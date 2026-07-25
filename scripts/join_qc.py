#!/usr/bin/env python
"""ONT join QC on a (cut) assembly.

For each large scaffold: find N-gaps (scaffolder joins) and count ONT reads that
SPAN each gap (primary alignment covering gap +/- flank). A gap with no spanning
reads is an unsupported join = candidate cut. Also scans mid-contig windows for
unexpectedly low spanning (a join the scaffolder didn't mark with an N-gap).

Usage: join_qc.py assembly.fa reads.bam [min_scaffold_mb=10] [flank_bp=3000] [min_span=5]
Requires pysam. BAM must be coordinate-sorted + indexed.
"""
import sys, re, os
import pysam

fa, bam = sys.argv[1], sys.argv[2]
MINSCAF = int(float(sys.argv[3]) * 1e6) if len(sys.argv) > 3 else 10_000_000
FLANK   = int(sys.argv[4]) if len(sys.argv) > 4 else 3000
MINSPAN = int(sys.argv[5]) if len(sys.argv) > 5 else 5
STEP = 500_000

def read_fa(fn):
    name = None; seq = []
    for l in open(fn):
        if l.startswith('>'):
            if name: yield name, ''.join(seq)
            name = l[1:].split()[0]; seq = []
        else:
            seq.append(l.strip())
    if name: yield name, ''.join(seq)

def spanning(bf, ctg, pos, flank=FLANK):
    lo, hi = pos - flank, pos + flank
    n = 0
    for r in bf.fetch(ctg, max(0, lo), hi):
        if r.is_unmapped or r.is_secondary or r.is_supplementary:
            continue
        if r.reference_start <= lo and r.reference_end >= hi:
            n += 1
    return n

# fetch()/count() need a .bai; map_ont's index is temp() and may be gone -> self-index.
if not (os.path.exists(bam + ".bai") or os.path.exists(re.sub(r"\.bam$", ".bai", bam))):
    sys.stderr.write("join_qc: no BAM index found -> pysam.index(bam)\n")
    pysam.index(bam)
bf = pysam.AlignmentFile(bam)
seqs = dict(read_fa(fa))
print(f"### ONT join QC (scaffolds >= {MINSCAF/1e6:.0f} Mb) ; span +/-{FLANK} bp ; min_span={MINSPAN} ###")
for name in sorted(seqs, key=lambda x: -len(seqs[x])):
    seq = seqs[name]; L = len(seq)
    if L < MINSCAF:
        continue
    gaps = [(m.start(), m.end()) for m in re.finditer('N+', seq)]
    print(f"\n== {name} ({L/1e6:.1f} Mb) ; {len(gaps)} N-gaps ==")
    for gs, ge in gaps:
        mid = (gs + ge) // 2
        dep = bf.count(name, max(0, mid - 1), mid + 1)
        sp = spanning(bf, name, mid)
        tag = "spanned/ok" if sp >= MINSPAN else "*** UNSUPPORTED (candidate cut)"
        print(f"   GAP @ {mid/1e6:6.2f}Mb (size {ge-gs}bp): depth {dep:4d}, spanning {sp:4d} -> {tag}")
    low = [f"{p/1e6:.1f}Mb" for p in range(STEP, L - STEP, STEP)
           if not any(gs <= p <= ge for gs, ge in gaps) and spanning(bf, name, p) < MINSPAN]
    print("   mid-contig low-spanning (no gap): " + (", ".join(low) if low else "none"))
