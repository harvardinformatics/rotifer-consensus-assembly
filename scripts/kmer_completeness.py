#!/usr/bin/env python
"""Copy-number-stratified k-mer completeness from a merqury spectra-cn histogram.

The headline merqury completeness counts ALL reliable read k-mers; this splits
them by read multiplicity (single-copy / 2-copy / repeat) and reports what
fraction of each stratum the assembly contains. For a purged haploid assembly of
a heterozygous genome the single-copy peak is ~half heterozygous (the dropped
allele reads as "missing" by design); the repeat strata measure true repeat
capture. Interpret with the spectra-cn plot.

Usage: kmer_completeness.py <sample>.spectra-cn.hist <filt_threshold>
"""
import sys

hist, filt = sys.argv[1], int(sys.argv[2])
present = {}; missing = {}; maxm = 0
for i, l in enumerate(open(hist)):
    if i == 0: continue
    f = l.split()
    if len(f) < 3: continue
    cls, m, n = f[0], int(f[1]), int(f[2]); maxm = max(maxm, m)
    (missing if cls == 'read-only' else present).setdefault(m, 0)
    (missing if cls == 'read-only' else present)[m] += n

def tot(d, lo, hi): return sum(v for k, v in d.items() if lo <= k < hi)

best, c = 0, filt
for m in range(filt, min(maxm, filt * 8)):
    t = present.get(m, 0) + missing.get(m, 0)
    if t > best: best, c = t, m

P, M = tot(present, filt, 10**18), tot(missing, filt, 10**18)
print(f"reliable threshold(filt)={filt} ; detected 1-copy peak multiplicity c={c}")
print(f"OVERALL reliable completeness = {100*P/(P+M):.2f}%   (present={P:,}  missing={M:,})")
bins = [(filt, round(1.5*c), 'single-copy (unique+het)'),
        (round(1.5*c), round(2.5*c), '2-copy (homozygous/dup)'),
        (round(2.5*c), round(6*c), '3-6x'),
        (round(6*c), 10**18, 'high-copy repeat (>6x)')]
print(f"{'stratum':30}{'present':>15}{'missing':>15}{'completeness':>14}")
for lo, hi, lab in bins:
    p, m = tot(present, lo, hi), tot(missing, lo, hi)
    comp = 100*p/(p+m) if p+m else 0
    print(f"{lab:30}{p:>15,}{m:>15,}{comp:>13.1f}%")
if M:
    print(f"\nOf ALL missing reliable k-mers ({M:,}):")
    print(f"  single-copy [{filt},{round(1.5*c)}): {100*tot(missing,filt,round(1.5*c))/M:.1f}%")
    print(f"  2-copy [{round(1.5*c)},{round(2.5*c)}): {100*tot(missing,round(1.5*c),round(2.5*c))/M:.1f}%")
    print(f"  >=3-copy (repeat): {100*tot(missing,round(2.5*c),10**18)/M:.1f}%")
