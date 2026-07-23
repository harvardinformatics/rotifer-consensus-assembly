#!/usr/bin/env python
"""Emit N evenly-spaced fixed-width window regions for a contig (samtools faidx
-r format: contig:start-end). Used by the contamination-screen BLAST step.

Usage: make_windows.py <contig> <length> <n_windows> <window_bp>
"""
import sys

contig, L, n, w = sys.argv[1], int(sys.argv[2]), int(sys.argv[3]), int(sys.argv[4])
margin = max(int(L * 0.02), 0)
usable = L - 2 * margin
if usable < n * w:
    margin = max(int(L * 0.01), 0)
    usable = L - 2 * margin
if usable <= 0 or n < 1:
    sys.exit(0)
step = usable // (n + 1)
for i in range(1, n + 1):
    s = margin + i * step
    e = s + w - 1
    if 0 < s and e < L:
        print(f"{contig}:{s}-{e}")
