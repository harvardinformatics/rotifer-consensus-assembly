#!/usr/bin/env python
"""Select FASTA records by name.

Usage:
  fasta_select.py in.fa --drop name1 name2 ...      # write all EXCEPT these
  fasta_select.py in.fa --keep name1 name2 ...       # write ONLY these
  fasta_select.py in.fa --drop-file names.txt        # drop names listed in file
  fasta_select.py in.fa --keep-file names.txt        # keep only names listed in file
Writes to stdout. A file may have comments (#) and extra columns (first token = name).
Empty drop set (or empty file) passes everything through.
"""
import sys

fa = sys.argv[1]
a = sys.argv[2:]
mode = None; names = set()

def inline(flag):
    i = a.index(flag)
    return set(x for x in a[i+1:] if not x.startswith("--"))

def fromfile(flag):
    f = a[a.index(flag) + 1]
    return set(l.split()[0] for l in open(f) if l.strip() and not l.startswith("#"))

if "--drop" in a:        mode, names = "drop", inline("--drop")
elif "--keep" in a:      mode, names = "keep", inline("--keep")
elif "--drop-file" in a: mode, names = "drop", fromfile("--drop-file")
elif "--keep-file" in a: mode, names = "keep", fromfile("--keep-file")

def read_fa(fn):
    name = None; buf = []
    for l in open(fn):
        if l.startswith('>'):
            if name: yield name, buf
            name = l[1:].split()[0]; buf = [l]
        else:
            buf.append(l)
    if name: yield name, buf

out = sys.stdout; dropped = 0
for name, lines in read_fa(fa):
    keep = True
    if mode == "drop" and name in names: keep = False
    if mode == "keep" and name not in names: keep = False
    if keep: out.writelines(lines)
    else: dropped += 1
sys.stderr.write(f"fasta_select: mode={mode} named={len(names)} dropped={dropped}\n")
