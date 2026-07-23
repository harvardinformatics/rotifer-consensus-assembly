#!/bin/bash
# filt_len.sh <in.fastq[.gz]> <out.fastq.gz> <min_len>
# Keep ONT reads with sequence length >= min_len. Dependency-free (zcat|awk|gzip).
set -eo pipefail
in="$1"; out="$2"; minlen="$3"
zcat -f "$in" | awk -v L="$minlen" '
  NR%4==1{h=$0} NR%4==2{s=$0} NR%4==3{p=$0}
  NR%4==0{ if(length(s)>=L) print h ORS s ORS p ORS $0 }' | gzip > "$out"
