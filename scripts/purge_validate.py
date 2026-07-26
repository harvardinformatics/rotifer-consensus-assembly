#!/usr/bin/env python
"""Validate a purge_dups run: is hap.fa (the sequence purge REMOVED) genuinely
redundant haplotypic duplication, or did purge cut UNIQUE content (over-purge)?

Deeply-diverged ohnologs are NOT a concern here -- they are unalignable, so
purge_dups' self-alignment never acts on them; whatever purge removed it removed
because it aligned to a retained copy, i.e. it should be allelic/haplotypic.
This script tests that expectation from three independent angles, per isolate:

  (1) compleasm markers  -- OVER-purge  = markers Complete in hap.fa but Missing
                            in purged.fa; UNDER-purge = Duplicated in purged.fa.
  (2) hap->purged alignment -- each removed contig's aligned fraction back to the
                            retained assembly (a real haplotypic dup aligns back).
  (3) coverage on the COMBINED (purged + hap) reference, from BOTH ONT and 10x
                            independently -- a redundant hap copy sits at ~HALF the
                            retained (diploid) depth (its reads split with the
                            purged homolog); UNIQUE/over-purged sequence sits at
                            ~FULL depth with no homolog.

Writes a text report + a multi-panel PNG. Best-effort: a failure in any one input
degrades that section rather than killing the report.

Usage: see argparse below (driven by the purge_validate Snakemake rule).
"""
import sys, argparse, gzip, statistics
from collections import defaultdict
import matplotlib; matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np

HAP_PREFIX = "HAP__"          # purgeqc_ref prefixes hap contig names in the combined ref
BIG = 5_000_000               # "big retained contig" length for estimating full/diploid depth


def read_labels(fn):
    """contig -> (cls in {purged,hap}, length)   (combined-ref names; hap = HAP__*)"""
    d = {}
    for l in open(fn):
        f = l.rstrip("\n").split("\t")
        if len(f) >= 3 and f[1] in ("purged", "hap"):
            d[f[0]] = (f[1], int(f[2]))
    return d


def read_mosdepth_summary(fn):
    """contig -> mean depth (mosdepth *.mosdepth.summary.txt; skip total/_region rows)"""
    d = {}
    try:
        for l in open(fn):
            f = l.rstrip("\n").split("\t")
            if len(f) < 4 or f[0] in ("chrom", "total") or f[0].endswith("_region"):
                continue
            try:
                d[f[0]] = float(f[3])
            except ValueError:
                pass
    except Exception as e:
        sys.stderr.write(f"[warn] mosdepth summary {fn}: {e}\n")
    return d


def read_regions(fn):
    """list of (contig, win_len, mean_depth) from mosdepth *.regions.bed.gz"""
    out = []
    try:
        op = gzip.open if fn.endswith(".gz") else open
        with op(fn, "rt") as fh:
            for l in fh:
                f = l.rstrip("\n").split("\t")
                if len(f) < 4:
                    continue
                try:
                    out.append((f[0], int(f[2]) - int(f[1]), float(f[3])))
                except ValueError:
                    pass
    except Exception as e:
        sys.stderr.write(f"[warn] regions {fn}: {e}\n")
    return out


def parse_compleasm_summary(fn):
    """category letter -> (pct, count); plus 'N' total.  Robust to S/D/F/I/M ordering."""
    cats = {}
    try:
        for l in open(fn):
            l = l.strip()
            if len(l) >= 2 and l[0] in "SDFIMN" and l[1] == ":":
                rest = l[2:]
                if l[0] == "N":
                    cats["N"] = (None, int(rest.split()[0].replace(",", "")))
                else:
                    pct = float(rest.split("%")[0])
                    cnt = int(rest.split(",")[1]) if "," in rest else 0
                    cats[l[0]] = (pct, cnt)
    except Exception as e:
        sys.stderr.write(f"[warn] compleasm summary {fn}: {e}\n")
    return cats


def parse_full_table(fn):
    """busco_id -> status  (compleasm full_table.tsv; col0=id, col1=Status)."""
    d = {}
    try:
        for l in open(fn):
            if l.startswith("#") or not l.strip():
                continue
            f = l.rstrip("\n").split("\t")
            if len(f) >= 2:
                # keep the "best" status if an id appears twice. compleasm's native
                # full_table.tsv uses "Single" for complete single-copy (not BUSCO's "Complete").
                prev = d.get(f[0])
                rank = {"Single": 3, "Complete": 3, "Duplicated": 3, "Fragmented": 2, "Incomplete": 1, "Missing": 0}
                if prev is None or rank.get(f[1], 0) > rank.get(prev, 0):
                    d[f[0]] = f[1]
    except Exception as e:
        sys.stderr.write(f"[warn] full_table {fn}: {e}\n")
    return d


def hap_redundancy(paf):
    """query(hap contig, stripped of HAP__) -> (aligned_frac, wtd_pct_id) vs purged.
    Merges query intervals so overlapping alignments don't inflate the fraction."""
    ivs = defaultdict(list)      # q -> [(qs,qe,de,alen)]
    qlen = {}
    try:
        for l in open(paf):
            f = l.rstrip("\n").split("\t")
            if len(f) < 12:
                continue
            q = f[0]; qlen[q] = int(f[1])
            de = next((float(x[5:]) for x in f[12:] if x.startswith("de:f:")), None)
            ivs[q].append((int(f[2]), int(f[3]), de, int(f[10])))
    except Exception as e:
        sys.stderr.write(f"[warn] hap paf {paf}: {e}\n")
    out = {}
    for q, lst in ivs.items():
        lst.sort()
        merged = 0; cs, ce = None, None
        for s, e, _, _ in lst:
            if cs is None:
                cs, ce = s, e
            elif s <= ce:
                ce = max(ce, e)
            else:
                merged += ce - cs; cs, ce = s, e
        if cs is not None:
            merged += ce - cs
        dew = sum((d if d is not None else 0) * a for _, _, d, a in lst)
        aln = sum(a for _, _, _, a in lst)
        pid = 100 * (1 - dew / aln) if aln else float("nan")
        out[q] = (merged / qlen[q] if qlen.get(q) else 0.0, pid)
    return out


def depth_mode(regions, contigs):
    """dominant depth (histogram peak) over a set of contigs, weighted by window bp."""
    vals = [(d, wl) for c, wl, d in regions if c in contigs and d > 0]
    if not vals:
        return float("nan")
    d = np.array([v[0] for v in vals]); w = np.array([v[1] for v in vals])
    hi = np.percentile(d, 99.5)
    bins = np.linspace(0, max(hi, 1), 200)
    h, edges = np.histogram(d, bins=bins, weights=w)
    return 0.5 * (edges[h.argmax()] + edges[h.argmax() + 1])


def main():
    p = argparse.ArgumentParser()
    p.add_argument("--iso", required=True)
    p.add_argument("--labels", required=True)
    p.add_argument("--cmpl-prepped"); p.add_argument("--cmpl-purged"); p.add_argument("--cmpl-hap")
    p.add_argument("--ft-prepped"); p.add_argument("--ft-purged"); p.add_argument("--ft-hap")
    p.add_argument("--hap-paf", required=True)
    p.add_argument("--ont-summary", required=True); p.add_argument("--ont-regions", required=True)
    p.add_argument("--tenx-summary", required=True); p.add_argument("--tenx-regions", required=True)
    p.add_argument("--genomescope")
    p.add_argument("--min-frac", type=float, default=0.5)
    p.add_argument("--seqkit-stats")
    p.add_argument("--out-report", required=True)
    p.add_argument("--out-plot", required=True)
    a = p.parse_args()

    labels = read_labels(a.labels)
    purged = {c for c, (cl, _) in labels.items() if cl == "purged"}
    hap    = {c for c, (cl, _) in labels.items() if cl == "hap"}
    ont = read_mosdepth_summary(a.ont_summary); tenx = read_mosdepth_summary(a.tenx_summary)
    ont_reg = read_regions(a.ont_regions);      tenx_reg = read_regions(a.tenx_regions)
    red = hap_redundancy(a.hap_paf)             # keyed by stripped hap name

    # retained "full/diploid" depth = mode over big purged contigs' windows
    bigpurged = {c for c in purged if labels[c][1] >= BIG}
    D_ont = depth_mode(ont_reg, bigpurged or purged)
    D_10x = depth_mode(tenx_reg, bigpurged or purged)

    # ---- classify each hap contig -------------------------------------------
    rows = []
    tot = defaultdict(lambda: defaultdict(int))     # signal -> class -> bp
    for c in sorted(hap, key=lambda x: -labels[x][1]):
        L = labels[c][1]
        short = c[len(HAP_PREFIX):] if c.startswith(HAP_PREFIX) else c
        frac, pid = red.get(short, (0.0, float("nan")))
        do = ont.get(c, 0.0); dt = tenx.get(c, 0.0)
        ro = do / D_ont if D_ont and D_ont == D_ont else float("nan")
        rt = dt / D_10x if D_10x and D_10x == D_10x else float("nan")

        def classify(ratio):
            if ratio != ratio:              # nan
                return "?"
            if ratio < 0.15:
                return "lowcov/junk"
            if frac >= a.min_frac and ratio <= 0.75:
                return "haplotypic-dup"
            if frac < a.min_frac and ratio >= 0.75:
                return "UNIQUE(over-purge)"
            return "ambiguous"
        co, ct = classify(ro), classify(rt)
        agree = "yes" if co == ct else "NO"
        tot["ont"][co] += L; tot["10x"][ct] += L
        rows.append((short, L, frac, pid, do, ro, co, dt, rt, ct, agree))

    # ---- compleasm over/under-purge -----------------------------------------
    cmpl = {k: parse_compleasm_summary(v) for k, v in
            (("prepped", a.cmpl_prepped), ("purged", a.cmpl_purged), ("hap", a.cmpl_hap)) if v}
    ft = {k: parse_full_table(v) for k, v in
          (("prepped", a.ft_prepped), ("purged", a.ft_purged), ("hap", a.ft_hap)) if v}
    overpurge_ids = []
    PRESENT = ("Single", "Complete", "Duplicated")
    if "hap" in ft and "purged" in ft:
        for bid, st in ft["hap"].items():
            if st in PRESENT and ft["purged"].get(bid, "Missing") not in PRESENT:
                overpurge_ids.append(bid)

    # ---- write report -------------------------------------------------------
    with open(a.out_report, "w") as o:
        o.write(f"################ purge_dups VALIDATION -- {a.iso} ################\n\n")
        o.write(f"retained (diploid) depth mode:  ONT={D_ont:.1f}x   10x={D_10x:.1f}x\n")
        o.write(f"  -> a genuine haplotypic dup should read ~half of these ({D_ont/2:.0f}x / {D_10x/2:.0f}x)\n")
        o.write(f"hap contigs: {len(hap)}   purged contigs: {len(purged)}\n\n")

        o.write("== (1) compleasm completeness (single/dup/frag/incomplete/missing) ==\n")
        for k in ("prepped", "purged", "hap"):
            if k in cmpl:
                c = cmpl[k]
                cell = lambda x: f"{c.get(x,(0,0))[1]}" + (f"({c[x][0]:.1f}%)" if x in c and c[x][0] is not None else "")
                o.write(f"  {k:8}: S={cell('S')} D={cell('D')} F={cell('F')} I={cell('I')} M={cell('M')} N={c.get('N',(None,'?'))[1]}\n")
        o.write(f"\n  OVER-PURGE markers (Complete in hap.fa but Missing/Fragmented in purged.fa): {len(overpurge_ids)}\n")
        if overpurge_ids:
            o.write("    " + " ".join(overpurge_ids[:40]) + (" ..." if len(overpurge_ids) > 40 else "") + "\n")
        if "purged" in cmpl:
            dpct = cmpl["purged"].get("D", (0, 0))[0]
            o.write(f"  UNDER-PURGE signal (Duplicated% still in purged.fa): {dpct:.2f}%\n")

        o.write("\n== (3) hap-contig classification (redundancy + independent ONT & 10x depth) ==\n")
        o.write(f"{'hap_contig':<22}{'Mb':>7}{'aln%':>7}{'id%':>7}{'ONTx':>7}{'r':>6}{'ONTcls':>18}"
                f"{'10xx':>7}{'r':>6}{'10xcls':>18}{'agree':>6}\n")
        for (c, L, frac, pid, do, ro, co, dt, rt, ct, ag) in rows:
            o.write(f"{c:<22}{L/1e6:7.2f}{100*frac:7.1f}{pid:7.1f}{do:7.1f}{ro:6.2f}{co:>18}"
                    f"{dt:7.1f}{rt:6.2f}{ct:>18}{ag:>6}\n")

        o.write("\n== hap.fa sequence budget by class (Mb) ==\n")
        classes = sorted(set(list(tot["ont"]) + list(tot["10x"])))
        o.write(f"{'class':<22}{'ONT_Mb':>10}{'10x_Mb':>10}\n")
        for cl in classes:
            o.write(f"{cl:<22}{tot['ont'].get(cl,0)/1e6:10.2f}{tot['10x'].get(cl,0)/1e6:10.2f}\n")

        gs = ""
        if a.genomescope:
            try:
                gs = open(a.genomescope).read()
            except Exception:
                gs = "(genomescope summary unavailable)"
        if gs:
            o.write("\n== GenomeScope2 (10x k-mers) genome-size expectation ==\n" + gs + "\n")
        if a.seqkit_stats:
            try:
                o.write("\n== seqkit stats (prepped -> purged) ==\n" + open(a.seqkit_stats).read() + "\n")
            except Exception:
                pass

        # ---- verdict ----
        uniq = tot["ont"].get("UNIQUE(over-purge)", 0)
        redn = tot["ont"].get("haplotypic-dup", 0)
        o.write("\n================ VERDICT ================\n")
        o.write(f"hap.fa redundant (haplotypic): {redn/1e6:.1f} Mb   |   UNIQUE (over-purge risk): {uniq/1e6:.1f} Mb   "
                f"(ONT); over-purge markers: {len(overpurge_ids)}\n")
        if uniq > 20e6 or len(overpurge_ids) > 20:
            o.write("  >>> OVER-PURGE LIKELY: substantial unique/full-depth sequence or complete markers were removed.\n")
        elif redn > 0 and uniq < 5e6 and len(overpurge_ids) <= 5:
            o.write("  >>> PURGE LOOKS VALID: removed sequence is redundant, half-depth, and adds ~no complete markers.\n")
        else:
            o.write("  >>> INCONCLUSIVE: inspect the ambiguous rows and both coverage signals by hand.\n")

    # ---- plots --------------------------------------------------------------
    try:
        fig, ax = plt.subplots(1, 3, figsize=(19, 5.2))
        # panel 0/1: depth histograms (purged vs hap), ONT then 10x
        for k, (name, reg, D) in enumerate([("ONT", ont_reg, D_ont), ("10x", tenx_reg, D_10x)]):
            axx = ax[k]
            pv = [d for c, w, d in reg if c in purged and d > 0]
            hv = [d for c, w, d in reg if c in hap and d > 0]
            hi = np.nanpercentile(np.array(pv + hv) if (pv or hv) else np.array([1.0]), 99)
            bins = np.linspace(0, max(hi, 1), 80)
            if pv: axx.hist(pv, bins=bins, color="#4c72b0", alpha=0.6, label="purged (retained)")
            if hv: axx.hist(hv, bins=bins, color="#dd8452", alpha=0.6, label="hap (removed)")
            if D == D:
                axx.axvline(D, color="k", ls="--", lw=1); axx.axvline(D/2, color="grey", ls=":", lw=1)
                axx.text(D, axx.get_ylim()[1]*0.9, f"full {D:.0f}x", rotation=90, va="top", fontsize=8)
                axx.text(D/2, axx.get_ylim()[1]*0.9, f"half {D/2:.0f}x", rotation=90, va="top", fontsize=8, color="grey")
            axx.set_title(f"{a.iso} {name} depth: purged vs hap windows")
            axx.set_xlabel(f"{name} depth"); axx.set_ylabel("windows"); axx.legend(fontsize=8)
        # panel 2: hap-contig scatter (redundancy vs ONT depth ratio)
        axx = ax[2]
        colmap = {"haplotypic-dup": "#55a868", "UNIQUE(over-purge)": "#c44e52",
                  "ambiguous": "#8172b3", "lowcov/junk": "#999999", "?": "#cccccc"}
        for (c, L, frac, pid, do, ro, co, dt, rt, ct, ag) in rows:
            if ro != ro:
                continue
            axx.scatter(100*frac, ro, s=max(8, min(400, L/2e5)), c=colmap.get(co, "#cccccc"),
                        alpha=0.7, edgecolor="k", linewidth=0.3)
        axx.axhline(0.75, color="grey", ls=":"); axx.axvline(100*a.min_frac, color="grey", ls=":")
        axx.set_xlabel("aligned fraction to purged (%)"); axx.set_ylabel("ONT depth / retained")
        axx.set_title(f"{a.iso} hap contigs (size = length)\ngreen=haplotypic-dup  red=UNIQUE  purple=ambiguous")
        axx.set_ylim(0, max(1.6, axx.get_ylim()[1]))
        plt.tight_layout(); plt.savefig(a.out_plot, dpi=140); plt.close()
    except Exception as e:
        fig, ax = plt.subplots(figsize=(8, 4))
        ax.text(0.5, 0.5, f"purge_validate plot failed:\n{e}", ha="center", va="center", wrap=True)
        ax.axis("off"); plt.savefig(a.out_plot, dpi=140); plt.close()
        sys.stderr.write(f"purge_validate plot ERROR: {e}\n")

    sys.stderr.write(f"purge_validate {a.iso}: D_ont={D_ont:.1f} D_10x={D_10x:.1f} "
                     f"hap={len(hap)} overpurge_markers={len(overpurge_ids)} -> {a.out_report}\n")


if __name__ == "__main__":
    main()
