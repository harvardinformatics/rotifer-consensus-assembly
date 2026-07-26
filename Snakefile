# ============================================================================
# autobuild -- reproducible M. quadricornifera pangenome-ready reference build.
#
# INPUT  : MA/MM primary assemblies + 10x linked reads + raw ONT fastq_pass
# OUTPUT : single cross-strain consensus reference (decontaminated) + QC
#
# FOUR PHASES, human break between each (see README "Running in pieces"):
#   snakemake phase0   # contamination screen on the primaries -> review, make {iso}.remove.txt
#   snakemake phaseA   # prep(decontam) -> purge -> ARKS scaffold -> QC on scaffolds
#   snakemake phaseB   # apply cuts -> re-QC -> cross-isolate correspondence
#   snakemake          # Phase C: build consensus -> contam net -> final QC
#
# Software is provisioned by Snakemake, NOT hard-coded paths:
#   conda envs   (envs/*.yaml)   -- add --sw-deployment-method conda
#   a container  (FCS-GX)        -- add --sw-deployment-method apptainer  (for phase0 FCS-GX)
# e.g.  snakemake --sw-deployment-method conda apptainer \
#         --executor slurm -j 200 \
#         --default-resources slurm_account=informatics slurm_partition=sapphire <target>
# All tunables live in config.yaml.
# ============================================================================
import os, glob

configfile: "config.yaml"

ISO = config["isolates"]
WD  = config["workdir"].rstrip("/")
T   = config["threads"]
SD  = os.path.join(workflow.basedir, "scripts")

# ---- post-purge QC (phaseA_qc) constants -----------------------------------
CMPL_LIN = config["purge_qc"]["compleasm_lineage"]           # e.g. metazoa_odb10 (compleasm run -l)
CMPL_DL  = CMPL_LIN.replace("_odb10", "")                    # e.g. metazoa (compleasm download)
CMPL_LIB = config["purge_qc"]["compleasm_libdir"].rstrip("/")
MK       = config["purge_qc"]["meryl_k"]

def purgeqc_fa(iso, which):
    return {"prepped": f"{WD}/{iso}/prep/{iso}.prepped.fa",
            "purged":  f"{WD}/{iso}/purge/{iso}.purged.fa",
            "hap":     f"{WD}/{iso}/purge/hap.fa"}[which]

# nt ships as many volumes; the Phase-0 blast is parallelized by DB VOLUME (search
# each volume-chunk against the full window set in parallel -> nt scanned once
# total). Volumes are discovered at parse time from the nt db path.
NTDIR     = os.path.dirname(config["contam"]["nt_db"])
_NTVOLS   = sorted(f[:-4] for f in glob.glob(config["contam"]["nt_db"] + ".*.nin"))
_VPJ      = config["contam"]["blast_volumes_per_job"]
NT_CHUNKS = [_NTVOLS[i:i + _VPJ] for i in range(0, len(_NTVOLS), _VPJ)]

wildcard_constraints:
    iso   = "|".join(ISO),
    stage = "scaffold|cut",
    which = "prepped|purged|hap",
    reads = "ont|tenx",

localrules: write_cuts, write_consensus_spec

def asm(stage):
    return f"{WD}/{{iso}}/{stage}/{{iso}}.{stage}.fa"

# ---------------------------------------------------------------------------
# Phase targets (run in pieces; edit config at each break)
# ---------------------------------------------------------------------------
rule phase0:
    input:
        expand(f"{WD}/{{iso}}/prep/{{iso}}.contam_report.tsv",     iso=ISO),
        expand(f"{WD}/{{iso}}/prep/{{iso}}.remove_candidates.txt", iso=ISO),
        expand(f"{WD}/{{iso}}/prep/{{iso}}.blob.png",              iso=ISO),

rule phaseA:   # decontam + purge_dups, then STOP -- review the purge before heavy ARKS scaffolding
    input:
        expand(f"{WD}/{{iso}}/purge/{{iso}}.purge_stats.txt", iso=ISO),
        expand(f"{WD}/{{iso}}/purge/{{iso}}.purge_hist.png",  iso=ISO),

rule phaseA_qc:   # post-purge VALIDATION -- is hap.fa really redundant haplotypic seq?
    input:
        expand(f"{WD}/{{iso}}/purge_qc/{{iso}}.purge_validation.txt", iso=ISO),
        expand(f"{WD}/{{iso}}/purge_qc/{{iso}}.purge_validation.png", iso=ISO),

rule phaseB:   # ARKS scaffold + scaffold QC + cross-isolate synteny
    input:
        expand(f"{WD}/{{iso}}/qc_scaffold/cov.{{iso}}.txt", iso=ISO),
        expand(f"{WD}/{{iso}}/qc_scaffold/join_qc.txt",     iso=ISO),
        f"{WD}/synteny_scaffold/correspondence.tsv",
        f"{WD}/synteny_scaffold/MAxMM_dotplot.png",

rule phaseC:   # de-chimerize (apply cut_sites) + re-QC + correspondence
    input:
        expand(f"{WD}/{{iso}}/qc_cut/cov.{{iso}}.txt", iso=ISO),
        expand(f"{WD}/{{iso}}/qc_cut/join_qc.txt",     iso=ISO),
        f"{WD}/synteny_cut/correspondence.tsv",
        f"{WD}/synteny_cut/MAxMM_dotplot.png",

rule all:   # Phase C
    input:
        f"{WD}/consensus/MAMM_final.fa",
        f"{WD}/consensus/contam_report.tsv",
        f"{WD}/consensus/cov_qc.txt",
        f"{WD}/consensus/asm_stats.txt",
        expand(f"{WD}/{{iso}}/merqury/{{iso}}.qv",                    iso=ISO),
        expand(f"{WD}/{{iso}}/merqury/{{iso}}.kmer_completeness.txt", iso=ISO),

# ===========================================================================
# INPUT PREP
# ===========================================================================
# ONT reads are a FINALIZED input (>=10 kb, both runs pooled per isolate) supplied
# via config `ont_reads`; exact provenance is in source-data/ont-2024-filt/README.md.
# The pipeline consumes them directly -- no in-workflow filtering.

# 10x: R1 16 bp barcode -> BX:Z tag, interleave ALL lanes (temp: only ARKS consumes it)
rule prep_10x:
    input:
        r1 = lambda wc: config["tenx_R1"][wc.iso],
        r2 = lambda wc: config["tenx_R2"][wc.iso],
    output: bx = temp(f"{WD}/{{iso}}/scaffold/reads.bx.fq.gz")
    params: bl = config["barcode_len"]
    conda: "envs/bioinf.yaml"
    threads: 16
    resources: mem_mb=8000, runtime=720   # pigz-parallel; single-thread gzip timed out at 5h
    shell:
        r"""
        mkdir -p $(dirname {output.bx})
        paste <(pigz -dc {input.r1} | paste - - - -) <(pigz -dc {input.r2} | paste - - - -) \
        | awk -F'\t' -v BL={params.bl} 'BEGIN{{OFS="\n"}} {{
              bx=substr($2,1,BL);
              h1=$1; sub(/^@/,"",h1); sub(/ .*/,"",h1);
              h2=$5; sub(/^@/,"",h2); sub(/ .*/,"",h2);
              print "@"h1" BX:Z:"bx"-1", substr($2,BL+1), "+", substr($4,BL+1);
              print "@"h2" BX:Z:"bx"-1", $6, "+", $8;
          }}' | pigz -p {threads} -c > {output.bx}
        """

# ===========================================================================
# PHASE 0: prep reference -- contamination screen on the PRIMARIES
# ===========================================================================
rule prep_cov:
    input:
        pri = lambda wc: config["primary"][wc.iso],
        ont = lambda wc: config["ont_reads"][wc.iso],
    output: cov = f"{WD}/{{iso}}/prep/cov.txt"
    params: pre = config["ont_map_preset"], d = f"{WD}/{{iso}}/prep"
    conda: "envs/bioinf.yaml"
    threads: T
    resources: mem_mb=128000, runtime=600   # ONT map|sort peaks ~38G; 40G OOM-killed both isolates (2026-07-24)
    shell:
        r"""
        mkdir -p {params.d}
        minimap2 -ax {params.pre} -t {threads} {input.pri} {input.ont} 2>{params.d}/cov.mm.log \
          | samtools sort -@8 -m3G -o {params.d}/{wildcards.iso}.primary.ont.bam -
        samtools index {params.d}/{wildcards.iso}.primary.ont.bam
        samtools coverage {params.d}/{wildcards.iso}.primary.ont.bam > {output.cov}
        # BAM kept (was rm -f'd) -- ONT vs raw primary, useful for inspection
        """

# FCS-GX via the fcs.py wrapper + Singularity image -- the NCBI-documented
# approach (FCS-GX quickstart) and what dkhost's screen.slurm used. fcs.py drives
# the container itself, so NO Snakemake `container:` directive. Compute nodes have
# internet here, so tool/db fetch run as ordinary jobs.
rule fcsgx_tools:
    output:
        py  = config["prep"]["fcsgx_tooldir"].rstrip("/") + "/fcs.py",
        sif = config["prep"]["fcsgx_tooldir"].rstrip("/") + "/fcs-gx.sif",
    params: py_url = config["prep"]["fcs_py_url"], sif_url = config["prep"]["fcsgx_sif_url"]
    resources: mem_mb=4000, runtime=120
    shell:
        r"""
        mkdir -p $(dirname {output.py})
        curl -Ls {params.py_url}  -o {output.py}
        curl -Ls {params.sif_url} -o {output.sif}
        """

rule fetch_gxdb:   # fcs.py db get --mft <S3 manifest> --dir <gxdb>  (~470 GB, one-time)
    input:  py = rules.fcsgx_tools.output.py, sif = rules.fcsgx_tools.output.sif
    output: touch(config["prep"]["fcsgx_db_dir"].rstrip("/") + "/.synced")
    params: db = config["prep"]["fcsgx_db_dir"], mft = config["prep"]["fcsgx_manifest"]
    conda: "envs/bioinf.yaml"
    resources: mem_mb=16000, runtime=1440
    shell:
        r"""
        mkdir -p {params.db}
        FCS_DEFAULT_IMAGE={input.sif} python3 {input.py} db get --mft "{params.mft}" --dir {params.db}
        """

def fcsgx_inputs(wc):
    d = {"pri": config["primary"][wc.iso]}
    if config["prep"]["run_fcsgx"]:
        d["gxdb"] = config["prep"]["fcsgx_db_dir"].rstrip("/") + "/.synced"
        d["sif"]  = rules.fcsgx_tools.output.sif
        d["py"]   = rules.fcsgx_tools.output.py
    return d

# FRESH FCS-GX with a CURRENT db (the original report used a 2023-01-24 db).
rule prep_fcsgx:
    input: unpack(fcsgx_inputs)
    output: rpt = f"{WD}/{{iso}}/prep/fcsgx_report.txt"
    params:
        run = int(bool(config["prep"]["run_fcsgx"])),
        db = config["prep"]["fcsgx_db_dir"], taxid = config["prep"]["fcsgx_taxid"],
        sif = config["prep"]["fcsgx_tooldir"].rstrip("/") + "/fcs-gx.sif",
        py  = config["prep"]["fcsgx_tooldir"].rstrip("/") + "/fcs.py",
        d = f"{WD}/{{iso}}/prep/fcsgx",
    conda: "envs/bioinf.yaml"
    threads: T
    resources:   # FCS-GX loads the db to RAM (NCBI recommends ~512 GB); trivial when off
        mem_mb  = lambda wc: 512000 if config["prep"]["run_fcsgx"] else 2000,
        runtime = lambda wc: 720 if config["prep"]["run_fcsgx"] else 10,
    shell:
        r"""
        mkdir -p {params.d}
        if [ "{params.run}" = "1" ]; then
            FCS_DEFAULT_IMAGE={params.sif} python3 {params.py} screen genome \
                --fasta {input.pri} --out-dir {params.d} --gx-db {params.db} --tax-id {params.taxid}
            cp {params.d}/*.{params.taxid}.fcs_gx_report.txt {output.rpt}
        else
            echo "#FCS-GX not run (prep.run_fcsgx=false)." > {output.rpt}
        fi
        """

# Blobtools-style taxonomy. Window the GC-flagged contigs (prep_windows), then
# blastn the SAME window set against each nt VOLUME-CHUNK in parallel (blast_chunk;
# DB split, not query split -> nt scanned once total, spread across jobs), and
# merge best hits (blast_merge). See README for why DB-split beats query-batching.
rule prep_windows:
    input:  pri = lambda wc: config["primary"][wc.iso]
    output:
        gc  = f"{WD}/{{iso}}/prep/gc.tsv",
        win = f"{WD}/{{iso}}/prep/win.fa",
    params:
        d = f"{WD}/{{iso}}/prep", gcflag = config["prep"]["blob_gc_flag"],
        nwin = config["prep"]["blob_windows_per_contig"], wbp = config["prep"]["blob_window_bp"],
    conda: "envs/bioinf.yaml"
    threads: 8
    resources: mem_mb=8000, runtime=120
    shell:
        r"""
        mkdir -p {params.d}
        samtools faidx {input.pri}
        seqkit fx2tab -j {threads} -n -l -g {input.pri} > {output.gc}    # name length GC%
        awk -v G={params.gcflag} '$3/100 > G {{print $1}}' {output.gc} > {params.d}/blast_targets.txt
        : > {params.d}/regions.txt
        while read c; do
            L=$(awk -v c="$c" '$1==c{{print $2}}' {input.pri}.fai)
            python {SD}/make_windows.py "$c" "$L" {params.nwin} {params.wbp} >> {params.d}/regions.txt
        done < {params.d}/blast_targets.txt
        if [ -s {params.d}/regions.txt ]; then samtools faidx -r {params.d}/regions.txt {input.pri} > {output.win}
        else : > {output.win}; fi
        """

rule blast_chunk:   # one nt volume-chunk; -dbsize keeps per-volume e-values on the full-nt scale
    input:  win = f"{WD}/{{iso}}/prep/win.fa"
    output: tsv = f"{WD}/{{iso}}/prep/blast/chunk_{{c}}.tsv"
    params:
        vols = lambda wc: " ".join(NT_CHUNKS[int(wc.c)]),
        ntdir = NTDIR, dbsize = config["contam"]["nt_dbsize"],
    conda: "envs/blast.yaml"
    threads: 8
    resources: mem_mb=16000, runtime=240
    shell:
        r"""
        mkdir -p $(dirname {output.tsv})
        FMT="6 qseqid pident length evalue bitscore staxids sscinames sblastnames sskingdoms stitle"
        if [ -s {input.win} ]; then
            BLASTDB={params.ntdir} blastn -task megablast -db "{params.vols}" -query {input.win} \
                -dbsize {params.dbsize} -num_threads {threads} -max_target_seqs 5 -evalue 1e-10 \
                -outfmt "$FMT" > {output.tsv}
        else : > {output.tsv}; fi
        """

rule blast_merge:   # gather chunks; sort so the best (max-bitscore) hit is first per window
    input: lambda wc: expand(f"{WD}/{wc.iso}/prep/blast/chunk_{{c}}.tsv", c=range(len(NT_CHUNKS)))
    output: tsv = f"{WD}/{{iso}}/prep/blastn.tsv"
    resources: mem_mb=8000, runtime=60
    shell: "cat {input} | sort -t$'\t' -k1,1 -k5,5nr > {output.tsv}"

rule prep_report:   # combine GC + coverage + taxonomy (+ FCS-GX) -> report + candidates + plots
    input:
        gc     = f"{WD}/{{iso}}/prep/gc.tsv",
        cov    = f"{WD}/{{iso}}/prep/cov.txt",
        fcs    = f"{WD}/{{iso}}/prep/fcsgx_report.txt",
        blastn = f"{WD}/{{iso}}/prep/blastn.tsv",
    output:
        report = f"{WD}/{{iso}}/prep/{{iso}}.contam_report.tsv",
        cand   = f"{WD}/{{iso}}/prep/{{iso}}.remove_candidates.txt",
        blob   = f"{WD}/{{iso}}/prep/{{iso}}.blob.png",
    params: d = f"{WD}/{{iso}}/prep", gcflag = config["prep"]["blob_gc_flag"]
    conda: "envs/blast.yaml"
    resources: mem_mb=16000, runtime=120
    shell:
        r"""
        python {SD}/blob_report.py --gc {input.gc} --cov {input.cov} --fcsgx {input.fcs} \
            --blastn {input.blastn} --gcflag {params.gcflag} \
            --report {output.report} --candidates {output.cand} --plotprefix {params.d}/{wildcards.iso}
        """

# start of Phase A. Two decontamination actions on the primary, each driven by a curated
# file the user materializes at the Phase-0 manual gate (both live in {iso}/prep/):
#   (1) DROP whole contigs listed in {iso}.remove.txt (confident FCS EXCLUDEs, one id/line).
#   (2) TRIM terminal foreign spans in {iso}.trim.rpt (a TRIM-only FCS action report) -- a
#       contaminant at a contig END is not reliably endogenous, so we cut it.
# FIX (internal) and REVIEW are absent from both files and LEFT WHOLE, never masked, because
# an internal foreign span in a bdelloid is plausibly real HGT. TRIM is applied by
# `fcs.py clean genome` (native -- no container/db -- per the NCBI quickstart); the rule
# re-filters trim.rpt to TRIM rows only, so a stray FIX/EXCLUDE line could never mask/drop
# sequence here. Empty files => primary copied through unchanged.
rule prep_reference:
    input:
        pri     = lambda wc: config["primary"][wc.iso],
        rmlist  = f"{WD}/{{iso}}/prep/{{iso}}.remove.txt",
        trimrpt = f"{WD}/{{iso}}/prep/{{iso}}.trim.rpt",
    output: fa = f"{WD}/{{iso}}/prep/{{iso}}.prepped.fa"
    params:
        py      = config["prep"]["fcsgx_tooldir"].rstrip("/") + "/fcs.py",
        eff     = f"{WD}/{{iso}}/prep/{{iso}}.trim_applied.rpt",
        kept    = f"{WD}/{{iso}}/prep/{{iso}}.kept.fa",
        trimmed = f"{WD}/{{iso}}/prep/{{iso}}.trimmed_ends.fa",
    conda: "envs/bioinf.yaml"
    resources: mem_mb=8000, runtime=120
    shell:
        r"""
        # (1) drop the curated whole-contig EXCLUDEs
        if [ -s {input.rmlist} ]; then seqkit grep -v -f {input.rmlist} {input.pri} -o {params.kept}
        else cp {input.pri} {params.kept}; fi
        # (2) apply ONLY the TRIM rows from the curated trim report (FIX/EXCLUDE/REVIEW ignored)
        grep -E '^#' {input.trimrpt} > {params.eff} || true
        awk -F'\t' '!/^#/ && $5=="TRIM"' {input.trimrpt} >> {params.eff} || true
        if [ "$(awk -F'\t' '!/^#/ && $5=="TRIM"' {params.eff} | wc -l)" -gt 0 ]; then
            cat {params.kept} | python3 {params.py} clean genome \
                --action-report {params.eff} --output {output.fa} \
                --contam-fasta-out {params.trimmed}
            rm -f {params.kept}
        else
            mv {params.kept} {output.fa}
        fi
        """

# ===========================================================================
# PHASE A: purge_dups -> ARKS scaffold  (on the DECONTAMINATED primary)
# ===========================================================================
rule purge_dups:
    input:
        pri = f"{WD}/{{iso}}/prep/{{iso}}.prepped.fa",
        ont = lambda wc: config["ont_reads"][wc.iso],
    output:
        purged = f"{WD}/{{iso}}/purge/{{iso}}.purged.fa",
        hap    = f"{WD}/{{iso}}/purge/hap.fa",      # REMOVED haplotigs -- tracked & kept for phaseA_qc
        stat   = f"{WD}/{{iso}}/purge/PB.stat",     # depth histogram (for purge_qc review)
        cuts   = f"{WD}/{{iso}}/purge/cutoffs",     # calcuts thresholds actually used
    params:
        d = f"{WD}/{{iso}}/purge", map_pre = config["ont_map_preset"],
        self_pre = config["purge"]["self_preset"], self_fl = config["purge"]["self_flags"],
        calcuts = lambda wc: config["purge"]["calcuts_opts"].get(wc.iso, ""),
        getseqs = config["purge"]["get_seqs_opts"],
    conda: "envs/purge_dups.yaml"
    threads: T
    resources: mem_mb=48000, runtime=480
    shell:
        r"""
        mkdir -p {params.d} && cd {params.d}
        minimap2 -x {params.map_pre} -t {threads} {input.pri} {input.ont} 2>mm.log | gzip -c > reads.paf.gz
        pbcstat reads.paf.gz
        calcuts {params.calcuts} PB.stat > cutoffs 2> calcuts.log
        split_fa {input.pri} > pri.split
        minimap2 -x {params.self_pre} {params.self_fl} -t {threads} pri.split pri.split 2>self.log | gzip -c > pri.split.self.paf.gz
        purge_dups -2 -T cutoffs -c PB.base.cov pri.split.self.paf.gz > dups.bed 2> pd.log
        get_seqs {params.getseqs} dups.bed {input.pri}
        test -s purged.fa && mv purged.fa {output.purged}
        rm -f reads.paf.gz pri.split pri.split.self.paf.gz    # large intermediates
        """

# Phase-A review gate: summarize the purge so a human can sanity-check calcuts BEFORE the
# heavy ARKS scaffolding. purge_hist.py overlays the calcuts cutoffs on the ONT depth
# histogram (the key check for MM's possibly-bimodal / degenerate-tetraploid depth); the
# stats file gives the cutoffs + before/after seqkit stats (how much sequence was purged).
rule purge_qc:
    input:
        prepped = f"{WD}/{{iso}}/prep/{{iso}}.prepped.fa",
        purged  = f"{WD}/{{iso}}/purge/{{iso}}.purged.fa",
        stat    = f"{WD}/{{iso}}/purge/PB.stat",
        cuts    = f"{WD}/{{iso}}/purge/cutoffs",
    output:
        stats = f"{WD}/{{iso}}/purge/{{iso}}.purge_stats.txt",
        hist  = f"{WD}/{{iso}}/purge/{{iso}}.purge_hist.png",
    conda: "envs/bioinf.yaml"
    resources: mem_mb=4000, runtime=30
    shell:
        r"""
        echo "## calcuts cutoffs (haploid/diploid depth thresholds chosen by calcuts):" > {output.stats}
        cat {input.cuts} >> {output.stats}
        printf '\n## seqkit stats -- prepped (pre-purge) then purged (post-purge):\n' >> {output.stats}
        seqkit stats -aT {input.prepped} {input.purged} >> {output.stats}
        python {SD}/purge_hist.py --stat {input.stat} --cuts {input.cuts} \
            --title "{wildcards.iso} purge_dups ONT depth" --out {output.hist}
        """

# ===========================================================================
# PHASE A QC: is the purge VALID?  (target: phaseA_qc)
# Three independent lines of evidence that hap.fa is redundant haplotypic seq:
#   (1) compleasm markers on prepped/purged/hap  (over- & under-purge)
#   (2) hap->purged alignment                    (does each removed contig map back?)
#   (3) coverage on the combined purged+hap ref, ONT AND 10x independently
#       (a real haplotypic dup reads ~half the retained/diploid depth)
# ===========================================================================
rule compleasm_db:                    # one-time lineage download (compute nodes have internet)
    output: touch(CMPL_LIB + "/.synced")
    params: dl=CMPL_DL, lib=CMPL_LIB
    conda: "envs/compleasm.yaml"
    resources: mem_mb=8000, runtime=180
    shell:
        r"""
        # compleasm 0.2.6 is broken against the CURRENT ezlab data two ways, both inside
        # Downloader -- whose __init__ runs UNCONDITIONALLY, including from `compleasm run`:
        #   (a) download_placement() parses names via `strain.split(".")` expecting 3 (or 4)
        #       fields, but names like 'x_odb10.DATE.tar.gz' now split into more -> ValueError.
        #       Patch the split to a bounded maxsplit (fixes the download AND .done-skip branch).
        #   (b) the placement-file URLs now 404. Placement files are only for --autolineage (we
        #       pass -l explicitly), so touch placement_files.done to force the no-download branch.
        # Together (+ clearing any stale .tmp lock) both `download` AND `run` succeed. (2026-07-26)
        CP=$(ls $CONDA_PREFIX/lib/python*/site-packages/compleasm.py 2>/dev/null | head -1 || true)
        if [ -n "$CP" ]; then
          sed -i \
            -e 's/prefix, version, sufix = strain\.split("\.")/prefix, version, sufix = strain.split(".", 2)/g' \
            -e 's/prefix, aln, version, sufix = strain\.split("\.")/prefix, aln, version, sufix = strain.split(".", 3)/g' \
            "$CP"
        fi
        mkdir -p {params.lib}/placement_files
        rm -f {params.lib}/placement_files.tmp
        touch {params.lib}/placement_files.done
        compleasm download {params.dl} --library_path {params.lib}
        """

rule compleasm:                       # {which} in {prepped,purged,hap}
    input:
        fa = lambda wc: purgeqc_fa(wc.iso, wc.which),
        db = CMPL_LIB + "/.synced",
    output:
        summ = f"{WD}/{{iso}}/purge_qc/compleasm.{{which}}.summary.txt",
        ft   = f"{WD}/{{iso}}/purge_qc/compleasm.{{which}}.full_table.tsv",
    params:
        d   = f"{WD}/{{iso}}/purge_qc/cmpl_{{which}}",
        lin = CMPL_LIN, lib = CMPL_LIB,
    conda: "envs/compleasm.yaml"
    threads: T
    resources: mem_mb=32000, runtime=300
    shell:
        r"""
        rm -rf {params.d} && mkdir -p $(dirname {output.summ})
        compleasm run -a {input.fa} -o {params.d} -t {threads} -l {params.lin} -L {params.lib}
        cp {params.d}/summary.txt {output.summ}
        cp "$(find {params.d} -name full_table.tsv | head -1)" {output.ft}
        """

rule purgeqc_ref:                     # combined purged + hap (hap contigs get a HAP__ prefix) + labels
    input:
        purged = f"{WD}/{{iso}}/purge/{{iso}}.purged.fa",
        hap    = f"{WD}/{{iso}}/purge/hap.fa",
    output:
        fa     = f"{WD}/{{iso}}/purge_qc/{{iso}}.combined.fa",
        labels = f"{WD}/{{iso}}/purge_qc/labels.tsv",
    conda: "envs/bioinf.yaml"
    resources: mem_mb=8000, runtime=60
    shell:
        r"""
        mkdir -p $(dirname {output.fa})
        seqkit replace -p '(.+)' -r 'HAP__$1' {input.hap} > {output.fa}.haptmp
        cat {input.purged} {output.fa}.haptmp > {output.fa}
        rm -f {output.fa}.haptmp
        samtools faidx {output.fa}
        seqkit fx2tab -n -l {input.purged} | awk 'BEGIN{{OFS="\t"}} {{print $1,"purged",$NF}}'   > {output.labels}
        seqkit fx2tab -n -l {input.hap}    | awk 'BEGIN{{OFS="\t"}} {{print "HAP__"$1,"hap",$NF}}' >> {output.labels}
        """

rule purgeqc_map:                     # {reads} in {ont,tenx}: map to combined ref, KEEP bam, mosdepth
    input:
        fa    = f"{WD}/{{iso}}/purge_qc/{{iso}}.combined.fa",
        reads = lambda wc: config["ont_reads"][wc.iso] if wc.reads == "ont" else config["tenx_R2"][wc.iso],
    output:
        bam  = f"{WD}/{{iso}}/purge_qc/{{iso}}.{{reads}}.bam",
        bai  = f"{WD}/{{iso}}/purge_qc/{{iso}}.{{reads}}.bam.bai",
        summ = f"{WD}/{{iso}}/purge_qc/{{iso}}.{{reads}}.mosdepth.summary.txt",
        reg  = f"{WD}/{{iso}}/purge_qc/{{iso}}.{{reads}}.regions.bed.gz",
    params:
        pre = lambda wc: "map-ont" if wc.reads == "ont" else "sr",
        w   = config["purge_qc"]["cov_window_bp"],
        pfx = f"{WD}/{{iso}}/purge_qc/{{iso}}.{{reads}}",
    conda: "envs/bioinf.yaml"
    threads: T
    resources: mem_mb=128000, runtime=720   # same ONT map|sort OOM headroom as map_ont
    shell:
        r"""
        minimap2 -ax {params.pre} -t {threads} {input.fa} {input.reads} 2>{params.pfx}.mm.log \
          | samtools sort -@8 -m3G -T {params.pfx}.srt -o {output.bam} -
        samtools index -@ {threads} {output.bam}
        mosdepth -t {threads} --by {params.w} --no-per-base {params.pfx} {output.bam}
        """

rule hap_aln:                         # hap (query) vs purged (target); -c --cs -> real %id + aligned bp
    input:
        purged = f"{WD}/{{iso}}/purge/{{iso}}.purged.fa",
        hap    = f"{WD}/{{iso}}/purge/hap.fa",
    output: paf = f"{WD}/{{iso}}/purge_qc/hap_vs_purged.paf"
    params: pre = config["purge_qc"]["hap_aln_preset"]
    conda: "envs/bioinf.yaml"
    threads: T
    resources: mem_mb=32000, runtime=180
    shell:
        r"""
        mkdir -p $(dirname {output.paf})
        minimap2 -cx {params.pre} --cs -t {threads} {input.purged} {input.hap} > {output.paf} 2>{output.paf}.log
        """

rule meryl_hist:                      # 10x R2 k-mer histogram for GenomeScope2
    input: r2 = lambda wc: config["tenx_R2"][wc.iso]
    output: hist = f"{WD}/{{iso}}/purge_qc/{{iso}}.k{MK}.hist"
    params: k = MK, mem_gb = 60, d = f"{WD}/{{iso}}/purge_qc/{{iso}}.k{MK}.meryl"
    conda: "envs/meryl.yaml"
    threads: T
    resources: mem_mb=64000, runtime=480
    shell:
        r"""
        mkdir -p $(dirname {output.hist})
        rm -rf {params.d}
        meryl k={params.k} threads={threads} memory={params.mem_gb} count {input.r2} output {params.d}
        meryl histogram {params.d} > {output.hist}
        rm -rf {params.d}
        """

rule genomescope:                     # genome-size expectation from 10x k-mers (best-effort; never blocks)
    input: hist = f"{WD}/{{iso}}/purge_qc/{{iso}}.k{MK}.hist"
    output: summ = f"{WD}/{{iso}}/purge_qc/{{iso}}.genomescope.summary.txt"
    params:
        k = MK, p = config["purge_qc"]["genomescope_ploidy"],
        d = f"{WD}/{{iso}}/purge_qc/{{iso}}.genomescope",
    conda: "envs/genomescope.yaml"
    resources: mem_mb=8000, runtime=60
    shell:
        r"""
        mkdir -p {params.d}
        if genomescope2 -i {input.hist} -o {params.d} -k {params.k} -p {params.p} -n {wildcards.iso} > {params.d}/gs.log 2>&1; then
            cp {params.d}/{wildcards.iso}_summary.txt {output.summ}
        else
            echo "GenomeScope2 FAILED -- see {params.d}/gs.log" > {output.summ}
        fi
        """

rule purge_validate:                  # synthesize all three evidence lines -> report + plot
    input:
        labels       = f"{WD}/{{iso}}/purge_qc/labels.tsv",
        cmpl_prepped = f"{WD}/{{iso}}/purge_qc/compleasm.prepped.summary.txt",
        cmpl_purged  = f"{WD}/{{iso}}/purge_qc/compleasm.purged.summary.txt",
        cmpl_hap     = f"{WD}/{{iso}}/purge_qc/compleasm.hap.summary.txt",
        ft_prepped   = f"{WD}/{{iso}}/purge_qc/compleasm.prepped.full_table.tsv",
        ft_purged    = f"{WD}/{{iso}}/purge_qc/compleasm.purged.full_table.tsv",
        ft_hap       = f"{WD}/{{iso}}/purge_qc/compleasm.hap.full_table.tsv",
        hap_paf      = f"{WD}/{{iso}}/purge_qc/hap_vs_purged.paf",
        ont_summ     = f"{WD}/{{iso}}/purge_qc/{{iso}}.ont.mosdepth.summary.txt",
        ont_reg      = f"{WD}/{{iso}}/purge_qc/{{iso}}.ont.regions.bed.gz",
        tenx_summ    = f"{WD}/{{iso}}/purge_qc/{{iso}}.tenx.mosdepth.summary.txt",
        tenx_reg     = f"{WD}/{{iso}}/purge_qc/{{iso}}.tenx.regions.bed.gz",
        gs           = f"{WD}/{{iso}}/purge_qc/{{iso}}.genomescope.summary.txt",
        stats        = f"{WD}/{{iso}}/purge/{{iso}}.purge_stats.txt",
    output:
        report = f"{WD}/{{iso}}/purge_qc/{{iso}}.purge_validation.txt",
        plot   = f"{WD}/{{iso}}/purge_qc/{{iso}}.purge_validation.png",
    params: minfrac = config["purge_qc"]["redundancy_min_frac"]
    conda: "envs/bioinf.yaml"
    resources: mem_mb=16000, runtime=60
    shell:
        r"""
        python {SD}/purge_validate.py --iso {wildcards.iso} --labels {input.labels} \
          --cmpl-prepped {input.cmpl_prepped} --cmpl-purged {input.cmpl_purged} --cmpl-hap {input.cmpl_hap} \
          --ft-prepped {input.ft_prepped} --ft-purged {input.ft_purged} --ft-hap {input.ft_hap} \
          --hap-paf {input.hap_paf} \
          --ont-summary {input.ont_summ} --ont-regions {input.ont_reg} \
          --tenx-summary {input.tenx_summ} --tenx-regions {input.tenx_reg} \
          --genomescope {input.gs} --seqkit-stats {input.stats} \
          --min-frac {params.minfrac} \
          --out-report {output.report} --out-plot {output.plot}
        """

rule scaffold_arks:
    input:
        purged = f"{WD}/{{iso}}/purge/{{iso}}.purged.fa",
        bx     = f"{WD}/{{iso}}/scaffold/reads.bx.fq.gz",
    output: scaf = f"{WD}/{{iso}}/scaffold/{{iso}}.scaffold.fa"
    params: d = f"{WD}/{{iso}}/scaffold", ap = config["arks"]["params"], glob = config["arks"]["output_glob"]
    conda: "envs/arcs.yaml"
    threads: T
    resources: mem_mb=100000, runtime=720     # ARKS/LINKS memory scales with 10x barcodes; tune if OOM
    shell:
        r"""
        cd {params.d}
        ln -sf {input.purged} draft.fa
        ln -sf {input.bx}     reads.fq.gz
        arcs-make arks draft=draft reads=reads t={threads} {params.ap}
        n=$(ls {params.glob} 2>/dev/null | wc -l)
        [ "$n" -eq 1 ] || {{ echo "ERROR: expected 1 ARKS output matching '{params.glob}', found $n"; ls {params.glob}; exit 1; }}
        cp "$(ls {params.glob})" {output.scaf}
        """

# ===========================================================================
# PHASE B: de-chimerization cuts (config cut_sites)
# ===========================================================================
rule write_cuts:
    output: f"{WD}/{{iso}}/cut/cuts.tsv"
    params: cuts = lambda wc: config["cut_sites"].get(wc.iso, [])
    run:
        os.makedirs(os.path.dirname(output[0]), exist_ok=True)
        with open(output[0], "w") as o:
            for c in params.cuts:
                o.write(f"{c['scaffold']}\t{c['pos_mb']}\n")

rule split_cuts:
    input:
        scaf = f"{WD}/{{iso}}/scaffold/{{iso}}.scaffold.fa",
        cuts = f"{WD}/{{iso}}/cut/cuts.tsv",
    output: cut = f"{WD}/{{iso}}/cut/{{iso}}.cut.fa"
    params: snap = config["cut_snap_bp"]
    conda: "envs/bioinf.yaml"
    resources: mem_mb=8000, runtime=60
    shell: "python {SD}/apply_cuts.py {input.scaf} {input.cuts} {output.cut} {params.snap}"

# ===========================================================================
# QC (stage-parametrized: identical on 'scaffold' and 'cut')
# ===========================================================================
rule map_ont:
    input:
        fa  = lambda wc: asm(wc.stage).format(iso=wc.iso),
        ont = lambda wc: config["ont_reads"][wc.iso],
    output:
        bam = f"{WD}/{{iso}}/qc_{{stage}}/{{iso}}.{{stage}}.ont.bam",       # KEEP (not temp) -- ONT->assembly BAM
        bai = f"{WD}/{{iso}}/qc_{{stage}}/{{iso}}.{{stage}}.ont.bam.bai",
        cov = f"{WD}/{{iso}}/qc_{{stage}}/cov.{{iso}}.txt",
    params: pre = config["ont_map_preset"]
    conda: "envs/bioinf.yaml"
    threads: T
    resources: mem_mb=128000, runtime=600   # same ONT map|sort OOM as prep_cov (2026-07-24)
    shell:
        r"""
        mkdir -p $(dirname {output.bam})
        minimap2 -ax {params.pre} -t {threads} {input.fa} {input.ont} 2>{output.bam}.log \
          | samtools sort -@8 -m3G -o {output.bam} -
        samtools index {output.bam}
        samtools coverage {output.bam} > {output.cov}
        """

rule join_qc:
    input:
        fa  = lambda wc: asm(wc.stage).format(iso=wc.iso),
        bam = f"{WD}/{{iso}}/qc_{{stage}}/{{iso}}.{{stage}}.ont.bam",
    output: rpt = f"{WD}/{{iso}}/qc_{{stage}}/join_qc.txt"
    params:
        minscaf = config["qc"]["join_qc_min_scaffold_mb"],
        flank   = config["qc"]["join_qc_flank_bp"],
        minspan = config["qc"]["join_qc_min_span"],
    conda: "envs/bioinf.yaml"
    resources: mem_mb=8000, runtime=120
    shell: "python {SD}/join_qc.py {input.fa} {input.bam} {params.minscaf} {params.flank} {params.minspan} > {output.rpt}"

rule synteny:
    input:
        ma = f"{WD}/MA/{{stage}}/MA.{{stage}}.fa",
        mm = f"{WD}/MM/{{stage}}/MM.{{stage}}.fa",
    output:
        paf   = f"{WD}/synteny_{{stage}}/MAxMM.paf",
        png   = f"{WD}/synteny_{{stage}}/MAxMM_dotplot.png",
        tsv   = f"{WD}/synteny_{{stage}}/correspondence.tsv",
        edges = f"{WD}/synteny_{{stage}}/edges.tsv",       # per-homolog aln_bp + real %id
    params:
        preset = config["synteny"]["preset"], minaln = config["synteny"]["min_aln_bp"],
        minlen = config["synteny"]["dotplot_min_len_bp"], minid = config["synteny"]["min_id"],
    conda: "envs/bioinf.yaml"
    threads: T
    resources: mem_mb=64000, runtime=300
    shell:
        r"""
        mkdir -p $(dirname {output.paf})
        # BASE-LEVEL alignment (-c --cs) so the PAF carries de:f: (gap-compressed divergence)
        # -> real %id, not the chaining estimate. Order is `minimap2 MM MA` so query=MA (=X axis).
        minimap2 -cx {params.preset} --cs -t {threads} {input.mm} {input.ma} > {output.paf} 2>{output.paf}.log
        python {SD}/paf_dotplot.py {output.paf} {output.png} {params.minlen}
        python {SD}/synteny_classify.py {output.paf} {output.tsv} {output.edges} {params.minaln} {params.minid}
        """

# ===========================================================================
# PHASE C: consensus reference -> contamination net -> final QC
# ===========================================================================
rule write_consensus_spec:
    output: f"{WD}/consensus/consensus_spec.json"
    run:
        import json
        os.makedirs(os.path.dirname(output[0]), exist_ok=True)
        spec = dict(config["consensus"])
        spec["cut"] = {i: f"{WD}/{i}/cut/{i}.cut.fa" for i in ISO}
        with open(output[0], "w") as o:
            json.dump(spec, o, indent=2)

rule build_consensus:
    input:
        spec = f"{WD}/consensus/consensus_spec.json",
        cut  = expand(f"{WD}/{{iso}}/cut/{{iso}}.cut.fa", iso=ISO),
    output: fa = f"{WD}/consensus/MAMM_merged.fa"
    conda: "envs/bioinf.yaml"
    resources: mem_mb=8000, runtime=60
    shell: "python {SD}/build_consensus.py {input.spec} {output.fa}"

rule asm_stats:
    input:  f"{WD}/consensus/MAMM_final.fa"
    output: f"{WD}/consensus/asm_stats.txt"
    conda: "envs/bioinf.yaml"
    resources: mem_mb=4000, runtime=30
    shell:  "seqkit stats -a -T {input} > {output}"

# end-of-pipeline contamination NET on the consensus (second line after Phase 0)
rule contam_screen:
    input:  fa = f"{WD}/consensus/MAMM_merged.fa"
    output: rpt = f"{WD}/consensus/contam_report.tsv"
    params:
        gc = config["contam"]["gc_flag_high"], do_blast = int(bool(config["contam"]["run_blast"])),
        nt = config["contam"]["nt_db"], nr = config["contam"]["nr_db"],
        nwin = config["contam"]["windows_per_contig"], wbp = config["contam"]["window_bp"],
        d = f"{WD}/consensus/contam",
    conda: "envs/blast.yaml"
    threads: T
    resources: mem_mb=32000, runtime=600
    shell:
        r"""
        mkdir -p {params.d}
        samtools faidx {input.fa}
        seqkit fx2tab -j {threads} -n -l -g {input.fa} | awk -v G={params.gc} '$3/100 > G {{print $1}}' > {params.d}/flagged.txt
        rm -f {params.d}/blastn.tsv {params.d}/blastx.tsv
        if [ "{params.do_blast}" = "1" ] && [ -s {params.d}/flagged.txt ]; then
            : > {params.d}/regions.txt
            while read c; do
                L=$(awk -v c="$c" '$1==c{{print $2}}' {input.fa}.fai)
                python {SD}/make_windows.py "$c" "$L" {params.nwin} {params.wbp} >> {params.d}/regions.txt
            done < {params.d}/flagged.txt
            samtools faidx -r {params.d}/regions.txt {input.fa} > {params.d}/win.fa
            FMT="6 qseqid pident length evalue bitscore staxids sscinames sblastnames sskingdoms stitle"
            BLASTDB={params.nt} blastn -task megablast -db {params.nt} -query {params.d}/win.fa \
                -num_threads {threads} -max_target_seqs 5 -evalue 1e-10 -outfmt "$FMT" > {params.d}/blastn.tsv || true
            BLASTDB={params.nr} blastx -db {params.nr} -query {params.d}/win.fa \
                -num_threads {threads} -max_target_seqs 5 -evalue 1e-5 -outfmt "$FMT" > {params.d}/blastx.tsv || true
        fi
        python {SD}/contam_report.py {input.fa} {params.gc} {params.d} > {output.rpt}
        """

rule finalize_ref:
    input:  fa = f"{WD}/consensus/MAMM_merged.fa", rpt = f"{WD}/consensus/contam_report.tsv"
    output: fa = f"{WD}/consensus/MAMM_final.fa"
    params: drop = lambda wc: " ".join(config["contam"]["drop_contigs"])
    conda: "envs/bioinf.yaml"
    resources: mem_mb=4000, runtime=60
    shell:
        r"""
        if [ -n "{params.drop}" ]; then printf '%s\n' {params.drop} | seqkit grep -v -f - {input.fa} -o {output.fa}
        else cp {input.fa} {output.fa}; fi
        """

rule qc_coverage:
    input:
        fa  = f"{WD}/consensus/MAMM_final.fa",
        ont = [config["ont_reads"][i] for i in ISO],
    output: cov = f"{WD}/consensus/cov_qc.txt"
    params:
        d = f"{WD}/consensus/covqc", isos = " ".join(ISO),
        collapse = config["qc"]["cov_collapse_ratio"], lowcov = config["qc"]["cov_lowcov_frac"],
        tk = config["qc"]["telo_topkmer_frac"], tc = config["qc"]["telo_canonical_frac"],
        pre = config["ont_map_preset"],
    conda: "envs/bioinf.yaml"
    threads: T
    resources: mem_mb=40000, runtime=600
    shell:
        r"""
        mkdir -p {params.d}
        isos=({params.isos}); reads=({input.ont})
        for k in ${{!isos[@]}}; do
            X=${{isos[$k]}}; R=${{reads[$k]}}
            minimap2 -ax {params.pre} -t {threads} {input.fa} $R 2>{params.d}/$X.log \
              | samtools sort -@8 -m3G -o {params.d}/$X.consensus.ont.bam -
            samtools index {params.d}/$X.consensus.ont.bam
            samtools coverage {params.d}/$X.consensus.ont.bam > {params.d}/cov_$X.txt
            # BAM kept (was rm -f'd) -- ONT vs final consensus, per isolate
        done
        python {SD}/cov_qc.py {input.fa} {params.d}/cov_MA.txt {params.d}/cov_MM.txt \
             {params.collapse} {params.lowcov} {params.tk} {params.tc} > {output.cov}
        """

rule qc_merqury:
    input:
        fa = f"{WD}/{{iso}}/cut/{{iso}}.cut.fa",
        r2 = lambda wc: config["tenx_R2"][wc.iso],
    output:
        qv = f"{WD}/{{iso}}/merqury/{{iso}}.qv",
        cn = f"{WD}/{{iso}}/merqury/{{iso}}.{{iso}}.cut.spectra-cn.hist",
    params: d = f"{WD}/{{iso}}/merqury", k = config["qc"]["meryl_k"], name = "{iso}"
    conda: "envs/merqury.yaml"
    threads: T
    resources: mem_mb=48000, runtime=600
    shell:
        r"""
        mkdir -p {params.d}; cd {params.d}
        meryl k={params.k} threads={threads} count {input.r2} output {params.name}.meryl
        merqury.sh {params.name}.meryl {input.fa} {params.name} 2>&1 | tail -5
        """

rule qc_kmer_completeness:
    input:  cn = f"{WD}/{{iso}}/merqury/{{iso}}.{{iso}}.cut.spectra-cn.hist"
    output: f"{WD}/{{iso}}/merqury/{{iso}}.kmer_completeness.txt"
    params: filt = f"{WD}/{{iso}}/merqury/{{iso}}.filt"
    conda: "envs/bioinf.yaml"
    resources: mem_mb=4000, runtime=30
    shell:  "python {SD}/kmer_completeness.py {input.cn} $(cat {params.filt}) > {output}"
