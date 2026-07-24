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

rule phaseA:
    input:
        expand(f"{WD}/{{iso}}/qc_scaffold/cov.{{iso}}.txt", iso=ISO),
        expand(f"{WD}/{{iso}}/qc_scaffold/join_qc.txt",     iso=ISO),
        f"{WD}/synteny_scaffold/correspondence.tsv",
        f"{WD}/synteny_scaffold/MAxMM_dotplot.png",

rule phaseB:
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
    threads: 4
    resources: mem_mb=8000, runtime=300
    shell:
        r"""
        mkdir -p $(dirname {output.bx})
        paste <(zcat {input.r1} | paste - - - -) <(zcat {input.r2} | paste - - - -) \
        | awk -F'\t' -v BL={params.bl} 'BEGIN{{OFS="\n"}} {{
              bx=substr($2,1,BL);
              h1=$1; sub(/^@/,"",h1); sub(/ .*/,"",h1);
              h2=$5; sub(/^@/,"",h2); sub(/ .*/,"",h2);
              print "@"h1" BX:Z:"bx"-1", substr($2,BL+1), "+", substr($4,BL+1);
              print "@"h2" BX:Z:"bx"-1", $6, "+", $8;
          }}' | gzip -c > {output.bx}
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
          | samtools sort -@8 -m3G -o {params.d}/cov.bam -
        samtools index {params.d}/cov.bam
        samtools coverage {params.d}/cov.bam > {output.cov}
        rm -f {params.d}/cov.bam {params.d}/cov.bam.bai
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

# start of Phase A: drop the VERIFIED removal list from the primary (seqkit grep -v).
# {iso}.remove.txt is created by the user from Phase-0 candidates (empty = keep all).
rule prep_reference:
    input:
        pri     = lambda wc: config["primary"][wc.iso],
        rmlist  = f"{WD}/{{iso}}/prep/{{iso}}.remove.txt",
    output: fa = f"{WD}/{{iso}}/prep/{{iso}}.prepped.fa"
    conda: "envs/bioinf.yaml"
    resources: mem_mb=4000, runtime=60
    shell:
        r"""
        if [ -s {input.rmlist} ]; then seqkit grep -v -f {input.rmlist} {input.pri} -o {output.fa}
        else cp {input.pri} {output.fa}; fi
        """

# ===========================================================================
# PHASE A: purge_dups -> ARKS scaffold  (on the DECONTAMINATED primary)
# ===========================================================================
rule purge_dups:
    input:
        pri = f"{WD}/{{iso}}/prep/{{iso}}.prepped.fa",
        ont = lambda wc: config["ont_reads"][wc.iso],
    output: purged = f"{WD}/{{iso}}/purge/{{iso}}.purged.fa"
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
        bam = temp(f"{WD}/{{iso}}/qc_{{stage}}/{{iso}}.{{stage}}.ont.bam"),
        bai = temp(f"{WD}/{{iso}}/qc_{{stage}}/{{iso}}.{{stage}}.ont.bam.bai"),
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
        paf = f"{WD}/synteny_{{stage}}/MAxMM.paf",
        png = f"{WD}/synteny_{{stage}}/MAxMM_dotplot.png",
        tsv = f"{WD}/synteny_{{stage}}/correspondence.tsv",
    params:
        preset = config["synteny"]["preset"], minaln = config["synteny"]["min_aln_bp"],
        minlen = config["synteny"]["dotplot_min_len_bp"],
    conda: "envs/bioinf.yaml"
    threads: T
    resources: mem_mb=32000, runtime=300
    shell:
        r"""
        mkdir -p $(dirname {output.paf})
        minimap2 -x {params.preset} -t {threads} {input.ma} {input.mm} > {output.paf} 2>{output.paf}.log
        python {SD}/paf_dotplot.py {output.paf} {output.png} {params.minlen}
        python {SD}/synteny_classify.py {output.paf} {output.tsv} {params.minaln}
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
              | samtools sort -@8 -m3G -o {params.d}/$X.bam -
            samtools index {params.d}/$X.bam
            samtools coverage {params.d}/$X.bam > {params.d}/cov_$X.txt
            rm -f {params.d}/$X.bam {params.d}/$X.bam.bai
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
