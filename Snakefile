# ============================================================================
# Snakefile -- reproducible M. quadricornifera pangenome-ready reference build.
#
# INPUT  : MA/MM primary assemblies + 10x linked reads + raw ONT fastq_pass
# OUTPUT : single cross-strain consensus reference (contamination-cleaned) + QC
#
# THREE PHASES with a human break between each (see README "Running in pieces"):
#   snakemake phaseA   # filter -> purge -> scaffold -> QC on scaffolds
#                      #   >>> inspect QC, set config cut_sites <<<
#   snakemake phaseB   # apply cuts -> re-QC -> cross-isolate correspondence
#                      #   >>> inspect correspondence, set config consensus <<<
#   snakemake          # build consensus -> contamination screen -> final QC
#
# All tunables live in config.yaml; nothing tunable is hard-coded here.
# Dependencies flow through files (no SLURM afterok). SLURM (Snakemake >=8):
#   snakemake --executor slurm -j 200 \
#     --default-resources slurm_account=informatics slurm_partition=sapphire
# ============================================================================

configfile: "config.yaml"

ISO = config["isolates"]
WD  = config["workdir"].rstrip("/")
T   = config["threads"]
SD  = config["scriptdir"].rstrip("/")
CB  = config["conda_base"].rstrip("/")

MM2 = f"{CB}/envs/{config['env_map']}/bin/minimap2"
ST  = f"{CB}/envs/{config['env_map']}/bin/samtools"
PY  = f"{CB}/envs/{config['env_map']}/bin/python"
PDB = f"{CB}/envs/{config['env_linked']}/bin"
ARCS= f"{CB}/envs/{config['env_linked']}/bin/arcs-make"
MQ  = f"{CB}/envs/{config['env_merqury']}/bin"
BL  = f"{CB}/envs/{config['env_blast']}/bin"

wildcard_constraints:
    iso   = "|".join(ISO),
    stage = "scaffold|cut",

# housekeeping rules run in the Snakemake process (not submitted to SLURM)
localrules: write_cuts, write_consensus_spec, versions

def asm(stage):   # per-isolate assembly at a given stage
    return f"{WD}/{{iso}}/{stage}/{{iso}}.{stage}.fa"

# ---------------------------------------------------------------------------
# Phase target rules (run the build in pieces; break to edit config between)
# ---------------------------------------------------------------------------
rule phase0:   # contamination screen on the primaries -> review, then create {iso}.remove.txt
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
        expand(f"{WD}/{{iso}}/merqury/{{iso}}.qv",           iso=ISO),
        expand(f"{WD}/{{iso}}/merqury/{{iso}}.kmer_completeness.txt", iso=ISO),
        f"{WD}/versions.txt",

# ===========================================================================
# INPUT PREP
# ===========================================================================
# ONT: raw fastq_pass -> >=ont_min_len filtered set (in-workflow reproducibility)
rule filter_ont:
    input:  lst = lambda wc: config["ont_fastq_list"][wc.iso]
    output: filt = f"{WD}/{{iso}}/reads/{{iso}}.ont.filt.fastq.gz"
    params: tmp = f"{WD}/{{iso}}/reads/tmp", minlen = config["ont_min_len"]
    threads: 16
    resources: mem_mb=16000, runtime=360
    shell:
        r"""
        rm -rf {params.tmp}; mkdir -p {params.tmp}
        nl -ba {input.lst} | while read i f; do
            printf '%s %s/part_%05d.fq.gz %s\n' "$f" "{params.tmp}" "$i" "{params.minlen}"
        done > {params.tmp}/jobs.txt
        cat {params.tmp}/jobs.txt | xargs -P {threads} -L1 bash {SD}/filt_len.sh
        cat {params.tmp}/part_*.fq.gz > {output.filt}
        rm -rf {params.tmp}
        """

# 10x: move the 16 bp R1 barcode into a BX:Z tag, interleave ALL lanes for ARKS
rule prep_10x:
    input:
        r1 = lambda wc: config["tenx_R1"][wc.iso],
        r2 = lambda wc: config["tenx_R2"][wc.iso],
    output: bx = f"{WD}/{{iso}}/scaffold/reads.bx.fq.gz"
    params: bl = config["barcode_len"]
    threads: 4
    resources: mem_mb=8000, runtime=300
    shell:
        r"""
        mkdir -p $(dirname {output.bx})
        paste <(zcat {input.r1} | paste - - - -) <(zcat {input.r2} | paste - - - -) \
        | awk -F'\t' -v BL={params.bl} 'BEGIN{{OFS="\n"}} {{
              bx=substr($2,1,BL);
              r1s=substr($2,BL+1); r1q=substr($4,BL+1);
              h=$1; sub(/^@/,"",h); sub(/ .*/,"",h);
              print "@"h" BX:Z:"bx"-1", r1s, "+", r1q;
              h2=$5; sub(/^@/,"",h2); sub(/ .*/,"",h2);
              print "@"h2" BX:Z:"bx"-1", $6, "+", $8;
          }}' | gzip -c > {output.bx}
        """

# ===========================================================================
# PHASE 0: prep reference -- contamination screen on the PRIMARIES
# ===========================================================================
# coverage of the primary (the blobtools coverage axis)
rule prep_cov:
    input:
        pri = lambda wc: config["primary"][wc.iso],
        ont = f"{WD}/{{iso}}/reads/{{iso}}.ont.filt.fastq.gz",
    output: cov = f"{WD}/{{iso}}/prep/cov.txt"
    params: pre = config["ont_map_preset"], d = f"{WD}/{{iso}}/prep"
    threads: T
    resources: mem_mb=96000, runtime=600
    shell:
        r"""
        mkdir -p {params.d}
        {MM2} -ax {params.pre} -t {threads} {input.pri} {input.ont} 2>{params.d}/cov.mm.log \
          | {ST} sort -@8 -m3G -o {params.d}/cov.bam -
        {ST} index {params.d}/cov.bam
        {ST} coverage {params.d}/cov.bam > {output.cov}
        rm -f {params.d}/cov.bam {params.d}/cov.bam.bai
        """

# FRESH FCS-GX with a CURRENT db (gated; the original used a 2023-01-24 db).
# Needs the gx container + ~470 GB db provisioned (see config prep). If off,
# emits a placeholder so the report step is uniform.
rule prep_fcsgx:
    input:  pri = lambda wc: config["primary"][wc.iso]
    output: rpt = f"{WD}/{{iso}}/prep/fcsgx_report.txt"
    params:
        run = int(bool(config["prep"]["run_fcsgx"])), py = config["prep"]["fcsgx_py"],
        sif = config["prep"]["fcsgx_sif"], db = config["prep"]["fcsgx_db"],
        taxid = config["prep"]["fcsgx_taxid"], d = f"{WD}/{{iso}}/prep/fcsgx",
    threads: T
    resources: mem_mb=512000, runtime=720
    shell:
        r"""
        mkdir -p {params.d}
        if [ "{params.run}" = "1" ]; then
            FCS_DEFAULT_IMAGE={params.sif} python3 {params.py} screen genome \
              --fasta {input.pri} --out-dir {params.d} --gx-db {params.db} --tax-id {params.taxid}
            cp {params.d}/*.fcs_gx_report.txt {output.rpt}
        else
            echo "#FCS-GX not run (prep.run_fcsgx=false); provision gx db + container to enable." > {output.rpt}
        fi
        """

# combine GC + coverage + taxonomy (+ FCS-GX) -> candidate report + list + plots
rule prep_report:
    input:
        pri = lambda wc: config["primary"][wc.iso],
        cov = f"{WD}/{{iso}}/prep/cov.txt",
        fcs = f"{WD}/{{iso}}/prep/fcsgx_report.txt",
    output:
        report = f"{WD}/{{iso}}/prep/{{iso}}.contam_report.tsv",
        cand   = f"{WD}/{{iso}}/prep/{{iso}}.remove_candidates.txt",
        blob   = f"{WD}/{{iso}}/prep/{{iso}}.blob.png",
    params:
        d = f"{WD}/{{iso}}/prep", gcflag = config["prep"]["blob_gc_flag"],
        nwin = config["prep"]["blob_windows_per_contig"], wbp = config["prep"]["blob_window_bp"],
        nt = config["contam"]["nt_db"], nr = config["contam"]["nr_db"],
        do_blastx = int(bool(config["prep"]["blob_run_blastx"])),
    threads: T
    resources: mem_mb=64000, runtime=720
    shell:
        r"""
        mkdir -p {params.d}
        {ST} faidx {input.pri}
        {PY} {SD}/asm_stats.py {input.pri} --per-contig > {params.d}/gc.tsv
        awk -v G={params.gcflag} 'NR>1 && $3>G{{print $1}}' {params.d}/gc.tsv > {params.d}/blast_targets.txt
        : > {params.d}/regions.txt
        while read c; do
            L=$(awk -v c="$c" '$1==c{{print $2}}' {input.pri}.fai)
            {PY} {SD}/make_windows.py "$c" "$L" {params.nwin} {params.wbp} >> {params.d}/regions.txt
        done < {params.d}/blast_targets.txt
        FMT="6 qseqid pident length evalue bitscore staxids sscinames sblastnames sskingdoms stitle"
        rm -f {params.d}/blastn.tsv {params.d}/blastx.tsv
        if [ -s {params.d}/regions.txt ]; then
            {ST} faidx -r {params.d}/regions.txt {input.pri} > {params.d}/win.fa
            BLASTDB={params.nt} {BL}/blastn -task megablast -db {params.nt} -query {params.d}/win.fa \
                -num_threads {threads} -max_target_seqs 5 -evalue 1e-10 -outfmt "$FMT" > {params.d}/blastn.tsv || true
            if [ "{params.do_blastx}" = "1" ]; then
                BLASTDB={params.nr} {BL}/blastx -db {params.nr} -query {params.d}/win.fa \
                    -num_threads {threads} -max_target_seqs 5 -evalue 1e-5 -outfmt "$FMT" > {params.d}/blastx.tsv || true
            fi
        fi
        {PY} {SD}/blob_report.py --gc {params.d}/gc.tsv --cov {input.cov} --fcsgx {input.fcs} \
            --blastn {params.d}/blastn.tsv --gcflag {params.gcflag} \
            --report {output.report} --candidates {output.cand} --plotprefix {params.d}/{wildcards.iso}
        """

# start of Phase A: drop the VERIFIED removal list from the primary. The user
# creates {iso}.remove.txt from Phase-0 candidates (empty file = remove nothing).
rule prep_reference:
    input:
        pri    = lambda wc: config["primary"][wc.iso],
        remove = f"{WD}/{{iso}}/prep/{{iso}}.remove.txt",
    output: fa = f"{WD}/{{iso}}/prep/{{iso}}.prepped.fa"
    shell:  "{PY} {SD}/fasta_select.py {input.pri} --drop-file {input.remove} > {output.fa}"

# ===========================================================================
# PHASE A: purge_dups -> ARKS scaffold  (on the DECONTAMINATED primary)
# ===========================================================================
rule purge_dups:
    input:
        pri = f"{WD}/{{iso}}/prep/{{iso}}.prepped.fa",
        ont = f"{WD}/{{iso}}/reads/{{iso}}.ont.filt.fastq.gz",
    output: purged = f"{WD}/{{iso}}/purge/{{iso}}.purged.fa"
    params:
        d        = f"{WD}/{{iso}}/purge",
        map_pre  = config["ont_map_preset"],
        self_pre = config["purge"]["self_preset"],
        self_fl  = config["purge"]["self_flags"],
        calcuts  = lambda wc: config["purge"]["calcuts_opts"].get(wc.iso, ""),
        getseqs  = config["purge"]["get_seqs_opts"],
    threads: T
    resources: mem_mb=96000, runtime=480
    shell:
        r"""
        mkdir -p {params.d} && cd {params.d}
        {MM2} -x {params.map_pre} -t {threads} {input.pri} {input.ont} 2>mm.log | gzip -c > reads.paf.gz
        {PDB}/pbcstat reads.paf.gz
        {PDB}/calcuts {params.calcuts} PB.stat > cutoffs 2> calcuts.log
        {PDB}/split_fa {input.pri} > pri.split
        {MM2} -x {params.self_pre} {params.self_fl} -t {threads} pri.split pri.split 2>self.log | gzip -c > pri.split.self.paf.gz
        {PDB}/purge_dups -2 -T cutoffs -c PB.base.cov pri.split.self.paf.gz > dups.bed 2> pd.log
        {PDB}/get_seqs {params.getseqs} dups.bed {input.pri}    # writes purged.fa + hap.fa in cwd
        test -s purged.fa                                       # fail loudly if get_seqs produced nothing
        mv purged.fa {output.purged}                            # -> {{iso}}.purged.fa
        """

rule scaffold_arks:
    input:
        purged = f"{WD}/{{iso}}/purge/{{iso}}.purged.fa",
        bx     = f"{WD}/{{iso}}/scaffold/reads.bx.fq.gz",
    output: scaf = f"{WD}/{{iso}}/scaffold/{{iso}}.scaffold.fa"
    params: d = f"{WD}/{{iso}}/scaffold", ap = config["arks"]["params"], glob = config["arks"]["output_glob"]
    threads: T
    resources: mem_mb=120000, runtime=720
    shell:
        r"""
        cd {params.d}
        ln -sf {input.purged} draft.fa
        ln -sf {input.bx}     reads.fq.gz
        {ARCS} arks draft=draft reads=reads t={threads} {params.ap}
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
        import os
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
    shell:  "{PY} {SD}/apply_cuts.py {input.scaf} {input.cuts} {output.cut} {params.snap}"

# ===========================================================================
# QC (stage-parametrized: runs identically on 'scaffold' and 'cut')
# ===========================================================================
rule map_ont:
    input:
        fa  = lambda wc: asm(wc.stage).format(iso=wc.iso),
        ont = f"{WD}/{{iso}}/reads/{{iso}}.ont.filt.fastq.gz",
    output:
        bam = f"{WD}/{{iso}}/qc_{{stage}}/{{iso}}.{{stage}}.ont.bam",
        cov = f"{WD}/{{iso}}/qc_{{stage}}/cov.{{iso}}.txt",
    params: pre = config["ont_map_preset"]
    threads: T
    resources: mem_mb=96000, runtime=600
    shell:
        r"""
        mkdir -p $(dirname {output.bam})
        {MM2} -ax {params.pre} -t {threads} {input.fa} {input.ont} 2>{output.bam}.log \
          | {ST} sort -@ 8 -m 3G -o {output.bam} -
        {ST} index {output.bam}
        {ST} coverage {output.bam} > {output.cov}
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
    shell: "{PY} {SD}/join_qc.py {input.fa} {input.bam} {params.minscaf} {params.flank} {params.minspan} > {output.rpt}"

rule synteny:
    input:
        ma = f"{WD}/MA/{{stage}}/MA.{{stage}}.fa",
        mm = f"{WD}/MM/{{stage}}/MM.{{stage}}.fa",
    output:
        paf = f"{WD}/synteny_{{stage}}/MAxMM.paf",
        png = f"{WD}/synteny_{{stage}}/MAxMM_dotplot.png",
        tsv = f"{WD}/synteny_{{stage}}/correspondence.tsv",
    params:
        preset = config["synteny"]["preset"],
        minaln = config["synteny"]["min_aln_bp"],
        minlen = config["synteny"]["dotplot_min_len_bp"],
    threads: T
    resources: mem_mb=64000, runtime=300
    shell:
        r"""
        mkdir -p $(dirname {output.paf})
        {MM2} -x {params.preset} -t {threads} {input.ma} {input.mm} > {output.paf} 2>{output.paf}.log
        {PY} {SD}/paf_dotplot.py {output.paf} {output.png} {params.minlen}
        {PY} {SD}/synteny_classify.py {output.paf} {output.tsv} {params.minaln}
        """

# ===========================================================================
# PHASE C: consensus reference -> contamination screen -> final QC
# ===========================================================================
# emit the consensus decisions (config) + cut-fasta paths as a JSON spec, so
# build_consensus.py needs no YAML parser in the job env.
rule write_consensus_spec:
    output: f"{WD}/consensus/consensus_spec.json"
    run:
        import os, json
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
    shell:  "{PY} {SD}/build_consensus.py {input.spec} {output.fa}"

rule asm_stats:
    input:  f"{WD}/consensus/MAMM_final.fa"
    output: f"{WD}/consensus/asm_stats.txt"
    shell:  "{PY} {SD}/asm_stats.py {input} > {output}"

# contamination screen: per-contig GC (+ optional windowed BLAST on GC-flagged)
rule contam_screen:
    input:  fa = f"{WD}/consensus/MAMM_merged.fa"
    output: rpt = f"{WD}/consensus/contam_report.tsv"
    params:
        gc      = config["contam"]["gc_flag_high"],
        do_blast= int(bool(config["contam"]["run_blast"])),
        nt      = config["contam"]["nt_db"],
        nr      = config["contam"]["nr_db"],
        nwin    = config["contam"]["windows_per_contig"],
        wbp     = config["contam"]["window_bp"],
        d       = f"{WD}/consensus/contam",
    threads: T
    resources: mem_mb=64000, runtime=600
    shell:
        r"""
        mkdir -p {params.d}
        {ST} faidx {input.fa}
        # flag high-GC contigs
        {PY} {SD}/asm_stats.py {input.fa} --per-contig \
          | awk -v G={params.gc} 'NR>1 && $3>G {{print $1}}' > {params.d}/flagged.txt
        if [ "{params.do_blast}" = "1" ] && [ -s {params.d}/flagged.txt ]; then
            # windowed blastn->nt + blastx->nr on flagged contigs only
            : > {params.d}/regions.txt
            while read c; do
                L=$(awk -v c=$c '$1==c{{print $2}}' {input.fa}.fai)
                {PY} {SD}/make_windows.py $c $L {params.nwin} {params.wbp} >> {params.d}/regions.txt
            done < {params.d}/flagged.txt
            {ST} faidx -r {params.d}/regions.txt {input.fa} > {params.d}/win.fa
            FMT="6 qseqid pident length evalue bitscore staxids sscinames sblastnames sskingdoms stitle"
            BLASTDB={params.nt} {BL}/blastn -task megablast -db {params.nt} -query {params.d}/win.fa \
                -num_threads {threads} -max_target_seqs 5 -evalue 1e-10 -outfmt "$FMT" > {params.d}/blastn.tsv || true
            BLASTDB={params.nr} {BL}/blastx -db {params.nr} -query {params.d}/win.fa \
                -num_threads {threads} -max_target_seqs 5 -evalue 1e-5 -outfmt "$FMT" > {params.d}/blastx.tsv || true
        fi
        {PY} {SD}/contam_report.py {input.fa} {params.gc} {params.d} > {output.rpt}
        """

# final reference = consensus minus the reviewed drop_contigs list
rule finalize_ref:
    input:  fa = f"{WD}/consensus/MAMM_merged.fa", rpt = f"{WD}/consensus/contam_report.tsv"
    output: fa = f"{WD}/consensus/MAMM_final.fa"
    params: drop = " ".join(config["contam"]["drop_contigs"]) if config["contam"]["drop_contigs"] else ""
    shell:  "{PY} {SD}/fasta_select.py {input.fa} --drop {params.drop} > {output.fa}"

# coverage QC: BOTH isolates' ONT -> final reference (collapse / telomere)
rule qc_coverage:
    input:
        fa  = f"{WD}/consensus/MAMM_final.fa",
        ont = expand(f"{WD}/{{iso}}/reads/{{iso}}.ont.filt.fastq.gz", iso=ISO),
    output: cov = f"{WD}/consensus/cov_qc.txt"
    params:
        d = f"{WD}/consensus/covqc", isos = " ".join(ISO),
        collapse = config["qc"]["cov_collapse_ratio"], lowcov = config["qc"]["cov_lowcov_frac"],
        tk = config["qc"]["telo_topkmer_frac"], tc = config["qc"]["telo_canonical_frac"],
        pre = config["ont_map_preset"],
    threads: T
    resources: mem_mb=96000, runtime=600
    shell:
        r"""
        mkdir -p {params.d}
        for X in {params.isos}; do
            {MM2} -ax {params.pre} -t {threads} {input.fa} {WD}/$X/reads/$X.ont.filt.fastq.gz 2>{params.d}/$X.log \
              | {ST} sort -@8 -m3G -o {params.d}/$X.bam -
            {ST} index {params.d}/$X.bam
            {ST} coverage {params.d}/$X.bam > {params.d}/cov_$X.txt
        done
        {PY} {SD}/cov_qc.py {input.fa} {params.d}/cov_MA.txt {params.d}/cov_MM.txt \
             {params.collapse} {params.lowcov} {params.tk} {params.tc} > {output.cov}
        """

# merqury per isolate (QV + completeness + spectra-cn) on the per-isolate cut asm
rule qc_merqury:
    input:
        fa = f"{WD}/{{iso}}/cut/{{iso}}.cut.fa",
        r2 = lambda wc: config["tenx_R2"][wc.iso],
    output:
        qv  = f"{WD}/{{iso}}/merqury/{{iso}}.qv",
        cn  = f"{WD}/{{iso}}/merqury/{{iso}}.{{iso}}.cut.spectra-cn.hist",
    params: d = f"{WD}/{{iso}}/merqury", k = config["qc"]["meryl_k"], name = "{iso}"
    threads: T
    resources: mem_mb=96000, runtime=600
    shell:
        r"""
        mkdir -p {params.d}; cd {params.d}
        {MQ}/meryl k={params.k} threads={threads} count {input.r2} output {params.name}.meryl
        {MQ}/merqury.sh {params.name}.meryl {input.fa} {params.name} 2>&1 | tail -5
        """

rule qc_kmer_completeness:
    input:  cn = f"{WD}/{{iso}}/merqury/{{iso}}.{{iso}}.cut.spectra-cn.hist"
    output: f"{WD}/{{iso}}/merqury/{{iso}}.kmer_completeness.txt"
    params: filt = f"{WD}/{{iso}}/merqury/{{iso}}.filt"
    shell:  "{PY} {SD}/kmer_completeness.py {input.cn} $(cat {params.filt}) > {output}"

rule versions:
    output: f"{WD}/versions.txt"
    shell:
        r"""
        {{
          echo "minimap2: $({MM2} --version)"
          echo "samtools: $({ST} --version | head -1)"
          echo "purge_dups: $({PDB}/purge_dups 2>&1 | head -1 || true)"
          echo "arcs-make: $({ARCS} --version 2>&1 | head -1 || true)"
          echo "meryl: $({MQ}/meryl --version 2>&1 | head -1 || true)"
          echo "merqury.sh: $({MQ}/merqury.sh 2>&1 | head -1 || true)"
          echo "blastn: $({BL}/blastn -version | head -1)"
          echo "dustmasker: $({BL}/dustmasker -version 2>&1 | head -1 || true)"
        }} > {output}
        """
