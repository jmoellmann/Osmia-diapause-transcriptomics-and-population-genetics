rule feature_bed_cds:
    input:
        gff=lambda wc: config["reference"][wc.species]["annotation_gff"],
    output:
        bed="{results}/{species}/beds/cds.bed",
    log:
        "{results}/{species}/logs/feature_bed_cds.log",
    conda:
        "../envs/coreutils.yaml"
    shell:
        r"""
        (awk -F'\t' -v OFS='\t' '$3 == "CDS" && !/^#/ {{
            print $1, $4-1, $5, $9, "NA", $7;
        }}' {input.gff} | sort -k1,1 -k2,2n > {output.bed}) > {log} 2>&1
        """


rule pixy_10kb:
    input:
        vcf="{results}/{species}/final_calling/allsites.filtered.vcf.gz",
        csi="{results}/{species}/final_calling/allsites.filtered.vcf.gz.csi",
        populations=lambda wc: config["reference"][wc.species]["sample2ploidy"],
    output:
        pi="{results}/{species}/pixy_10kb/10kb_pi.txt",
        fst="{results}/{species}/pixy_10kb/10kb_fst.txt",
        dxy="{results}/{species}/pixy_10kb/10kb_dxy.txt",
        theta="{results}/{species}/pixy_10kb/10kb_watterson_theta.txt",
        tajima_d="{results}/{species}/pixy_10kb/10kb_tajima_d.txt",
    params:
        outdir=lambda wc, output: os.path.dirname(output.pi),
        window=10000,
    log:
        "{results}/{species}/logs/pixy_10kb.log",
    threads: 30
    conda:
        "../envs/pixy.yaml"
    shell:
        "pixy --stats pi fst dxy watterson_theta tajima_d "
        "--window_size {params.window} --n_cores {threads} "
        "--populations {input.populations} --vcf {input.vcf} "
        "--output_folder {params.outdir} --output_prefix 10kb > {log} 2>&1"


rule pixy_cds:
    input:
        vcf="{results}/{species}/final_calling/allsites.filtered.vcf.gz",
        csi="{results}/{species}/final_calling/allsites.filtered.vcf.gz.csi",
        populations=lambda wc: config["reference"][wc.species]["sample2ploidy"],
        bed="{results}/{species}/beds/cds.bed",
    output:
        pi="{results}/{species}/pixy_cds/cds_pi.txt",
        fst="{results}/{species}/pixy_cds/cds_fst.txt",
        dxy="{results}/{species}/pixy_cds/cds_dxy.txt",
        theta="{results}/{species}/pixy_cds/cds_watterson_theta.txt",
        tajima_d="{results}/{species}/pixy_cds/cds_tajima_d.txt",
    params:
        outdir=lambda wc, output: os.path.dirname(output.pi),
    log:
        "{results}/{species}/logs/pixy_cds.log",
    threads: 30
    conda:
        "../envs/pixy.yaml"
    shell:
        "pixy --stats pi fst dxy watterson_theta tajima_d "
        "--n_cores {threads} --populations {input.populations} --vcf {input.vcf} "
        "--bed_file {input.bed} --output_folder {params.outdir} --output_prefix cds "
        "> {log} 2>&1"


rule plink_prune:
    input:
        vcf="{results}/{species}/final_calling/variants.filtered.snps-only.haploids.vcf.gz",
    output:
        multiext("{results}/{species}/admixture/pruned_dataset", ".bed", ".bim", ".fam"),
        prune_in="{results}/{species}/admixture/pruned_dataset.prune.in",
    params:
        out=lambda wc, output: output[0][: -len(".bed")],
    log:
        "{results}/{species}/logs/plink_prune.log",
    threads: 32
    conda:
        "../envs/plink.yaml"
    shell:
        "plink --vcf {input.vcf} --allow-extra-chr --double-id "
        "--indep-pairwise 50 10 0.2 --make-bed --threads {threads} "
        "--out {params.out} > {log} 2>&1"


rule plink_extract_pruned:
    input:
        bed="{results}/{species}/admixture/pruned_dataset.bed",
        prune_in="{results}/{species}/admixture/pruned_dataset.prune.in",
    output:
        multiext("{results}/{species}/admixture/final_pruned", ".bed", ".bim", ".fam"),
    params:
        bfile=lambda wc, input: input.bed[: -len(".bed")],
        out=lambda wc, output: output[0][: -len(".bed")],
    log:
        "{results}/{species}/logs/plink_extract_pruned.log",
    threads: 32
    conda:
        "../envs/plink.yaml"
    shell:
        "plink --bfile {params.bfile} --allow-extra-chr --double-id "
        "--extract {input.prune_in} --threads {threads} --make-bed "
        "--out {params.out} > {log} 2>&1"


rule admixture_k:
    input:
        bed="{results}/{species}/admixture/final_pruned.bed",
    output:
        log_out="{results}/{species}/admixture/log_{k}.out",
    params:
        outdir=lambda wc, output: os.path.dirname(output.log_out),
        bed_name=lambda wc, input: os.path.basename(input.bed),
    log:
        "{results}/{species}/logs/admixture_k/{k}.log",
    threads: 4
    conda:
        "../envs/admixture.yaml"
    shell:
        "(cd {params.outdir} && admixture --cv -j{threads} {params.bed_name} {wildcards.k} "
        "| tee log_{wildcards.k}.out) > {log} 2>&1"


rule admixture_summary:
    input:
        logs=expand(
            "{{results}}/{{species}}/admixture/log_{k}.out", k=range(1, 9),
        ),
    output:
        summary="{results}/{species}/admixture/summary.txt",
    log:
        "{results}/{species}/logs/admixture_summary.log",
    conda:
        "../envs/coreutils.yaml"
    shell:
        "grep -h CV {input.logs} | sort -k3 -n > {output.summary} 2> {log}"


rule ld_decay:
    input:
        vcf="{results}/{species}/final_calling/variants.filtered.vcf.gz",
    output:
        stat="{results}/{species}/lddecay/out.stat.gz",
    params:
        outprefix=lambda wc, output: output.stat[: -len(".stat.gz")],
    log:
        "{results}/{species}/logs/ld_decay.log",
    conda:
        "../envs/popldecay.yaml"
    shell:
        "PopLDdecay -InVCF {input.vcf} -OutStat {params.outprefix} -OutType 2 > {log} 2>&1"


rule general_vcf_stats:
    input:
        vcf="{results}/{species}/final_calling/allsites.filtered.vcf.gz",
    output:
        frq="{results}/{species}/general_stats/stats.frq",
        idepth="{results}/{species}/general_stats/stats.idepth",
        ldepth_mean="{results}/{species}/general_stats/stats.ldepth.mean",
        lqual="{results}/{species}/general_stats/stats.lqual",
        imiss="{results}/{species}/general_stats/stats.imiss",
        lmiss="{results}/{species}/general_stats/stats.lmiss",
        het="{results}/{species}/general_stats/stats.het",
    params:
        out=lambda wc, output: output.frq[: -len(".frq")],
    log:
        "{results}/{species}/logs/general_vcf_stats.log",
    conda:
        "../envs/vcftools.yaml"
    shell:
        """
        (vcftools --gzvcf {input.vcf} --freq2 --out {params.out} --max-alleles 2
        vcftools --gzvcf {input.vcf} --depth --out {params.out}
        vcftools --gzvcf {input.vcf} --site-mean-depth --out {params.out}
        vcftools --gzvcf {input.vcf} --site-quality --out {params.out}
        vcftools --gzvcf {input.vcf} --missing-indv --out {params.out}
        vcftools --gzvcf {input.vcf} --missing-site --out {params.out}
        vcftools --gzvcf {input.vcf} --het --out {params.out}) > {log} 2>&1
        """
