rule call_variants_chrom:
    input:
        unpack(bams_of_species),
        genome=genome_fasta,
        sample2ploidy=lambda wc: config["reference"][wc.species]["sample2ploidy"],
    output:
        vcf=temp("{results}/{species}/final_calling/tmp/raw_{chromosome}.vcf.gz"),
    params:
        q=config["mpileup"]["min_mq"],
        Q=config["mpileup"]["min_bq"],
        d=config["mpileup"]["max_depth"],
        L=config["mpileup"]["max_depth_per_file"],
    log:
        "{results}/{species}/logs/call_variants/{chromosome}.log",
    conda:
        "../envs/bcftools.yaml"
    shell:
        "(bcftools mpileup -q {params.q} -Q {params.Q} -d {params.d} -L {params.L} -Ou "
        "-a AD,ADF,ADR,DP,SP,SCR,INFO/AD -f {input.genome} -r {wildcards.chromosome} {input.bam} "
        "| bcftools call -m -Oz -a GQ --samples-file {input.sample2ploidy} "
        "-o {output.vcf}) > {log} 2>&1"


rule concat_variants:
    input:
        vcf=lambda wc: expand(
            "{{results}}/{{species}}/final_calling/tmp/raw_{chromosome}.vcf.gz",
            chromosome=chromosome_list(wc.species),
        ),
    output:
        vcf="{results}/{species}/final_calling/reheadered.vcf.gz",
    log:
        "{results}/{species}/logs/concat_variants.log",
    conda:
        "../envs/bcftools.yaml"
    shell:
        "bcftools concat {input.vcf} -O z -o {output.vcf} > {log} 2>&1"


rule index_reheadered:
    input:
        vcf="{results}/{species}/final_calling/reheadered.vcf.gz",
    output:
        csi="{results}/{species}/final_calling/reheadered.vcf.gz.csi",
    log:
        "{results}/{species}/logs/index_reheadered.log",
    conda:
        "../envs/bcftools.yaml"
    shell:
        "bcftools index --threads 1 {input.vcf} > {log} 2>&1"


rule filter_variants_chrom:
    input:
        vcf="{results}/{species}/final_calling/reheadered.vcf.gz",
        csi="{results}/{species}/final_calling/reheadered.vcf.gz.csi",
    output:
        vcf=temp("{results}/{species}/final_calling/tmp/chr_{chromosome}.filtered.vcf.gz"),
    params:
        expr=(
            f"QUAL>={config['variant_filter']['min_qual']} && "
            f"MEAN(FORMAT/DP)>={config['variant_filter']['min_mean_dp']} && "
            f"MEAN(FORMAT/DP)<={config['variant_filter']['max_mean_dp']} && "
            f"F_MISSING<={config['variant_filter']['max_f_missing']}"
        ),
    log:
        "{results}/{species}/logs/filter_variants/{chromosome}.log",
    conda:
        "../envs/bcftools.yaml"
    shell:
        "bcftools view -r {wildcards.chromosome} -i '{params.expr}' {input.vcf} "
        "-O z -o {output.vcf} > {log} 2>&1"


rule index_filtered_chrom:
    input:
        vcf="{results}/{species}/final_calling/tmp/chr_{chromosome}.filtered.vcf.gz",
    output:
        csi=temp("{results}/{species}/final_calling/tmp/chr_{chromosome}.filtered.vcf.gz.csi"),
    log:
        "{results}/{species}/logs/index_filtered_chrom/{chromosome}.log",
    conda:
        "../envs/bcftools.yaml"
    shell:
        "bcftools index --threads 1 {input.vcf} > {log} 2>&1"


rule concat_filtered:
    input:
        vcf=lambda wc: expand(
            "{{results}}/{{species}}/final_calling/tmp/chr_{chromosome}.filtered.vcf.gz",
            chromosome=chromosome_list(wc.species),
        ),
        csi=lambda wc: expand(
            "{{results}}/{{species}}/final_calling/tmp/chr_{chromosome}.filtered.vcf.gz.csi",
            chromosome=chromosome_list(wc.species),
        ),
    output:
        vcf="{results}/{species}/final_calling/allsites.filtered.vcf.gz",
    log:
        "{results}/{species}/logs/concat_filtered.log",
    conda:
        "../envs/bcftools.yaml"
    shell:
        "bcftools concat {input.vcf} -O z -o {output.vcf} > {log} 2>&1"


rule index_allsites_filtered:
    input:
        vcf="{results}/{species}/final_calling/allsites.filtered.vcf.gz",
    output:
        csi="{results}/{species}/final_calling/allsites.filtered.vcf.gz.csi",
    log:
        "{results}/{species}/logs/index_allsites_filtered.log",
    conda:
        "../envs/bcftools.yaml"
    shell:
        "bcftools index --threads 1 {input.vcf} > {log} 2>&1"


rule subset_variant_type:
    input:
        vcf="{results}/{species}/final_calling/allsites.filtered.vcf.gz",
        csi="{results}/{species}/final_calling/allsites.filtered.vcf.gz.csi",
    output:
        vcf="{results}/{species}/final_calling/variants.filtered{suffix}.vcf.gz",
    wildcard_constraints:
        suffix="|.snps-only|.indels-only",
    params:
        v=lambda wc: {
            "": "snps,indels", ".snps-only": "snps", ".indels-only": "indels",
        }[wc.suffix],
    log:
        "{results}/{species}/logs/subset_variant_type/variants.filtered{suffix}.log",
    conda:
        "../envs/bcftools.yaml"
    shell:
        "bcftools view --threads 1 -v {params.v} {input.vcf} -O z -o {output.vcf} > {log} 2>&1"


rule subset_haploid_snps:
    input:
        vcf="{results}/{species}/final_calling/variants.filtered.snps-only.vcf.gz",
        haploids=lambda wc: config["reference"][wc.species]["haploids"],
    output:
        vcf="{results}/{species}/final_calling/variants.filtered.snps-only.haploids.vcf.gz",
    log:
        "{results}/{species}/logs/subset_haploid_snps.log",
    conda:
        "../envs/bcftools.yaml"
    shell:
        "bcftools view -S {input.haploids} {input.vcf} -O z -o {output.vcf} > {log} 2>&1"
