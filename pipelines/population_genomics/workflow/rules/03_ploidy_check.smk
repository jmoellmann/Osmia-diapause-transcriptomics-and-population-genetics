rule ploidy_check_call_chrom:
    input:
        unpack(bams_of_species),
        genome=genome_fasta,
    output:
        vcf=temp("{results}/{species}/ploidy_check/tmp/{chromosome}.vcf.gz"),
    log:
        "{results}/{species}/logs/ploidy_check_call/{chromosome}.log",
    conda:
        "../envs/bcftools.yaml"
    shell:
        "(bcftools mpileup -q 20 -Q 20 -d 250 -L 250 -Ou "
        "-a AD,ADF,ADR,DP,SP,SCR,INFO/AD -f {input.genome} -r {wildcards.chromosome} {input.bam} "
        "| bcftools call -mv -Oz -a GQ --ploidy 2 -o {output.vcf}) > {log} 2>&1"


rule ploidy_check_concat:
    input:
        vcf=lambda wc: expand(
            "{{results}}/{{species}}/ploidy_check/tmp/{chromosome}.vcf.gz",
            chromosome=chromosome_list(wc.species),
        ),
    output:
        vcf="{results}/{species}/ploidy_check/reheadered.vcf.gz",
    log:
        "{results}/{species}/logs/ploidy_check_concat.log",
    conda:
        "../envs/bcftools.yaml"
    shell:
        "bcftools concat {input.vcf} -O z -o {output.vcf} > {log} 2>&1"


rule ploidy_check_filter:
    input:
        vcf="{results}/{species}/ploidy_check/reheadered.vcf.gz",
    output:
        vcf="{results}/{species}/ploidy_check/reheadered.filtered.vcf.gz",
    params:
        expr="QUAL>20 && INFO/DP>5",
    log:
        "{results}/{species}/logs/ploidy_check_filter.log",
    conda:
        "../envs/bcftools.yaml"
    shell:
        "bcftools view -i '{params.expr}' {input.vcf} -O z -o {output.vcf} > {log} 2>&1"


rule ploidy_check_het:
    input:
        vcf="{results}/{species}/ploidy_check/reheadered.filtered.vcf.gz",
    output:
        het="{results}/{species}/ploidy_check/het_check.het",
    params:
        out=lambda wc, output: output.het[: -len(".het")],
    log:
        "{results}/{species}/logs/ploidy_check_het.log",
    conda:
        "../envs/vcftools.yaml"
    shell:
        "vcftools --gzvcf {input.vcf} --het --out {params.out} > {log} 2>&1"


rule ploidy_check_depth:
    input:
        vcf="{results}/{species}/ploidy_check/reheadered.filtered.vcf.gz",
    output:
        idepth="{results}/{species}/ploidy_check/depth_check.idepth",
    params:
        out=lambda wc, output: output.idepth[: -len(".idepth")],
    log:
        "{results}/{species}/logs/ploidy_check_depth.log",
    conda:
        "../envs/vcftools.yaml"
    shell:
        "vcftools --gzvcf {input.vcf} --depth --out {params.out} > {log} 2>&1"
