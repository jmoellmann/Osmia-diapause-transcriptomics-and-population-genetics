rule nsl_chrom_rename_map:
    input:
        chromosomes=chromosomes_file,
    output:
        map="{results}/{species}/nsl/chrom_rename_map.tsv",
    log:
        "{results}/{species}/logs/nsl_chrom_rename_map.log",
    conda:
        "../envs/python.yaml"
    script:
        "../scripts/nsl_chrom_rename_map.py"


rule nsl_filter:
    input:
        vcf="{results}/{species}/final_calling/variants.filtered.snps-only.haploids.vcf.gz",
    output:
        vcf=temp("{results}/{species}/nsl/selscan-ready.vcf.gz"),
    log:
        "{results}/{species}/logs/nsl_filter/{species}.log",
    conda:
        "../envs/bcftools.yaml"
    shell:
        r"""
        (bcftools view {input.vcf} -m2 -M2 -v snps -i 'F_MISSING=0' -Ou \
            | bcftools annotate -x INFO,^FORMAT/GT -Oz -o {output.vcf}) > {log} 2>&1
        """


rule nsl_pseudophase:
    input:
        vcf="{results}/{species}/nsl/selscan-ready.vcf.gz",
    output:
        vcf=temp("{results}/{species}/nsl/selscan-ready2.vcf.gz"),
        csi=temp("{results}/{species}/nsl/selscan-ready2.vcf.gz.csi"),
    log:
        "{results}/{species}/logs/nsl_pseudophase/{species}.log",
    conda:
        "../envs/bcftools.yaml"
    shell:
        r"""
        (zcat {input.vcf} | awk 'BEGIN{{OFS="\t"}}
            /^#/ {{print; next}}
            {{
                for(i=10; i<=NF; i++) {{
                    if($i == "0") $i = "0|0"
                    else if($i == "1") $i = "1|1"
                }}
                print
            }}' | bcftools view -Oz -o {output.vcf}
        bcftools index {output.vcf}) > {log} 2>&1
        """


rule nsl_rename_chroms:
    input:
        vcf="{results}/{species}/nsl/selscan-ready2.vcf.gz",
        csi="{results}/{species}/nsl/selscan-ready2.vcf.gz.csi",
        map="{results}/{species}/nsl/chrom_rename_map.tsv",
    output:
        vcf=temp("{results}/{species}/nsl/selscan-ready.renamed.vcf.gz"),
        csi=temp("{results}/{species}/nsl/selscan-ready.renamed.vcf.gz.csi"),
    log:
        "{results}/{species}/logs/nsl_rename_chroms/{species}.log",
    conda:
        "../envs/bcftools.yaml"
    shell:
        r"""
        (bcftools annotate --rename-chrs {input.map} {input.vcf} -Oz -o {output.vcf}
        bcftools index {output.vcf}) > {log} 2>&1
        """


rule nsl_subset_chrom:
    input:
        vcf="{results}/{species}/nsl/selscan-ready.renamed.vcf.gz",
        csi="{results}/{species}/nsl/selscan-ready.renamed.vcf.gz.csi",
    output:
        vcf=temp("{results}/{species}/nsl/chr{chrom_idx}.vcf.gz"),
    wildcard_constraints:
        chrom_idx=r"\d+",
    log:
        "{results}/{species}/logs/nsl_subset_chrom/{species}.chr{chrom_idx}.log",
    conda:
        "../envs/bcftools.yaml"
    shell:
        "bcftools view {input.vcf} -r {wildcards.chrom_idx} -Oz -o {output.vcf} > {log} 2>&1"


rule nsl_selscan:
    input:
        vcf="{results}/{species}/nsl/chr{chrom_idx}.vcf.gz",
    output:
        nsl_out=temp("{results}/{species}/nsl/chr{chrom_idx}.selscan.nsl.nsl.out"),
    params:
        out_prefix="{results}/{species}/nsl/chr{chrom_idx}.selscan.nsl",
    wildcard_constraints:
        chrom_idx=r"\d+",
    log:
        "{results}/{species}/logs/nsl_selscan/{species}.chr{chrom_idx}.log",
    conda:
        "../envs/selscan.yaml"
    shell:
        "selscan --nsl --vcf {input.vcf} --unphased --out {params.out_prefix} > {log} 2>&1"


rule nsl_norm:
    input:
        nsl_out=lambda wc: expand(
            "{results}/{species}/nsl/chr{chrom_idx}.selscan.nsl.nsl.out",
            results=wc.results, species=wc.species,
            chrom_idx=range(1, len(chromosome_list(wc.species)) + 1),
        ),
    output:
        norm=temp(expand(
            "{{results}}/{{species}}/nsl/chr{chrom_idx}.selscan.nsl.nsl.out.100bins.norm",
            chrom_idx=range(1, 17),
        )),
    log:
        "{results}/{species}/logs/nsl_norm.log",
    conda:
        "../envs/selscan.yaml"
    shell:
        "norm --nsl --files {input.nsl_out} > {log} 2>&1"


rule nsl_combine:
    input:
        norm=lambda wc: expand(
            "{results}/{species}/nsl/chr{chrom_idx}.selscan.nsl.nsl.out.100bins.norm",
            results=wc.results, species=wc.species,
            chrom_idx=range(1, len(chromosome_list(wc.species)) + 1),
        ),
    output:
        tsv="{results}/{species}/nsl/combined.nsl.norm.tsv",
    log:
        "{results}/{species}/logs/nsl_combine.log",
    conda:
        "../envs/coreutils.yaml"
    shell:
        r"""
        (rm -f {output.tsv}
        chr=1
        for f in {input.norm}; do
            awk -v chr="$chr" 'BEGIN{{OFS="\t"}} {{$1=chr; print}}' "$f" >> {output.tsv}
            chr=$((chr + 1))
        done) > {log} 2>&1
        """
