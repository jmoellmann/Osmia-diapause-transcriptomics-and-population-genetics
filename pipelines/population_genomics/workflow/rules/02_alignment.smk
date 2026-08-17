rule bowtie2_index:
    input:
        genome="{results}/resources/genomes/{species}/genome.fa",
    output:
        multiext(
            "{results}/resources/genomes/{species}/genome.fa.bt2idx",
            ".1.bt2", ".2.bt2", ".3.bt2", ".4.bt2", ".rev.1.bt2", ".rev.2.bt2",
        ),
    params:
        prefix=lambda wc, output: output[0][: -len(".1.bt2")],
    log:
        "{results}/{species}/logs/bowtie2_index.log",
    threads: 20
    conda:
        "../envs/bowtie2.yaml"
    shell:
        "bowtie2-build --threads {threads} {input.genome} {params.prefix} > {log} 2>&1"


rule align_bowtie2:
    input:
        r1="{results}/{species}/fastq/trimmed/{sample}_1.fq.gz",
        r2="{results}/{species}/fastq/trimmed/{sample}_2.fq.gz",
        index=multiext(
            "{results}/resources/genomes/{species}/genome.fa.bt2idx",
            ".1.bt2", ".2.bt2", ".3.bt2", ".4.bt2", ".rev.1.bt2", ".rev.2.bt2",
        ),
    output:
        bam=temp("{results}/{species}/bam/raw/{sample}.bam"),
    params:
        prefix=lambda wc, input: input.index[0][: -len(".1.bt2")],
    log:
        "{results}/{species}/logs/align/{sample}.log",
    threads: 20
    conda:
        "../envs/bowtie2.yaml"
    shell:
        "(bowtie2 -p {threads} -x {params.prefix} -1 {input.r1} -2 {input.r2} "
        "| samtools view -@ {threads} -o {output.bam}) > {log} 2>&1"


rule multiqc_bowtie2:
    input:
        logs=lambda wc: [
            f"{wc.results}/{sp}/logs/align/{s}.log"
            for sp in config["species"] for s in samples_of(sp)
        ],
    output:
        html=report(
            "{results}/qc/multiqc_bowtie2_report.html",
            caption="MultiQC summary of bowtie2 alignment rates, both species.",
            category="Alignment QC",
        ),
    params:
        log_dirs=lambda wc: sorted({
            f"{wc.results}/{sp}/logs/align" for sp in config["species"]
        }),
        out_dir=lambda wc, output: os.path.dirname(output.html),
    log:
        "{results}/logs/multiqc_bowtie2.log",
    conda:
        "../envs/multiqc.yaml"
    shell:
        "multiqc --force -o {params.out_dir} {params.log_dirs} > {log} 2>&1"


rule reheader_bam:
    input:
        bam="{results}/{species}/bam/raw/{sample}.bam",
    output:
        bam=temp("{results}/{species}/bam/reheadered/{sample}.bam"),
    log:
        "{results}/{species}/logs/reheader_bam/{sample}.log",
    conda:
        "../envs/samtools.yaml"
    shell:
        "samtools addreplacerg -r ID:{wildcards.sample} -r SM:{wildcards.sample} "
        "-o {output.bam} {input.bam} > {log} 2>&1"


rule sort_bam:
    input:
        bam="{results}/{species}/bam/reheadered/{sample}.bam",
    output:
        bam="{results}/{species}/bam/sorted/{sample}.bam",
    log:
        "{results}/{species}/logs/sort_bam/{sample}.log",
    threads: 3
    resources:
        mem_mb=6000,
    conda:
        "../envs/samtools.yaml"
    shell:
        "samtools sort --threads {threads} -m 2000M -o {output.bam} {input.bam} > {log} 2>&1"


rule index_bam:
    input:
        bam="{results}/{species}/bam/sorted/{sample}.bam",
    output:
        bai="{results}/{species}/bam/sorted/{sample}.bam.bai",
    log:
        "{results}/{species}/logs/index_bam/{sample}.log",
    conda:
        "../envs/samtools.yaml"
    shell:
        "samtools index {input.bam} > {log} 2>&1"
