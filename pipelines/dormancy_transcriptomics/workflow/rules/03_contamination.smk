rule kraken2_classify:
    input:
        r1="{results}/{species}/fastq/trimmed/{sample}_1.fq.gz",
        r2="{results}/{species}/fastq/trimmed/{sample}_2.fq.gz",
        db=config["kraken"]["db"],
    output:
        report="{results}/{species}/contamination/kraken2/{sample}.kraken.report.out",
        out=temp("{results}/{species}/contamination/kraken2/{sample}.kraken.out"),
    log:
        "{results}/{species}/logs/kraken2_classify/{sample}.log",
    threads: 16
    conda:
        "../envs/kraken2.yaml"
    shell:
        "kraken2 --paired --threads {threads} --db {input.db} --use-names "
        "--report {output.report} --output {output.out} {input.r1} {input.r2} "
        "> {log} 2>&1"
