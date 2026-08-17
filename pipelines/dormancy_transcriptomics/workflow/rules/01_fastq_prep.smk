rule merge_lanes:
    input:
        r1=lambda wc: raw_lane_files(wc, "_1"),
        r2=lambda wc: raw_lane_files(wc, "_2"),
    output:
        r1=temp("{results}/{species}/fastq/merged/{sample}_1.fq.gz"),
        r2=temp("{results}/{species}/fastq/merged/{sample}_2.fq.gz"),
    log:
        "{results}/{species}/logs/merge_lanes/{sample}.log",
    conda:
        "../envs/coreutils.yaml"
    shell:
        "cat {input.r1} > {output.r1} 2> {log} && cat {input.r2} > {output.r2} 2>> {log}"


rule fastqc_raw:
    input:
        r1="{results}/{species}/fastq/merged/{sample}_1.fq.gz",
        r2="{results}/{species}/fastq/merged/{sample}_2.fq.gz",
    output:
        r1_zip="{results}/{species}/qc/fastqc_raw/{sample}_1_fastqc.zip",
        r1_html="{results}/{species}/qc/fastqc_raw/{sample}_1_fastqc.html",
        r2_zip="{results}/{species}/qc/fastqc_raw/{sample}_2_fastqc.zip",
        r2_html="{results}/{species}/qc/fastqc_raw/{sample}_2_fastqc.html",
    params:
        outdir=lambda wc, output: os.path.dirname(output.r1_zip),
    log:
        "{results}/{species}/logs/fastqc_raw/{sample}.log",
    threads: 2
    conda:
        "../envs/fastqc.yaml"
    shell:
        "fastqc -t {threads} -o {params.outdir} {input.r1} {input.r2} > {log} 2>&1"


rule multiqc_raw:
    input:
        zips=lambda wc: [
            f"{wc.results}/{sp}/qc/fastqc_raw/{s}_{mate}_fastqc.zip"
            for sp in config["species"] for s in samples_of(sp) for mate in ("1", "2")
        ],
    output:
        html=report(
            "{results}/qc/multiqc_before_trimming/multiqc_report.html",
            caption="MultiQC summary of raw (pre-trim) FastQC reports, both species.",
            category="Raw-read QC",
        ),
    params:
        fastqc_dirs=lambda wc: sorted({
            f"{wc.results}/{sp}/qc/fastqc_raw" for sp in config["species"]
        }),
        out_dir=lambda wc, output: os.path.dirname(output.html),
    log:
        "{results}/logs/multiqc_raw.log",
    conda:
        "../envs/multiqc.yaml"
    shell:
        "multiqc --force -o {params.out_dir} {params.fastqc_dirs} > {log} 2>&1"


rule trim_fastp:
    input:
        r1="{results}/{species}/fastq/merged/{sample}_1.fq.gz",
        r2="{results}/{species}/fastq/merged/{sample}_2.fq.gz",
    output:
        r1="{results}/{species}/fastq/trimmed/{sample}_1.fq.gz",
        r2="{results}/{species}/fastq/trimmed/{sample}_2.fq.gz",
        json="{results}/{species}/qc/fastp/{sample}_fastp.json",
        html="{results}/{species}/qc/fastp/{sample}_fastp.html",
    log:
        "{results}/{species}/logs/trim_fastp/{sample}.log",
    threads: 2
    conda:
        "../envs/fastp.yaml"
    shell:
        "fastp -w {threads} "
        "--in1 {input.r1} --in2 {input.r2} "
        "--out1 {output.r1} --out2 {output.r2} "
        "--json {output.json} --html {output.html} "
        "> {log} 2>&1"
