rule star_index:
    input:
        genome_fasta=lambda wc: config["reference"][wc.species]["genome_fasta"],
        gtf=lambda wc: config["reference"][wc.species]["gtf"],
    output:
        index_dir=directory("{results}/{species}/star_index"),
    log:
        "{results}/{species}/logs/star_index.log",
    threads: 8
    resources:
        mem_mb=32000,
    conda:
        "../envs/star.yaml"
    shell:
        "mkdir -p {output.index_dir} && "
        "STAR --runThreadN {threads} --runMode genomeGenerate --genomeSAindexNbases 12 "
        "--genomeDir {output.index_dir} --genomeFastaFiles {input.genome_fasta} "
        "--sjdbGTFfile {input.gtf} > {log} 2>&1"


rule star_align:
    input:
        r1="{results}/{species}/fastq/trimmed/{sample}_1.fq.gz",
        r2="{results}/{species}/fastq/trimmed/{sample}_2.fq.gz",
        index_dir="{results}/{species}/star_index",
    output:
        bam="{results}/{species}/star/{sample}/Aligned.sortedByCoord.out.bam",
        transcriptome_bam="{results}/{species}/star/{sample}/Aligned.toTranscriptome.out.bam",
        gene_counts="{results}/{species}/star/{sample}/ReadsPerGene.out.tab",
        log_final="{results}/{species}/star/{sample}/Log.final.out",
        unmapped1="{results}/{species}/star/{sample}/Unmapped.out.mate1",
        unmapped2="{results}/{species}/star/{sample}/Unmapped.out.mate2",
    params:
        out_prefix=lambda wc, output: os.path.dirname(output.bam) + "/",
    log:
        "{results}/{species}/logs/star_align/{sample}.log",
    threads: 5
    resources:
        mem_mb=16000,
    conda:
        "../envs/star.yaml"
    shell:
        "STAR --runThreadN {threads} --readFilesCommand gunzip -c "
        "--outReadsUnmapped Fastx --outSAMtype BAM SortedByCoordinate "
        "--limitBAMsortRAM 12000000000 --twopassMode Basic "
        "--quantMode TranscriptomeSAM GeneCounts "
        "--genomeDir {input.index_dir} --readFilesIn {input.r1} {input.r2} "
        "--outFileNamePrefix {params.out_prefix} > {log} 2>&1"


rule kallisto_index:
    input:
        cds_fasta=lambda wc: config["reference"][wc.species]["cds_fasta"],
    output:
        index="{results}/{species}/kallisto_index/{species}.index",
    log:
        "{results}/{species}/logs/kallisto_index.log",
    conda:
        "../envs/kallisto.yaml"
    shell:
        "kallisto index -i {output.index} {input.cds_fasta} > {log} 2>&1"


rule kallisto_quant:
    input:
        r1="{results}/{species}/fastq/trimmed/{sample}_1.fq.gz",
        r2="{results}/{species}/fastq/trimmed/{sample}_2.fq.gz",
        index="{results}/{species}/kallisto_index/{species}.index",
    output:
        h5="{results}/{species}/kallisto/{sample}/abundance.h5",
        tsv="{results}/{species}/kallisto/{sample}/abundance.tsv",
        run_info="{results}/{species}/kallisto/{sample}/run_info.json",
    params:
        out_dir=lambda wc, output: os.path.dirname(output.h5),
    log:
        "{results}/{species}/logs/kallisto_quant/{sample}.log",
    threads: 16
    conda:
        "../envs/kallisto.yaml"
    shell:
        "kallisto quant --threads {threads} -i {input.index} -o {params.out_dir} "
        "{input.r1} {input.r2} > {log} 2>&1"


rule multiqc_kallisto:
    input:
        run_info=lambda wc: expand(
            "{results}/{species}/kallisto/{sample}/run_info.json",
            results=wc.results, species=wc.species, sample=samples_of(wc.species),
        ),
    output:
        html=report(
            "{results}/{species}/qc/multiqc_report.html",
            caption="MultiQC summary of kallisto pseudoalignment rates, per sample.",
            category="Kallisto QC",
        ),
        txt="{results}/{species}/qc/multiqc_data/multiqc_kallisto.txt",
    params:
        kallisto_dir=lambda wc, input: os.path.dirname(os.path.dirname(input.run_info[0])),
        out_dir=lambda wc, output: os.path.dirname(output.html),
    log:
        "{results}/{species}/logs/multiqc_kallisto.log",
    conda:
        "../envs/multiqc.yaml"
    shell:
        "multiqc --force -o {params.out_dir} {params.kallisto_dir} > {log} 2>&1"
