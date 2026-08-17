rule download_genome:
    output:
        fa="{results}/resources/genomes/{species}/genome.fa",
    params:
        url=lambda wc: config["reference"][wc.species]["genome_fasta_url"],
    log:
        "{results}/logs/download_genome/{species}.log",
    conda:
        "../envs/download.yaml"
    shell:
        "(wget -qO- {params.url} | gunzip -c > {output.fa}) > {log} 2>&1"


rule download_gtf:
    output:
        gtf="{results}/resources/genomes/{species}/annotation.gtf",
    params:
        url=lambda wc: config["reference"][wc.species]["gtf_url"],
    log:
        "{results}/logs/download_gtf/{species}.log",
    conda:
        "../envs/download.yaml"
    shell:
        "(wget -qO- {params.url} | gunzip -c > {output.gtf}) > {log} 2>&1"


rule download_cds_fasta:
    output:
        fa="{results}/resources/genomes/{species}/cds.fa",
    params:
        url=lambda wc: config["reference"][wc.species]["cds_fasta_url"],
    log:
        "{results}/logs/download_cds_fasta/{species}.log",
    conda:
        "../envs/download.yaml"
    shell:
        "(wget -qO- {params.url} | gunzip -c > {output.fa}) > {log} 2>&1"


rule download_annotation_gff:
    output:
        gff="{results}/resources/genomes/{species}/annotation.gff",
    params:
        url=lambda wc: config["reference"][wc.species]["annotation_gff_url"],
    log:
        "{results}/logs/download_annotation_gff/{species}.log",
    conda:
        "../envs/download.yaml"
    shell:
        "(wget -qO- {params.url} | gunzip -c > {output.gff}) > {log} 2>&1"


rule download_kraken_db:
    output:
        directory("{results}/resources/kraken_db"),
    params:
        url=config["kraken"]["db_url"],
    log:
        "{results}/logs/download_kraken_db.log",
    conda:
        "../envs/download.yaml"
    shell:
        r"""
        (mkdir -p {output}
        wget -qO- {params.url} | tar -xz -C {output} --strip-components=1
        ) > {log} 2>&1
        """
