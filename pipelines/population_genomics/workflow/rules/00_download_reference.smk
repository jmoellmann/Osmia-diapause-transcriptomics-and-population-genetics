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


rule download_proteome:
    output:
        fa="{results}/resources/orthofinder/full_proteomes/{proteome}.faa",
    params:
        url=lambda wc: config["orthofinder"]["proteomes"][f"{wc.proteome}_url"],
    log:
        "{results}/logs/download_proteome/{proteome}.log",
    conda:
        "../envs/download.yaml"
    shell:
        "(wget -qO- {params.url} | gunzip -c > {output.fa}) > {log} 2>&1"


rule orthofinder_primary_transcripts:
    input:
        fa="{results}/resources/orthofinder/full_proteomes/{proteome}.faa",
        script=workflow.source_path("../scripts/vendor/primary_transcript.py"),
    output:
        fa="{results}/resources/orthofinder/primary_transcripts/{proteome}.faa",
    params:
        out_dir=lambda wc, input: os.path.join(os.path.dirname(input.fa), "primary_transcripts"),
    log:
        "{results}/logs/orthofinder_primary_transcripts/{proteome}.log",
    conda:
        "../envs/python.yaml"
    shell:
        r"""
        (python {input.script} {input.fa}
        mv {params.out_dir}/{wildcards.proteome}.faa {output.fa}
        ) > {log} 2>&1
        """


rule prepare_orthofinder_proteomes_dir:
    input:
        fa=lambda wc: expand(
            "{results}/resources/orthofinder/primary_transcripts/{proteome}.faa",
            results=wc.results, proteome=ORTHOFINDER_PROTEOMES,
        ),
    output:
        directory("{results}/resources/orthofinder/proteomes_dir"),
    log:
        "{results}/logs/prepare_orthofinder_proteomes_dir.log",
    conda:
        "../envs/coreutils.yaml"
    shell:
        r"""
        (mkdir -p {output}
        cp {input.fa} {output}/
        ) > {log} 2>&1
        """
