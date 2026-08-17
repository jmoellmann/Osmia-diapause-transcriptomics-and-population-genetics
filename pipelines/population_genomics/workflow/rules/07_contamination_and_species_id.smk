rule extract_unmapped_reads:
    input:
        bam="{results}/{species}/bam/sorted/{sample}.bam",
        bai="{results}/{species}/bam/sorted/{sample}.bam.bai",
    output:
        fq=temp("{results}/{species}/contamination/unmapped/{sample}.unmapped.fq"),
    log:
        "{results}/{species}/logs/extract_unmapped_reads/{sample}.log",
    conda:
        "../envs/samtools.yaml"
    shell:
        "samtools fastq -f 4 -@ 4 {input.bam} > {output.fq} 2> {log}"


rule kraken2_classify:
    input:
        fq="{results}/{species}/contamination/unmapped/{sample}.unmapped.fq",
        db=config["kraken"]["db"],
    output:
        report="{results}/{species}/contamination/kraken2/{sample}.kraken.report.out",
        out=temp("{results}/{species}/contamination/kraken2/{sample}.kraken.out"),
    log:
        "{results}/{species}/logs/kraken2_classify/{sample}.log",
    threads: 8
    conda:
        "../envs/kraken2.yaml"
    shell:
        "kraken2 --threads {threads} --db {input.db} --use-names "
        "--report {output.report} --output {output.out} {input.fq} > {log} 2>&1"


rule gather_kraken_reports:
    input:
        reports=lambda wc: [
            f"{wc.results}/{sp}/contamination/kraken2/{s}.kraken.report.out"
            for sp in config["species"] for s in samples_of(sp)
        ],
    output:
        csv="{results}/contamination/kraken2_report_summary.csv",
    log:
        "{results}/logs/gather_kraken_reports.log",
    conda:
        "../envs/coreutils.yaml"
    shell:
        r"""
        (rm -f {output.csv}
        for report in {input.reports}; do
            sample_name=$(basename "$report" | cut -d "_" -f 1,2)
            n=$(wc -l < "$report")
            cat "$report" | tr -s " " | sed 's/\t/,/g' | sed 's/ //1' \
                | paste -d, <(printf "%${{n}}s" | sed "s/ /${{sample_name}}\n/g") - \
                >> {output.csv}
        done) > {log} 2>&1
        """


rule kraken2_contamination_summary:
    input:
        kraken_summary="{results}/contamination/kraken2_report_summary.csv",
        samples=config["samples"],
    output:
        summary_tsv="{results}/contamination/genus_summary.tsv",
        count_png=report(
            "{results}/contamination/count_clade_by_genus.png",
            caption="Kraken2 count_clade by genus, filtered to genera with >=5 samples.",
        ),
        perc_png=report(
            "{results}/contamination/percent_by_genus.png",
            caption="Kraken2 read percentage by genus, filtered to genera with >=5 samples.",
        ),
    log:
        "{results}/logs/kraken2_contamination_summary.log",
    conda:
        "../envs/r-popgen.yaml"
    script:
        "../scripts/kraken2_report_summary.R"


rule co1_consensus:
    input:
        bam="{results}/{species}/bam/sorted/{sample}.bam",
        bai="{results}/{species}/bam/sorted/{sample}.bam.bai",
    output:
        fa="{results}/{species}/species_id/consensus/{sample}.region{region_idx}.consensus.fa",
    params:
        region=lambda wc: region_coord(wc.species, wc.region_idx),
    log:
        "{results}/{species}/logs/co1_consensus/{sample}.region{region_idx}.log",
    conda:
        "../envs/samtools.yaml"
    shell:
        "samtools consensus -r {params.region} {input.bam} > {output.fa} 2> {log}"


rule co1_consensus_revcomp:
    input:
        fa="{results}/{species}/species_id/consensus/{sample}.region{region_idx}.consensus.fa",
    output:
        fa="{results}/{species}/species_id/consensus/{sample}.region{region_idx}.revcomp.fa",
    params:
        script=workflow.source_path("../scripts/reverse_complement.py"),
    log:
        "{results}/{species}/logs/co1_consensus_revcomp/{sample}.region{region_idx}.log",
    conda:
        "../envs/biopython.yaml"
    shell:
        "python {params.script} {input.fa} {output.fa} > {log} 2>&1"


rule query_bold:
    input:
        fa="{results}/{species}/species_id/consensus/{sample}.region{region_idx}.{orientation}.fa",
    output:
        xml="{results}/{species}/species_id/bold_xml/{sample}.region{region_idx}.{orientation}.out.xml",
    wildcard_constraints:
        orientation="consensus|revcomp",
    log:
        "{results}/{species}/logs/query_bold/{sample}.region{region_idx}.{orientation}.log",
    conda:
        "../envs/wget.yaml"
    shell:
        r"""
        (seq=$(grep -v ">" {input.fa} | xargs echo -n | tr -d '[:blank:]')
        wget -O {output.xml} "http://www.boldsystems.org/index.php/Ids_xml?db=COX1_SPECIES_PUBLIC&sequence=${{seq}}"
        ) > {log} 2>&1
        """


rule summarise_bold_output:
    input:
        xml=lambda wc: [
            f"{wc.results}/{sp}/species_id/bold_xml/{s}.region{ri}.{orient}.out.xml"
            for sp in config["species"] for s in samples_of(sp)
            for ri in range(len(config["species_id"]["co1_regions"][sp]))
            for orient in ["consensus", "revcomp"]
        ],
    output:
        tsv="{results}/species_id/CO1_regions.tsv",
    params:
        min_size=40,
        script=workflow.source_path("../scripts/parse_bold_output_xml.py"),
    log:
        "{results}/logs/summarise_bold_output.log",
    conda:
        "../envs/biopython.yaml"
    shell:
        r"""
        (rm -f {output.tsv}
        for xml in {input.xml}; do
            file_size=$(wc -c < "$xml")
            if [ "$file_size" -gt {params.min_size} ]; then
                fname=$(basename "$xml")
                co1_region=$(echo "$fname" | cut -d "." -f 2)
                sample=$(echo "$fname" | cut -d "." -f 1)
                best=$(python {params.script} "$xml" | head -n 1)
                echo -e "${{sample}}\t${{co1_region}}\t${{best}}" >> {output.tsv}
            fi
        done) > {log} 2>&1
        """


rule species_id_final_summary:
    input:
        tsv="{results}/species_id/CO1_regions.tsv",
    output:
        summary="{results}/species_id/summary.txt",
    log:
        "{results}/logs/species_id_final_summary.log",
    conda:
        "../envs/coreutils.yaml"
    shell:
        r"""grep "region0" {input.tsv} | cut -f 3 | sort | uniq -c > {output.summary} 2> {log}"""
