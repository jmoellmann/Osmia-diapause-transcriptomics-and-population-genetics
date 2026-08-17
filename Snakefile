"""
Top-level composed workflow — dormancy-popgen-combined.

Wires together the two self-contained pipelines under pipelines/ (population_genomics,
dormancy_transcriptomics) for the cross-pipeline dependencies:

  - population_genomics' `orthofinder` rule (comparative-genomics orthogroups, bicornis
    + cornuta + 4 outgroup species) feeds dormancy_transcriptomics'
    `expression_analysis_{bicornis,cornuta}` rules (GO/KEGG annotation via Drosophila
    orthology).
  - dormancy_transcriptomics' `expression_analysis_{bicornis,cornuta}` rules (DESeq2
    results + GO term lists) feed back into population_genomics' `popgen_analysis`
    rule (nSL-vs-differential-expression synthesis figure).

"""


module population_genomics:
    snakefile: "pipelines/population_genomics/Snakefile"
    prefix: "pipelines/population_genomics"

use rule * from population_genomics as population_genomics_*


module dormancy_transcriptomics:
    snakefile: "pipelines/dormancy_transcriptomics/Snakefile"
    prefix: "pipelines/dormancy_transcriptomics"

use rule * from dormancy_transcriptomics as dormancy_transcriptomics_*


# --- Cross-pipeline overrides -------------------------------------------------------


use rule expression_analysis_bicornis from dormancy_transcriptomics as \
        dormancy_transcriptomics_expression_analysis_bicornis with:
    input:
        kallisto=lambda wc: expand(
            "pipelines/dormancy_transcriptomics/{results}/bicornis/kallisto/{sample}/abundance.h5",
            results=dormancy_transcriptomics.config["results_dir"],
            sample=dormancy_transcriptomics.samples_of("bicornis"),
        ),
        tx2gene=dormancy_transcriptomics.config["expression_analysis"]["tx2gene"]["bicornis"],
        orthogroups="pipelines/population_genomics/{results}/orthofinder/N0.tsv".format(
            results=population_genomics.config["results_dir"]
        ),
        dmel_flybase2ncbi=dormancy_transcriptomics.config["expression_analysis"]["dmel_flybase2ncbi"],
        dmel_gene_summaries=dormancy_transcriptomics.config["expression_analysis"]["dmel_gene_summaries"],
        dme2keggpw=dormancy_transcriptomics.config["expression_analysis"]["dme2keggpw"],
        ncbi2dme=dormancy_transcriptomics.config["expression_analysis"]["ncbi2dme"],
        dm_ncbi_to_kegg=dormancy_transcriptomics.config["expression_analysis"]["dm_ncbi_to_kegg"],
        go_basic=dormancy_transcriptomics.config["expression_analysis"]["go_basic"],
        gene2go=dormancy_transcriptomics.config["expression_analysis"]["gene2go"],
        dmel_gene2go=dormancy_transcriptomics.config["expression_analysis"]["dmel_gene2go"],
        annotations="pipelines/dormancy_transcriptomics/{results}/bicornis/expression_analysis/annotations.tsv",


use rule expression_analysis_cornuta from dormancy_transcriptomics as \
        dormancy_transcriptomics_expression_analysis_cornuta with:
    input:
        kallisto=lambda wc: expand(
            "pipelines/dormancy_transcriptomics/{results}/cornuta/kallisto/{sample}/abundance.h5",
            results=dormancy_transcriptomics.config["results_dir"],
            sample=dormancy_transcriptomics.samples_of("cornuta"),
        ),
        tx2gene=dormancy_transcriptomics.config["expression_analysis"]["tx2gene"]["cornuta"],
        orthogroups="pipelines/population_genomics/{results}/orthofinder/N0.tsv".format(
            results=population_genomics.config["results_dir"]
        ),
        dmel_flybase2ncbi=dormancy_transcriptomics.config["expression_analysis"]["dmel_flybase2ncbi"],
        dmel_gene_summaries=dormancy_transcriptomics.config["expression_analysis"]["dmel_gene_summaries"],
        dme2keggpw=dormancy_transcriptomics.config["expression_analysis"]["dme2keggpw"],
        ncbi2dme=dormancy_transcriptomics.config["expression_analysis"]["ncbi2dme"],
        dm_ncbi_to_kegg=dormancy_transcriptomics.config["expression_analysis"]["dm_ncbi_to_kegg"],
        go_basic=dormancy_transcriptomics.config["expression_analysis"]["go_basic"],
        gene2go=dormancy_transcriptomics.config["expression_analysis"]["gene2go"],
        dmel_gene2go=dormancy_transcriptomics.config["expression_analysis"]["dmel_gene2go"],
        annotations="pipelines/dormancy_transcriptomics/{results}/cornuta/expression_analysis/annotations.tsv",


use rule popgen_analysis from population_genomics as population_genomics_popgen_analysis with:
    input:
        sequenced_samples_metadata=population_genomics.config["popgen_analysis"]["sequenced_samples_metadata"],
        raw_read_counts="pipelines/population_genomics/{results}/qc/multiqc_before_trimming/multiqc_data/fastqc_sequence_counts_plot.txt".format(
            results=population_genomics.config["results_dir"]
        ),
        bicornis_10kb_pi="pipelines/population_genomics/{results}/bicornis/pixy_10kb/10kb_pi.txt",
        bicornis_10kb_theta="pipelines/population_genomics/{results}/bicornis/pixy_10kb/10kb_watterson_theta.txt",
        bicornis_annotation_gff="pipelines/population_genomics/{results}/resources/genomes/bicornis/annotation.gff",
        bicornis_cds_pi="pipelines/population_genomics/{results}/bicornis/pixy_cds/cds_pi.txt",
        cornuta_10kb_pi="pipelines/population_genomics/{results}/cornuta/pixy_10kb/10kb_pi.txt",
        cornuta_10kb_theta="pipelines/population_genomics/{results}/cornuta/pixy_10kb/10kb_watterson_theta.txt",
        cornuta_annotation_gff="pipelines/population_genomics/{results}/resources/genomes/cornuta/annotation.gff",
        cornuta_cds_pi="pipelines/population_genomics/{results}/cornuta/pixy_cds/cds_pi.txt",
        bicornis_nsl_chr_lengths=population_genomics.config["reference"]["bicornis"]["chromosome_lengths"],
        bicornis_nsl="pipelines/population_genomics/{results}/bicornis/nsl/combined.nsl.norm.tsv",
        cornuta_nsl_chr_lengths=population_genomics.config["reference"]["cornuta"]["chromosome_lengths"],
        cornuta_nsl="pipelines/population_genomics/{results}/cornuta/nsl/combined.nsl.norm.tsv",
        expression_results_bicornis="pipelines/dormancy_transcriptomics/{results}/bicornis/expression_analysis/expression-results-bicornis.Rdata",
        expression_results_cornuta="pipelines/dormancy_transcriptomics/{results}/cornuta/expression_analysis/expression-results-cornuta.Rdata",
        go_bicornis="pipelines/dormancy_transcriptomics/{results}/bicornis/expression_analysis/bicornis-go.RData",
        go_cornuta="pipelines/dormancy_transcriptomics/{results}/cornuta/expression_analysis/cornuta-go.RData",


rule all:
    input:
        rules.population_genomics_all.input,
        rules.dormancy_transcriptomics_all.input,
