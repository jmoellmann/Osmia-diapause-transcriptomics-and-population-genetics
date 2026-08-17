rule gff_annotations:
    input:
        gff=lambda wc: config["reference"][wc.species]["annotation_gff"],
    output:
        tsv="{results}/{species}/expression_analysis/annotations.tsv",
    params:
        script=workflow.source_path("../scripts/annotations_from_gff.sh"),
    log:
        "{results}/{species}/logs/gff_annotations.log",
    conda:
        "../envs/coreutils.yaml"
    shell:
        "bash {params.script} {input.gff} {output.tsv} > {log} 2>&1"


rule expression_analysis_bicornis:
    input:
        kallisto=lambda wc: expand(
            "{results}/bicornis/kallisto/{sample}/abundance.h5",
            results=wc.results, sample=samples_of("bicornis"),
        ),
        tx2gene=config["expression_analysis"]["tx2gene"]["bicornis"],
        orthogroups="../population_genomics/results/orthofinder/N0.tsv",
        dmel_flybase2ncbi=config["expression_analysis"]["dmel_flybase2ncbi"],
        dmel_gene_summaries=config["expression_analysis"]["dmel_gene_summaries"],
        dme2keggpw=config["expression_analysis"]["dme2keggpw"],
        ncbi2dme=config["expression_analysis"]["ncbi2dme"],
        dm_ncbi_to_kegg=config["expression_analysis"]["dm_ncbi_to_kegg"],
        go_basic=config["expression_analysis"]["go_basic"],
        gene2go=config["expression_analysis"]["gene2go"],
        dmel_gene2go=config["expression_analysis"]["dmel_gene2go"],
        annotations="{results}/bicornis/expression_analysis/annotations.tsv",
    output:
        go_rdata="{results}/bicornis/expression_analysis/bicornis-go.RData",
        combined_pca_png=report(
            "{results}/bicornis/expression_analysis/combined-pca-bicornis.png",
            caption="Full-data PCA (temperature, time), O. bicornis.",
        ),
        combined_gsea_pdf=report(
            "{results}/bicornis/expression_analysis/combined-gsea-plus6.pdf",
            caption="Combined GSEA plot, +6C timepoint, O. bicornis.",
        ),
        heatmap_combined_png=report(
            "{results}/bicornis/expression_analysis/heatmap_combined.png",
            caption="Combined DEG heatmap, O. bicornis.",
        ),
        expression_results_rdata="{results}/bicornis/expression_analysis/expression-results-bicornis.Rdata",
        res_t1_tsv="{results}/bicornis/expression_analysis/output_tables/res_t1.tsv",
        res_t2_tsv="{results}/bicornis/expression_analysis/output_tables/res_t2.tsv",
        res_t3_tsv="{results}/bicornis/expression_analysis/output_tables/res_t3.tsv",
        res_sexes_tsv="{results}/bicornis/expression_analysis/output_tables/res_sexes.tsv",
        res_sexes_interaction_tsv="{results}/bicornis/expression_analysis/output_tables/res_sexes_interaction.tsv",
        degs_xlsx="{results}/bicornis/expression_analysis/output_tables/bicornis-degs.xlsx",
    params:
        kallisto_dir=lambda wc: f"{wc.results}/bicornis/kallisto",
        metadata_path=config["metadata_tsv"],
        excluded_samples=config["expression_analysis"]["bicornis_excluded_samples"],
    log:
        "{results}/bicornis/logs/expression_analysis.log",
    conda:
        "../envs/r-expression-analysis.yaml"
    script:
        "../scripts/expression_analysis_bicornis.R"


rule expression_analysis_cornuta:
    input:
        kallisto=lambda wc: expand(
            "{results}/cornuta/kallisto/{sample}/abundance.h5",
            results=wc.results, sample=samples_of("cornuta"),
        ),
        tx2gene=config["expression_analysis"]["tx2gene"]["cornuta"],
        orthogroups="../population_genomics/results/orthofinder/N0.tsv",
        dmel_flybase2ncbi=config["expression_analysis"]["dmel_flybase2ncbi"],
        dmel_gene_summaries=config["expression_analysis"]["dmel_gene_summaries"],
        dme2keggpw=config["expression_analysis"]["dme2keggpw"],
        ncbi2dme=config["expression_analysis"]["ncbi2dme"],
        dm_ncbi_to_kegg=config["expression_analysis"]["dm_ncbi_to_kegg"],
        go_basic=config["expression_analysis"]["go_basic"],
        gene2go=config["expression_analysis"]["gene2go"],
        dmel_gene2go=config["expression_analysis"]["dmel_gene2go"],
        annotations="{results}/cornuta/expression_analysis/annotations.tsv",
    output:
        go_rdata="{results}/cornuta/expression_analysis/cornuta-go.RData",
        combined_pca_png=report(
            "{results}/cornuta/expression_analysis/combined-pca-cornuta.png",
            caption="Full-data PCA (temperature, time), O. cornuta.",
        ),
        combined_gsea_pdf=report(
            "{results}/cornuta/expression_analysis/combined-gsea-plus6-cornuta.pdf",
            caption="Combined GSEA plot, +6C timepoint, O. cornuta.",
        ),
        heatmap_combined_png=report(
            "{results}/cornuta/expression_analysis/heatmap_combined_cornuta.png",
            caption="Combined DEG heatmap, O. cornuta.",
        ),
        expression_results_rdata="{results}/cornuta/expression_analysis/expression-results-cornuta.Rdata",
        res_t1_tsv="{results}/cornuta/expression_analysis/output_tables/res_t1.tsv",
        res_t2_tsv="{results}/cornuta/expression_analysis/output_tables/res_t2.tsv",
        res_sexes_tsv="{results}/cornuta/expression_analysis/output_tables/res_sexes.tsv",
        res_sexes_interaction_tsv="{results}/cornuta/expression_analysis/output_tables/res_sexes_interaction.tsv",
        degs_xlsx="{results}/cornuta/expression_analysis/output_tables/cornuta-degs.xlsx",
    params:
        kallisto_dir=lambda wc: f"{wc.results}/cornuta/kallisto",
        metadata_path=config["metadata_tsv"],
    log:
        "{results}/cornuta/logs/expression_analysis.log",
    conda:
        "../envs/r-expression-analysis.yaml"
    script:
        "../scripts/expression_analysis_cornuta.R"
