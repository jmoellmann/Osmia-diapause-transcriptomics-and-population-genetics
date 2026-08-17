rule popgen_analysis:
    input:
        sequenced_samples_metadata=config["popgen_analysis"]["sequenced_samples_metadata"],
        raw_read_counts="{results}/qc/multiqc_before_trimming/multiqc_data/fastqc_sequence_counts_plot.txt",
        bicornis_10kb_pi="{results}/bicornis/pixy_10kb/10kb_pi.txt",
        bicornis_10kb_theta="{results}/bicornis/pixy_10kb/10kb_watterson_theta.txt",
        bicornis_annotation_gff=config["reference"]["bicornis"]["annotation_gff"],
        bicornis_cds_pi="{results}/bicornis/pixy_cds/cds_pi.txt",
        cornuta_10kb_pi="{results}/cornuta/pixy_10kb/10kb_pi.txt",
        cornuta_10kb_theta="{results}/cornuta/pixy_10kb/10kb_watterson_theta.txt",
        cornuta_annotation_gff=config["reference"]["cornuta"]["annotation_gff"],
        cornuta_cds_pi="{results}/cornuta/pixy_cds/cds_pi.txt",
        bicornis_nsl_chr_lengths=config["reference"]["bicornis"]["chromosome_lengths"],
        bicornis_nsl="{results}/bicornis/nsl/combined.nsl.norm.tsv",
        cornuta_nsl_chr_lengths=config["reference"]["cornuta"]["chromosome_lengths"],
        cornuta_nsl="{results}/cornuta/nsl/combined.nsl.norm.tsv",
        expression_results_bicornis="../dormancy_transcriptomics/results/bicornis/expression_analysis/expression-results-bicornis.Rdata",
        expression_results_cornuta="../dormancy_transcriptomics/results/cornuta/expression_analysis/expression-results-cornuta.Rdata",
        go_bicornis="../dormancy_transcriptomics/results/bicornis/expression_analysis/bicornis-go.RData",
        go_cornuta="../dormancy_transcriptomics/results/cornuta/expression_analysis/cornuta-go.RData",
    output:
        bicornis_nsl_manhattan_png=report(
            "{results}/popgen_analysis/bicornis-nsl-manhattan.png",
            caption="nSL selection-scan Manhattan plot, O. bicornis.",
        ),
        cornuta_nsl_manhattan_png=report(
            "{results}/popgen_analysis/cornuta-nsl-manhattan.png",
            caption="nSL selection-scan Manhattan plot, O. cornuta.",
        ),
        combined_nsl_manhattan_png=report(
            "{results}/popgen_analysis/combined_nsl_manhattan.png",
            caption="Combined nSL Manhattan plot, both species.",
        ),
        combined_nsl_vs_lfc_png=report(
            "{results}/popgen_analysis/combined_nsl_vs_lfc.png",
            caption="nSL vs. dormancy-experiment log2FoldChange, both species.",
        ),
        combined_nsl_all_png=report(
            "{results}/popgen_analysis/combined_nsl_all.png",
            caption="Combined nSL Manhattan + nSL-vs-LFC figure, both species.",
        ),
        figure5_png=report(
            "{results}/popgen_analysis/figure5.png",
            caption="GO-enrichment synthesis figure (manuscript Figure 5): nSL-selected "
            "genes vs. dormancy-experiment DE results, both species.",
        ),
    log:
        "{results}/logs/popgen_analysis.log",
    conda:
        "../envs/r-popgen-analysis.yaml"
    script:
        "../scripts/popgen_analysis.R"
