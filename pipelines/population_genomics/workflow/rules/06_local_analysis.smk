rule ploidy_check_plot:
    input:
        het="{results}/{species}/ploidy_check/het_check.het",
        depth=lambda wc: config["manual_inputs"]["depth_check"][wc.species],
        samples=config["samples"],
    output:
        het_png=report(
            "{results}/{species}/ploidy_check/het.png",
            caption="Observed vs. expected homozygosity, colored by recorded sex.",
        ),
        depth_png=report(
            "{results}/{species}/ploidy_check/depth_by_sex.png",
            caption="Mean per-individual depth density, colored by recorded sex.",
        ),
    log:
        "{results}/{species}/logs/ploidy_check_plot.log",
    conda:
        "../envs/r-popgen.yaml"
    script:
        "../scripts/ploidy_check_plot.R"


rule admixture_plot_cv:
    input:
        bicornis_summary="{results}/bicornis/admixture/summary.txt",
        cornuta_summary="{results}/cornuta/admixture/summary.txt",
    output:
        png=report(
            "{results}/admixture_cv.png",
            caption="ADMIXTURE cross-validation error by K, both species.",
        ),
    log:
        "{results}/logs/admixture_plot_cv.log",
    conda:
        "../envs/r-popgen.yaml"
    script:
        "../scripts/admixture_plot_cv.R"


rule general_vcf_stats_plots:
    input:
        frq="{results}/{species}/general_stats/stats.frq",
        idepth="{results}/{species}/general_stats/stats.idepth",
        ldepth_mean="{results}/{species}/general_stats/stats.ldepth.mean",
        lqual="{results}/{species}/general_stats/stats.lqual",
        imiss="{results}/{species}/general_stats/stats.imiss",
        lmiss="{results}/{species}/general_stats/stats.lmiss",
        het="{results}/{species}/general_stats/stats.het",
    output:
        lqual_png="{results}/{species}/general_stats/lqual.png",
        ldepth_png="{results}/{species}/general_stats/ldepth.mean.png",
        lmiss_png="{results}/{species}/general_stats/lmiss.png",
        frq_png="{results}/{species}/general_stats/frq.png",
        idepth_png="{results}/{species}/general_stats/idepth.png",
        imiss_png="{results}/{species}/general_stats/imiss.png",
        het_png="{results}/{species}/general_stats/het.png",
    params:
        stats_dir=lambda wc, input: os.path.dirname(input.frq),
        out_dir=lambda wc, output: os.path.dirname(output.lqual_png),
        script=workflow.source_path("../scripts/general_stats.R"),
    log:
        "{results}/{species}/logs/general_vcf_stats_plots.log",
    conda:
        "../envs/r-popgen.yaml"
    shell:
        "Rscript {params.script} {params.stats_dir} {params.out_dir} > {log} 2>&1"
