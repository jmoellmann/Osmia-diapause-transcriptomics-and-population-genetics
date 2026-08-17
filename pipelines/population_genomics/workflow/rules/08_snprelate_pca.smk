rule snprelate_pca:
    input:
        vcf_cornuta="{results}/cornuta/final_calling/variants.filtered.snps-only.vcf.gz",
        vcf_bicornis="{results}/bicornis/final_calling/variants.filtered.snps-only.vcf.gz",
        samples=config["samples"],
    output:
        gds_cornuta="{results}/cornuta/snprelate/variants.filtered.snps-only.gds",
        gds_bicornis="{results}/bicornis/snprelate/variants.filtered.snps-only.gds",
        combined_png=report(
            "{results}/snprelate/combined-pca.png",
            caption="PCA (LD-pruned, male samples) for both species.",
        ),
        per_chr_dir=directory("{results}/snprelate/per_chromosome"),
        supplemental_png=report(
            "{results}/snprelate/supplemental_figure3.png",
            caption="Per-chromosome PCA (2 representative chromosomes), bicornis, pruned vs. unpruned.",
        ),
    log:
        "{results}/logs/snprelate_pca.log",
    conda:
        "../envs/r-snprelate.yaml"
    script:
        "../scripts/snprelate_pca.R"
