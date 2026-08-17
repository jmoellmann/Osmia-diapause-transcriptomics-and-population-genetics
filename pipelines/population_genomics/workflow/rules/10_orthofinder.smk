rule orthofinder:
    input:
        proteomes_dir=config["orthofinder"]["proteomes_dir"],
    output:
        orthogroups=report(
            "{results}/orthofinder/N0.tsv",
            caption="OrthoFinder phylogenetic hierarchical orthogroups (bicornis, "
            "cornuta, + 4 outgroup species).",
        ),
    params:
        out_dir=lambda wc, output: os.path.dirname(output.orthogroups),
        run_dir=lambda wc: f"{wc.results}/orthofinder/_run",
    log:
        "{results}/logs/orthofinder.log",
    threads: 16
    resources:
        mem_mb=32000,
        runtime=240,
    conda:
        "../envs/orthofinder.yaml"
    shell:
        r"""
        (rm -rf {params.run_dir}
        orthofinder -t {threads} -f {input.proteomes_dir} -o {params.run_dir}
        cp {params.run_dir}/Results_*/Phylogenetic_Hierarchical_Orthogroups/N0.tsv {output.orthogroups}
        ) > {log} 2>&1
        """
