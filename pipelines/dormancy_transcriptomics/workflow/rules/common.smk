import os
import pandas as pd


def _require_config(value, key):
    if not value:
        raise WorkflowError(
            f'config value "{key}" is not set. Supply it with --config, e.g.:\n'
            f"  snakemake --config {key}=..."
        )
    return value


_require_config(config.get("raw_fastq_dir"), "raw_fastq_dir")
for _sp in config["species"]:
    for _key in ("genome_fasta_url", "gtf_url", "cds_fasta_url", "annotation_gff_url"):
        _require_config(config["reference"][_sp].get(_key), f"reference.{_sp}.{_key}")

for _sp in config["reference"]:
    config["reference"][_sp]["genome_fasta"] = f"{{results}}/resources/genomes/{_sp}/genome.fa"
    config["reference"][_sp]["gtf"] = f"{{results}}/resources/genomes/{_sp}/annotation.gtf"
    config["reference"][_sp]["cds_fasta"] = f"{{results}}/resources/genomes/{_sp}/cds.fa"
    config["reference"][_sp]["annotation_gff"] = (
        f"{{results}}/resources/genomes/{_sp}/annotation.gff"
    )

KRAKEN_ENABLED = bool(config["kraken"].get("db_url"))
config["kraken"]["db"] = "{results}/resources/kraken_db"

_ea_resources_dir = config["expression_analysis"]["resources_dir"]
for _key, _val in config["expression_analysis"].items():
    if isinstance(_val, str) and "{expression_analysis[resources_dir]}" in _val:
        config["expression_analysis"][_key] = _val.replace(
            "{expression_analysis[resources_dir]}", _ea_resources_dir
        )
    elif isinstance(_val, dict):
        for _sp, _v in _val.items():
            if isinstance(_v, str) and "{expression_analysis[resources_dir]}" in _v:
                config["expression_analysis"][_key][_sp] = _v.replace(
                    "{expression_analysis[resources_dir]}", _ea_resources_dir
                )


samples_df = pd.read_csv(config["samples"], sep="\t", dtype=str).set_index("id", drop=False)
samples_df = samples_df[samples_df["species"].isin(config["species"])]


def all_target_samples():
    return list(samples_df["id"])


def samples_of(species):
    return list(samples_df[samples_df["species"] == species]["id"])


def species_of(sample):
    return samples_df.loc[sample, "species"]


def raw_lane_files(wildcards, read):
    import glob
    pattern = f"{config['raw_fastq_dir']}/{wildcards.sample}/*{read}.fq.gz"
    files = sorted(glob.glob(pattern))
    if not files:
        raise WorkflowError(f"No raw fastq files matching {pattern}")
    return files


def region_coord(species, idx):
    return config["species_id"]["co1_regions"][species][int(idx)]
