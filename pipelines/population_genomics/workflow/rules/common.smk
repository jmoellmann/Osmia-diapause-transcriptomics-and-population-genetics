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
    _require_config(
        config["reference"][_sp].get("genome_fasta_url"), f"reference.{_sp}.genome_fasta_url"
    )
    _require_config(
        config["reference"][_sp].get("annotation_gff_url"),
        f"reference.{_sp}.annotation_gff_url",
    )

for _sp in config["reference"]:
    config["reference"][_sp]["genome_fasta"] = f"{{results}}/resources/genomes/{_sp}/genome.fa"
    config["reference"][_sp]["annotation_gff"] = (
        f"{{results}}/resources/genomes/{_sp}/annotation.gff"
    )

KRAKEN_ENABLED = bool(config["kraken"].get("db_url"))
config["kraken"]["db"] = "{results}/resources/kraken_db"

ORTHOFINDER_PROTEOMES = [
    "bicornis", "cornuta",
    "bombus_terrestris", "drosophila_melanogaster", "megachile_rotundata", "osmia_lignaria",
]
ORTHOFINDER_ENABLED = all(
    config["orthofinder"]["proteomes"].get(f"{_p}_url") for _p in ORTHOFINDER_PROTEOMES
)
config["orthofinder"]["proteomes_dir"] = "{results}/resources/orthofinder/proteomes_dir"


def all_target_samples():
    return list(samples_df["id"])

samples_df = pd.read_csv(config["samples"], sep="\t", dtype=str).set_index("id", drop=False)

samples_df = samples_df[samples_df["species"].isin(config["species"])]


def samples_of(species):
    return list(samples_df[samples_df["species"] == species]["id"])


def species_of(sample):
    return samples_df.loc[sample, "species"]


def aligner_of(sample):
    return config["reference"][species_of(sample)]["aligner"]


def genome_fasta(wildcards):
    return config["reference"][wildcards.species]["genome_fasta"]


def chromosomes_file(wildcards):
    return config["reference"][wildcards.species]["chromosomes"]


def chromosome_list(species):
    with open(config["reference"][species]["chromosomes"]) as f:
        return [line.strip() for line in f if line.strip()]


def raw_lane_files(wildcards, read):
    import glob
    pattern = f"{config['raw_fastq_dir']}/{wildcards.sample}/*{read}.fq.gz"
    files = sorted(glob.glob(pattern))
    if not files:
        raise WorkflowError(f"No raw fastq files matching {pattern}")
    return files


def bams_of_species(wildcards):
    return {
        "bam": expand(
            "{results}/{species}/bam/sorted/{sample}.bam",
            results=wildcards.results, species=wildcards.species,
            sample=samples_of(wildcards.species),
        ),
        "bai": expand(
            "{results}/{species}/bam/sorted/{sample}.bam.bai",
            results=wildcards.results, species=wildcards.species,
            sample=samples_of(wildcards.species),
        ),
    }


def region_coord(species, idx):
    return config["species_id"]["co1_regions"][species][int(idx)]
