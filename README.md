# Dormancy-Popgen-Study

Snakemake pipelines for a combined population-genomics /
dormancy-transcriptomics study of *Osmia bicornis* and *O. cornuta*

## Structure

- **`pipelines/population_genomics/`** — WGS-based population genomics: alignment,
  variant calling, pixy/ADMIXTURE/PopLDdecay/nSL population-genetic statistics,
  contamination screening, species-ID verification, SNPRelate PCA, OrthoFinder
  orthogroup inference, and a cross-pipeline synthesis analysis. Self-contained; see
  its own `README.md`.
- **`pipelines/dormancy_transcriptomics/`** — RNA-seq differential expression for a
  temperature x time-point dormancy experiment: STAR/kallisto alignment and
  quantification, contamination screening, species-ID verification, DESeq2/GSEA
  differential expression. Self-contained; see its own `README.md`.


## How to run

**Either pipeline standalone** 

```bash
cd pipelines/population_genomics  # or pipelines/dormancy_transcriptomics
conda activate snakemake
snakemake --use-conda --cores <N> --config raw_fastq_dir=/path/to/fastq
```

**Everything, composed** 

```bash
conda activate snakemake
snakemake all --dry-run --printshellcmds
snakemake all --use-conda --cores <N>
```

On a SLURM cluster, you can use a Snakemake executor plugin/profile for `sbatch`
(`snakemake-executor-plugin-slurm`) with `cluster.account` from each pipeline's
`config/config.yaml`
