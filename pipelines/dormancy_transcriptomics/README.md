# Dormancy transcriptomics pipeline

RNA-seq differential expression for *Osmia bicornis* and *O. cornuta* under a
temperature x time-point dormancy experiment: read trimming, per-species alignment
(STAR) + pseudo-alignment/quantification (kallisto), contamination screening, species-
ID verification, and DESeq2/GSEA differential expression analysis.

## Inputs

- **Raw paired-end RNA-seq reads**, one directory per sample. `raw_fastq_dir` needs to 
be supplied on every run via `snakemake --config raw_fastq_dir=/path/to/fastq` and 
linked to the raw fastq data top directory. IDs should match `config/samples.tsv`.
- **Reference genome/annotation**: same two assemblies as `pipelines/population_genomics`
  (bicornis NCBI RefSeq GCF_907164935.1_iOsmBic2.1; cornuta long-read based hifiasm assembly,
  see Möllmann et al. 2026, https://doi.org/10.1093/gbe/evaf224 and 
  https://github.com/LIB-insect-comparative-genomics/Osmia.cornuta.genome). Downloaded by
  `00_download_reference.smk`; see `config/config.yaml` → `reference` for the URLs.
- **`metadata.tsv`**: Some analysis branches require additional metadata. See the example 
metadata.tsv for format requirements.
- **Expression-analysis reference/annotation resources**: GO term definitions,
  Drosophila-to-Osmia orthology/KEGG mapping tables, NCBI gene2go dumps — external
  reference data , tracked in this pipeline's own `resources/` directory
  (`config["expression_analysis"]["resources_dir"]`).

## Workflow

1. **Fastq prep**: merge lanes per sample, FastQC on the raw merged reads
   (`fastqc_raw`), adapter/quality trim with fastp.
2. **Alignment**: STAR (genomic, feeds species-IDing) and kallisto (transcript-level,
   feeds DESeq2) run independently over the same trimmed reads per sample.
3. **Contamination screening**: Kraken2 on trimmed reads — per-sample
   reports (see `workflow/rules/03_contamination.smk`).
4. **Species-ID verification**: CO1 barcode consensus (2 regions/species, both
   orientations) from STAR's genomic BAM -> BOLD Systems query -> best-match summary.
   Same live-API caveat as population_genomics (needs internet connection for calls 
   against the BOLD DB).
5. **Differential expression**: per-species DESeq2 LRT tests (temperature x
   time-point, full model + per-timepoint + sex/interaction), GO/KEGG annotation via
   Drosophila orthology (from population_genomics' OrthoFinder output), fgsea, DEG
   heatmaps, results tables export.
6. **QC reporting**: MultiQC summaries of raw (pre-trim) read quality (`multiqc_raw`,
   both species combined) and
   kallisto pseudoalignment rates (`multiqc_kallisto`, per species). 


## Running

```bash
cd pipelines/dormancy_transcriptomics
conda activate snakemake     # or: mamba env create for the pipeline's own envs via --use-conda
snakemake --use-conda --cores <N>
```

On the SLURM cluster, you can use a Snakemake executor plugin/profile for `sbatch`
(`snakemake-executor-plugin-slurm`) with `cluster.account` from `config/config.yaml`.
