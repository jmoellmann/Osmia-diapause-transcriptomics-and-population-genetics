# Population genomics pipeline

WGS-based population genomics for *Osmia bicornis* and *O. cornuta*: read trimming,
per-species bowtie2 alignment (RefSeq assembly for bicornis, 
long-read Möllmann et al. 2026 assembly for cornuta), BAM preparation, 
ploidy checking, joint variant calling, population-genetic statistics
(pixy π/Fst/Dxy/Watterson's θ/Tajima's D, ADMIXTURE ancestry, PopLDdecay linkage
decay, nSL selection scan), contamination screening, species-ID verification,
SNPRelate PCA, OrthoFinder orthogroup inference.

This pipeline only works in full once the RNAseq sibling pipeline has been run.


## Inputs

- **Raw paired-end WGS reads**, one directory per sample. `raw_fastq_dir` has no
  default — raw sequencing data, supplied every run via 
  `snakemake --config raw_fastq_dir=/path/to/fastq`. Samples may be split across
  multiple sequencing lanes/runs; all matching files per sample/read-mate are
  concatenated.
- **Reference genomes**: bicornis (NCBI RefSeq GCF_907164935.1_iOsmBic2.1) and 
cornuta (Möllmann et al. 2026, see https://doi.org/10.1093/gbe/evaf224 and 
https://github.com/LIB-insect-comparative-genomics/Osmia.cornuta.genome). Both are 
downloaded by `00_download_reference.smk`.
- **Sample metadata**: `config/samples.tsv`

## Workflow

1. **Fastq prep**: merge lanes per sample, adapter/quality trim with fastp.
2. **Alignment**: bowtie2 (both species), read-group tagging, sort, index.
3. **Ploidy check**: joint diploid calling per chromosome across all samples,
   QUAL/DP filtering, heterozygosity and depth stats, plot.
4. **Final variant calling**: joint calling per chromosome using each sample's actual
   ploidy, filtering (QUAL/mean-depth/missingness), subsetting to SNPs/indels/both.
5. **Population genetics**: pixy in 10kb windows and CDS regions; ADMIXTURE (K=1–8)
   on LD-pruned haploid-sample SNPs; PopLDdecay; general VCF QC stats; nSL selection
   scan (per-chromosome selscan, both species).
6. **Contamination screening**: Kraken2 on unmapped reads, aggregated genus-level
   summary + plots.
7. **Species-ID verification**: CO1 barcode consensus (2 regions/species, both
   orientations) → BOLD Systems query → best-match summary.
8. **SNPRelate PCA**: population-structure PCA (both species combined + bicornis
   per-chromosome breakdown).
9. **OrthoFinder**: comparative-genomics orthogroup inference across bicornis, cornuta,
   and 4 outgroup species (feeds `../dormancy_transcriptomics`' DE scripts).
10. **`popgen_analysis`**: nucleotide-diversity + nSL Manhattan plots, and a GO-enrichment
    synthesis figure (`figure5.png`) combining nSL-selected genes with differential
    expression results from `../dormancy_transcriptomics`.
11. **QC reporting**: MultiQC summaries for raw (pre-trim) read quality
    (`multiqc_raw`) and bowtie2 mapping rates (`multiqc_bowtie2`), both species
    combined

## Running

```bash
cd pipelines/population_genomics
conda activate snakemake     # or: mamba env create for the pipeline's own envs via --use-conda
snakemake --use-conda --cores <N>
```

On the SLURM cluster, you can use a Snakemake executor plugin/profile for `sbatch`
(`snakemake-executor-plugin-slurm`) with `cluster.account` from `config/config.yaml`.
