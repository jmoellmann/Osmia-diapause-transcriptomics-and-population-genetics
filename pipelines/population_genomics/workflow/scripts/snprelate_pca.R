# SNPRelate PCA (population structure), both species combined into one figure, plus a
# bicornis-only per-chromosome PCA breakdown. Only LD-pruned, male samples are used;
# a couple of known outlier samples per species are excluded explicitly below.
library(SNPRelate)
library(dplyr)
library(stringr)
library(ggplot2)
library(ggrepel)
library(cowplot)

samples <- read.delim(snakemake@input[["samples"]], sep = "\t") %>%
  rename(ID = metadata_id) %>%
  mutate(ID = str_replace_all(ID, "-", "_"))

metadata_cornuta <- samples %>% filter(species == "cornuta") %>%
  dplyr::select(ID, location, sex, latitude, longitude, elevation, habitat, avg_temp, bee_hotel)
metadata_bicornis <- samples %>% filter(species == "bicornis") %>%
  dplyr::select(ID, location, sex, latitude, longitude, elevation, habitat, avg_temp, bee_hotel)

genofile_path_cornuta <- snakemake@output[["gds_cornuta"]]
genofile_path_bicornis <- snakemake@output[["gds_bicornis"]]

snpgdsVCF2GDS(snakemake@input[["vcf_cornuta"]], genofile_path_cornuta, method = "biallelic.only")
snpgdsVCF2GDS(snakemake@input[["vcf_bicornis"]], genofile_path_bicornis, method = "biallelic.only")

snpgdsSummary(genofile_path_cornuta)
genofile_cornuta <- snpgdsOpen(genofile_path_cornuta)

outliers <- c("KLN1_01", "MUNH1_01")
males_cornuta <- metadata_cornuta %>% filter(sex == "male", !(ID %in% outliers)) %>% pull(ID)

set.seed(1000)
selected_snps_cornuta <- snpgdsSelectSNP(genofile_cornuta, maf = 0.05,
                                          sample.id = males_cornuta, autosome.only = FALSE)
snpset_hardpruned_cornuta <- snpgdsLDpruning(genofile_cornuta,
                                              snp.id = selected_snps_cornuta,
                                              sample.id = males_cornuta,
                                              ld.threshold = 0.2, autosome.only = FALSE)
hardpruned_snps_cornuta <- unlist(snpset_hardpruned_cornuta, use.names = FALSE)
cat("cornuta SNPs after LD pruning:", length(hardpruned_snps_cornuta), "\n")

pca_pruned_cornuta <- snpgdsPCA(genofile_cornuta, num.thread = 2,
                                 snp.id = hardpruned_snps_cornuta,
                                 sample.id = males_cornuta, autosome.only = FALSE)
pc.percent_cornuta <- pca_pruned_cornuta$varprop * 100

tab_cornuta <- tibble(sample.id = as.factor(pca_pruned_cornuta$sample.id),
                       EV1 = pca_pruned_cornuta$eigenvect[, 1],
                       EV2 = pca_pruned_cornuta$eigenvect[, 2])
joined_cornuta <- tab_cornuta %>% left_join(metadata_cornuta, by = c("sample.id" = "ID"))

pca_cornuta <- ggplot(joined_cornuta, aes(x = EV1 * 100, y = EV2 * 100, color = latitude)) +
  geom_point(pch = 21, size = 1.5, stroke = 0.8, fill = "#7D8AED") +
  scale_x_continuous(name = paste0("PC1 (", round(pc.percent_cornuta[1], 2), "% variance)")) +
  scale_y_continuous(name = paste0("PC2 (", round(pc.percent_cornuta[2], 2), "% variance)")) +
  scale_color_gradientn(colours = c("grey70", "black"), name = "Latitude [°N]: ") +
  theme_bw() +
  theme(legend.position = "right",
        legend.text = element_text(angle = 90, vjust = 2, hjust = 0.5),
        legend.title.position = "bottom",
        legend.title = element_text(angle = 90),
        legend.key.height = unit(5, "mm"),
        legend.key.width = unit(3, "mm"),
        panel.border = element_rect(color = "grey25", fill = NA, linewidth = 0.5))

snpgdsSummary(genofile_path_bicornis)
genofile_bicornis <- snpgdsOpen(genofile_path_bicornis)

sample_ids_bicornis <- read.gdsn(index.gdsn(genofile_bicornis, "sample.id"))
snp_ids_bicornis <- read.gdsn(index.gdsn(genofile_bicornis, "snp.id"))
chromosomes_bicornis <- read.gdsn(index.gdsn(genofile_bicornis, "snp.chromosome"))
unique_chrs_bicornis <- sort(unique(chromosomes_bicornis))

outliers_bicornis <- c("XNB_68")
males_bicornis <- metadata_bicornis %>% filter(sex == "male", !(ID %in% outliers_bicornis)) %>% pull(ID)

set.seed(1000)
selected_snps_bicornis <- snpgdsSelectSNP(genofile_bicornis, maf = 0.05,
                                           sample.id = males_bicornis, autosome.only = FALSE)
snpset_hardpruned_bicornis <- snpgdsLDpruning(genofile_bicornis,
                                               snp.id = selected_snps_bicornis,
                                               sample.id = males_bicornis,
                                               ld.threshold = 0.2, autosome.only = FALSE)
hardpruned_snps_bicornis <- unlist(snpset_hardpruned_bicornis, use.names = FALSE)
cat("bicornis SNPs after LD pruning:", length(hardpruned_snps_bicornis), "\n")

pca_pruned_bicornis <- snpgdsPCA(genofile_bicornis, num.thread = 2,
                                  snp.id = hardpruned_snps_bicornis,
                                  sample.id = males_bicornis, autosome.only = FALSE)
pc.percent_bicornis <- pca_pruned_bicornis$varprop * 100

tab_bicornis <- tibble(sample.id = as.factor(pca_pruned_bicornis$sample.id),
                        EV1 = pca_pruned_bicornis$eigenvect[, 1],
                        EV2 = pca_pruned_bicornis$eigenvect[, 2])
joined_bicornis <- tab_bicornis %>% left_join(metadata_bicornis, by = c("sample.id" = "ID"))

pca_bicornis <- ggplot(joined_bicornis, aes(x = EV1 * 100, y = EV2 * 100, color = latitude)) +
  geom_point(shape = 21, size = 1.5, stroke = 0.8, fill = "#F5B829") +
  scale_x_continuous(name = paste0("PC1 (", round(pc.percent_bicornis[1], 2), "% variance)")) +
  scale_y_continuous(name = paste0("PC2 (", round(pc.percent_bicornis[2], 2), "% variance)")) +
  scale_color_gradientn(colours = c("grey70", "black"), name = "Latitude [°N]: ") +
  theme_bw() +
  theme(legend.position = "none",
        panel.border = element_rect(color = "grey25", fill = NA, linewidth = 0.5))

combined <- plot_grid(pca_bicornis, pca_cornuta, ncol = 2, labels = c("C", "D"),
                      rel_widths = c(1, 1.15))
ggsave(filename = snakemake@output[["combined_png"]], width = 183, height = 70,
       units = "mm", plot = combined, dpi = 600)

# Per-chromosome PCA breakdown, bicornis only.
dir.create(snakemake@output[["per_chr_dir"]], showWarnings = FALSE, recursive = TRUE)

pcas_unpruned <- list()
for (chromosome in unique_chrs_bicornis) {
  chr_snp_ids <- snp_ids_bicornis[chromosomes_bicornis == chromosome]
  pca_unpruned <- snpgdsPCA(genofile_bicornis, num.thread = 2, snp.id = chr_snp_ids,
                             sample.id = males_bicornis, autosome.only = FALSE)
  pc.percent <- pca_unpruned$varprop * 100
  tab <- tibble(sample.id = as.factor(pca_unpruned$sample.id),
                EV1 = pca_unpruned$eigenvect[, 1], EV2 = pca_unpruned$eigenvect[, 2])
  joined <- tab %>% left_join(metadata_bicornis, by = c("sample.id" = "ID"))
  chromosome_newname <- as.numeric(str_extract(str_remove_all(chromosome, "NC_0602"), "\\d+")) - 15
  pcas_unpruned[[chromosome]] <- joined %>%
    ggplot(aes(x = EV1 * 100, y = EV2 * 100, label = sample.id, color = latitude)) +
    geom_point(size = 2, pch = 21, stroke = 0.5, fill = "#F5B829") +
    geom_text_repel(aes(color = latitude), size = 3, force = 1.5, max.overlaps = 600) +
    xlab(paste0("PC1 (", round(pc.percent[1], 2), "% variance)")) +
    ylab(paste0("PC2 (", round(pc.percent[2], 2), "% variance)")) +
    ggtitle(paste0("Chromosome ", chromosome_newname)) +
    theme_bw() +
    theme(plot.title = element_text(hjust = 0.5, face = "bold"),
          panel.border = element_rect(color = "black", fill = NA, linewidth = 0.5))
}
for (chromosome in names(pcas_unpruned)) {
  ggsave(plot = pcas_unpruned[[chromosome]],
         filename = file.path(snakemake@output[["per_chr_dir"]], paste0(chromosome, "_unpruned_pca.png")))
}

pcas_pruned <- list()
for (chromosome in unique_chrs_bicornis) {
  chr_snp_ids <- snp_ids_bicornis[chromosomes_bicornis == chromosome]
  pca_pruned <- snpgdsPCA(genofile_bicornis, num.thread = 2,
                           snp.id = intersect(chr_snp_ids, hardpruned_snps_bicornis),
                           sample.id = males_bicornis, autosome.only = FALSE)
  pc.percent <- pca_pruned$varprop * 100
  tab <- tibble(sample.id = as.factor(pca_pruned$sample.id),
                EV1 = pca_pruned$eigenvect[, 1], EV2 = pca_pruned$eigenvect[, 2])
  joined <- tab %>% left_join(metadata_bicornis, by = c("sample.id" = "ID"))
  chromosome_newname <- as.numeric(str_extract(str_remove_all(chromosome, "NC_0602"), "\\d+")) - 15
  pcas_pruned[[chromosome]] <- joined %>%
    ggplot(aes(x = EV1 * 100, y = EV2 * 100, label = sample.id, color = latitude)) +
    geom_point(size = 2, pch = 21, stroke = 0.5, fill = "#F5B829") +
    geom_text_repel(aes(color = latitude), size = 3, force = 1.5, max.overlaps = 600) +
    xlab(paste0("PC1 (", round(pc.percent[1], 2), "% variance)")) +
    ylab(paste0("PC2 (", round(pc.percent[2], 2), "% variance)")) +
    ggtitle(paste0("Chromosome ", chromosome_newname, " (LD-pruned)")) +
    scale_x_reverse() +
    theme_bw() +
    theme(plot.title = element_text(hjust = 0.5, face = "bold"),
          panel.border = element_rect(color = "black", fill = NA, linewidth = 0.5))
}
for (chromosome in names(pcas_pruned)) {
  ggsave(plot = pcas_pruned[[chromosome]],
         filename = file.path(snakemake@output[["per_chr_dir"]], paste0(chromosome, "_pruned_pca.png")))
}

# Supplemental figure: two representative chromosomes (4 and 10), pruned vs. unpruned.
ggsave(plot = plot_grid(
    pcas_unpruned[[4]] + guides(color = "none"), pcas_pruned[[4]],
    pcas_unpruned[[10]] + guides(color = "none"), pcas_pruned[[10]],
    ncol = 2, rel_widths = c(1, 1.25), labels = c("A", "B", "C", "D")),
  filename = snakemake@output[["supplemental_png"]], width = 183, height = 160,
  unit = "mm", dpi = 600, scale = 1.2)

snpgdsClose(genofile_cornuta)
snpgdsClose(genofile_bicornis)
