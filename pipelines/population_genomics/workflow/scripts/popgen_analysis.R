library(tidyverse)
library(ggplot2)
library(fgsea)
library(data.table)
library(patchwork)
library(cowplot)

# Nucleotide-diversity and nSL Manhattan plots, nSL-vs-differential-expression
# regressions, and a GO-enrichment synthesis figure combining nSL selection signal
# with dormancy-experiment expression results.

metadata_path <- snakemake@input[["sequenced_samples_metadata"]]
read_counts_path <- snakemake@input[["raw_read_counts"]]
bicornis_10kb_pi_path <- snakemake@input[["bicornis_10kb_pi"]]
bicornis_10kb_theta_path <- snakemake@input[["bicornis_10kb_theta"]]
bicornis_annotation_path <- snakemake@input[["bicornis_annotation_gff"]]
bicornis_cds_pi_path <- snakemake@input[["bicornis_cds_pi"]]
cornuta_10kb_pi_path <- snakemake@input[["cornuta_10kb_pi"]]
cornuta_10kb_theta_path <- snakemake@input[["cornuta_10kb_theta"]]
cornuta_annotation_path <- snakemake@input[["cornuta_annotation_gff"]]
cornuta_cds_pi_path <- snakemake@input[["cornuta_cds_pi"]]
bicornis_chr_lengths_path <- snakemake@input[["bicornis_nsl_chr_lengths"]]
bicornis_nsl_path <- snakemake@input[["bicornis_nsl"]]
cornuta_chr_lengths_path <- snakemake@input[["cornuta_nsl_chr_lengths"]]
cornuta_nsl_path <- snakemake@input[["cornuta_nsl"]]
expression_results_bicornis_path <- snakemake@input[["expression_results_bicornis"]]
expression_results_cornuta_path <- snakemake@input[["expression_results_cornuta"]]
go_bicornis_path <- snakemake@input[["go_bicornis"]]
go_cornuta_path <- snakemake@input[["go_cornuta"]]
bicornis_nsl_manhattan_out <- snakemake@output[["bicornis_nsl_manhattan_png"]]
cornuta_nsl_manhattan_out <- snakemake@output[["cornuta_nsl_manhattan_png"]]
combined_nsl_manhattan_out <- snakemake@output[["combined_nsl_manhattan_png"]]
combined_nsl_vs_lfc_out <- snakemake@output[["combined_nsl_vs_lfc_png"]]
combined_nsl_all_out <- snakemake@output[["combined_nsl_all_png"]]
figure5_out <- snakemake@output[["figure5_png"]]

# Read count summary ------------------------------------------------------
metadata <- read_csv(metadata_path) %>% 
  dplyr::mutate(Sample = str_replace(ID, "-", "_")) %>% 
  dplyr::mutate(Sample= str_replace(Sample, "Ö", "O")) %>% 
  dplyr::mutate(Sample= str_replace(Sample, "Ü", "U"))

read_counts <- read_tsv(read_counts_path) %>% 
  dplyr::mutate(Sample = str_extract(Sample, ".+?_\\d+")) %>% 
  dplyr::left_join(metadata) %>%  
  group_by(species, sex, Sample) %>% summarise(
    total_reads = sum(`Unique Reads`) + sum(`Duplicate Reads`)) %>% 
  arrange(total_reads)

read_counts %>% dplyr::filter(species == "bicornis", sex == "male") %>% summary()
read_counts %>% dplyr::filter(species == "cornuta", sex == "male") %>% summary()

cornuta_len <- 256000000
bicornis_len <- 182000000

read_counts %>% dplyr::filter(species == "bicornis", sex == "male") %>% 
  pull(total_reads) %>% mean() * 150 / 2 / bicornis_len


# Bicornis ----------------------------------------------------------------

bicornis_10kb_pi <- read_tsv(bicornis_10kb_pi_path) %>% 
  dplyr::filter(pop == "haploid") %>% arrange(chromosome, window_pos_1, window_pos_2) %>% 
  mutate(chromosome = as.numeric(
    str_extract(str_remove(chromosome, "NC_0602"), "\\d+")) - 15) %>% 
  mutate(chromosome = factor(chromosome, levels = c(1:16))) %>% 
  dplyr::rename(start = window_pos_1, end = window_pos_2)

bicornis_10kb_pi_nona <- bicornis_10kb_pi %>% dplyr::filter(!is.na(avg_pi))


for (chr in unique(bicornis_10kb_pi_nona$chromosome)) {
  test <- wilcox.test(pull(filter(bicornis_10kb_pi_nona, chromosome == chr), "avg_pi"), 
                      pull(filter(bicornis_10kb_pi_nona, chromosome != chr), "avg_pi"),
                      alternative = "greater")
  print(paste0("Chromosome ", chr, ": ", test$p.value))
}



bicornis_10kb_theta <- read_tsv(bicornis_10kb_theta_path) %>% 
  dplyr::filter(pop == "haploid") %>% arrange(chromosome, window_pos_1, window_pos_2) %>% 
  mutate(chromosome = as.numeric(
    str_extract(str_remove(chromosome, "Ocor"), "\\d+"))) %>% 
  mutate(chromosome = factor(chromosome, levels = c(1:16))) %>% 
  dplyr::rename(start = window_pos_1, end = window_pos_2)

bicornis_10kb_theta %>% pull(avg_watterson_theta) %>% median(na.rm = TRUE)

bicornis_annotation <- read_tsv(
  bicornis_annotation_path, comment = "#", col_names = c(
    "chromosome", "x1", "type", "start", "end", "x2", "strand", "x3", "info")) %>% 
  dplyr::arrange(chromosome, type, start, end) %>% 
  dplyr::mutate(geneID = str_extract(info, "(?<=gene=)[^;]+")) %>% 
  dplyr::select(chromosome, geneID, type, start, end, strand) %>% 
  dplyr::filter(!(type %in% c("CDS", "cDNA_match", "region"))) %>% 
  mutate(chromosome = as.numeric(
    str_extract(str_remove(chromosome, "NC_0602"), "\\d+")) - 15) %>% 
  mutate(chromosome = factor(chromosome, levels = c(1:16)))

# Genome-wide mean of nucleotide diversity
bicornis_10kb_pi %>% pull(avg_pi) %>% mean(na.rm = TRUE)

bicornis_gene_pi <- read_tsv(bicornis_cds_pi_path) %>% 
  dplyr::filter(pop == "haploid") %>% 
  arrange(chromosome, window_pos_1, window_pos_2) %>% 
  mutate(chromosome = as.numeric(
    str_extract(str_remove(chromosome, "NC_0602"), "\\d+")) - 15) %>% 
  mutate(chromosome = factor(chromosome, levels = c(1:16))) %>% 
  dplyr::rename(start = window_pos_1, end = window_pos_2) %>% 
  left_join(bicornis_annotation, by = c("chromosome", "start", "end")) %>% 
  dplyr::filter(type == "exon") %>% dplyr::group_by(geneID) %>% 
  reframe(
    chromosome,
    strand,
    start = min(start),
    end = max(end),
    avg_pi = sum(count_diffs, na.rm = TRUE) / sum(count_comparisons, na.rm = TRUE),
    no_sites = sum(no_sites, na.rm = TRUE),
    count_diffs = sum(count_diffs, na.rm = TRUE),
    count_comparisons = sum(count_comparisons, na.rm = TRUE),
    count_missing = sum(count_missing, na.rm = TRUE),
    n_windows = n()) %>% arrange(chromosome, start, end) %>% 
  dplyr::filter(no_sites > 1000) %>% distinct()

load(expression_results_bicornis_path)

bicornis_gene_pi_res <- bicornis_gene_pi %>% 
  dplyr::full_join(res, by = "geneID") %>% 
  dplyr::filter(!is.na(padj))

wilcox.test(bicornis_gene_pi_res %>% dplyr::filter(pvalue < 0.05) %>% pull(avg_pi), 
            bicornis_gene_pi_res %>% dplyr::filter(pvalue >= 0.05) %>% pull(avg_pi))


bicornis_gene_coords_res <- bicornis_annotation %>% 
  dplyr::filter(type == "gene") %>%
  dplyr::left_join(res, by = c("geneID")) %>% 
  mutate(
    # Calculate the midpoint of the gene feature
    centroid = (start + end) / 2
  ) %>% drop_na()


pixy_labeller <- as_labeller(c(avg_pi = "pi"),
                             default = label_parsed)

bicornis_10kb_pi_plt_rdy <- bicornis_10kb_pi %>% 
  full_join(bicornis_gene_coords_res, suffix = c("", ".y"), by = join_by(
    chromosome, start <= centroid, end > centroid)) %>% 
  group_by(chromosome, start, end) %>% 
  reframe(avg_pi, n = sum(padj < 0.05, na.rm = TRUE)) %>% 
  distinct() %>% 
  mutate(chromosome = as.numeric(factor(chromosome))) %>% 
  mutate(chrom_color_group = case_when(
    n > 0 ~ "window with temperature-responsive genes",
    chromosome %% 2 != 0 ~ "even", 
    chromosome %% 2 == 0 ~ "odd")) %>%
  # Add a new variable for high values above 95% quantile
  mutate(above_quantile = avg_pi > quantile(
    bicornis_10kb_pi$avg_pi, 0.95, na.rm = TRUE),
         below_quantile = avg_pi < quantile(
           bicornis_10kb_pi$avg_pi, 0.05, na.rm = TRUE)) %>%
  mutate(chromosome = factor(chromosome, levels = c(1:16))) %>%
  mutate(id = row_number()) %>% group_by(chromosome) %>% 
  dplyr::filter(!is.na(start), !is.na(end)) %>% 
  dplyr::mutate(chromosome_centroid = (
    min(start) + max(end)) / 2 - 5) %>% 
  ungroup() %>% rowwise() %>% dplyr::mutate(dist_to_centroid = min(
    abs(start - chromosome_centroid), abs(end - chromosome_centroid))) %>% 
  group_by(chromosome) %>% 
  dplyr::mutate(has_centroid = c(dist_to_centroid == min(dist_to_centroid))) %>% 
  ungroup()

bicornis_10kb_pi_plt_rdy %>% dplyr::count(
  chrom_color_group == "window with temperature-responsive genes",
  below_quantile) %>% drop_na() %>% pull(n) %>% 
  matrix(nrow = 2, ncol = 2) %>% 
  fisher.test(alternative = "greater")

bicornis_10kb_pi_plt_rdy %>% dplyr::count(
  chrom_color_group == "window with temperature-responsive genes",
  above_quantile) %>% drop_na() %>% pull(n) %>% 
  matrix(nrow = 2, ncol = 2) %>% 
  fisher.test(alternative = "greater")

bicornis_mean_pi <-  bicornis_10kb_pi %>% dplyr::pull(avg_pi) %>% mean(na.rm = TRUE)

bicornis_pi_manhattan <- bicornis_10kb_pi_plt_rdy %>% drop_na() %>% 
  ggplot(aes(x = id, 
             y = avg_pi)) +
  # Plot points with conditional coloring
  geom_point(data = . %>% filter(chrom_color_group != "window with temperature-responsive genes"),
             aes(color = chrom_color_group), 
             size = 0.5, alpha = 0.8, stroke = 0.5) +
  geom_point(data = . %>% filter(chrom_color_group == "window with temperature-responsive genes"), 
             aes(color = chrom_color_group), 
             size = 0.5, alpha = 0.8, stroke = 0.5) +
  # Add horizontal line at 95% quantile
  geom_hline(yintercept = bicornis_mean_pi, color = "grey10", linetype = "dashed", alpha = 0.5) +
  ylab("Nucl. diversity \u03c0") +
  scale_x_continuous(
    expand = c(0.01, 0.01),
    breaks = bicornis_10kb_pi_plt_rdy$id[bicornis_10kb_pi_plt_rdy$has_centroid],
    labels = unique(bicornis_10kb_pi_plt_rdy$chromosome)) + 
  # Update color scale to include red for high values
  scale_color_manual(values = c(
    "even" = "#EDC24C", "odd" = "#D6A424",
    "window with temperature-responsive genes"  = "#C73A3AF8"), 
    breaks = c("window with temperature-responsive genes")) + 
  guides(color = guide_legend(override.aes = list(size = 2))) +
  theme_bw() +  theme(
    axis.title.x = element_blank(),
    legend.position = "none",
    panel.grid.major.x = element_blank(),
    panel.grid.minor.x = element_line(linewidth = 0.5),
    panel.grid.minor.y = element_blank(),
    panel.grid.major.y = element_blank(),
    panel.border = element_rect(
      color = "black", fill = NA, linewidth = 0.5)) + 
  scale_y_continuous(expand = c(0, 0), limits = c(NA, 0.010), breaks = seq(0, 0.010, 0.002))


# ggsave(filename = "bicornis-pi-manhattan.png", width = 183, height = 65, units = "mm",
#        plot = bicornis_pi_manhattan, dpi = 600)


load(expression_results_bicornis_path)

degs <- res %>% dplyr::filter(padj < 0.05) %>% pull(geneID)
degs_t1 <- res_t1 %>% dplyr::filter(padj < 0.05) %>% pull(geneID)
degs_t2 <- res_t2 %>% dplyr::filter(padj < 0.05) %>% pull(geneID)
degs_t2 <- res_t3 %>% dplyr::filter(padj < 0.05) %>% pull(geneID)


bicornis_gene_pi_res %>%
  dplyr::mutate(bin = floor(start / 10000)) %>% 
  group_by(chromosome, type, bin) %>% 
  summarise(geneID, avg_pi = mean(avg_pi, na.rm = TRUE)) %>% ungroup() %>% 
  dplyr::mutate(tmp = ifelse(geneID %in% degs, "DEG", "background")) %>% 
  ggplot(aes(x = avg_pi, y = tmp)) + geom_violin() +
  stat_summary(fun = median, geom = "point")


bicornis_gene_pi_res %>%
  dplyr::mutate(tmp = ifelse(geneID %in% degs, "DEG", "background")) %>% 
  ggplot(aes(x = avg_pi, y = tmp)) + geom_violin() +
  stat_summary(fun = median, geom = "point")

# wilcox.test(
#   bicornis_max_nsl %>% dplyr::filter(geneID %in% degs) %>% pull(max_nsl),
#   bicornis_max_nsl %>% dplyr::filter(!(geneID %in% degs)) %>% pull(max_nsl))

bicornis_pi_mod_t1 <- bicornis_gene_pi %>% 
  dplyr::left_join(res_t1, by = "geneID") %>% 
  dplyr::filter(pvalue < 0.05) %>% 
  lm(formula = avg_pi ~ log2FoldChange.y) %>% summary()

bicornis_pi_mod_t2 <- bicornis_gene_pi %>% 
  dplyr::left_join(res_t2, by = "geneID") %>% 
  dplyr::filter(pvalue < 0.05) %>% 
  lm(formula = avg_pi ~ log2FoldChange.y) %>% summary()

bicornis_pi_mod_t3 <- bicornis_gene_pi %>% 
  dplyr::left_join(res_t3, by = "geneID") %>% 
  dplyr::filter(pvalue < 0.05) %>% 
  lm(formula = avg_pi ~ log2FoldChange.y) %>% summary()

bicornis_pi_vs_lfc_t1 <- bicornis_gene_pi %>% 
  dplyr::left_join(res_t1, by = "geneID") %>% 
  dplyr::mutate(tmp = factor(sign(log2FoldChange.y))) %>% 
  dplyr::filter(pvalue < 0.05) %>% 
  ggplot(aes(x = log2FoldChange.y, y = avg_pi)) + 
  geom_point(color = "#F5B829", alpha = 0.75) + 
  geom_smooth(method = "lm", color = "grey50", alpha = 0.5) + 
  annotate("text", x = 1.5, y = 0.005, size = 4, label = paste0(
    "P < ", signif(bicornis_pi_mod_t1$coefficients[8], digits = 2), "\n", 
    "R² = ", signif(bicornis_pi_mod_t1$adj.r.squared, digits = 2)), 
  ) +
  theme_bw() + ylab("Nucleotide diversity pi") + 
  scale_x_continuous(expand = c(0.01, 0.01), limits = c(-2, 2)) + 
  scale_y_continuous(expand = c(0.0001, 0.0001), limits = c(0, 0.006)) +
  theme(axis.title.x = element_blank())


bicornis_pi_vs_lfc_t2 <- bicornis_gene_pi %>% 
  dplyr::left_join(res_t2, by = "geneID") %>% 
  dplyr::mutate(tmp = factor(sign(log2FoldChange.y))) %>% 
  dplyr::filter(pvalue < 0.05) %>% 
  ggplot(aes(x = log2FoldChange.y, y = avg_pi)) + 
  geom_point(color = "#F5B829", alpha = 0.75) + 
  geom_smooth(method = "lm", color = "grey50", alpha = 0.5) + 
  annotate("text", x = 1.5, y = 0.005, size = 4, label = paste0(
    "P < ", signif(bicornis_pi_mod_t2$coefficients[8], digits = 2), "\n", 
    "R² = ", signif(bicornis_pi_mod_t2$adj.r.squared, digits = 2)), 
  ) +
  theme_bw() + xlab("log\U2082(Fold change)") + ylab("Nucleotide diversity") + 
  scale_x_continuous(expand = c(0.01, 0.01), limits = c(-2, 2)) + 
  scale_y_continuous(expand = c(0.0001, 0.0001), limits = c(0, 0.006)) +
  theme(axis.title.x = element_blank())


bicornis_pi_vs_lfc_t3 <- bicornis_gene_pi %>% 
  dplyr::left_join(res_t3, by = "geneID") %>% 
  dplyr::mutate(tmp = factor(sign(log2FoldChange.y))) %>% 
  dplyr::filter(pvalue < 0.05) %>% 
  ggplot(aes(x = log2FoldChange.y, y = avg_pi)) + 
  geom_point(color = "#F5B829", alpha = 0.75) + 
  geom_smooth(method = "lm", color = "grey50", alpha = 0.5) + 
  annotate("text", x = 1.5, y = 0.005, size = 4, label = paste0(
    "P < ", signif(bicornis_pi_mod_t3$coefficients[8], digits = 2), "\n", 
    "R² = ", signif(bicornis_pi_mod_t3$adj.r.squared, digits = 2)), 
  ) +
  theme_bw() + xlab("log\U2082(Fold change)") + ylab("Nucleotide diversity") + 
  scale_x_continuous(expand = c(0.01, 0.01), limits = c(-2, 2)) + 
  scale_y_continuous(expand = c(0.0001, 0.0001), limits = c(0, 0.006))


bicornis_combined_pi_vs_lfc <- cowplot::plot_grid(
  bicornis_pi_vs_lfc_t1, bicornis_pi_vs_lfc_t2, 
  bicornis_pi_vs_lfc_t3, nrow = 3, rel_heights = c(1, 1, 1.075),
  labels = c("A", "B", "C"))

# Cornuta -----------------------------------------------------------------

cornuta_10kb_pi <- read_tsv(cornuta_10kb_pi_path) %>% 
  dplyr::filter(pop == "haploid") %>% arrange(chromosome, window_pos_1, window_pos_2) %>% 
  mutate(chromosome = as.numeric(
    str_extract(str_remove(chromosome, "Ocor"), "\\d+"))) %>% 
  mutate(chromosome = factor(chromosome, levels = c(1:16))) %>% 
  dplyr::rename(start = window_pos_1, end = window_pos_2)

cornuta_10kb_pi %>% pull(avg_pi) %>% median(na.rm = TRUE)

cornuta_10kb_theta <- read_tsv(cornuta_10kb_theta_path) %>% 
  dplyr::filter(pop == "haploid") %>% arrange(chromosome, window_pos_1, window_pos_2) %>% 
  mutate(chromosome = as.numeric(
    str_extract(str_remove(chromosome, "Ocor"), "\\d+"))) %>% 
  mutate(chromosome = factor(chromosome, levels = c(1:16))) %>% 
  dplyr::rename(start = window_pos_1, end = window_pos_2)

cornuta_10kb_theta %>% pull(avg_watterson_theta) %>% median(na.rm = TRUE)

cornuta_annotation <- read_tsv(
  cornuta_annotation_path, comment = "#", col_names = c(
    "chromosome", "x1", "type", "start", "end", "x2", "strand", "x3", "info")) %>% 
  dplyr::arrange(chromosome, type, start, end) %>% 
  dplyr::mutate(geneID = str_extract(info, "gene\\d+")) %>% 
  dplyr::select(chromosome, geneID, type, start, end, strand) %>% 
  dplyr::filter(!(type %in% c("CDS", "cDNA_match", "region"))) %>% 
  dplyr::mutate(chromosome = as.integer(str_extract(chromosome, "\\d+"))) %>% 
  mutate(chromosome = factor(chromosome, levels = c(1:16)))


cornuta_gene_pi <- read_tsv(cornuta_cds_pi_path) %>% 
  dplyr::filter(pop == "haploid") %>% 
  arrange(chromosome, window_pos_1, window_pos_2) %>% 
  mutate(chromosome = as.numeric(
    str_extract(str_remove(chromosome, "Ocor"), "\\d+"))) %>% 
  mutate(chromosome = factor(chromosome, levels = c(1:16))) %>% 
  dplyr::rename(start = window_pos_1, end = window_pos_2) %>% 
  left_join(cornuta_annotation, by = c("chromosome", "start", "end")) %>% 
  dplyr::filter(type == "exon") %>% dplyr::group_by(geneID) %>% 
  reframe(
    chromosome,
    strand,
    start = min(start),
    end = max(end),
    avg_pi = sum(count_diffs, na.rm = TRUE) / sum(count_comparisons, na.rm = TRUE),
    no_sites = sum(no_sites, na.rm = TRUE),
    count_diffs = sum(count_diffs, na.rm = TRUE),
    count_comparisons = sum(count_comparisons, na.rm = TRUE),
    count_missing = sum(count_missing, na.rm = TRUE),
    n_windows = n()) %>% arrange(chromosome, start, end) %>% 
  dplyr::filter(no_sites > 1000) %>% distinct()

load(expression_results_cornuta_path)

cornuta_gene_pi_res <- cornuta_gene_pi %>% 
  dplyr::full_join(res, by = "geneID") %>% 
  dplyr::filter(!is.na(padj))

wilcox.test(cornuta_gene_pi_res %>% dplyr::filter(pvalue < 0.05) %>% pull(avg_pi), 
            cornuta_gene_pi_res %>% dplyr::filter(pvalue >= 0.05) %>% pull(avg_pi))


cornuta_gene_coords_res <- cornuta_annotation %>% 
  dplyr::filter(type == "gene") %>%
  dplyr::left_join(res, by = c("geneID")) %>% 
  mutate(
    # Calculate the midpoint of the gene feature
    centroid = (start + end) / 2
  ) %>% drop_na()


pixy_labeller <- as_labeller(c(avg_pi = "pi"),
                             default = label_parsed)

cornuta_10kb_pi_plt_rdy <- cornuta_10kb_pi %>% 
  full_join(cornuta_gene_coords_res, suffix = c("", ".y"), by = join_by(
    chromosome, start <= centroid, end > centroid)) %>% 
  group_by(chromosome, start, end) %>% 
  reframe(avg_pi, n = sum(padj < 0.05, na.rm = TRUE)) %>% 
  distinct() %>% 
  mutate(chromosome = as.numeric(factor(chromosome))) %>% 
  mutate(chrom_color_group = case_when(
    n > 0 ~ "window with temperature-responsive genes",
    chromosome %% 2 != 0 ~ "even", 
    chromosome %% 2 == 0 ~ "odd")) %>%
  # Add a new variable for high values above 95% quantile
  mutate(above_quantile = avg_pi > quantile(
    bicornis_10kb_pi$avg_pi, 0.95, na.rm = TRUE),
    below_quantile = avg_pi < quantile(
      bicornis_10kb_pi$avg_pi, 0.05, na.rm = TRUE)) %>% 
  mutate(chromosome = factor(chromosome, levels = c(1:16))) %>%
  mutate(id = row_number()) %>% group_by(chromosome) %>% 
  dplyr::filter(!is.na(start), !is.na(end)) %>% 
  dplyr::mutate(chromosome_centroid = (
    min(start) + max(end)) / 2 - 5) %>% 
  ungroup() %>% rowwise() %>% dplyr::mutate(dist_to_centroid = min(
    abs(start - chromosome_centroid), abs(end - chromosome_centroid))) %>% 
  group_by(chromosome) %>% 
  dplyr::mutate(has_centroid = c(dist_to_centroid == min(dist_to_centroid))) %>% 
  ungroup()

cornuta_10kb_pi_plt_rdy %>% dplyr::count(
  chrom_color_group == "window with temperature-responsive genes",
  below_quantile) %>% drop_na() %>% pull(n) %>% 
  matrix(nrow = 2, ncol = 2) %>% 
  fisher.test(alternative = "greater")

cornuta_10kb_pi_plt_rdy %>% dplyr::count(
  chrom_color_group == "window with temperature-responsive genes",
  above_quantile) %>% drop_na() %>% pull(n) %>% 
  matrix(nrow = 2, ncol = 2) %>% 
  fisher.test(alternative = "greater")

cornuta_mean_pi <- cornuta_10kb_pi %>% dplyr::pull(avg_pi) %>% mean(na.rm = TRUE)

cornuta_pi_manhattan <- cornuta_10kb_pi_plt_rdy %>% 
  drop_na() %>% 
  ggplot(aes(x = id, 
             y = avg_pi)) +
  # Plot points with conditional coloring
  geom_point(data = . %>% filter(chrom_color_group != "window with temperature-responsive genes"),
             aes(color = chrom_color_group), 
             size = 0.5, alpha = 0.8, stroke = 0.5) +
  geom_point(data = . %>% filter(chrom_color_group == "window with temperature-responsive genes"), 
             aes(color = chrom_color_group), 
             size = 0.5, alpha = 0.8, stroke = 0.5) +
  # Add horizontal line at 95% quantile
  geom_hline(yintercept = cornuta_mean_pi, color = "grey10", linetype = "dashed", alpha = 0.5) +
  xlab("Chromosome")+
  ylab("Nucl. diversity \u03c0") +
  scale_x_continuous(
    expand = c(0.01, 0.01),
    breaks = cornuta_10kb_pi_plt_rdy$id[cornuta_10kb_pi_plt_rdy$has_centroid],
    labels = unique(cornuta_10kb_pi_plt_rdy$chromosome)) + 
  # Update color scale to include red for high values
  scale_color_manual(values = c(
    "even" = "#989FE6", "odd" = "#717AC7",
    "window with temperature-responsive genes"  = "#C73A3AF8"), 
    breaks = c("window with temperature-responsive genes")) + 
  guides(color = guide_legend(override.aes = list(size = 2))) +
  theme_bw() +  theme(
    legend.margin = margin(t = 0, b = 0),
    legend.key.height = unit(1, "mm"),
    legend.title = element_blank(),
    panel.grid.major.x = element_blank(),
    panel.grid.minor.x = element_line(linewidth = 0.5),
    panel.grid.minor.y = element_blank(),
    panel.grid.major.y = element_blank(),
    legend.position = "bottom",
    panel.border = element_rect(
      color = "black", fill = NA, linewidth = 0.5)) + 
  scale_y_continuous(expand = c(0, 0), limits = c(NA, 0.010), breaks = seq(0, 0.010, 0.002))


# ggsave(filename = "cornuta-pi-manhattan.png", width = 183, height = 65, units = "mm",
#        plot = cornuta_pi_manhattan, dpi = 600)


combined_pi_manhattan <- cowplot::plot_grid(bicornis_pi_manhattan, cornuta_pi_manhattan, nrow = 2,
                   labels = c("A", "B"), rel_heights = c(1, 1.22))

# ggsave(filename = "combined-pi-manhattan.png", width = 183, height = 110, units = "mm",
#        plot = combined_pi_manhattan, dpi = 600)



load(expression_results_cornuta_path)

degs <- res %>% dplyr::filter(padj < 0.05) %>% pull(geneID)

cornuta_gene_pi_res %>%
  dplyr::mutate(bin = floor(start / 10000)) %>% 
  group_by(chromosome, type, bin) %>% 
  summarise(geneID, avg_pi = mean(avg_pi, na.rm = TRUE)) %>% ungroup() %>% 
  dplyr::mutate(tmp = ifelse(geneID %in% degs, "DEG", "background")) %>% 
  ggplot(aes(x = avg_pi, y = tmp)) + geom_violin() +
  stat_summary(fun = median, geom = "point")


cornuta_gene_pi_res %>%
  dplyr::mutate(tmp = ifelse(geneID %in% degs, "DEG", "background")) %>% 
  ggplot(aes(x = avg_pi, y = tmp)) + geom_violin() +
  stat_summary(fun = median, geom = "point")

# wilcox.test(
#   cornuta_max_nsl %>% dplyr::filter(geneID %in% degs) %>% pull(max_nsl),
#   cornuta_max_nsl %>% dplyr::filter(!(geneID %in% degs)) %>% pull(max_nsl))

cornuta_pi_mod_t1 <- cornuta_gene_pi %>% 
  dplyr::left_join(res_t1, by = "geneID") %>% 
  dplyr::filter(pvalue < 0.05) %>% 
  lm(formula = avg_pi ~ log2FoldChange.y) %>% summary()

cornuta_pi_mod_t2 <- cornuta_gene_pi %>% 
  dplyr::left_join(res_t2, by = "geneID") %>% 
  dplyr::filter(pvalue < 0.05) %>% 
  lm(formula = avg_pi ~ log2FoldChange.y) %>% summary()

cornuta_pi_vs_lfc_t1 <- cornuta_gene_pi %>% 
  dplyr::left_join(res_t1, by = "geneID") %>% 
  dplyr::mutate(tmp = factor(sign(log2FoldChange.y))) %>% 
  dplyr::filter(pvalue < 0.05) %>% 
  ggplot(aes(x = log2FoldChange.y, y = avg_pi)) + 
  geom_point(color = "#F5B829", alpha = 0.75) + 
  geom_smooth(method = "lm", color = "grey50", alpha = 0.5) + 
  annotate("text", x = 1.5, y = 0.005, size = 4, label = paste0(
    "P < ", signif(cornuta_pi_mod_t1$coefficients[8], digits = 2), "\n", 
    "R² = ", signif(cornuta_pi_mod_t1$adj.r.squared, digits = 2)), 
  ) +
  theme_bw() + ylab("Nucleotide diversity pi") + 
  scale_x_continuous(expand = c(0.01, 0.01), limits = c(-2, 2)) + 
  scale_y_continuous(expand = c(0.0001, 0.0001), limits = c(0, 0.006)) +
  theme(axis.title.x = element_blank())


cornuta_pi_vs_lfc_t2 <- cornuta_gene_pi %>% 
  dplyr::left_join(res_t2, by = "geneID") %>% 
  dplyr::mutate(tmp = factor(sign(log2FoldChange.y))) %>% 
  dplyr::filter(pvalue < 0.05) %>% 
  ggplot(aes(x = log2FoldChange.y, y = avg_pi)) + 
  geom_point(color = "#F5B829", alpha = 0.75) + 
  geom_smooth(method = "lm", color = "grey50", alpha = 0.5) + 
  annotate("text", x = 1.5, y = 0.005, size = 4, label = paste0(
    "P < ", signif(cornuta_pi_mod_t2$coefficients[8], digits = 2), "\n", 
    "R² = ", signif(cornuta_pi_mod_t2$adj.r.squared, digits = 2)), 
  ) +
  theme_bw() + xlab("log\U2082(Fold change)") + ylab("Nucleotide diversity") + 
  scale_x_continuous(expand = c(0.01, 0.01), limits = c(-2, 2)) + 
  scale_y_continuous(expand = c(0.0001, 0.0001), limits = c(0, 0.006)) +
  theme(axis.title.x = element_blank())


cornuta_combined_pi_vs_lfc <- cowplot::plot_grid(
  cornuta_pi_vs_lfc_t1, cornuta_pi_vs_lfc_t2, 
  nrow = 2, rel_heights = c(1, 1.075),
  labels = c("D", "E"))

# NSL bicornis ------------------------------------------------------------

bicornis_chr_lengths <- read_tsv(bicornis_chr_lengths_path,
                                col_names = c("chromosome", "length")) %>% 
  dplyr::mutate(chromosome = as.factor(as.integer(str_extract(chromosome, "\\d+"))))

bicornis_nsl <- read_tsv(
  bicornis_nsl_path,
  col_names = c(
    "chromosome", "position", "maf", "x1", "x2", "nsl", "nsl_norm", "x3")) %>% 
  dplyr::mutate(chromosome = as.factor(as.integer(str_extract(chromosome, "\\d+"))))

# bicornis_nsl_windows <- read_tsv(
#   "~/PhD/projects/pop-genomics/nsl/bicornis/combined.nsl.norm.10kb.windows.tsv",
#   col_names = c("chromosome", "window_pos_2", "n_snp", "maf", "percentile", "nsl_norm")) %>%
#   dplyr::mutate(window_pos_2 = window_pos_2 - 1) %>% 
#   dplyr::mutate(nsl_norm = ifelse(is.na(nsl_norm), 0, nsl_norm))

bicornis_nsl_manual_windows <- bicornis_nsl %>% 
  dplyr::mutate(
    end = (position %/% 10000) * 10000,
    start = end - 10000) %>%
  group_by(chromosome, end) %>%
  summarise(n_snp = length(nsl_norm),
            prop = sum(abs(nsl_norm) > 2, na.rm = TRUE)/length(nsl_norm),
            nsl_norm2 = mean(nsl_norm, na.rm = TRUE),
            nsl_norm = mean(abs(nsl_norm), na.rm = TRUE)) %>%
  ungroup() %>% dplyr::mutate(start = end - 10000)

bicornis_nsl_manual_windows_imputed <- bicornis_chr_lengths %>%
  mutate(end = map2(chromosome, length, ~seq(10000, .y, by = 10000))) %>%
  unnest(end) %>%
  left_join(bicornis_nsl_manual_windows, by = c("chromosome", "end")) %>%
  replace_na(list(nsl_norm = 0, nsl_norm2 = 0, n_snp = 0, prop = 0)) %>% 
  dplyr::mutate(id = row_number()) %>% 
  dplyr::mutate(start = end - 10000)


bicornis_10kb_nsl_plt_rdy <- bicornis_nsl_manual_windows %>% 
  full_join(bicornis_gene_coords_res, suffix = c("", ".y"), by = join_by(
    chromosome, start <= centroid, end > centroid)) %>% 
  group_by(chromosome, start, end) %>% 
  reframe(nsl_norm, n = sum(padj < 0.05, na.rm = TRUE)) %>% distinct() %>% 
  mutate(chromosome = as.numeric(factor(chromosome))) %>% 
  mutate(chrom_color_group = case_when(
    n > 0 ~ "window with temperature-responsive genes",
    chromosome %% 2 != 0 ~ "even", 
    chromosome %% 2 == 0 ~ "odd")) %>%
  # Add a new variable for high values above 95% quantile
  mutate(above_quantile = nsl_norm > quantile(
    bicornis_nsl_manual_windows$nsl_norm, 0.95, na.rm = TRUE)) %>% 
  mutate(chromosome = factor(chromosome, levels = c(1:16))) %>%
  mutate(id = row_number()) %>% group_by(chromosome) %>% 
  dplyr::filter(!is.na(start), !is.na(end)) %>% 
  dplyr::mutate(chromosome_centroid = (
    min(start) + max(end)) / 2 - 5) %>% 
  ungroup() %>% rowwise() %>% dplyr::mutate(dist_to_centroid = min(
    abs(start - chromosome_centroid))) %>% 
  group_by(chromosome) %>% 
  dplyr::mutate(has_centroid = c(dist_to_centroid == min(dist_to_centroid))) %>% 
  ungroup()

  
bicornis_10kb_nsl_plt_rdy %>% dplyr::count(
  chrom_color_group == "window with temperature-responsive genes",
  above_quantile) %>% drop_na() %>% pull(n) %>% 
  matrix(nrow = 2, ncol = 2) %>% 
  fisher.test(alternative = "greater")

bicornis_mean_nsl <- bicornis_nsl_manual_windows %>% dplyr::pull(nsl_norm) %>% mean(na.rm = TRUE)

bicornis_nsl_manhattan <- bicornis_10kb_nsl_plt_rdy %>% 
  drop_na() %>% 
  ggplot(aes(x = id, 
             y = nsl_norm)) +
  # Plot points with conditional coloring
  geom_point(data = . %>% filter(chrom_color_group != "window with temperature-responsive genes"),
             aes(color = chrom_color_group), 
             size = 0.5, alpha = 0.8, stroke = 0.5) +
  geom_point(data = . %>% filter(chrom_color_group == "window with temperature-responsive genes"), 
             aes(color = chrom_color_group), 
             size = 0.5, alpha = 0.8, stroke = 0.5) +
  # Add horizontal line at 95% quantile
  geom_hline(yintercept = bicornis_mean_nsl, color = "grey10", linetype = "dashed", alpha = 0.5) +
  xlab("Chromosome")+
  ylab("|nSL| score") +
  scale_x_continuous(
    expand = c(0.01, 0.01),
    breaks = bicornis_10kb_nsl_plt_rdy$id[bicornis_10kb_nsl_plt_rdy$has_centroid],
    labels = unique(bicornis_10kb_nsl_plt_rdy$chromosome)) + 
  # Update color scale to include red for high values
  scale_color_manual(values = c(
    "even" = "#EDC24C", "odd" = "#D6A424",
    "window with temperature-responsive genes"  = "#C73A3AF8"), 
    breaks = c("window with temperature-responsive genes")) + 
  guides(color = guide_legend(override.aes = list(size = 2))) +
  theme_bw() +  theme(
    axis.title.x = element_blank(),
    legend.position = "none",
    panel.grid.major.x = element_blank(),
    panel.grid.minor.x = element_line(linewidth = 0.5),
    panel.grid.minor.y = element_blank(),
    panel.grid.major.y = element_blank(),
    panel.border = element_rect(
      color = "black", fill = NA, linewidth = 0.5)) + 
  scale_y_continuous(expand = c(0, 0), limits = c(0, 4), breaks = 1:3)

ggsave(filename = bicornis_nsl_manhattan_out, width = 183, height = 65, units = "mm",
       plot = bicornis_nsl_manhattan, dpi = 600)



# Convert to data.tables if not already
setDT(bicornis_annotation)
setDT(bicornis_nsl)

# Create temporary position column
bicornis_nsl[, pos_temp := position]

# Set keys for efficient joining
setkey(bicornis_nsl, chromosome, position)
bicornis_annotation2 <- bicornis_annotation %>% 
  dplyr::mutate(start = start, end = end)
setkey(bicornis_annotation2, chromosome)

# Perform the calculation
bicornis_max_nsl <- bicornis_nsl[bicornis_annotation2,
                               {
                                 abs_nsl <- abs(nsl_norm)
                                 list(geneID = geneID,
                                      type = type,
                                      max_nsl = max(abs_nsl, na.rm = TRUE),
                                      n = .N,
                                      max_position = pos_temp[which.max(abs_nsl)])
                               },
                               by = .EACHI,
                               on = .(chromosome,
                                      position >= start,
                                      position <= end),
                               nomatch = NA]

# Rename using column positions (2nd and 3rd columns)
setnames(bicornis_max_nsl,
         old = c(2, 3),
         new = c("start", "end"))

# Clean up temporary column
bicornis_nsl[, pos_temp := NULL]

# Convert to tibble
bicornis_max_nsl <- as_tibble(bicornis_max_nsl) %>%
  dplyr::mutate(max_nsl = ifelse(n == 0, 0, max_nsl))

summary1 <- bicornis_max_nsl %>% group_by(chromosome) %>%
  dplyr::count(max_position) %>%
  dplyr::filter(n > 1, !is.na(max_position)) %>%
  arrange(desc(n)) %>% ungroup() %>% dplyr::arrange(desc(n))

load(expression_results_bicornis_path)

degs <- res %>% dplyr::filter(padj < 0.05) %>% pull(geneID)
degs_t1 <- res_t1 %>% dplyr::filter(padj < 0.05) %>% pull(geneID)
degs_t2 <- res_t2 %>% dplyr::filter(padj < 0.05) %>% pull(geneID)
degs_t2 <- res_t3 %>% dplyr::filter(padj < 0.05) %>% pull(geneID)


bicornis_max_nsl %>% dplyr::filter(max_nsl > 0) %>% 
  dplyr::mutate(bin = floor(start / 10000)) %>% 
  group_by(chromosome, type, bin) %>% 
  summarise(geneID, max_nsl = mean(max_nsl, na.rm = TRUE)) %>% ungroup() %>% 
  dplyr::filter(type %in% c("gene", "exon", "intron", "intergenic_region")) %>% 
  dplyr::mutate(tmp = ifelse(geneID %in% degs, "DEG", "background")) %>% 
  ggplot(aes(x = max_nsl, y = tmp)) + geom_violin() + facet_wrap(vars(type)) + 
  stat_summary(fun = median, geom = "point")

bicornis_max_nsl %>% dplyr::filter(max_nsl > 0) %>% 
  dplyr::filter(type %in% c("gene", "exon", "intron", "intergenic_region")) %>% 
  dplyr::mutate(tmp = ifelse(geneID %in% degs, "DEG", "background")) %>% 
  ggplot(aes(x = max_nsl, y = tmp)) + geom_violin() + facet_wrap(vars(type)) + 
  stat_summary(fun = median, geom = "point")

# wilcox.test(
#   bicornis_max_nsl %>% dplyr::filter(geneID %in% degs) %>% pull(max_nsl),
#   bicornis_max_nsl %>% dplyr::filter(!(geneID %in% degs)) %>% pull(max_nsl))

bicornis_nsl_mod_t1 <- bicornis_max_nsl %>% 
  dplyr::filter(type == "gene", max_nsl > 0) %>% 
  dplyr::left_join(res_t1, by = "geneID") %>% 
  dplyr::filter(pvalue < 0.05) %>% 
  lm(formula = max_nsl ~ log2FoldChange.y) %>% summary()

bicornis_nsl_mod_t2 <- bicornis_max_nsl %>% 
  dplyr::filter(type == "gene", max_nsl > 0) %>% 
  dplyr::left_join(res_t2, by = "geneID") %>% 
  dplyr::filter(pvalue < 0.05) %>% 
  lm(formula = max_nsl ~ log2FoldChange.y) %>% summary()

bicornis_nsl_mod_t3 <- bicornis_max_nsl %>% 
  dplyr::filter(type == "gene", max_nsl > 0) %>% 
  dplyr::left_join(res_t3, by = "geneID") %>% 
  dplyr::filter(pvalue < 0.05) %>% 
  lm(formula = max_nsl ~ log2FoldChange.y) %>% summary()

bicornis_nsl_vs_lfc_t1 <- bicornis_max_nsl %>% 
  dplyr::filter(max_nsl > 0, type == "gene") %>% 
  dplyr::left_join(res_t1, by = "geneID") %>% 
  dplyr::mutate(tmp = factor(sign(log2FoldChange.y))) %>% 
  dplyr::filter(pvalue < 0.05) %>% 
  ggplot(aes(x = log2FoldChange.y, y = max_nsl)) + 
  geom_point(color = "#F5B829", alpha = 0.75) + 
  geom_smooth(method = "lm", color = "grey50", alpha = 0.5) + 
  annotate("text", x = 1.5, y = 4.25, size = 4, label = paste0(
    "P < ", signif(bicornis_nsl_mod_t1$coefficients[8], digits = 2), "\n", 
    "R² = ", signif(bicornis_nsl_mod_t1$adj.r.squared, digits = 2)), 
  ) +
  theme_bw() + ylab("|nSL| score") + 
  scale_x_continuous(expand = c(0.01, 0.01), limits = c(-2, 2)) + 
  scale_y_continuous(expand = c(0.01, 0.01), limits = c(0, 5)) +
  theme(axis.title.x = element_blank(),
        plot.title = element_text(face = "bold", hjust = 0.5)) + 
  ggtitle("Week 1")


bicornis_nsl_vs_lfc_t2 <- bicornis_max_nsl %>% 
  dplyr::filter(max_nsl > 0, type == "gene") %>% 
  dplyr::left_join(res_t2, by = "geneID") %>% 
  dplyr::mutate(tmp = factor(sign(log2FoldChange.y))) %>% 
  dplyr::filter(pvalue < 0.05) %>% 
  ggplot(aes(x = log2FoldChange.y, y = max_nsl)) + 
  geom_point(color = "#F5B829", alpha = 0.75) + 
  geom_smooth(method = "lm", color = "grey50", alpha = 0.5) + 
  annotate("text", x = 1.5, y = 4.25, size = 4, label = paste0(
    "P < ", signif(bicornis_nsl_mod_t2$coefficients[8], digits = 2), "\n", 
    "R² = ", signif(bicornis_nsl_mod_t2$adj.r.squared, digits = 2)), 
  ) +
  theme_bw() + xlab("log\U2082(Fold change)") + ylab("|nSL| score") + 
  scale_x_continuous(expand = c(0.01, 0.01), limits = c(-2, 2)) + 
  scale_y_continuous(expand = c(0.01, 0.01), limits = c(0, 5)) +
  theme(axis.title.x = element_blank(),
        plot.title = element_text(face = "bold", hjust = 0.5)) + 
  ggtitle("Week 7")


bicornis_nsl_vs_lfc_t3 <- bicornis_max_nsl %>% 
  dplyr::filter(max_nsl > 0, type == "gene") %>% 
  dplyr::left_join(res_t3, by = "geneID") %>% 
  dplyr::mutate(tmp = factor(sign(log2FoldChange.y))) %>% 
  dplyr::filter(pvalue < 0.05) %>% 
  ggplot(aes(x = log2FoldChange.y, y = max_nsl)) + 
  geom_point(color = "#F5B829", alpha = 0.75) + 
  geom_smooth(method = "lm", color = "grey50", alpha = 0.5) + 
  annotate("text", x = 1.5, y = 4.25, size = 4, label = paste0(
    "P < ", signif(bicornis_nsl_mod_t3$coefficients[8], digits = 2), "\n", 
    "R² = ", signif(bicornis_nsl_mod_t3$adj.r.squared, digits = 2)), 
  ) +
  theme_bw() + xlab("log\U2082(Fold change)") + ylab("|nSL| score") + 
  scale_x_continuous(expand = c(0.01, 0.01), limits = c(-2, 2)) + 
  scale_y_continuous(expand = c(0.01, 0.01), limits = c(0, 5)) + 
  theme(plot.title = element_text(face = "bold", hjust = 0.5)) + 
  ggtitle("Week 14")


bicornis_combined_nsl_vs_lfc <- cowplot::plot_grid(
  bicornis_nsl_vs_lfc_t1, bicornis_nsl_vs_lfc_t2, 
  bicornis_nsl_vs_lfc_t3, nrow = 3, rel_heights = c(1, 1, 1.075),
  labels = c("C", "D", "E"))



# NSL cornuta ------------------------------------------------------------

cornuta_chr_lengths <- read_tsv(cornuta_chr_lengths_path,
                                col_names = c("chromosome", "length")) %>% 
  dplyr::mutate(chromosome = as.factor(as.integer(str_extract(chromosome, "\\d+"))))

cornuta_nsl <- read_tsv(
  cornuta_nsl_path,
  col_names = c(
    "chromosome", "position", "maf", "x1", "x2", "nsl", "nsl_norm", "x3")) %>% 
  dplyr::mutate(chromosome = as.factor(as.integer(str_extract(chromosome, "\\d+"))))

# cornuta_nsl_windows <- read_tsv(
#   "~/PhD/projects/pop-genomics/nsl/cornuta/combined.nsl.norm.10kb.windows.tsv",
#   col_names = c("chromosome", "window_pos_2", "n_snp", "maf", "percentile", "nsl_norm")) %>%
#   dplyr::mutate(window_pos_2 = window_pos_2 - 1) %>% 
#   dplyr::mutate(nsl_norm = ifelse(is.na(nsl_norm), 0, nsl_norm))

cornuta_nsl_manual_windows <- cornuta_nsl %>% 
  dplyr::mutate(
    end = (position %/% 10000) * 10000,
    start = end - 10000) %>%
  group_by(chromosome, end) %>%
  summarise(n_snp = length(nsl_norm),
            prop = sum(abs(nsl_norm) > 2, na.rm = TRUE)/length(nsl_norm),
            nsl_norm2 = mean(nsl_norm, na.rm = TRUE),
            nsl_norm = mean(abs(nsl_norm), na.rm = TRUE)) %>%
  ungroup() %>% dplyr::mutate(start = end - 10000)

cornuta_nsl_manual_windows_imputed <- cornuta_chr_lengths %>%
  mutate(end = map2(chromosome, length, ~seq(10000, .y, by = 10000))) %>%
  unnest(end) %>%
  left_join(cornuta_nsl_manual_windows, by = c("chromosome", "end")) %>%
  replace_na(list(nsl_norm = 0, nsl_norm2 = 0, n_snp = 0, prop = 0)) %>% 
  dplyr::mutate(id = row_number()) %>% 
  dplyr::mutate(start = end - 10000)




cornuta_10kb_nsl_plt_rdy <- cornuta_nsl_manual_windows %>% 
  full_join(cornuta_gene_coords_res, suffix = c("", ".y"), by = join_by(
    chromosome, start <= centroid, end > centroid)) %>% 
  group_by(chromosome, start, end) %>% 
  reframe(nsl_norm, n = sum(padj < 0.05, na.rm = TRUE)) %>% distinct() %>% 
  mutate(chromosome = as.numeric(factor(chromosome))) %>% 
  mutate(chrom_color_group = case_when(
    n > 0 ~ "window with temperature-responsive genes",
    chromosome %% 2 != 0 ~ "even", 
    chromosome %% 2 == 0 ~ "odd")) %>%
  # Add a new variable for high values above 95% quantile
  mutate(above_quantile = nsl_norm > quantile_95) %>%
  mutate(chromosome = factor(chromosome, levels = c(1:16))) %>%
  mutate(id = row_number()) %>% group_by(chromosome) %>% 
  dplyr::filter(!is.na(start), !is.na(end)) %>% 
  dplyr::mutate(chromosome_centroid = (
    min(start) + max(end)) / 2 - 5) %>% 
  ungroup() %>% rowwise() %>% dplyr::mutate(dist_to_centroid = min(
    abs(start - chromosome_centroid))) %>% 
  group_by(chromosome) %>% 
  dplyr::mutate(has_centroid = c(dist_to_centroid == min(dist_to_centroid))) %>% 
  ungroup()

cornuta_mean_nsl <- cornuta_nsl_manual_windows %>% dplyr::pull(nsl_norm) %>% mean(na.rm = TRUE)

cornuta_nsl_manhattan <- cornuta_10kb_nsl_plt_rdy %>% 
  drop_na() %>% 
  ggplot(aes(x = id, 
             y = nsl_norm)) +
  # Plot points with conditional coloring
  geom_point(data = . %>% filter(chrom_color_group != "window with temperature-responsive genes"),
             aes(color = chrom_color_group), 
             size = 0.5, alpha = 0.8, stroke = 0.5) +
  geom_point(data = . %>% filter(chrom_color_group == "window with temperature-responsive genes"), 
             aes(color = chrom_color_group), 
             size = 0.5, alpha = 0.8, stroke = 0.5) +
  # Add horizontal line at 95% quantile
  geom_hline(yintercept = cornuta_mean_nsl, color = "grey10", linetype = "dashed", alpha = 0.5) +
  xlab("Chromosome")+
  ylab("|nSL| score") +
  scale_x_continuous(
    expand = c(0.01, 0.01),
    breaks = cornuta_10kb_nsl_plt_rdy$id[cornuta_10kb_nsl_plt_rdy$has_centroid],
    labels = unique(cornuta_10kb_nsl_plt_rdy$chromosome)) + 
  # Update color scale to include red for high values
  scale_color_manual(values = c(
    "even" = "#989FE6", "odd" = "#717AC7",
    "window with temperature-responsive genes"  = "#C73A3AF8"), 
    breaks = c("window with temperature-responsive genes")) + 
  guides(color = guide_legend(override.aes = list(size = 2))) +
  theme_bw() +  theme(
    legend.margin = margin(t = 0, b = 0),
    legend.key.height = unit(1, "mm"),
    legend.title = element_blank(),
    panel.grid.major.x = element_blank(),
    panel.grid.minor.x = element_line(linewidth = 0.5),
    panel.grid.minor.y = element_blank(),
    panel.grid.major.y = element_blank(),
    legend.position = "bottom",
    panel.border = element_rect(
      color = "black", fill = NA, linewidth = 0.5)) + 
  scale_y_continuous(expand = c(0, 0), limits = c(0, 4), breaks = 1:3)


ggsave(filename = cornuta_nsl_manhattan_out, width = 183, height = 65, units = "mm",
       plot = cornuta_nsl_manhattan, dpi = 600)


combined_nsl_manhattan <- cowplot::plot_grid(bicornis_nsl_manhattan, cornuta_nsl_manhattan, nrow = 2,
                                            labels = c("A", "B"), rel_heights = c(1, 1.22))


ggsave(filename = combined_nsl_manhattan_out, width = 183, height = 110, units = "mm",
       plot = combined_nsl_manhattan, dpi = 600)


# cornuta_gene_coords <- read_tsv(
#   "nucleotide_diversity/cornuta/cds.bed",
#   col_names = c("chromosome", "window_pos_1", "window_pos_2", "transcript_id")) %>%
#   dplyr::arrange(chromosome, window_pos_1, window_pos_2) %>%
#   dplyr::mutate(geneID = str_extract(transcript_id, "gene\\d+")) %>%
#   group_by(chromosome, geneID) %>%
#   dplyr::summarise(window_pos_1 = min(window_pos_1),
#                    window_pos_2 = max(window_pos_2)) %>% ungroup() %>%
#   dplyr::mutate(chromosome = dense_rank(chromosome)) %>%
#   dplyr::filter(chromosome <= 16)

# cornuta_annotation <- read_tsv(
#   "data/cornuta/annotation-with-intron-with-intergenic.gff",
#   comment = "#", col_names = c(
#     "chromosome", "x1", "type", "start", "end", "x2", "strand", "x3", "info")) %>%
#   dplyr::arrange(chromosome, type, start, end) %>%
#   dplyr::mutate(geneID = str_extract(info, "gene\\d+")) %>%
#   dplyr::select(chromosome, geneID, type, start, end, strand) %>%
#   dplyr::filter(type != "CDS") %>%
#   dplyr::mutate(chromosome = as.integer(str_extract(chromosome, "\\d+")))

# cornuta_mean_nsl <- cornuta_gene_coords %>%
#   rowwise() %>%
#   mutate(
#     avg_nsl_norm = max(
#       cornuta_nsl$nsl_norm[
#         cornuta_nsl$chromosome == chromosome &
#           cornuta_nsl$position >= window_pos_1  &
#           cornuta_nsl$position <= window_pos_2],
#       na.rm = TRUE),
#     n = length(cornuta_nsl$nsl_norm[
#       cornuta_nsl$chromosome == chromosome &
#         cornuta_nsl$position >= window_pos_1 &
#         cornuta_nsl$position <= window_pos_2])) %>% ungroup()
#
#
# cornuta_mean_nsl <- cornuta_gene_coords %>%
#   left_join(cornuta_nsl, by = "chromosome", relationship = "many-to-many") %>%
#   filter(position >= window_pos_1 & position <= window_pos_2) %>%
#   group_by(geneID) %>%
#   summarise(
#     avg_nsl_norm = min(nsl_norm, na.rm = TRUE),
#     n = n(),
#     .groups = "drop"
#   )


# Convert to data.tables if not already
setDT(cornuta_annotation)
setDT(cornuta_nsl)

# Create temporary position column
cornuta_nsl[, pos_temp := position]

# Set keys for efficient joining
setkey(cornuta_nsl, chromosome, position)
setkey(cornuta_annotation, chromosome)

# Perform the calculation
cornuta_max_nsl <- cornuta_nsl[cornuta_annotation,
                                {
                                  abs_nsl <- abs(nsl_norm)
                                  list(geneID = geneID,
                                       type = type,
                                       max_nsl = max(abs_nsl, na.rm = TRUE),
                                       n = .N,
                                       max_position = pos_temp[which.max(abs_nsl)])
                                },
                                by = .EACHI,
                                on = .(chromosome,
                                       position >= start,
                                       position <= end),
                                nomatch = NA]

# Rename using column positions (2nd and 3rd columns)
setnames(cornuta_max_nsl,
         old = c(2, 3),
         new = c("start", "end"))

# Clean up temporary column
cornuta_nsl[, pos_temp := NULL]

# Convert to tibble
cornuta_max_nsl <- as_tibble(cornuta_max_nsl) %>%
  dplyr::mutate(max_nsl = ifelse(n == 0, 0, max_nsl))

summary1 <- cornuta_max_nsl %>% group_by(chromosome) %>%
  dplyr::count(max_position) %>%
  dplyr::filter(n > 1, !is.na(max_position)) %>%
  arrange(desc(n)) %>% ungroup() %>% dplyr::arrange(desc(n))

load(expression_results_cornuta_path)

degs <- res %>% dplyr::filter(padj < 0.05) %>% pull(geneID)
degs_t1 <- res_t1 %>% dplyr::filter(padj < 0.05) %>% pull(geneID)
degs_t2 <- res_t2 %>% dplyr::filter(padj < 0.05) %>% pull(geneID)


cornuta_max_nsl %>% 
  dplyr::filter(max_nsl > 0) %>% 
  dplyr::mutate(bin = floor(start / 10000)) %>% 
  group_by(chromosome, type, bin) %>% 
  summarise(geneID, max_nsl = mean(max_nsl, na.rm = TRUE)) %>% ungroup() %>% 
  dplyr::filter(type %in% c("gene", "exon", "intron", "intergenic_region")) %>% 
  dplyr::mutate(tmp = ifelse(geneID %in% degs, "DEG", "background")) %>% 
  ggplot(aes(x = max_nsl, y = tmp)) + geom_violin() +
  facet_wrap(vars(type)) + 
  stat_summary(fun = median, geom = "point")

cornuta_max_nsl %>%   dplyr::filter(max_nsl > 0) %>% 
  dplyr::filter(type %in% c("gene", "exon", "intron", "intergenic_region")) %>% 
  dplyr::mutate(tmp = ifelse(geneID %in% degs, "DEG", "background")) %>% 
  ggplot(aes(x = max_nsl, y = tmp)) + geom_violin() +
  facet_wrap(vars(type)) + 
  stat_summary(fun = median, geom = "point")


# wilcox.test(
#   cornuta_max_nsl %>% dplyr::filter(geneID %in% degs) %>% pull(max_nsl),
#   cornuta_max_nsl %>% dplyr::filter(!(geneID %in% degs)) %>% pull(max_nsl))

cornuta_nsl_mod_t1 <- cornuta_max_nsl %>% 
  dplyr::filter(type == "gene", max_nsl > 0) %>% 
  dplyr::left_join(res_t1, by = "geneID") %>% 
  dplyr::filter(padj < 0.05) %>% 
  lm(formula = max_nsl ~ log2FoldChange.y) %>% summary()

cornuta_nsl_mod_t2 <- cornuta_max_nsl %>% 
  dplyr::filter(type == "gene", max_nsl > 0) %>% 
  dplyr::left_join(res_t2, by = "geneID") %>% 
  dplyr::filter(padj < 0.05) %>% 
  lm(formula = max_nsl ~ log2FoldChange.y) %>% summary()

cornuta_nsl_vs_lfc_t1 <- cornuta_max_nsl %>% 
  dplyr::filter(max_nsl > 0, type == "gene") %>% 
  dplyr::left_join(res_t1, by = "geneID") %>% 
  dplyr::mutate(tmp = factor(sign(log2FoldChange.y))) %>% 
  dplyr::filter(pvalue < 0.05) %>% 
  ggplot(aes(x = log2FoldChange.y, y = max_nsl)) + 
  geom_point(color = "#7D8AED", alpha = 0.75) + 
  geom_smooth(method = "lm", color = "grey50", alpha = 0.5) + 
  annotate("text", x = 1.5, y = 4.25, size = 4, label = paste0(
    "P < ", signif(cornuta_nsl_mod_t1$coefficients[8], digits = 2), "\n", 
    "R² = ", signif(cornuta_nsl_mod_t1$adj.r.squared, digits = 2)), 
           ) + 
  theme_bw() + ylab("|nSL| score") + 
  scale_x_continuous(expand = c(0.01, 0.01), limits = c(-2, 2)) + 
  scale_y_continuous(expand = c(0.01, 0.01), limits = c(0, 5)) + 
  ggtitle("Week 1") + theme(axis.title.x = element_blank(),
        plot.title = element_text(face = "bold", hjust = 0.5))

cornuta_nsl_vs_lfc_t2 <- cornuta_max_nsl %>% 
  dplyr::filter(max_nsl > 0, type == "gene") %>% 
  dplyr::left_join(res_t2, by = "geneID") %>% 
  dplyr::mutate(tmp = factor(sign(log2FoldChange.y))) %>% 
  dplyr::filter(pvalue < 0.05) %>% 
  ggplot(aes(x = log2FoldChange.y, y = max_nsl)) + 
  geom_point(color = "#7D8AED", alpha = 0.75) + 
  geom_smooth(method = "lm", color = "grey50", alpha = 0.5) + 
  annotate("text", x = 1.5, y = 4.25, size = 4, label = paste0(
    "P < ", signif(cornuta_nsl_mod_t2$coefficients[8], digits = 2), "\n", 
    "R² = ", signif(cornuta_nsl_mod_t2$adj.r.squared, digits = 2)), 
  ) +
  theme_bw() + xlab("log\U2082(Fold change)") + ylab("|nSL| score") + 
  scale_x_continuous(expand = c(0.01, 0.01), limits = c(-2, 2)) + 
  scale_y_continuous(expand = c(0.01, 0.01), limits = c(0, 5)) + 
  ggtitle("Week 7") + theme(plot.title = element_text(face = "bold", hjust = 0.5))

cornuta_combined_nsl_vs_lfc <- cowplot::plot_grid(
  cornuta_nsl_vs_lfc_t1, cornuta_nsl_vs_lfc_t2, nrow = 2,
  rel_heights = c(1, 1.075), labels = c("F", "G"))

cornuta_combined_nsl_vs_lfc2 <- cowplot::plot_grid(
  cornuta_combined_nsl_vs_lfc, NULL, nrow = 2, rel_heights = c(2, 1))

combined_nsl_vs_lfc <- cowplot::plot_grid(bicornis_combined_nsl_vs_lfc, 
                   cornuta_combined_nsl_vs_lfc2,
                   ncol = 2)

ggsave(filename = combined_nsl_vs_lfc_out, width = 183, height = 120, units = "mm",
       plot = combined_nsl_vs_lfc, dpi = 600, scale = 1.25)

combined_nsl_all <- cowplot::plot_grid(combined_nsl_manhattan, 
                                       combined_nsl_vs_lfc, nrow = 2)

ggsave(filename = combined_nsl_all_out, width = 183, height = 230, units = "mm",
       plot = combined_nsl_all, dpi = 600, scale = 1.25)


# cornuta_gene_coords <- read_tsv(
#   "nucleotide_diversity/cornuta/cds.bed",  
#   col_names = c("chromosome", "window_pos_1", "window_pos_2", "transcript_id")) %>% 
#   dplyr::arrange(chromosome, window_pos_1, window_pos_2) %>% 
#   dplyr::mutate(geneID = str_extract(transcript_id, "gene\\d+")) %>% 
#   group_by(chromosome, geneID) %>% 
#   dplyr::summarise(window_pos_1 = min(window_pos_1), 
#                    window_pos_2 = max(window_pos_2)) %>% ungroup() %>% 
#   dplyr::mutate(chromosome = dense_rank(chromosome)) %>% 
#   dplyr::filter(chromosome <= 16)
# 
# cornuta_annotation <- read_tsv(
#   "data/cornuta/annotation-with-intron-with-intergenic.gff", 
#   comment = "#", col_names = c(
#     "chromosome", "x1", "type", "start", "end", "x2", "strand", "x3", "info")) %>% 
#   dplyr::arrange(chromosome, type, start, end) %>% 
#   dplyr::mutate(geneID = str_extract(info, "gene\\d+")) %>% 
#   dplyr::select(chromosome, geneID, type, start, end, strand) %>% 
#   dplyr::filter(type != "CDS") %>% 
#   dplyr::mutate(chromosome = as.integer(str_extract(chromosome, "\\d+")))




plot_gsea <- function(gsea, max_y, title = "none", y_title = "none", 
                      show_guide = FALSE, low = "#FFC847", high = "#EBA50E"){
  plot <- gsea %>% dplyr::filter(NES > 0) %>%
    dplyr::mutate(type = factor(type, levels = c("K", "BP", "MF", "CC"))) %>%
    dplyr::group_by(type) %>% 
    dplyr::arrange(desc(NES), .by_group = TRUE) %>%
    dplyr::mutate(
      pathway = ifelse(padj >= 0.05 & padj < 0.1,
                       paste0("<i>", pathway, "</i>"), 
                       pathway),
      pathway = ifelse(padj <= 0.01,
                       paste0("<b>", pathway, "</b>"), 
                       pathway),
      pathway = str_wrap(pathway, 50),
      pathway = str_replace_all(pathway, "\n", "<br>"),
      pathway = factor(pathway, levels = unique(pathway))
    ) %>% 
    dplyr::ungroup() %>% 
    ggplot() + aes(y = NES, x = pathway, group = type, fill = size) + geom_col() + 
    scale_x_discrete() + 
    scale_y_continuous(limits = c(0, max_y)) +
    scale_fill_gradient(low = low, high = high, name = "Set size: ", 
                        limits = c(10,200), breaks = c(50, 100, 150)) +
    theme_bw() + 
    ggtitle(title) + 
    facet_grid(cols = vars(type), scales = "free_x", space = "free_x",
               switch = "x") +
    labs(y = NULL, x = NULL) +
    theme(
      plot.title = element_text(face = "bold", hjust = 0.5),
      plot.margin = margin(r = 0),
      axis.title.y = element_text(size = 8),
      axis.text.x = ggtext::element_markdown(angle = 90, hjust = 1, vjust = 0.5),
      panel.border = element_rect(colour = "grey25", fill = NA, linewidth = 0.5)
    )
  if(!show_guide){
    plot <- plot + guides(fill = "none")
  }
  if(y_title != "none"){
    plot <- plot + ylab(y_title)
  }
  
  
  return(plot)
}






jaccard_similarity <- function(set1, set2) {
  intersection <- length(intersect(set1, set2))
  union <- length(union(set1, set2))
  intersection / union
}

# Filter redundant GO terms based on gene overlap
filter_redundant <- function(df, similarity_threshold = 0.75, size = FALSE) {
  
  n <- nrow(df)
  if(n < 6){
    return(df)
  }
  keep <- rep(TRUE, n)
  
  
  df <- df %>% arrange(desc(abs(NES)))
  
  if(size){
    # Sort by size (keep shortest)
    df <- df %>% arrange(nchar(pathway))
  }
  for(i in 1:(n-1)) {
    if(!keep[i]) next
    
    for(j in (i+1):n) {
      if(!keep[j]) next
      
      # Calculate gene overlap
      genes1 <- df$leadingEdge[[i]]
      genes2 <- df$leadingEdge[[j]]
      similarity <- jaccard_similarity(genes1, genes2)
      
      # If too similar, remove the less significant one
      if(sign(df$NES[[i]]) == sign(df$NES[[j]])){
        if(similarity > similarity_threshold) {
          keep[j] <- FALSE
          print(paste0("Removing ", df$pathway[[j]]))
        }
      }
    }
  }
  
  df[keep, ]
}


# Bicornis GO enrichment --------------------------------------------------

load(expression_results_bicornis_path)
load(go_bicornis_path)

MINGSSIZE = 10
MAXGSSIZE = 200

bicornis_bp_t1 <- bicornis_max_nsl %>% dplyr::filter(max_nsl > 0, type == "gene") %>%
  left_join(res_t1, by = "geneID") %>% dplyr::filter(padj < 0.05) %>%
  dplyr::pull(max_nsl, geneID) %>%
  fgsea::fgsea(pathways = gene2bpGO, minSize = MINGSSIZE, maxSize = MAXGSSIZE,
               gseaParam = 0, scoreType = "pos") %>% as_tibble() %>%
  dplyr::mutate(type = "BP")

bicornis_mf_t1 <- bicornis_max_nsl %>% dplyr::filter(max_nsl > 0, type == "gene") %>%
  left_join(res_t1, by = "geneID") %>% dplyr::filter(padj < 0.05) %>%
  dplyr::pull(max_nsl, geneID) %>%
  fgsea::fgsea(pathways = gene2mfGO, minSize = MINGSSIZE, maxSize = MAXGSSIZE,
               gseaParam = 0, scoreType = "pos") %>% as_tibble() %>%
  dplyr::mutate(type = "MF")

bicornis_cc_t1 <- bicornis_max_nsl %>% dplyr::filter(max_nsl > 0, type == "gene") %>%
  left_join(res_t1, by = "geneID") %>% dplyr::filter(padj < 0.05) %>%
  dplyr::pull(max_nsl, geneID) %>%
  fgsea::fgsea(pathways = gene2ccGO, minSize = MINGSSIZE, maxSize = MAXGSSIZE,
               gseaParam = 0, scoreType = "pos") %>% as_tibble() %>%
  dplyr::mutate(type = "CC")

bicornis_kegg_t1 <- bicornis_max_nsl %>% dplyr::filter(max_nsl > 0, type == "gene") %>%
  left_join(res_t1, by = "geneID") %>% dplyr::filter(padj < 0.05) %>%
  dplyr::pull(max_nsl, geneID) %>%
  fgsea::fgsea(pathways = gene2KEGG, minSize = MINGSSIZE,
               gseaParam = 0, scoreType = "pos") %>% as_tibble() %>%
  dplyr::mutate(type = "K")

bicornis_combined_t1 <- bicornis_kegg_t1 %>% bind_rows(bicornis_bp_t1) %>%
  bind_rows(bicornis_mf_t1) %>% bind_rows(bicornis_cc_t1) %>%
  dplyr::filter(padj < 0.1)



bicornis_bp_t2 <- bicornis_max_nsl %>% dplyr::filter(max_nsl > 0, type == "gene") %>% 
  left_join(res_t2, by = "geneID") %>% dplyr::filter(padj < 0.05) %>%
  dplyr::pull(max_nsl, geneID) %>% 
  fgsea::fgsea(pathways = gene2bpGO, minSize = MINGSSIZE, maxSize = MAXGSSIZE, 
               gseaParam = 0, scoreType = "pos") %>% as_tibble() %>% 
  dplyr::mutate(type = "BP") %>% 
  dplyr::filter(padj < 0.1) %>%  
  filter_redundant() %>% 
  dplyr::arrange(desc(NES)) %>%
  dplyr::filter(NES > 0 & rank(NES) > n() - 5 | NES < 0 & rank(NES) <= 5)

bicornis_mf_t2 <- bicornis_max_nsl %>% dplyr::filter(max_nsl > 0, type == "gene") %>% 
  left_join(res_t2, by = "geneID") %>% dplyr::filter(padj < 0.05) %>%
  dplyr::pull(max_nsl, geneID) %>% 
  fgsea::fgsea(pathways = gene2mfGO, minSize = MINGSSIZE, maxSize = MAXGSSIZE, 
               gseaParam = 0, scoreType = "pos") %>% as_tibble() %>% 
  dplyr::mutate(type = "MF") %>% 
  dplyr::filter(padj < 0.1) %>%  
  filter_redundant() %>% 
  dplyr::arrange(desc(NES)) %>%
  dplyr::filter(NES > 0 & rank(NES) > n() - 3 | NES < 0 & rank(NES) <= 3) %>% 
  dplyr::mutate(pathway = ifelse(
    pathway == "DNA-binding transcription factor activity, RNA polymerase II-specific",
    "RNA polymerase II-specific TF activity", pathway)) %>% 
  dplyr::mutate(pathway = ifelse(
    pathway == "RNA polymerase II cis-regulatory region sequence-specific DNA binding",
    "RNA polymerase II-specific DNA binding", pathway))

bicornis_cc_t2 <- bicornis_max_nsl %>% dplyr::filter(max_nsl > 0, type == "gene") %>% 
  left_join(res_t2, by = "geneID") %>% dplyr::filter(padj < 0.05) %>%
  dplyr::pull(max_nsl, geneID) %>% 
  fgsea::fgsea(pathways = gene2ccGO, minSize = MINGSSIZE, maxSize = MAXGSSIZE, 
               gseaParam = 0, scoreType = "pos") %>% as_tibble() %>% 
  dplyr::mutate(type = "CC") %>% 
  dplyr::filter(padj < 0.1) %>%  
  filter_redundant() %>% 
  dplyr::arrange(desc(NES)) %>%
  dplyr::filter(NES > 0 & rank(NES) > n() - 3 | NES < 0 & rank(NES) <= 3)

bicornis_kegg_t2 <- bicornis_max_nsl %>% dplyr::filter(max_nsl > 0, type == "gene") %>% 
  left_join(res_t2, by = "geneID") %>% dplyr::filter(padj < 0.05) %>%
  dplyr::pull(max_nsl, geneID) %>% 
  fgsea::fgsea(pathways = gene2KEGG, minSize = MINGSSIZE,
               gseaParam = 0, scoreType = "pos") %>% as_tibble() %>% 
  dplyr::mutate(type = "K") %>% 
  dplyr::filter(padj < 0.1) %>%  
  filter_redundant() %>% 
  dplyr::arrange(desc(NES)) %>%
  dplyr::filter(NES > 0 & rank(NES) > n() - 3 | NES < 0 & rank(NES) <= 3)

bicornis_combined_t2 <- bicornis_kegg_t2 %>% bind_rows(bicornis_bp_t2) %>% 
  bind_rows(bicornis_mf_t2) %>% bind_rows(bicornis_cc_t2)


bicornis_bp_t3 <- bicornis_max_nsl %>% dplyr::filter(max_nsl > 0, type == "gene") %>% 
  left_join(res_t3, by = "geneID") %>% dplyr::filter(padj < 0.05) %>%
  dplyr::pull(max_nsl, geneID) %>% 
  fgsea::fgsea(pathways = gene2bpGO, minSize = MINGSSIZE, maxSize = MAXGSSIZE, 
               gseaParam = 0, scoreType = "pos") %>% as_tibble() %>% 
  dplyr::mutate(type = "BP") %>% 
  dplyr::filter(padj < 0.1) %>%  
  filter_redundant() %>% 
  dplyr::arrange(desc(NES)) %>%
  dplyr::filter(NES > 0 & rank(NES) > n() - 5 | NES < 0 & rank(NES) <= 5)

bicornis_mf_t3 <- bicornis_max_nsl %>% dplyr::filter(max_nsl > 0, type == "gene") %>% 
  left_join(res_t3, by = "geneID") %>% dplyr::filter(padj < 0.05) %>%
  dplyr::pull(max_nsl, geneID) %>% 
  fgsea::fgsea(pathways = gene2mfGO, minSize = MINGSSIZE, maxSize = MAXGSSIZE, 
               gseaParam = 0, scoreType = "pos") %>% as_tibble() %>% 
  dplyr::mutate(type = "MF") %>% 
  dplyr::filter(padj < 0.1) %>%  
  dplyr::arrange(desc(NES)) %>%
  filter_redundant() %>% 
  dplyr::filter(NES > 0 & rank(NES) > n() - 3 | NES < 0 & rank(NES) <= 3)

bicornis_cc_t3 <- bicornis_max_nsl %>% dplyr::filter(max_nsl > 0, type == "gene") %>% 
  left_join(res_t3, by = "geneID") %>% dplyr::filter(padj < 0.05) %>%
  dplyr::pull(max_nsl, geneID) %>% 
  fgsea::fgsea(pathways = gene2ccGO, minSize = MINGSSIZE, maxSize = MAXGSSIZE, 
               gseaParam = 0, scoreType = "pos") %>% as_tibble() %>% 
  dplyr::mutate(type = "CC") %>% 
  dplyr::filter(padj < 0.1) %>%  
  dplyr::arrange(desc(NES)) %>%
  filter_redundant() %>% 
  dplyr::filter(NES > 0 & rank(NES) > n() - 3 | NES < 0 & rank(NES) <= 3)

bicornis_kegg_t3 <- bicornis_max_nsl %>% dplyr::filter(max_nsl > 0, type == "gene") %>% 
  left_join(res_t3, by = "geneID") %>% dplyr::filter(padj < 0.05) %>%
  dplyr::pull(max_nsl, geneID) %>% 
  fgsea::fgsea(pathways = gene2KEGG, minSize = MINGSSIZE,
               gseaParam = 0, scoreType = "pos") %>% as_tibble() %>% 
  dplyr::mutate(type = "K") %>% 
  dplyr::filter(padj < 0.1) %>%  
  dplyr::arrange(desc(NES)) %>%
  filter_redundant() %>% 
  dplyr::filter(NES > 0 & rank(NES) > n() - 3 | NES < 0 & rank(NES) <= 3)

bicornis_combined_t3 <- bicornis_kegg_t3 %>% bind_rows(bicornis_bp_t3) %>% 
  bind_rows(bicornis_mf_t3) %>% bind_rows(bicornis_cc_t3)

max_y <- max(bicornis_combined_t2$NES, bicornis_combined_t3$NES)

bicornis_t2_plot <- plot_gsea(bicornis_combined_t2, max_y = max_y, title = "Week 7", 
                              y_title = "Normalized enrichment score")
bicornis_t3_plot <- plot_gsea(bicornis_combined_t3, show_guide = TRUE, max_y = max_y, title = "Week 14")
bicornis_t3_plot2 <- plot_grid(bicornis_t3_plot, NULL, nrow = 2, rel_heights = c(1, 0.02))

bicornis_combined_go_plot <- plot_grid(
  bicornis_t2_plot, bicornis_t3_plot2, ncol = 2, rel_widths = c(1, 1.125),
          labels = c("A", "B"))

# Cornuta GO enrichment ---------------------------------------------------



load(expression_results_cornuta_path)
load(go_cornuta_path)

MINGSSIZE = 10
MAXGSSIZE = 200

cornuta_bp_t1 <- cornuta_max_nsl %>% dplyr::filter(max_nsl > 0, type == "gene") %>%
  left_join(res_t1, by = "geneID") %>% dplyr::filter(padj < 0.05) %>%
  dplyr::pull(max_nsl, geneID) %>%
  fgsea::fgsea(pathways = gene2bpGO, minSize = MINGSSIZE, maxSize = MAXGSSIZE,
               gseaParam = 0, scoreType = "pos") %>% as_tibble() %>%
  dplyr::mutate(type = "BP")

cornuta_mf_t1 <- cornuta_max_nsl %>% dplyr::filter(max_nsl > 0, type == "gene") %>%
  left_join(res_t1, by = "geneID") %>% dplyr::filter(padj < 0.05) %>%
  dplyr::pull(max_nsl, geneID) %>%
  fgsea::fgsea(pathways = gene2mfGO, minSize = MINGSSIZE, maxSize = MAXGSSIZE,
               gseaParam = 0, scoreType = "pos") %>% as_tibble() %>%
  dplyr::mutate(type = "MF")

cornuta_cc_t1 <- cornuta_max_nsl %>% dplyr::filter(max_nsl > 0, type == "gene") %>%
  left_join(res_t1, by = "geneID") %>% dplyr::filter(padj < 0.05) %>%
  dplyr::pull(max_nsl, geneID) %>%
  fgsea::fgsea(pathways = gene2ccGO, minSize = MINGSSIZE, maxSize = MAXGSSIZE,
               gseaParam = 0, scoreType = "pos") %>% as_tibble() %>%
  dplyr::mutate(type = "CC")

cornuta_kegg_t1 <- cornuta_max_nsl %>% dplyr::filter(max_nsl > 0, type == "gene") %>%
  left_join(res_t1, by = "geneID") %>% dplyr::filter(padj < 0.05) %>%
  dplyr::pull(max_nsl, geneID) %>%
  fgsea::fgsea(pathways = gene2KEGG, minSize = MINGSSIZE,
               gseaParam = 0, scoreType = "pos") %>% as_tibble() %>%
  dplyr::mutate(type = "K")

cornuta_combined_t1 <- cornuta_kegg_t1 %>% bind_rows(cornuta_bp_t1) %>%
  bind_rows(cornuta_mf_t1) %>% bind_rows(cornuta_cc_t1) %>%
  dplyr::filter(padj < 0.1)



cornuta_bp_t2 <- cornuta_max_nsl %>% dplyr::filter(max_nsl > 0, type == "gene") %>% 
  left_join(res_t2, by = "geneID") %>% dplyr::filter(padj < 0.05) %>%
  dplyr::pull(max_nsl, geneID) %>% 
  fgsea::fgsea(pathways = gene2bpGO, minSize = MINGSSIZE, maxSize = MAXGSSIZE, 
               gseaParam = 0, scoreType = "pos") %>% as_tibble() %>% 
  dplyr::mutate(type = "BP") %>% 
  dplyr::filter(padj < 0.1) %>%  
  filter_redundant() %>% 
  dplyr::arrange(desc(NES)) %>%
  dplyr::filter(NES > 0 & rank(NES) > n() - 5 | NES < 0 & rank(NES) <= 5)

cornuta_mf_t2 <- cornuta_max_nsl %>% dplyr::filter(max_nsl > 0, type == "gene") %>% 
  left_join(res_t2, by = "geneID") %>% dplyr::filter(padj < 0.05) %>%
  dplyr::pull(max_nsl, geneID) %>% 
  fgsea::fgsea(pathways = gene2mfGO, minSize = MINGSSIZE, maxSize = MAXGSSIZE, 
               gseaParam = 0, scoreType = "pos") %>% as_tibble() %>% 
  dplyr::mutate(type = "MF") %>% 
  dplyr::filter(padj < 0.1) %>%  
  filter_redundant() %>% 
  dplyr::arrange(desc(NES)) %>%
  dplyr::filter(NES > 0 & rank(NES) > n() - 3 | NES < 0 & rank(NES) <= 3) %>% 
  dplyr::mutate(pathway = ifelse(
    pathway == "DNA-binding transcription factor activity, RNA polymerase II-specific",
    "RNA polymerase II-specific TF activity", pathway)) %>% 
  dplyr::mutate(pathway = ifelse(
    pathway == "RNA polymerase II cis-regulatory region sequence-specific DNA binding",
    "RNA polymerase II-specific DNA binding", pathway))

cornuta_cc_t2 <- cornuta_max_nsl %>% dplyr::filter(max_nsl > 0, type == "gene") %>% 
  left_join(res_t2, by = "geneID") %>% dplyr::filter(padj < 0.05) %>%
  dplyr::pull(max_nsl, geneID) %>% 
  fgsea::fgsea(pathways = gene2ccGO, minSize = MINGSSIZE, maxSize = MAXGSSIZE, 
               gseaParam = 0, scoreType = "pos") %>% as_tibble() %>% 
  dplyr::mutate(type = "CC") %>% 
  dplyr::filter(padj < 0.1) %>%  
  filter_redundant() %>% 
  dplyr::arrange(desc(NES)) %>%
  dplyr::filter(NES > 0 & rank(NES) > n() - 3 | NES < 0 & rank(NES) <= 3)

cornuta_kegg_t2 <- cornuta_max_nsl %>% dplyr::filter(max_nsl > 0, type == "gene") %>% 
  left_join(res_t2, by = "geneID") %>% dplyr::filter(padj < 0.05) %>%
  dplyr::pull(max_nsl, geneID) %>% 
  fgsea::fgsea(pathways = gene2KEGG, minSize = MINGSSIZE,
               gseaParam = 0, scoreType = "pos") %>% as_tibble() %>% 
  dplyr::mutate(type = "K") %>% 
  dplyr::filter(padj < 0.1) %>%  
  filter_redundant() %>% 
  dplyr::arrange(desc(NES)) %>%
  dplyr::filter(NES > 0 & rank(NES) > n() - 3 | NES < 0 & rank(NES) <= 3)

cornuta_combined_t2 <- cornuta_kegg_t2 %>% bind_rows(cornuta_bp_t2) %>% 
  bind_rows(cornuta_mf_t2) %>% bind_rows(cornuta_cc_t2)

max_y <- max(cornuta_combined_t1$NES, cornuta_combined_t2$NES)

cornuta_t1_plot <- plot_gsea(cornuta_combined_t1, max_y = max_y, 
                             title = "Week 1", low = "#A8A8FA", high = "#6A5CFF",
                             y_title = "Normalized enrichment score")
cornuta_t2_plot <- plot_gsea(cornuta_combined_t2, show_guide = TRUE, max_y = max_y, 
                             title = "Week 7", low = "#A8A8FA", high = "#6A5CFF")
cornuta_t2_plot2 <- plot_grid(cornuta_t2_plot, NULL, nrow = 2, rel_heights = c(1, 0.005))

cornuta_combined_go_plot <- plot_grid(
  cornuta_t1_plot, cornuta_t2_plot2, ncol = 2, rel_widths = c(1, 1.6),
          labels = c("C", "D"))

ggsave(plot_grid(bicornis_combined_go_plot, cornuta_combined_go_plot,
                              nrow = 2), filename = figure5_out, 
       width = 183, height = 183, unit = "mm", dpi = 600, scale = 1.25)
