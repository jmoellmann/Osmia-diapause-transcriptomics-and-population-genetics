library(tidyverse)
library(tximport)
library(DESeq2)
library(topGO)
library(fgsea)
library(ComplexHeatmap)
library(clusterProfiler)
library(limma)
library(RColorBrewer)
library(patchwork)
library(cowplot)

#library(factoextra)

# Differential expression (DESeq2, temperature x time-point design) and GSEA for
# O. cornuta — same structure as expression_analysis_bicornis.R, without a
# sample-exclusion filter (none needed for this species).
#
# kallisto_dir points at this pipeline's own kallisto_quant output.
#
# Reads osmia_bicornis_gene2go_refseq.tsv — a bicornis-specific NCBI gene2go dump,
# mapped onto cornuta genes via orthology since no cornuta-specific dump exists.

kallisto_dir <- snakemake@params[["kallisto_dir"]]
metadata_path <- snakemake@params[["metadata_path"]]
tx2gene_path <- snakemake@input[["tx2gene"]]
orthogroups_path <- snakemake@input[["orthogroups"]]
dmel_flybase2ncbi_path <- snakemake@input[["dmel_flybase2ncbi"]]
dmel_gene_summaries_path <- snakemake@input[["dmel_gene_summaries"]]
dme2keggpw_path <- snakemake@input[["dme2keggpw"]]
ncbi2dme_path <- snakemake@input[["ncbi2dme"]]
dm_ncbi_to_kegg_path <- snakemake@input[["dm_ncbi_to_kegg"]]
go_basic_path <- snakemake@input[["go_basic"]]
gene2go_path <- snakemake@input[["gene2go"]]
dmel_gene2go_path <- snakemake@input[["dmel_gene2go"]]
annotations_path <- snakemake@input[["annotations"]]
go_rdata_out <- snakemake@output[["go_rdata"]]
combined_pca_out <- snakemake@output[["combined_pca_png"]]
combined_gsea_out <- snakemake@output[["combined_gsea_pdf"]]
heatmap_combined_out <- snakemake@output[["heatmap_combined_png"]]
expression_results_out <- snakemake@output[["expression_results_rdata"]]
res_t1_out <- snakemake@output[["res_t1_tsv"]]
res_t2_out <- snakemake@output[["res_t2_tsv"]]
res_sexes_out <- snakemake@output[["res_sexes_tsv"]]
res_sexes_interaction_out <- snakemake@output[["res_sexes_interaction_tsv"]]
degs_xlsx_out <- snakemake@output[["degs_xlsx"]]
dir.create(dirname(res_t1_out), recursive = TRUE, showWarnings = FALSE)

term_categories <- c(
  "cilium movement involved in cell motility" = "olfaction",
  "dolichol-linked oligosaccharide biosynthetic process" = "lipid metabolism",
  "phospholipid biosynthetic process" = "lipid metabolism",
  "flagellated sperm motility" = "reproduction",
  "salivary gland boundary specification" = "development",
  "mesodermal cell migration" = "development",
  "dynein light intermediate chain binding" = "olfaction",
  "ATPase-coupled transmembrane transporter activity" = "metabolism",
  "translation elongation factor activity" = "gene regulation",
  "proton-transporting ATPase activity, rotational mechanism" = "metabolism",
  "Cytoskeleton in muscle cells" = "muscle function",
  "Toll-like receptor signaling pathway" = "immunity",
  "small molecule metabolic process" = "metabolism",
  "proton motive force-driven ATP synthesis" = "metabolism",
  "ribosomal small subunit assembly" = "gene regulation",
  "mitochondrial electron transport, cytochrome c to oxygen" = "metabolism",
  "habituation" = "neuronal function",
  "maintenance of presynaptic active zone structure" = "neuronal function",
  "inorganic cation transmembrane transport" = "other",
  "germarium-derived female germ-line cyst formation" = "reproduction",
  "synaptic transmission, cholinergic" = "neuronal function",
  "estradiol 17-beta-dehydrogenase [NAD(P)+] activity" = "reproduction",
  "proton-transporting ATP synthase activity, rotational mechanism" = "metabolism",
  "fatty acid binding" = "lipid metabolism",
  "cell adhesion molecule binding" = "signaling",
  "non-membrane spanning protein tyrosine kinase activity" = "signaling",
  "acetylcholine-gated monoatomic cation-selective channel activity" = "neuronal function",
  "Fatty acid degradation" = "lipid metabolism",
  "Dorso-ventral axis formation" = "development",
  "follicle cell of egg chamber migration" = "other",
  "chromatin organization" = "gene regulation",                                 
  "neurotransmitter:sodium symporter activity" = "neuronal function",
  "ligand-gated monoatomic ion channel activity" = "signaling",           
  "transmembrane receptor protein tyrosine kinase activity" = "signaling",
  "nuclear receptor activity" = "signaling",     
  "muscle contraction" = "muscle function",
  "cilium movement" = "olfaction",
  "detection of chemical stimulus involved in sensory perception of smell" = "olfaction",
  "calcium ion transmembrane transport" = "signaling",
  "sensory perception of smell" = "olfaction",
  "salivary gland development" = "development",
  "neuroblast development" = "development",
  "carbohydrate homeostasis" = "metabolism",
  "entrainment of circadian clock by photoperiod" = "other",
  "glial cell development" = "development",
  "ionotropic olfactory receptor activity" = "olfaction",
  "odorant binding" = "olfaction",
  "olfactory receptor activity" = "olfaction",
  "structural constituent of chromatin" = "gene regulation",
  "ubiquitin-like ligase-substrate adaptor activity" = "gene regulation",
  "hormone activity" = "signaling",
  "Butanoate metabolism" = "metabolism",
  "Valine, leucine and isoleucine degradation" = "metabolism",
  "Glycosylphosphatidylinositol (GPI)-anchor biosynthesis" = "lipid metabolism",
  "Pyrimidine metabolism" = "metabolism",
  "TGF-beta signaling pathway" = "development",
  "Hedgehog signaling pathway" = "development",
  "Notch signaling pathway" = "development",
  "Circadian rhythm" = "other",
  "myofibril assembly" = "muscle function",
  "fatty acid elongation" = "lipid metabolism",
  "very long-chain fatty acid biosynthetic process" = "lipid metabolism",
  "sphingolipid biosynthetic process" = "lipid metabolism",
  "compound eye development" = "development",
  "mRNA transport" = "gene regulation",
  "cytoplasmic translation" = "gene regulation",
  "regulation of membrane potential" = "other",
  "intracellular receptor signaling pathway" = "signaling",
  "L-amino acid transmembrane transporter activity" = "metabolism",
  "transaminase activity" = "metabolism",
  "amino acid transmembrane transporter activity" = "metabolism",
  "microtubule binding" = "other",
  "transcription coregulator activity" = "gene regulation",
  "Metabolism of xenobiotics by cytochrome P450" = "immunity",
  "Biosynthesis of unsaturated fatty acids" = "lipid metabolism",
  "Tyrosine metabolism" = "metabolism",
  "Folate biosynthesis" = "metabolism",
  "Hippo signaling pathway" = "development",
  "mRNA surveillance pathway" = "gene regulation",
  "Nucleocytoplasmic transport" = "gene regulation",
  "Viral life cycle" = "immunity"
)

term_categories_df <- tibble(pathway = names(term_categories), 
                             category = term_categories)



set.seed(0)

# Function to plot individual genes for manual diagnostics ----------------------------

plot_gene <- function(dds, gene){
  plotCounts(dds,
             gene = gene,
             intgroup = "temperature",
             returnData = TRUE) %>%
    {
      medians <- group_by(., temperature) %>%
        summarise(median_count = mean(count), .groups = "drop")
      
      ggplot(., aes(x = temperature, y = count)) +
        geom_violin(aes(fill = temperature), alpha = 0.3, trim = FALSE) +
        geom_point(position = position_jitter(height = 0, width = 0.05),
                   size = 3, alpha = 0.6) +
        geom_point(data = medians, aes(y = median_count),
                   size = 6, shape = 18, color = "red") +
        geom_line(data = medians, aes(y = median_count, group = 1),
                  color = "red", size = 1.5) +
        scale_fill_manual(values = c("Normal" = "gray70", "Hot" = "orange",
                                     "VeryHot" = "red3")) +
        theme_minimal() +
        theme(legend.position = "none") + 
        ggtitle(gene)
    }
}


plot_gene_time <- function(dds, gene){
  plotCounts(dds,
             gene = gene,
             intgroup = "time_point",
             returnData = TRUE) %>%
    {
      medians <- group_by(., time_point) %>%
        summarise(median_count = mean(count), .groups = "drop")
      
      ggplot(., aes(x = time_point, y = count)) +
        geom_violin(aes(fill = time_point), alpha = 0.3, trim = FALSE) +
        geom_point(position = position_jitter(height = 0, width = 0.05),
                   size = 3, alpha = 0.6) +
        geom_point(data = medians, aes(y = median_count),
                   size = 6, shape = 18, color = "red") +
        geom_line(data = medians, aes(y = median_count, group = 1),
                  color = "red", size = 1.5) +
        scale_fill_manual(values = c("Normal" = "gray70", "Hot" = "orange",
                                     "VeryHot" = "red3")) +
        theme_minimal() +
        theme(legend.position = "none") + 
        ggtitle(gene)
    }
}

load_files_and_make_dds <- function(metadata, design = ~ temperature){
  # Load files
  samples <- file.path(kallisto_dir, rownames(metadata))
  files <- paste0(samples, "/abundance.h5")
  txi <- tximport(files, type = "kallisto", tx2gene = tx2gene)
  
  #Generate dds and filter out low read count genes
  dds <- DESeqDataSetFromTximport(
    txi, colData = metadata, design = design)
  smallestGroupSize <- metadata %>% group_by(temperature) %>% 
    summarise(n = n()) %>% pull(n) %>% min()
  keep <- rowSums(counts(dds) >= 10) >= smallestGroupSize
  dds <- dds[keep,]
  
  return(dds)
}


plot_heatmap <- function(vsd, top_genes, metadata, column = "temperature",
                         colors = c("baseline" = "#7FB0F0", "intermediate" = "#7D8AED", 
                                    "highest" = "#B380F2"), anno = NULL, 
                         panel_height = NULL, row_gap = NULL,
                         height_per_gene = 0.005){
  # Count genes
  n_genes <- length(top_genes)
  
  # Calculate height based on number of genes
  heatmap_height <- unit(n_genes * height_per_gene, "cm")
  
  # Subset to genes of interest
  expr_data <- assay(vsd)[top_genes, ]
  
  # Z-score normalize (by row/gene)
  expr_scaled <- t(scale(t(expr_data)))
  
  # Column annotation (for samples)
  col_ann <- HeatmapAnnotation(
    border = TRUE,
    column = metadata$temperature,
    col = list(column = colors),
    annotation_name_gp = gpar(fontsize = 9),
    show_legend = FALSE,
    show_annotation_name = FALSE
  )
  
  # Create heatmap
  ht <- Heatmap(
    use_raster = FALSE,
    row_km = 2,
    row_km_repeats = 10,
    expr_scaled,
    border = TRUE,
    col = viridis::inferno(100),
    show_row_dend = FALSE,
    row_title = NULL,
    row_gap = if(!is.null(row_gap)) unit(row_gap, "cm") else unit(1, "mm"),
    #heatmap_height = heatmap_height,
    
    # Clustering
    cluster_rows = TRUE,
    cluster_columns = FALSE,
    
    # Annotations
    top_annotation = col_ann,
    
    # Row options (genes)
    show_row_names = FALSE,
    row_names_gp = gpar(fontsize = 8),
    
    # Column options (samples)
    column_split = factor(metadata$temperature,
                          levels = names(colors)),
    column_title_gp = gpar(fontsize = 10, fontface = "bold"),
    show_column_names = FALSE,
    show_heatmap_legend = FALSE
    #column_title = paste0("**", anno, " (n = ", n_genes, ")**")
  )
  # 
  # # Legend
  # heatmap_legend_param = list(
  #   title = "Expression\n(Z-score)",
  #   direction = "vertical",
  #   title_gp = gpar(fontsize = 10, fontface = "bold"),
  #   labels_gp = gpar(fontsize = 9)
  # )
  
  
  
  ht2 <- ggplotify::as.ggplot(
    grid.grabExpr(draw(ht, padding = unit(c(2, 2, 2, 2), "mm")))
  ) +
    ggtitle(paste0("**", anno, "**  (n = ", n_genes, ")")) +
    theme(
      plot.title = ggtext::element_markdown(hjust = 0.5, size = 14, margin = margin(b = 2)),
      plot.margin = margin(5, 5, 5, 5)
    )
  
  return(list(plot = ht2, n_genes = n_genes))  
}


get_results <- function(dds, 
                        contrast1 = c("temperature", "plus6", "normal"),
                        contrast2 = c("temperature", "plus3", "normal")){
  res <- results(dds, contrast = contrast1) %>% 
    as_tibble(rownames = "geneID")
  res2 <- results(dds, contrast = contrast2) %>% 
    as_tibble(rownames = "geneID") %>% 
    dplyr::select(-c(stat, pvalue, padj)) %>% 
    left_join(res, by = "geneID") %>% 
    dplyr::mutate(type = case_when(
      log2FoldChange.x > 0 & log2FoldChange.y > 0 ~ "up",
      log2FoldChange.x < 0 & log2FoldChange.y < 0 ~ "down",
      TRUE ~ "complex"))
  return(res2) 
}

fix_names <- function(metadata){
  metadata_fixed <- metadata %>% 
    dplyr::mutate(
      time_point = case_when(
        time_point == "pre_emerg_t1" ~ "week 1",
        time_point == "pre_emerg_t2" ~ "week 7"
      ),
      temperature = case_when(
        temperature == "normal" ~ "baseline",
        temperature == "plus3" ~ "intermediate",
        temperature == "plus6" ~ "highest"
      )
    )
  return(metadata_fixed)
}


create_gene2GO <- function(df) {
  gene2GO <- df %>% drop_na() %>% 
    group_by(!!sym(colnames(df)[1])) %>%
    summarise(GO_terms = list(unique(!!sym(colnames(df)[2])))) %>%
    deframe()
  gene2GO <- lapply(gene2GO, function(x) as.character(unlist(x)))
  return(gene2GO)
}

plot_time_pca <- function(vsd, seed = 1){
  set.seed(seed)
  design <- model.matrix(~time_point, data = colData(vsd))
  mat <- limma::removeBatchEffect(assay(vsd), batch = colData(vsd)$temperature, design = design)
  pca_corrected <- prcomp(t(mat))
  
  # Calculate variance explained
  timepoint_var_explained <- round(100 * pca_corrected$sdev^2 / sum(pca_corrected$sdev^2), 1)
  
  # Plot
  timepoint_pca_data <- data.frame(
    PC1 = pca_corrected$x[,1],
    PC2 = pca_corrected$x[,2],
    temperature = colData(vsd)$temperature,
    timepoint = colData(vsd)$time_point
  ) %>% as_tibble(rownames = "sampleID") %>%
    dplyr::mutate(timepoint = factor(case_when(
      timepoint == "pre_emerg_t1" ~ "week 1",
      timepoint == "pre_emerg_t2" ~ "week 7"
    ), levels = c("week 1", "week 7")))

 
  
  palette1 <- c("#7FB0F0", "#7D8AED")
  time_pca <- ggplot(timepoint_pca_data, aes(x = PC1, y = PC2, color = timepoint, shape = timepoint)) +
    geom_point(size = 1.5, stroke = 0.8, fill = "#7D8AED") + 
    #geom_text_repel(aes(label = sampleID)) +
    scale_color_manual(
      labels = levels(timepoint_pca_data$timepoint),
      values = c("grey50", "black"),
      name = ""
    ) +
    scale_shape_manual(
      label = levels(timepoint_pca_data$timepoint),
      values = c(23, 25),
      name = ""
    ) +
    scale_x_continuous(
      name = paste0("PC1 (", timepoint_var_explained[1], "% variance)"),
      breaks = seq(-20, 40, by = 20), limits = c(
        floor(min(timepoint_pca_data$PC1)) - 2, ceiling(max(timepoint_pca_data$PC1)) + 2)
    ) +
    scale_y_continuous(
      name = paste0("PC2 (", timepoint_var_explained[2], "% variance)"),
      breaks = seq(-40, 20, by = 20), limits = c(
        floor(min(timepoint_pca_data$PC2)) - 2, ceiling(max(timepoint_pca_data$PC2)) + 2)
    ) +
    theme_bw() + 
    theme(
      legend.margin = margin(t = 0, b = 0),
      legend.position = "bottom",
      legend.justification = "center",
      panel.border = element_rect(
        color = "grey25", fill = NA, linewidth = 0.5)    
    )

    return(time_pca)
}


plot_temp_pca <- function(vsd, seed = 1){
  set.seed(seed)
  
  design <- model.matrix(~temperature, data = colData(vsd))
  mat <- limma::removeBatchEffect(assay(vsd), batch = colData(vsd)$time_point, design = design)
  pca_corrected <- prcomp(t(mat))
  
  
  
  # Calculate variance explained
  temperature_var_explained <- round(100 * pca_corrected$sdev^2 / sum(pca_corrected$sdev^2), 1)
  
  # Plot
  temperature_pca_data <- data.frame(
    PC1 = pca_corrected$x[,1],
    PC2 = pca_corrected$x[,2],
    temperature = colData(vsd)$temperature,
    timepoint = colData(vsd)$time_point
  ) %>% as_tibble(rownames = "sampleID") %>%
    dplyr::mutate(temperature = factor(case_when(
      temperature == "normal" ~ "baseline",
      temperature == "plus3" ~ "intermediate",
      temperature == "plus6" ~ "highest"
    ), levels = c("baseline", "intermediate", "highest")))
  

  #palette2 <- c("#E2F513", "#FF9914", "#E00000")
  #palette2 <- c("#7FB0F0", "#7D8AED", "#B380F2")
  #customPalette1 <- c("#EDE100", "#E69F00", "#D55E00") # palette from organismal paper
  temp_pca <- ggplot(temperature_pca_data, aes(x = PC1, y = PC2, color = temperature, shape = temperature)) +
    geom_jitter(width = 2, height = 2, size = 1.5, stroke = 0.8, fill = "#7D8AED") + #color = "black") +
    #geom_text_repel(aes(label = sampleID)) +
    scale_color_manual(
      labels = levels(temperature_pca_data$temperature),
      values = c("grey66", "grey42", "black"),
      name = ""
    ) +
    scale_shape_manual(
      label = levels(temperature_pca_data$temperature),
      values = c(23, 25, 22),
      name = ""
    ) +
    scale_x_continuous(
      name = paste0("PC1 (", temperature_var_explained[1], "% variance)"),
      breaks = seq(-20, 40, by = 20), limits = c(
        floor(min(temperature_pca_data$PC1)) - 2, ceiling(max(temperature_pca_data$PC1)) + 2)
    ) +
    scale_y_continuous(
      name = paste0("PC2 (", temperature_var_explained[2], "% variance)"),
      breaks = seq(-40, 20, by = 20), limits = c(
        floor(min(temperature_pca_data$PC2)) - 2, ceiling(max(temperature_pca_data$PC2)) + 2)
    ) +
    theme_bw() + 
    theme(
      legend.margin = margin(t = 0, b = 0),
      legend.position = "bottom",
      legend.justification = "center",
      panel.border = element_rect(
        color = "grey25", fill = NA, linewidth = 0.5)
    )
  
  return(temp_pca)
}


plot_gseas <- function(gsea0, gsea1){
  
  # Calculate max rows for each week
  calc_max_rows <- function(gsea_data) {
    down_rows <- gsea_data %>% 
      dplyr::filter(NES < 0) %>% 
      nrow()
    up_rows <- gsea_data %>% 
      dplyr::filter(NES > 0) %>% 
      nrow()
    return(max(down_rows, up_rows, na.rm = TRUE))
  }
  
  xlim <- max(abs(c(gsea0$NES, gsea1$NES))) * 1.05
  print(xlim)
  max_rows_t1 <- calc_max_rows(gsea0)
  max_rows_t2 <- calc_max_rows(gsea1)

  # Function to create downregulated plot (NES < 0)
  create_down_plot <- function(gsea_data, show_x_axis = FALSE, show_guide = FALSE) {
    p <- gsea_data %>% 
      dplyr::filter(NES < 0) %>%
      
      dplyr::mutate(type = factor(type, levels = c("KEGG", "BP", "MF", "CC"))) %>%
      dplyr::group_by(type) %>% 
      dplyr::arrange(NES, .by_group = TRUE) %>%
      dplyr::mutate(
        pathway = ifelse(padj >= 0.05 & padj < 0.1,
                         paste0("<i>", pathway, "</i>"), 
                         pathway),
        pathway = ifelse(padj <= 0.01,
                         paste0("<b>", pathway, "</b>"), 
                         pathway),
        pathway = str_wrap(pathway, 60),
        pathway = str_replace_all(pathway, "\n", "<br>"),
        pathway = factor(pathway, levels = unique(pathway))
      ) %>% 
      dplyr::ungroup() %>% 
      ggplot() + 
      aes(x = NES, y = pathway, fill = size) + 
      geom_col() +
      scale_y_discrete(limits = rev) + 
      scale_x_continuous(limits = c(xlim * (-1), 0)) +
      scale_fill_gradient(low = "#A8A8FA", high = "#6A5CFF", name = "Set size: ") +
      cowplot::theme_half_open() + 
      cowplot::background_grid() + 
      facet_grid(rows = vars(type), scales = "free_y", space = "free_y") +
      labs(y = NULL, x = NULL) +
      theme(
        plot.margin = margin(r = 0),
        axis.text.y = ggtext::element_markdown(),
        panel.border = element_rect(colour = "black", fill = NA, linewidth = 1)
      )
    
    if (show_guide){
      p <- p + guides(fill = guide_colorbar(
        barwidth = 10,      # width of the bar
        barheight = 0.5,    # height of the bar
        title.position = "left",
        title.hjust = 0.5
      ))  
    } else {
      p <- p + guides(fill = "none")
    }
    
    if (!show_x_axis) {
      p <- p + theme(axis.text.x = element_blank(),
                     axis.ticks.x = element_blank())
    }
    
    return(p)
  }
  
  # Function to create upregulated plot (NES > 0)
  create_up_plot <- function(gsea_data, show_x_axis = FALSE) {
    p <- gsea_data %>% 
      dplyr::filter(NES > 0) %>%
      dplyr::mutate(type = factor(type, levels = c("KEGG", "BP", "MF", "CC"))) %>%
      dplyr::group_by(type) %>% 
      dplyr::arrange(desc(NES), .by_group = TRUE) %>%
      dplyr::mutate(
        pathway = ifelse(padj >= 0.05 & padj < 0.1,
                         paste0("<i>", pathway, "</i>"), 
                         pathway),
        pathway = ifelse(padj <= 0.01,
                         paste0("<b>", pathway, "</b>"), 
                         pathway),
        pathway = str_wrap(pathway, 60),
        pathway = str_replace_all(pathway, "\n", "<br>"),
        pathway = factor(pathway, levels = unique(pathway))
      ) %>% 
      dplyr::ungroup() %>% 
      ggplot() + 
      aes(x = NES, y = pathway, fill = size) + 
      geom_col() +
      scale_y_discrete(limits = rev, position = "right") + 
      scale_x_continuous(limits = c(0, xlim)) +
      scale_fill_gradient(low = "#A8A8FA", high = "#6A5CFF", name = "Set size: ") +
      cowplot::theme_half_open() + 
      cowplot::background_grid() + 
      facet_grid(rows = vars(type), scales = "free_y", space = "free_y", switch = "y") +
      labs(y = NULL, x = NULL) +
      theme(
        plot.margin = margin(l = 0),
        axis.text.y.right = ggtext::element_markdown(),
        strip.text.y.left = element_blank(),
        strip.background = element_blank(),
        panel.border = element_rect(colour = "black", fill = NA, linewidth = 1),
        panel.spacing.x = unit(0, "lines"),
      ) + 
      guides(fill = "none")
    
    if (!show_x_axis) {
      p <- p + theme(axis.text.x = element_blank(),
                     axis.ticks.x = element_blank())
    }
    
    return(p)
  }
  
  # Create plots for each timepoint (all without legends)
  p_t1_down <- create_down_plot(gsea0)
  p_t1_up <- create_up_plot(gsea0)
  
  p_t2_down <- create_down_plot(gsea1, show_x_axis = TRUE, show_guide = TRUE)
  p_t2_up <- create_up_plot(gsea1, show_x_axis = TRUE)

  
  # Combine each week's plots with title (using plot_spacer for proper centering)
  week1 <- (p_t1_down | p_t1_up) + 
    plot_layout(axes = "collect", widths = c(1, 1)) +
    plot_annotation(theme = theme(
      #plot.background = element_rect(colour = "black", fill = NA, linewidth = 1),
      plot.margin = margin(10, 5, 5, 5)
    ))
  
  week7 <- (p_t2_down | p_t2_up) + 
    plot_layout(axes = "collect", widths = c(1, 1)) +
    plot_annotation(theme = theme(
      #plot.background = element_rect(colour = "black", fill = NA, linewidth = 1),
      plot.margin = margin(10, 5, 5, 5)
    )) & labs(x = "Normalized Enrichment Score") & 
    theme(legend.position = "bottom",
          panel.border = element_rect(colour = "black", fill=NA, linewidth=1))
  
  
  
  # Combine all plots vertically with proportional heights
  combined <- plot_grid(
    week1, 
    week7, 
    ncol = 1,
    rel_heights = c(max_rows_t1, max_rows_t2 + 5.5),
    labels = NULL,
    align = "v",
    axis = "lr"
  )
  
  return(combined)
}

jaccard_similarity <- function(set1, set2) {
  intersection <- length(intersect(set1, set2))
  union <- length(union(set1, set2))
  intersection / union
}

# Filter redundant GO terms based on gene overlap
filter_redundant <- function(df, similarity_threshold = 0.75, size = FALSE) {
  
  n <- nrow(df)
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


# Load tables -------------------------------------------------------------

metadata <- read.delim(metadata_path, sep = "\t")

metadata2 <- metadata %>% dplyr::filter(species == "O. cornuta") %>% 
  dplyr::mutate(population = as.factor(population), 
                temperature = as.factor(temperature),
                time_point = as.factor(time_point),
                inferred_sex = as.factor(inferred_sex)) %>% 
  dplyr::filter(!(sampleID %in% c(""))) %>% 
  dplyr::arrange(time_point, temperature)

rownames(metadata2) <- metadata2$sampleID

tx2gene <- read_tsv(tx2gene_path, col_names = c("transcript", "gene")) %>% 
  mutate(transcript = str_remove(transcript, "rna-")) %>% 
  mutate(gene = str_remove(gene, "gene-"))

# Prepare GO and KEGG annotations

orthogroups <- read_tsv(orthogroups_path)

dmel_fb2ncbi <- read_tsv(
  dmel_flybase2ncbi_path, skip = 1, col_names = c("flybase", "ncbi")) %>% 
  distinct()

dmel2ocor <- orthogroups %>% 
  dplyr::select("Dmelanogaster", "Ocornuta") %>% 
  dplyr::filter(!is.na(Dmelanogaster) | !is.na(Ocornuta)) %>% 
  separate_longer_delim(c("Dmelanogaster"), delim = ",") %>% 
  mutate(Dmelanogaster = str_remove_all(Dmelanogaster, " ")) %>% 
  separate_longer_delim("Ocornuta", delim = ",") %>% 
  mutate(Ocornuta = str_remove_all(Ocornuta, " ")) %>% 
  distinct()

dmel_gene_summaries <- read_tsv(
  dmel_gene_summaries_path, comment = "#", 
  col_names = c("Dmelanogaster", "symbol", "source", "summary"))

dmel_ncbi2ocor <- dmel2ocor %>% 
  left_join(dmel_fb2ncbi, by = c("Dmelanogaster" = "flybase")) %>%
  drop_na() %>% dplyr::select(ncbi, Ocornuta) %>% 
  dplyr::rename(Dmelanogaster = ncbi) %>% distinct()
dm_all_genes <- dmel_ncbi2ocor %>% pull(Dmelanogaster)

dme2keggpw <- read_tsv(dme2keggpw_path, col_names = c("dme", "keggpw")) %>% 
  dplyr::mutate(keggpw = str_split_i(keggpw, " - ", 1))
ncbi2dme <- read_tsv(ncbi2dme_path, col_names = c("dme", "kegg")) %>% 
  dplyr::mutate(dme = str_split_i(dme, ":", 2), 
                kegg = str_split_i(kegg, ":", 2))

ncbi_to_kegg <- read_tsv(dm_ncbi_to_kegg_path)

dme_kegg_combined <- ncbi2dme %>% full_join(ncbi_to_kegg) %>% 
  full_join(dme2keggpw) %>% dplyr::mutate(ncbi = as.numeric(ncbi))

gene2KEGG_df <- dmel_ncbi2ocor %>% 
  left_join(dme_kegg_combined, by = c("Dmelanogaster" = "ncbi"), 
            relationship = "many-to-many") %>%
  drop_na() %>% dplyr::select(keggpw, Ocornuta) %>% 
  dplyr::mutate(keggpw = ifelse(
    keggpw == "Glycosylphosphatidylinositol (GPI)-anchor biosynthesis", 
    "GPI-anchor biosynthesis", keggpw))
gene2KEGG <- create_gene2GO(gene2KEGG_df)




go_basic <- read_tsv(go_basic_path)

obic_gene2go <- read_tsv(gene2go_path, skip = 8) %>% 
  dplyr::left_join(go_basic, by = "GO_ID") %>% 
  dplyr::mutate(GeneID = Symbol) %>% distinct()

obic2ocor <- orthogroups %>% 
  dplyr::select("Obicornis", "Ocornuta") %>% 
  dplyr::filter(!is.na(Obicornis) | !is.na(Ocornuta)) %>% 
  separate_longer_delim(c("Obicornis"), delim = ",") %>% 
  mutate(Obicornis = str_remove_all(Obicornis, " ")) %>% 
  separate_longer_delim("Ocornuta", delim = ",") %>% 
  mutate(Ocornuta = str_remove_all(Ocornuta, " ")) %>% 
  distinct()


ocor_gene2go <- obic_gene2go %>% inner_join(obic2ocor, by = c("GeneID" = "Obicornis")) %>% 
  dplyr::select(-GeneID, -Symbol) %>% dplyr::rename(GeneID = "Ocornuta") %>% 
  dplyr::mutate(Symbol = GeneID) %>% 
  dplyr::select(`!#DB`, GeneID, Symbol, everything()) %>% distinct()


dmel_gene2go <- read_tsv(dmel_gene2go_path, comment = "!",
                         col_names = colnames(ocor_gene2go)) %>% 
  right_join(dmel2ocor, by = c("GeneID" = "Dmelanogaster"), 
             relationship = "many-to-many") %>% 
  dplyr::filter(!is.na(Ocornuta)) %>% dplyr::select(-GeneID) %>% 
  dplyr::rename(GeneID = Ocornuta) %>%
  dplyr::left_join(go_basic, by = "GO_ID") %>% 
  dplyr::select(GeneID, GO_ID, Aspect, Name, Definition) %>% 
  distinct()

ocor_gene2go <- ocor_gene2go %>% 
  dplyr::select(GeneID, GO_ID, Aspect, Name, Definition) %>% 
  distinct()

combined_gene2go <- ocor_gene2go %>% bind_rows(dmel_gene2go) %>% 
  dplyr::select(GeneID, GO_ID, Aspect, Name, Definition) %>% distinct() %>% 
  dplyr::filter(!(str_detect(Name, "positive regulation")) & 
                  !(str_detect(Name, "negative regulation")))

gene2bpGO_df <- combined_gene2go %>% dplyr::filter(Aspect == "P") %>% 
  dplyr::select(Name, GeneID)

topGO_gene2bpGO_df <- combined_gene2go %>% dplyr::filter(Aspect == "P") %>% 
  dplyr::select(GeneID, GO_ID)

gene2mfGO_df <- combined_gene2go %>% dplyr::filter(Aspect == "F") %>% 
  dplyr::select(Name, GeneID)

gene2ccGO_df <- combined_gene2go %>% dplyr::filter(Aspect == "C") %>% 
  dplyr::select(Name, GeneID)


gene2bpGO <- create_gene2GO(gene2bpGO_df)
gene2mfGO <- create_gene2GO(gene2mfGO_df)
gene2ccGO <- create_gene2GO(gene2ccGO_df)

#Obic: 33150
#Dmel: 97754
#Combined: 107248



save(gene2bpGO, gene2mfGO, gene2ccGO, gene2KEGG, file = go_rdata_out)

# Full data model including interaction ------------------------------------

# Testing
metadata_females <- metadata2 %>% 
  dplyr::filter(inferred_sex == "female")

dds_pre_testing <- load_files_and_make_dds(
  metadata_females, design = ~ temperature + time_point + temperature:time_point)
dds <- DESeq(dds_pre_testing, test = "LRT", reduced = ~ time_point)
dds_interaction <- DESeq(dds_pre_testing, test = "LRT", reduced = ~ temperature + time_point)

# Get results
res <- get_results(dds)
degs <- res %>% dplyr::arrange(padj) %>% 
  dplyr::filter(padj < 0.05) %>% dplyr::pull(geneID)
length(degs)

res_interaction <- get_results(dds_interaction)
degs_interaction <- res_interaction %>% dplyr::arrange(padj) %>% 
  dplyr::filter(padj < 0.05) %>% dplyr::pull(geneID)
length(degs_interaction)

length(intersect(degs, degs_interaction))

# Full-data PCAs ----------------------------------------------------------

vsd <- vst(dds)

time_pca <- plot_time_pca(vsd)
temp_pca <- plot_temp_pca(vsd, seed = 12)

combined_pca <- plot_grid(temp_pca, time_pca, labels = c("B", "B"))

ggsave(plot = combined_pca, filename = combined_pca_out,
       width = 183, height = 95, units = "mm",
       dpi = 600)

# Time point 1 (after one week) -------------------------------------------

# Load data subset --------------------------------------------------------

metadata_t1 <- metadata_females %>% 
  dplyr::filter(time_point == "pre_emerg_t1")

# Testing
dds_t1 <- load_files_and_make_dds(metadata_t1)
dds_t1 <- DESeq(dds_t1, test = "LRT", reduced = ~ 1)

# Get results
res_t1 <- get_results(dds_t1)
res_t1_summary <- res_t1 %>% dplyr::left_join(dmel2ocor, by = c("geneID" = "Ocornuta")) %>% 
  left_join(dmel_gene_summaries)
degs_t1 <- res_t1 %>% dplyr::arrange(padj) %>% 
  dplyr::filter(padj < 0.05) %>% dplyr::pull(geneID)
lfc_t1 <- res_t1 %>% dplyr::arrange(log2FoldChange.y) %>% 
  pull(log2FoldChange.y, geneID)


# ORA

# Time point 2 (after 7 weeks) --------------------------------------------

# Load data subset --------------------------------------------------------

metadata_t2 <- metadata_females %>% 
  dplyr::filter(time_point == "pre_emerg_t2")

# Testing
dds_t2 <- load_files_and_make_dds(metadata_t2)
dds_t2 <- DESeq(dds_t2, test = "LRT", reduced = ~ 1)

# Get results
res_t2 <- get_results(dds_t2)
res_t2_summary <- res_t2 %>% dplyr::left_join(dmel2ocor, by = c("geneID" = "Ocornuta")) %>% 
  left_join(dmel_gene_summaries)
degs_t2 <- res_t2 %>% dplyr::arrange(padj) %>% 
  dplyr::filter(padj < 0.05) %>% dplyr::pull(geneID)
lfc_t2 <- res_t2 %>% dplyr::arrange(desc(log2FoldChange.y)) %>% 
  pull(log2FoldChange.y, geneID)


degs_t2_up <- res_t2 %>% dplyr::arrange(padj) %>% 
  dplyr::filter(padj < 0.05, log2FoldChange.y > 0) %>% dplyr::pull(geneID)
degs_t2_down <- res_t2 %>% dplyr::arrange(padj) %>% 
  dplyr::filter(padj < 0.05, log2FoldChange.y < 0) %>% dplyr::pull(geneID)


# Combined GSEA -----------------------------------------------------------

MINGSSIZE = 10
MAXGSSIZE = 200

# GSEA
gsea_t1_bp <- fgsea(gene2bpGO, lfc_t1, minSize = 10, maxSize = MAXGSSIZE, 
                    gseaParam = 0) %>% 
  as_tibble()  %>% 
  dplyr::mutate(type = "BP") #%>% 
# dplyr::left_join(reducedTermsPostprocessed, by = c("pathway" = "term"))

gsea_t1_bp_filtered <- gsea_t1_bp %>% dplyr::filter(pval < 0.1) %>%
  filter_redundant() %>% 
  dplyr::arrange(desc(NES)) %>%
  dplyr::filter(NES > 0 & rank(NES) > n() - 5 | NES < 0 & rank(NES) <= 5)

gsea_t1_mf <- fgsea(gene2mfGO, lfc_t1, minSize = MINGSSIZE, maxSize = MAXGSSIZE, 
                    gseaParam = 0) %>% as_tibble()  %>% 
  dplyr::mutate(type = "MF")# %>% 
#  dplyr::left_join(reducedTermsPostprocessed, by = c("pathway" = "term"))

gsea_t1_mf_filtered <- gsea_t1_mf %>%   
  dplyr::filter(padj < 0.1) %>%
  filter_redundant() %>% 
  dplyr::arrange(desc(NES)) %>%
  dplyr::filter(NES > 0 & rank(NES) > n() - 3 | NES < 0 & rank(NES) <= 3)

gsea_t1_cc <- fgsea(gene2ccGO, lfc_t1, minSize = MINGSSIZE, maxSize = MAXGSSIZE, 
                    gseaParam = 0) %>% as_tibble() %>% 
  dplyr::mutate(type = "CC")# %>% 

gsea_t1_cc_filtered <- gsea_t1_cc %>% 
  dplyr::filter(padj < 0.1) %>%
  filter_redundant() %>% 
  dplyr::arrange(desc(NES)) %>%
  dplyr::filter(NES > 0 & rank(NES) > n() - 3 | NES < 0 & rank(NES) <= 3)

gsea_t1_kegg <- fgsea(gene2KEGG, lfc_t1, minSize = MINGSSIZE, maxSize = MAXGSSIZE, 
                      gseaParam = 0) %>% as_tibble() %>% 
  dplyr::mutate(type = "KEGG")

gsea_t1_kegg_filtered <- gsea_t1_kegg %>% 
  dplyr::filter(pval < 0.1) %>%
  filter_redundant() %>% 
  dplyr::arrange(desc(NES)) %>%
  dplyr::filter(NES > 0 & rank(NES) > n() - 4 | NES < 0 & rank(NES) <= 4)


gsea_t1_combined <- gsea_t1_kegg_filtered %>% 
  bind_rows(gsea_t1_bp_filtered) %>% 
  bind_rows(gsea_t1_mf_filtered) %>%
  bind_rows(gsea_t1_cc_filtered) %>% 
  dplyr::mutate(
    pathway = ifelse(
      pathway == "adenylate cyclase-activating adrenergic receptor signaling pathway",
      "adrenergic adenylate cyclase-activating pathway",
      pathway),
    pathway = ifelse(
      pathway == "proton-transporting ATP synthase activity, rotational mechanism",
      "ATP synthase activity",
      pathway))



# GSEA
gsea_t2_bp <- fgsea(gene2bpGO, lfc_t2, minSize = MINGSSIZE, maxSize = MAXGSSIZE, 
                    gseaParam = 0) %>% 
  as_tibble()  %>% 
  dplyr::mutate(type = "BP") #%>% 
# dplyr::left_join(reducedTermsPostprocessed, by = c("pathway" = "term"))

gsea_t2_bp_filtered <- gsea_t2_bp %>% dplyr::filter(padj < 0.1) %>%
  filter_redundant() %>% 
  dplyr::arrange(desc(NES)) %>%
  dplyr::filter(NES > 0 & rank(NES) > n() - 5 | NES < 0 & rank(NES) <= 5)

gsea_t2_mf <- fgsea(gene2mfGO, lfc_t2, minSize = MINGSSIZE, maxSize = MAXGSSIZE, 
                    gseaParam = 0) %>% as_tibble()  %>% 
  dplyr::mutate(type = "MF")# %>% 
#  dplyr::left_join(reducedTermsPostprocessed, by = c("pathway" = "term"))

gsea_t2_mf_filtered <- gsea_t2_mf %>%  
  dplyr::filter(padj < 0.1) %>%
  filter_redundant() %>% 
  dplyr::arrange(desc(NES)) %>%
  dplyr::filter(NES > 0 & rank(NES) > n() - 3 | NES < 0 & rank(NES) <= 3)

gsea_t2_cc <- fgsea(gene2ccGO, lfc_t2, minSize = MINGSSIZE, maxSize = MAXGSSIZE, 
                    gseaParam = 0) %>% as_tibble() %>% 
  dplyr::mutate(type = "CC")# %>% 

gsea_t2_cc_filtered <- gsea_t2_cc %>% 
  dplyr::filter(padj < 0.1) %>%
  filter_redundant() %>% 
  dplyr::arrange(desc(NES)) %>%
  dplyr::filter(NES > 0 & rank(NES) > n() - 3 | NES < 0 & rank(NES) <= 3)

gsea_t2_kegg <- fgsea(gene2KEGG, lfc_t2, minSize = MINGSSIZE, maxSize = MAXGSSIZE, 
                      gseaParam = 0) %>% as_tibble() %>% 
  dplyr::mutate(type = "KEGG")

gsea_t2_kegg_filtered <- gsea_t2_kegg %>% 
  dplyr::filter(padj < 0.1) %>%
  filter_redundant() %>% 
  dplyr::arrange(desc(NES)) %>%
  dplyr::filter(NES > 0 & rank(NES) > n() - 4 | NES < 0 & rank(NES) <= 4)


gsea_t2_combined <- gsea_t2_kegg_filtered %>% 
  bind_rows(gsea_t2_bp_filtered) %>% 
  bind_rows(gsea_t2_mf_filtered) %>%
  bind_rows(gsea_t2_cc_filtered) %>% 
  dplyr::mutate(
    pathway = ifelse(
      pathway == "detection of chemical stimulus involved in sensory perception of smell",
      "detection of chemical stimulus / smell",
      pathway
    ),
    pathway = ifelse(
      pathway == "proteasome-mediated ubiquitin-dependent protein catabolic process",
      "Ubiquitin-proteasome protein degradation",
      pathway
    ),
    pathway = ifelse(
      pathway == "mitochondrial electron transport, NADH to ubiquinone",
      "mitochondrial electron transport",
      pathway),
    pathway = ifelse(
      pathway == "mitochondrial respiratory chain complex I assembly",
      "mitochondrial respiratory chain assembly",
      pathway
    ))

combined_gsea <- plot_gseas(gsea_t1_combined, gsea_t2_combined)

make_label <- function(df, angle = 0){
  p <- ggplot(df) + geom_rect(aes(xmin = xmin, xmax = xmax, ymin = ymin, ymax = ymax),
                              fill = "grey70", color = "black")
  if(angle == 0){
    p <- p +
      geom_text(aes(x = (xmin + xmax)/2, y = 1.5, label = label), angle = angle, size = 5) +
      theme_void()
    return(p)
  } else if(angle == 90){
    p <- p + 
      geom_text(aes(y = (ymin + ymax)/2, x = 1.5, label = label), angle = angle, size = 5) +
      theme_void()
  }
  return(p)
}


week1_label <- make_label(data.frame(xmin = 1, xmax = 2, ymin = 1, ymax = 10, label = "Week 1"), angle = 90)
week7_label <- make_label(data.frame(xmin = 1, xmax = 2, ymin = 1, ymax = 10, label = "Week 7"), angle = 90)
side_labels <- plot_grid(NULL, week1_label, week7_label, NULL, nrow = 4, 
                         rel_heights = c(0.6,5.1,5.2,2))

withlabel <- plot_grid(side_labels, NULL, combined_gsea, NULL, ncol = 4, rel_widths = c(1, 1.4, 32, 0.5))

ggsave(plot = withlabel, filename = combined_gsea_out,
       width = 11.56, height = 6, scale = 1)

# Combined heatmap --------------------------------------------------------

# Heatmap
vsd_t1 <- vst(dds_t1, blind=FALSE)
metadata_t1_fixed <- fix_names(metadata_t1)
heatmap_t1 <- plot_heatmap(vsd_t1, degs_t1, metadata_t1_fixed, anno = "Week 1")
#ggsave(plot = heatmap_t1, filename = "cornuta/heatmap_t1.png")

# Heatmap
vsd_t2 <- vst(dds_t2, blind=FALSE)
metadata_t2_fixed <- fix_names(metadata_t2)
heatmap_t2 <- plot_heatmap(vsd_t2, degs_t2, metadata_t2_fixed, anno = "Week 7")
#ggsave(plot = heatmap_t2, filename = "cornuta/heatmap_t2.png")

first_col <- cowplot::plot_grid(heatmap_t1$plot, NULL, nrow = 2, 
                                rel_heights = c(heatmap_t1$n_genes + 750, 
                                                heatmap_t2$n_genes - heatmap_t1$n_genes + 750))

combined_heatmap <- cowplot::plot_grid(first_col, heatmap_t2$plot, nrow = 1,
                                       labels = c("F", "G"))

#first_plot <- plot_grid(heatmap_t1$plot, NULL, ncol = 2, rel_widths = c(10, 2))

ggsave(plot = combined_heatmap, filename = heatmap_combined_out,
       width = 183, height = 200, units = "mm", scale = 0.8,
       dpi = 600)

# Sex-specific analysis ---------------------------------------------------

metadata_sexes <- metadata2 %>% dplyr::filter(time_point == "pre_emerg_t2")

# Testing
dds_sexes_pre_testing <- load_files_and_make_dds(
  metadata_sexes, design = 
    ~ temperature + inferred_sex + temperature:inferred_sex)
dds_sexes <- DESeq(dds_sexes_pre_testing, test = "LRT", 
                   reduced = ~ temperature)
dds_sexes_interaction <- DESeq(dds_sexes_pre_testing, test = "LRT", 
                               reduced = ~ temperature + inferred_sex)

# Get results
res_sexes <- get_results(dds_sexes, contrast2 = c("inferred_sex", "female", "male"))
degs_sexes <- res_sexes %>% dplyr::arrange(padj) %>% 
  dplyr::filter(padj < 0.05) %>% dplyr::pull(geneID)
length(degs_sexes)
lfc_sexes <- res_sexes %>% dplyr::arrange(desc(log2FoldChange.y)) %>% 
  pull(log2FoldChange.y, geneID)

res_sexes_interaction <- get_results(dds_sexes_interaction)
degs_sexes_interaction <- res_sexes_interaction %>% dplyr::arrange(padj) %>% 
  dplyr::filter(padj < 0.05) %>% dplyr::pull(geneID)
length(degs_sexes_interaction)

length(intersect(degs, degs_sexes_interaction))



save(res, res_t1, res_t2, res_sexes, res_interaction, 
     res_sexes_interaction, metadata, file = expression_results_out)



# Add annotation ----------------------------------------------------------
Ocor_anno <- read_tsv(annotations_path)

res_t1_anno <- res_t1 %>% dplyr::filter(padj < 0.05) %>% 
  left_join(Ocor_anno) %>% distinct() %>% 
  dplyr::filter(feature != "transcript") %>% dplyr::select(-baseMean.y, -stat) %>% 
  dplyr::rename(LFC.normal.plus3 = log2FoldChange.x, 
                LFC.normal.plus3.SE = lfcSE.x, 
                LFC.normal.plus6 = log2FoldChange.y, 
                LFC.normal.plus6.SE = lfcSE.y, 
                baseMean = baseMean.x)
res_t1_anno %>% 
  write_tsv(res_t1_out)

res_t2_anno <- res_t2 %>% dplyr::filter(padj < 0.05) %>% 
  left_join(Ocor_anno) %>% distinct() %>% 
  dplyr::filter(feature != "transcript") %>% dplyr::select(-baseMean.y, -stat) %>% 
  dplyr::rename(LFC.normal.plus3 = log2FoldChange.x, 
                LFC.normal.plus3.SE = lfcSE.x, 
                LFC.normal.plus6 = log2FoldChange.y, 
                LFC.normal.plus6.SE = lfcSE.y, 
                baseMean = baseMean.x)
res_t2_anno %>% 
  write_tsv(res_t2_out)

res_sexes_anno <- res_sexes %>% dplyr::filter(padj < 0.05) %>% 
  left_join(Ocor_anno) %>% distinct() %>% 
  dplyr::filter(feature != "transcript") %>% dplyr::select(
    -baseMean.y, -stat, -type, -log2FoldChange.y, -lfcSE.y) %>% 
  dplyr::rename(LFC.male.female = log2FoldChange.x, 
                LFC.male.female.SE = lfcSE.x, 
                baseMean = baseMean.x)
res_sexes_anno %>% 
  write_tsv(res_sexes_out)


res_sexes_interaction_anno <- res_sexes_interaction %>% 
  dplyr::filter(padj < 0.05) %>% left_join(Ocor_anno) %>% distinct() %>% 
  dplyr::filter(feature != "transcript") %>% dplyr::select(
    -baseMean.y, -stat, -type, -log2FoldChange.y, -lfcSE.y) %>% 
  dplyr::rename(LFC.male.female = log2FoldChange.x, 
                LFC.male.female.SE = lfcSE.x, 
                baseMean = baseMean.x)
res_sexes_interaction_anno %>% 
  write_tsv(res_sexes_interaction_out)

list(
  temperature.week1 = res_t1_anno, temperature.week7 = res_t2_anno, 
  temperatureXsex.week7 = res_sexes_interaction_anno, 
  sex.week7 = res_sexes_anno) %>% 
  write_xlsx(degs_xlsx_out)

