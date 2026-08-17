# Plots ADMIXTURE cross-validation error against K for both species, to help pick the
# most likely number of ancestral populations.
library(tidyverse)
library(cowplot)

bicornis_summary <- snakemake@input[["bicornis_summary"]]
cornuta_summary <- snakemake@input[["cornuta_summary"]]
out_png <- snakemake@output[["png"]]

bicornis_cv <- read_delim(bicornis_summary, delim = ": ",
                          col_names = c("K", "Err")) %>%
  dplyr::mutate(K = as.integer(str_extract(K, "\\d+"))) %>%
  dplyr::arrange(K) %>% dplyr::mutate(species = "bicornis")


cornuta_cv <- read_delim(cornuta_summary, delim = ": ",
                          col_names = c("K", "Err")) %>%
  dplyr::mutate(K = as.integer(str_extract(K, "\\d+"))) %>%
  dplyr::arrange(K) %>% dplyr::mutate(species = "cornuta")

combined <- bicornis_cv %>% dplyr::bind_rows(cornuta_cv)

combined_plot <- ggplot(combined) + aes(x = K, y = Err, color = species, fill = species) + 
  geom_line() + geom_point(color = "black", shape = 21) + 
  scale_x_continuous(breaks = 1:10) + 
  scale_fill_manual(values = c("#F5B829", "#7D8AED")) + 
  scale_color_manual(values = c("#F5B829", "#7D8AED")) + 
  theme_minimal() + 
  ylab("CV error") + 
  xlab("Putative number of ancestral populations (K)") + 
  theme(
    legend.title = element_blank(),
    legend.position = "right",
    legend.direction = "vertical",
    legend.key.spacing.y = unit(3, "mm"),
    legend.text = element_text(angle = 90, margin = margin(0.1, unit = "mm")),
    legend.text.position = "bottom",
    legend.key.width = unit(1, "mm"),
    panel.border = element_rect(
      color = "black", fill = NA, linewidth = 0.5))

combined_plot2 <- plot_grid(combined_plot, nrow = 1, labels = "E")

ggsave(plot = combined_plot2, filename = out_png, dpi = 600,
       width = 183, height = 45, units = "mm")
