# Plots heterozygosity (to find diploid males) and depth-by-sex, used to sanity-check
# the ploidy assignments in sample2ploidy-final.tsv. One script shared by both species.
library(tidyverse)

het <- read_tsv(snakemake@input[["het"]]) %>%
  mutate(INDV = str_replace(INDV, "_", "-"))

# samples.tsv (not a separate metadata.csv copy): its `metadata_id` column already
# carries the original hyphenated IDs this script's join key expects.
meta <- read_tsv(snakemake@input[["samples"]]) %>%
  rename(ID = metadata_id)

het2 <- het %>% left_join(meta, by = c("INDV" = "ID"))

p_het <- het2 %>% ggplot(aes(x = `O(HOM)`, y = `E(HOM)`, fill = sex)) +
  geom_point(alpha = 0.5, size = 3) +
  ggrepel::geom_label_repel(aes(label = INDV), force = 10, max.overlaps = 50)

ggsave(snakemake@output[["het_png"]], p_het)

depth <- read_tsv(snakemake@input[["depth"]]) %>%
  mutate(INDV = str_replace(INDV, "_", "-"))

depth2 <- depth %>% left_join(meta, by = c("INDV" = "ID"))

p_depth <- depth2 %>% ggplot(aes(x = MEAN_DEPTH, fill = sex)) + geom_density(alpha = 0.5)

ggsave(snakemake@output[["depth_png"]], p_depth)
