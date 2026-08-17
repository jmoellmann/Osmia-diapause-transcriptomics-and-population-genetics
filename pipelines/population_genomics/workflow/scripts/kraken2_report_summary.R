# Genus-level summary of Kraken2 contamination screening across all samples: filters
# to genera present in at least 5 samples, then joins against sample metadata for the
# species-level breakdown.
library(tidyverse)

df <- read_csv(
  snakemake@input[["kraken_summary"]],
  col_names = c("sample", "perc", "count_clade", "count_taxon",
                "rank", "ncbi_id", "name")) %>%
  # sample names in the report may carry umlauts from raw filenames; normalize before
  # extracting the clean sample ID so the join against samples.tsv matches.
  mutate(sample = stringi::stri_replace_all_fixed(
    sample,
    c("Ä", "Ö", "Ü"),
    c("A", "O", "U"),
    vectorize_all = FALSE
  )) %>%
  mutate(sample = str_extract(sample, '[A-Z]+[0-9]?_[0-9]{2}'))

sample2species <- read_tsv(snakemake@input[["samples"]]) %>%
  select(sample = id, species)

df <- df %>% left_join(sample2species, by = "sample")

filtered_df <- df %>% filter(count_clade > 10, rank == "G") %>% group_by(name) %>% 
  mutate(n = n()) %>% filter(n >= 5) %>% ungroup() %>% 
  filter(species %in% c("bicornis", "cornuta")) %>% 
  mutate(name_n = str_c(name, "(", n, ")")) %>% 
  arrange(desc(n))

summary_tbl <- filtered_df %>% group_by(name) %>%
  summarise(mean_count = mean(count_clade), max_count = max(count_clade),
            mean_perc = mean(perc)) %>%
  arrange(desc(mean_count))

write_tsv(summary_tbl, snakemake@output[["summary_tsv"]])

plot1 <- filtered_df %>%
  arrange(desc(perc)) %>%
  ggplot(aes(y = name_n, x = count_clade, color = species)) + geom_boxplot()
ggsave(snakemake@output[["count_png"]], plot1)

plot2 <- filtered_df %>%
  arrange(desc(perc)) %>%
  ggplot(aes(y = name_n, x = perc, color = species)) + geom_boxplot()
ggsave(snakemake@output[["perc_png"]], plot2)
