library(tidyverse)
library(ggrepel)

args <- commandArgs(trailingOnly = TRUE)
stats_dir <- args[1]
out_dir <- args[2]


# Variant quality ---------------------------------------------------------

var_qual <- read_delim(paste0(stats_dir, "/stats.lqual"), delim = "\t",
                       col_names = c("chr", "pos", "qual"), skip = 1)

a <- ggplot(var_qual, aes(qual)) +
  geom_density(fill = "dodgerblue1", colour = "black", alpha = 0.3) +
  theme_light() + xlim(c(0, 500))

ggsave(paste0(out_dir, "/lqual.png"), a)


# Variant mean depth ------------------------------------------------------

var_depth <- read_delim(paste0(stats_dir, "/stats.ldepth.mean"), delim = "\t",
                        col_names = c("chr", "pos", "mean_depth", "var_depth"),
                        skip = 1)

b <- ggplot(var_depth, aes(mean_depth)) +
  geom_density(fill = "dodgerblue1", colour = "black", alpha = 0.3) +
  theme_light() + xlim(c(0, 50))

ggsave(paste0(out_dir, "/ldepth.mean.png"), b)


# Variant missingness -----------------------------------------------------

var_miss <- read_delim(paste0(stats_dir, "/stats.lmiss"), delim = "\t",
                       col_names = c("chr", "pos", "nchr",
                                     "nfiltered", "nmiss", "fmiss"), skip = 1)
c <- ggplot(var_miss, aes(fmiss)) +
  geom_density(fill = "dodgerblue1", colour = "black", alpha = 0.3) +
  theme_light()

ggsave(paste0(out_dir, "/lmiss.png"), c)


# Minor allele frequency --------------------------------------------------

var_freq <- read_delim(paste0(stats_dir, "/stats.frq"), delim = "\t",
                       col_names = c("chr", "pos", "nalleles",
                                     "nchr", "a1", "a2"), skip = 1)

var_freq$maf <- var_freq %>% select(a1, a2) %>% apply(1, function(z) min(z))

d <- ggplot(var_freq, aes(maf)) +
  geom_density(fill = "dodgerblue1", colour = "black", alpha = 0.3) +
  theme_light()

ggsave(paste0(out_dir, "/frq.png"), d)


# Per individual depth ----------------------------------------------------


ind_depth <- read_delim(paste0(stats_dir, "/stats.idepth"), delim = "\t",
                        col_names = c("ind", "nsites", "depth"), skip = 1)

e <- ggplot(ind_depth, aes(depth)) +
  geom_histogram(fill = "dodgerblue1", colour = "black", alpha = 0.3) +
    geom_text_repel(aes(y = 1, label = ind),
                position = "dodge", max.time = 20, max.overlaps = 30) +
  theme_light()

ggsave(paste0(out_dir, "/idepth.png"), e)


# Per individual missingness ----------------------------------------------

ind_miss  <- read_delim(paste0(stats_dir, "/stats.imiss"), delim = "\t",
                        col_names = c("ind", "ndata", "nfiltered",
                                      "nmiss", "fmiss"), skip = 1)

f <- ggplot(ind_miss, aes(x = fmiss)) +
  geom_histogram(aes(y = after_stat(count)), fill = "dodgerblue1", colour = "black", alpha = 0.3) +
  geom_text_repel(aes(y = 1, label = ind), position = "dodge", max.time = 20, max.overlaps = 30) +
  theme_light()

ggsave(paste0(out_dir, "/imiss.png"), f)


# Per individual heterozygosity / inbreeding (f) --------------------------

ind_het <- read_delim(paste0(stats_dir, "/stats.het"), delim = "\t",
                      col_names = c("ind","ho", "he", "nsites", "f"), skip = 1)

g <- ggplot(ind_het, aes(f)) + 
  geom_histogram(fill = "dodgerblue1", colour = "black", alpha = 0.3) + 
  geom_text_repel(aes(y = 1, label = ind), 
                  position = "dodge", max.time = 20, max.overlaps = 30) + 
  theme_light()

ggsave(paste0(out_dir, "/het.png"), g)
