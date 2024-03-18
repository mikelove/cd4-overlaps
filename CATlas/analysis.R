library(here)
library(readr)
library(dplyr)

genome_tab_file <- "bindata.1000.hg38.tsv.gz"
gwas_tab_file <- "gwas-t1d-031824-EFO_0009756.tsv.gz"
ccre_bed_file <- "cCREs.bed.gz"

genome_tab <- read_delim(here("CATlas",genome_tab_file))
genome_tab <- genome_tab |>
  select(chr, start=start0, end, N, GC, tss, ldscore.sum) |>
  mutate(start = start + 1)
genome_tab <- genome_tab |>
  filter(N == 0, !is.na(ldscore.sum))

gwas_tab <- read_delim(here("CATlas",gwas_tab_file))
