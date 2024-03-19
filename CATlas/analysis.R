library(here)
library(readr)
library(dplyr)

bindata_file <- "bindata.1000.hg38.tsv.gz"
gwas_tab_file <- "gwas-t1d-031824-EFO_0009756.tsv.gz"
ccre_bed_file <- "T-lymphocyte-2-CD4.bed.gz"

bindata <- read_delim(here("CATlas",bindata_file))
bindata <- bindata |>
  select(seqnames=chr, start=start0, end,
         N, GC, tss, ldscore.sum) |>
  mutate(start = start + 1)
bindata <- bindata |>
  filter(N == 0, !is.na(ldscore.sum))

gwas_tab <- read_delim(here("CATlas",gwas_tab_file))
gwas_tab <- gwas_tab |>
  select(seqnames=CHR_ID, start=CHR_POS,
         gene=`REPORTED GENE(S)`,
         snp=SNPS, pvalue=`P-VALUE`) |>
  filter(!is.na(start)) |>
  mutate(width=1)

library(plyranges)
ccre <- read_narrowpeaks(here("CATlas",ccre_bed_file))
genome(ccre) <- "hg38"

bindata <- bindata |>
  as_granges()
genome(bindata) <- "hg38"

gwas <- gwas_tab |>
  as_granges()
genome(gwas) <- "GRCh38"

library(GenomeInfoDb)
seqlevelsStyle(gwas) <- "UCSC"
gwas <- gwas |>
  sortSeqlevels() |>
  sort()

ccre_small <- ccre |> select(score)

mhc <- data.frame(seqnames="chr6", start=29e6, end=33e6) |>
  as_granges()

gwas |>
  filter_by_overlaps(mhc)

gwas |>
  join_overlap_left(ccre_small, maxgap=1e4) |>
  filter(!duplicated(snp)) |>
  summarize(nhit = sum(is.na(score)))

library(plotgardener)
gwas_for_plot <- gwas |>
  join_overlap_left(ccre_small, maxgap=1e4) |>
  filter(!duplicated(snp))

gwas_for_plot <- gwas_for_plot |>
  as_tibble() |>
  dplyr::rename(pos=start) |>
  mutate(chrom = as.character(seqnames)) |>
  mutate(p = pmax(pvalue, 1e-50))

pageCreate(width = 5, height = 3)
manh <- plotManhattan(
  data = gwas_for_plot, assembly = "hg38",
  fill = c("grey30","grey70"), sigLine = TRUE,
  x = 0.5, y = 0.5, width = 4, height = 2,
)
annoGenomeLabel(plot = manh, x = 0.5,
                y = 2.5, fontsize=8)
annoYaxis(plot = manh, at = 0:10 * 5,
          axisLine = TRUE, fontsize = 8)
