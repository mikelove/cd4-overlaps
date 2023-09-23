library(SummarizedExperiment)
se <- readRDS("pseudobulk_for_michael_cd4_naive.rds")

library(tidySummarizedExperiment)
library(tidybulk)
library(tidyr)

gene_tab <- se |>
  pivot_transcript() |>
  mutate(log10_exprs = log10(rowMeans(assay(se, "counts_scaled")) + 1)) |>
  select(gene=.feature, padj = P_sex_adjusted___cd4.naive, log10_exprs) |>
  drop_na()

# GWAS data is hg19
library(TxDb.Hsapiens.UCSC.hg19.knownGene)
g <- genes(TxDb.Hsapiens.UCSC.hg19.knownGene)

library(plyranges)
library(org.Hs.eg.db)

g <- g |>
  mutate(gene = mapIds(org.Hs.eg.db,
                       keys=gene_id,
                       column="SYMBOL",
                       keytype="ENTREZID")) |>
  select(gene, entrez=gene_id)

g <- g |>
  filter(gene %in% gene_tab$gene)

mcols(g) <- mcols(g) |>
  as_tibble() |>
  left_join(gene_tab)

