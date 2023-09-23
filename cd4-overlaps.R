#####################################
## read in expression / DE results ##
#####################################

library(SummarizedExperiment)
se <- readRDS("pseudobulk_for_michael_cd4_naive.rds")

library(tidySummarizedExperiment)
library(tidybulk)
library(tidyr)

gene_tab <- se |>
  pivot_transcript() |>
  mutate(pb_ave_count = rowMeans(assay(se, "counts_scaled"))) |>
  dplyr::select(gene=.feature, padj = P_sex_adjusted___cd4.naive, pb_ave_count) |>
  drop_na()

######################
## load gene ranges ##
######################

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
  plyranges::select(gene, entrez=gene_id)

g <- g |>
  filter(gene %in% gene_tab$gene)

mcols(g) <- mcols(g) |>
  as_tibble() |>
  left_join(gene_tab)

g <- keepStandardChromosomes(g, pruning.mode = "coarse")

##########################
## read in GWAS results ##
##########################

library(readr)
diseases <- c("MS","RA","SLE")
csvs <- lapply(diseases, \(d) read_csv(paste0(d,"_hg19_gwas.csv")))
names(csvs) <- diseases

library(GenomeInfoDb)
ms <- csvs[["MS"]] |>
  plyranges::select(seqnames=Chromosome,
         start=`Position (hg19)`,
         gene=`Proximal Gene(s)`,
         type=Type,
         rsid=Effect) |>
  mutate(gene = sub("\\(.*","",gene), width=1) |>
  as_granges()

ra <- csvs[["RA"]] |>
  plyranges::select(seqnames=`Chr.`,
         start=Position,
         gene=`Gene name`,
         type=Gene,
         rsid=`Rs ID`) |>
  mutate(gene = sub("\\(.*","",gene), width=1) |>
  as_granges()

sle <- csvs[["SLE"]] |>
  plyranges::select(seqnames=Chr,
         start=Pos,
         gene=Gene,
         type=Annotation,
         rsid=rsid) |>
  mutate(gene = sub(";.*","",gene), width=1) |>
  as_granges()

gwas <- bind_ranges(ms=ms, ra=ra, sle=sle, .id="disease")
seqlevelsStyle(gwas) <- "UCSC"
seqlevels(gwas) <- seqlevels(g)
seqinfo(gwas) <- seqinfo(g)

table(gwas$type, gwas$disease)

library(forcats)
gwas <- gwas |>
  mutate(
    type = fct_collapse(
      type,
      intronic = c("intronic","intron"),
      exonic = c("exonic","nonsynonymous","synonymous","synonymous SNV","missense")
      ))

table(gwas$type, gwas$disease)

gwas <- gwas |>
  mutate(pos = start, gwas_gene=gene) |>
  plyranges::select(-gene)

# save(g, gwas, file="intermediate.rda")

####################################
## overlap DE genes and GWAS data ##
####################################

library(plyranges)
library(dplyr)

res <- g |>
  filter(padj < .1) |>
  join_overlap_inner(gwas, maxgap=5e4) |>
  mutate(tss_dist = ifelse(strand == "+", pos - start, end - pos)) |>
  as_tibble() |>
  dplyr::select(disease, chr=seqnames, de_gene=gene, gwas_gene, padj,
         pb_ave_count, rsid, type, tss_dist) |>
  arrange(disease, chr) 

print(res, n=100)

write.csv(res, file="results.csv", quote=FALSE, row.names=FALSE)

library(ggplot2)
library(ggrepel)
res |>
  filter(de_gene == gwas_gene) |>
  mutate(label = paste0(de_gene, ", ", round(padj,3))) |>
  ggplot(aes(tss_dist, pb_ave_count, col=disease, shape=type, label=label)) +
  geom_point(size=2) +
  geom_text_repel(show.legend=FALSE) +
  xlab("distance to DE gene TSS") +
  ylab("pseudobulk average count") +
  scale_y_log10(limits=c(1,1e5))

