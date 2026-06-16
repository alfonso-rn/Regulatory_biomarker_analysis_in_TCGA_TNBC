# 4. Differential expression analysis: TNBC vs non-TNBC ----
#
# Differential expression analysis is performed between clinically defined TNBC
# and non-TNBC primary tumor samples.
#
# DEGs are identified using FDR < 0.01 and |log2FC| >= 1, corresponding to
# statistically significant genes with at least a two-fold expression change.
#
# Since TNBC is provided as the second condition, positive logFC values are
# interpreted as higher expression in TNBC compared with non-TNBC.
#
# The resulting DEGs are summarized by regulation direction and gene type.

library(TCGAbiolinks)
library(dplyr)
library(ggplot2)

dataFilt <- readRDS("data/processed/dataFilt.rds")
samplesTNBC <- readRDS("data/processed/samplesTNBC.rds")
samplesNonTNBC <- readRDS("data/processed/samplesNonTNBC.rds")

# 4.1. Differential expression analysis (DEA) ----

dataDEA <- TCGAanalyze_DEA(
  mat1 = dataFilt[, samplesNonTNBC],
  mat2 = dataFilt[, samplesTNBC],
  Cond1type = "Non_TNBC",
  Cond2type = "TNBC",
  fdr.cut = 0.01, 
  logFC.cut = 1,
  method = "glmLRT"
)

cat("Total DEGs:", nrow(dataDEA), "\n",
    "Upregulated DEGs in TNBC:", sum(dataDEA$logFC > 0), "\n",
    "Downregulated DEGs in TNBC:", sum(dataDEA$logFC < 0))

# 4.2. DEGs table with expression values ----

dataDEGsLevel <- TCGAanalyze_LevelTab(
  FC_FDR_table_mRNA = dataDEA,
  typeCond1 = "Non_TNBC",
  typeCond2 = "TNBC",
  TableCond1 = dataFilt[, samplesNonTNBC],
  TableCond2 = dataFilt[, samplesTNBC]
)

# Add expression TNBC status column
dataDEA$expression_status <- ifelse(
  dataDEA$logFC > 0,
  "Up regulated in TNBC",
  "Down regulated in TNBC"
)

dataDEGsLevel$expression_status <- ifelse(
  dataDEGsLevel$logFC > 0,
  "Up regulated in TNBC",
  "Down regulated in TNBC"
)

saveRDS(dataDEA, "data/processed/dataDEA.rds")
saveRDS(dataDEGsLevel, "data/processed/dataDEGsLevel.rds")

# 4.3. Volcano plot representation ----

# Keep only one row per gene_name
dataDEA <- dataDEA[order(dataDEA$gene_name, -abs(dataDEA$logFC), dataDEA$FDR),]
dataDEA <- dataDEA[!duplicated(dataDEA$gene_name), ]

# Select top genes with |logFC| > 6
top_genes <- dataDEA$gene_name[abs(dataDEA$logFC) > 6]
cat("Genes with |logFC| > 6:", length(top_genes))
print(top_genes)

options(ggrepel.max.overlaps = Inf)

TCGAVisualize_volcano(
  x = dataDEA$logFC,
  y = dataDEA$FDR,
  filename = "results/figures/volcano_TNBC_vs_NonTNBC.png",
  x.cut = 1,
  y.cut = 0.01,
  names = dataDEA$gene_name,
  highlight = top_genes,
  show.names = "highlighted",
  color = c("grey", "red", "blue"),
  names.size = 1.5,
  xlab = "Gene expression fold change (Log2)",
  legend = "State",
  title = "Differential expression analysis: TNBC vs Non-TNBC",
  width = 10
)

# 4.4. Horizontal bar plot representation ----

barplot_gene_type <- ggplot(
  data = dataDEA,
  mapping = aes(
    x = factor(
      ifelse(logFC > 0, "Upregulated in TNBC", "Downregulated in TNBC"),
      levels = c("Downregulated in TNBC", "Upregulated in TNBC")
    ),
    fill = gene_type
  )
) +
  geom_bar(position = "fill") +
  coord_flip() +
  scale_y_continuous(
    labels = function(x) paste0(round(x * 100), "%")
  ) +
  labs(
    title = "Distribution of gene types among DEGs",
    subtitle = "TNBC vs non-TNBC",
    x = NULL,
    y = "Percentage of DEGs",
    fill = "Gene type"
  ) +
  theme_classic() + theme(legend.position = "bottom")

ggsave(filename = "results/figures/barplot_DEGs_type.png",
       plot = barplot_gene_type, width = 10, height = 6, dpi = 300 )

# 4.5. Gene type distribution summary table ----

smallRNA <- c("miRNA", "misc_RNA", "snRNA", "snoRNA", "scaRNA", "rRNA", "Mt_tRNA")

pseudogenes <- c("processed_pseudogene", "unprocessed_pseudogene", "unitary_pseudogene", "polymorphic_pseudogene",
  "transcribed_processed_pseudogene", "transcribed_unprocessed_pseudogene", "transcribed_unitary_pseudogene",
  "rRNA_pseudogene", "IG_V_pseudogene", "IG_C_pseudogene", "IG_J_pseudogene", "TR_V_pseudogene")

immune <- c("IG_V_gene", "IG_C_gene", "IG_J_gene", "IG_D_gene", "TR_V_gene", "TR_J_gene", "TR_C_gene")

gene_group_map <- c(
  protein_coding = "Protein-coding genes",
  lncRNA = "Long non-coding RNAs",
  setNames(rep("Small non-coding RNAs", length(smallRNA)), smallRNA),
  setNames(rep("Pseudogenes", length(pseudogenes)), pseudogenes),
  setNames(rep("Immune receptor genes", length(immune)), immune),
  TEC = "Uncertain"
)

dataDEA$gene_group <- unname(gene_group_map[dataDEA$gene_type])

gene_group_summary <- dataDEA %>%
  count(gene_group, gene_type, expression_status, name = "n") %>%
  as.data.frame()

gene_group_summary <- reshape(
  gene_group_summary,
  idvar = c("gene_group", "gene_type"),
  timevar = "expression_status",
  direction = "wide"
)

gene_group_summary[is.na(gene_group_summary)] <- 0

names(gene_group_summary) <- sub("n.Down regulated in TNBC", "n downregulated",
  names(gene_group_summary), fixed = TRUE)

names(gene_group_summary) <- sub("n.Up regulated in TNBC", "n upregulated",
  names(gene_group_summary), fixed = TRUE)

gene_group_summary <- gene_group_summary[ , c("gene_group", "gene_type", "n downregulated", "n upregulated")]

rownames(gene_group_summary) <- NULL

write.csv(gene_group_summary, "results/tables/gene_group_summary.csv")
