# 7. Functional enrichment analysis of CRM-associated genes ----
#
# Functional enrichment analysis (FEA) was performed on differentially expressed
# genes (DEGs) in TNBC using g:Profiler. Upregulated and downregulated genes were
# analysed separately to identify enriched Gene Ontology terms and biological pathways.
#
# Enriched terms containing genes associated with prioritized intragenic CRMs 
# were then identified to explore functional categories potentially linked to 
# regulatory regions affected by intronic mutations in TNBC.
#
# Enrichment sources:
# GO:BP = Biological Process / GO:MF = Molecular Function / GO:CC = Cellular Component
# REAC = Reactome pathways
# KEGG = Kyoto Encyclopedia of Genes and Genomes pathways
# WP = WikiPathways

library(dplyr)
library(gprofiler2)
library(ggplot2)

dataDEA <- readRDS("data/processed/dataDEA.rds")
intragenicCRMs <- readRDS("data/processed/intragenicCRMs.rds")

# 7.1. Functional enrichment analysis (FEA) ----

genes_up <- dataDEA %>%
  filter(expression_status == "Up regulated in TNBC") %>%
  pull(gene_name)

genes_down <- dataDEA %>%
  filter(expression_status == "Down regulated in TNBC") %>%
  pull(gene_name)

FEA <- gost(
  query = list(Upregulated_TNBC = genes_up, Downregulated_TNBC = genes_down),
  organism = "hsapiens",
  sources = c("GO:BP", "GO:MF", "GO:CC", "REAC", "KEGG", "WP"),
  correction_method = "fdr",
  evcodes = TRUE
)

# Results summary tables
FEA_DEGs <- FEA$result

FEA_DEGs_terms <- FEA_DEGs %>% count(query, source, name = "significant_terms")

write.csv(FEA_DEGs_terms, "results/tables/FEA_DEGs_terms.csv", row.names = FALSE)

FEA_top_abundant <- FEA_DEGs %>%
  arrange(desc(intersection_size)) %>%
  slice_head(n = 10) %>%
  mutate(parents = as.character(parents))

write.csv(FEA_top_abundant, "results/tables/FEA_top_abundant.csv", row.names = FALSE)

# Manhattan-like-plot
DEGs_FEAplot1 <- gostplot(FEA, capped = TRUE, interactive = FALSE)

DEGs_FEAplot2 <- publish_gostplot(DEGs_FEAplot1, 
                                  highlight_terms = FEA_top_abundant$term_id, 
                                  width = NA, height = NA, filename = NULL )

ggsave(filename = "results/figures/DEGs_FEAplot.png", plot = DEGs_FEAplot2,
       width = 10, height = 7, dpi = 300, bg = "white")

# 7.2. Identification of enriched terms containing self-targeting CRM genes ----

genesCRM <- intragenicCRMs$common_genes

FEA_DEGs$CRM_genes_in_term <- NA_character_

for (i in seq_len(nrow(FEA_DEGs))) {
  
  genes_intersection <- unlist(strsplit(FEA_DEGs$intersection[i], ","))
  genes_intersection <- trimws(genes_intersection)
  
  CRM_genes <- intersect(genes_intersection, genesCRM)
  
  if (length(CRM_genes) > 0) {
    FEA_DEGs$CRM_genes_in_term[i] <- paste(CRM_genes, collapse = "; ")
  }
}

FEA_intragenicCRMs <- FEA_DEGs[!is.na(FEA_DEGs$CRM_genes_in_term), ]
FEA_intragenicCRMs$parents <- as.character(FEA_intragenicCRMs$parents)

write.csv(FEA_intragenicCRMs, "results/tables/FEA_intragenicCRMs.csv", row.names = FALSE)

# 7.3. Top 5 enriched terms associated with self-targeting CRM genes ----

FEA_intragenicCRMs_top <- FEA_intragenicCRMs %>%
  mutate(
    CRM_gene_count = lengths(strsplit(CRM_genes_in_term, "; ")),
    minus_log10_p = -log10(p_value),
    term_label = paste(term_name, source, CRM_genes_in_term, sep = " | ")
  ) %>%
  arrange(query, source, p_value) %>%
  group_by(query, source) %>%
  slice_head(n = 5) %>%
  ungroup()

FEA_intragenicCRMs_dotplot <- ggplot(
  FEA_intragenicCRMs_top,
  aes(
    x = minus_log10_p,
    y = reorder(term_label, minus_log10_p),
    size = CRM_gene_count
  )
) +
  geom_point() +
  facet_wrap(~ query, scales = "free_y") +
  labs(
    x = "-log10(FDR-adjusted p-value)",
    y = NULL,
    size = "CRM genes\nin term"
  ) +
  theme_bw() +
  theme(
    axis.text.y = element_text(size = 6),
    strip.text = element_text(face = "bold")
  )

ggsave("results/figures/FEA_intragenicCRMs_dotplot.png", FEA_intragenicCRMs_dotplot,
       width = 13, height = 9, dpi = 300, bg = "white")
