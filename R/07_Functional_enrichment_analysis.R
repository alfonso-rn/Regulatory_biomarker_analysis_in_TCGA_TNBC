# 7. Functional enrichment analysis of CRM-associated genes ----
#
# Functional enrichment analysis (FEA) was performed on differentially expressed
# genes (DEGs) in TNBC using g:Profiler. Upregulated and downregulated genes were
# analysed separately to identify enriched Gene Ontology terms and biological pathways.
#
# Enriched terms containing genes associated with prioritized self-targeting CRMs 
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
selfTargetingCRMs <- readRDS("data/processed/selfTargetingCRMs.rds")

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
  slice_head(n = 10)

write.csv(FEA_top_abundant, "results/tables/FEA_top_abundant.csv", row.names = FALSE)

# Manhattan-like-plot
DEGs_FEAplot1 <- gostplot(FEA, capped = TRUE, interactive = FALSE)

DEGs_FEAplot2 <- publish_gostplot(DEGs_FEAplot1, 
                                  highlight_terms = FEA_top_abundant$term_id, 
                                  width = NA, height = NA, filename = NULL )

ggsave(filename = "results/figures/DEGs_FEAplot.png", plot = DEGs_FEAplot2,
       width = 10, height = 7, dpi = 300, bg = "white")

# 7.2. Identification of enriched terms containing self-targeting CRM genes ----

genesCRM <- selfTargetingCRMs$common_genes

FEA_DEGs$CRM_genes_in_term <- NA

for (i in seq_len(nrow(FEA_DEGs))) {
  
  genes_intersection <- unlist(strsplit(FEA_DEGs$intersection[i], ","))
  genes_intersection <- trimws(genes_intersection)
  
  genes_overlap <- genes_intersection[genes_intersection %in% genesCRM]
  
  if (length(genes_overlap) > 0) {
    FEA_DEGs$CRM_genes_in_term[i] <- paste(genes_overlap, collapse = ", ")
  }
}

FEA_selfTargetingCRMs <- FEA_DEGs[!is.na(FEA_DEGs$CRM_genes_in_term), ]

write.csv(FEA_selfTargetingCRMs, "results/tables/FEA_selfTargetingCRMs.csv", row.names = FALSE)

# 7.3. Top 5 enriched terms associated with self-targeting CRM genes ----

FEA_selfTargetingCRMs_top <- FEA_selfTargetingCRMs %>%
  mutate(
    CRM_gene_count = lengths(strsplit(CRM_genes_in_term, ",\\s*")),
    minus_log10_p = -log10(p_value),
    term_label = paste0(term_name, " [", source, "]",
                        "\nCRM genes: ", CRM_genes_in_term)
  ) %>%
  group_by(query, source) %>%
  arrange(p_value, desc(CRM_gene_count), desc(intersection_size), .by_group = TRUE) %>%
  slice_head(n = 5) %>%
  ungroup()

FEA_selfTargetingCRMs_dotplot <- ggplot(
  FEA_selfTargetingCRMs_top,
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

ggsave("results/figures/FEA_selfTargetingCRMs_dotplot.png", FEA_selfTargetingCRMs_dotplot,
       width = 13, height = 9, dpi = 300, bg = "white")
