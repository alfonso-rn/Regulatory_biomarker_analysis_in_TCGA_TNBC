# 7. GO annotation and enrichment analysis of CRM-associated genes ----
#
# Functional annotation is performed for genes associated with prioritized
# self-targeting CRMs potentially affected by intronic mutations in TNBC.
#
# CRM-associated genes are mapped to their encoded human protein products using
# RBioGateway. Protein-level annotations are then retrieved for the three Gene
# Ontology domains: cellular component, molecular function and biological process.
# 
# In parallel, CRM-associated genes are evaluated through TCGAbiolinks to identify
# enriched GO terms and pathways, providing a functional overview of the genes
# whose regulatory regions may be altered by mutations.

library(RBioGateway)
library(TCGAbiolinks)

selfTargetingCRMs <- readRDS("data/processed/selfTargetingCRMs.rds")

# 7.1. Identification of protein products encoded by CRM-associated genes ----

crm_genes <- unique(unlist(strsplit(selfTargetingCRMs$common_genes, ";\\s*")))

cat("Genes associated with selected CRMs:", length(crm_genes), "\n")
cat("Genes:", paste(crm_genes, collapse = "; "), "\n")

gene2prot_list <- list()

for (gene in crm_genes) {
  
  prots <- gene2prot(gene, taxon = "Homo sapiens")
  
  if (length(prots) > 0) {
    
    # Keep only human proteins
    prots_human <- prots[grepl("_HUMAN$", prots)]
    
    if (length(prots_human) > 0) {
      gene2prot_list[[gene]] <- data.frame(
        gene = gene,
        protein = prots_human )}
  }
}

gene2prot_human <- bind_rows(gene2prot_list)

# 7.2. Retrieval of cellular compartment annotations from human proteins ----

cc_list <- list()

for (i in seq_len(nrow(gene2prot_human))) {
  
  gene <- gene2prot_human$gene[i]
  protein <- gene2prot_human$protein[i]
  
  cc <- prot2cc(protein)
  
  # Skip proteins without available annotations
  if (is.null(cc) || !is.data.frame(cc) || nrow(cc) == 0) { next }
  
  cc$gene <- gene
  cc$protein <- protein
  
  cc_list[[i]] <- cc
}

crm_cc <- bind_rows(cc_list)

cat("Cellular compartment annotations retrieved:", nrow(crm_cc), "\n",
    "Proteins with cellular compartment annotations:", length(unique(crm_cc$protein)), "\n",
    "Genes with cellular compartment annotations:", length(unique(crm_cc$gene)), "\n")

# 7.3. Retrieval of molecular function annotations from human proteins ----

mf_list <- list()

for (i in seq_len(nrow(gene2prot_human))) {
  
  gene <- gene2prot_human$gene[i]
  protein <- gene2prot_human$protein[i]
  
  mf <- prot2mf(protein)
  
  # Skip proteins without available annotations
  if (is.null(mf) || !is.data.frame(mf) || nrow(mf) == 0) { next }
  
  mf$gene <- gene
  mf$protein <- protein
  
  mf_list[[i]] <- mf
}

crm_mf <- bind_rows(mf_list)

cat("Molecular function annotations retrieved:", nrow(crm_mf), "\n",
    "Proteins with molecular function annotations:", length(unique(crm_mf$protein)), "\n",
    "Genes with molecular function annotations:", length(unique(crm_mf$gene)))

# 7.4. Retrieval of biological process annotations from human proteins ----

bp_list <- list()

for (i in seq_len(nrow(gene2prot_human))) {
  
  gene <- gene2prot_human$gene[i]
  protein <- gene2prot_human$protein[i]
  
  bp <- prot2bp(protein)
  
  # Skip proteins without available annotations
  if (is.null(bp) || !is.data.frame(bp) || nrow(bp) == 0) { next }
  
  bp$gene <- gene
  bp$protein <- protein
  
  bp_list[[i]] <- bp
}

crm_bp <- bind_rows(bp_list)

cat("Biological process annotations retrieved:", nrow(crm_bp), "\n",
    "Proteins with biological process annotations:", length(unique(crm_bp$protein)), "\n",
    "Genes with biological process annotations:", length(unique(crm_bp$gene)))

write.csv(crm_cc, "results/tables/crm_GO_CC_RBioGateway.csv", row.names = FALSE)
write.csv(crm_mf, "results/tables/crm_GO_MF_RBioGateway.csv", row.names = FALSE)
write.csv(crm_bp, "results/tables/crm_GO_BP_RBioGateway.csv", row.names = FALSE)

# 7.5. Gene Ontology (GO) and Pathway enrichment in TCGAbiolinks

ansEA <- TCGAanalyze_EAcomplete(
  TFname = "CRM-associated DEGs in TNBC vs non-TNBC",
  RegulonList = crm_genes
)

TCGAvisualize_EAbarplot(
  tf = rownames(ansEA$ResBP), 
  GOBPTab = ansEA$ResBP,
  GOCCTab = ansEA$ResCC,
  GOMFTab = ansEA$ResMF,
  PathTab = ansEA$ResPat,
  nRGTab = crm_genes, 
  nBar = 10,
  filename = "results/figures/GO_EAbarplot.pdf"
)
