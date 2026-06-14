# 5. Identification of CRMs affected by mutations in DEGs ----
#
# Candidate cis-regulatory modules (CRMs) overlapping intronic mutations in
# differentially expressed genes (DEGs) from TNBC samples are identified using
# genomic coordinates from the mutation dataset.
#
# Identified CRMs are annotated with regulatory information from BioGateway and
# filtered to retain those active in breast tissue. Breast-associated CRMs are
# further evaluated according to their disease-related phenotypes, predicted
# target genes, and genomic overlap with annotated genes.
#
# CRMs are prioritized when the mutated DEG, the CRM-located gene, and the
# predicted CRM target gene correspond to the same gene, suggesting a potential
# self-regulatory effect of intronic mutations on DEG expression.
#
# Transcription factors associated with the prioritized CRMs are retrieved and
# filtered for breast tissue relevance.

library(RBioGateway)
library(dplyr)
library(GenomicRanges)
library(ensembldb)
library(EnsDb.Hsapiens.v86)

maf_intronic_TNBC <- readRDS("data/processed/maf_intronic_TNBC.rds")

# 5.1. Identification of CRMs overlapping with mutation coordinates ----

maf_intronic_TNBC$Chromosome <- sub("^chr", "chr-", maf_intronic_TNBC$Chromosome)

overlapCRMs_list <- vector("list", nrow(maf_intronic_TNBC))

for (i in seq_len(nrow(maf_intronic_TNBC))) {
  
  crms <- as.data.frame(
    getCRMs_by_overlap(
      maf_intronic_TNBC$Chromosome[i],
      maf_intronic_TNBC$Start_Position[i],
      maf_intronic_TNBC$End_Position[i]
    )
  )
  
  # Keep only cases where CRMs were found
  if (nrow(crms) > 0) {
    
    crms$tumor_sample_barcode <- maf_intronic_TNBC$Tumor_Sample_Barcode[i]
    crms$mut_DEGs <- maf_intronic_TNBC$Hugo_Symbol[i]
    crms$mut_start <- maf_intronic_TNBC$Start_Position[i]
    crms$mut_end <- maf_intronic_TNBC$End_Position[i]
    crms$ref_allele <- maf_intronic_TNBC$Reference_Allele[i]
    crms$mut_allele <- maf_intronic_TNBC$Tumor_Seq_Allele2[i]
    crms$variant_type <- maf_intronic_TNBC$Variant_Type[i]
    crms$expression_status <- maf_intronic_TNBC$expression_status[i]
    
    overlapCRMs_list[[i]] <- crms
  }
}

overlapCRMs <- bind_rows(overlapCRMs_list)

cat("CRMs overlapping intronic mutations:", nrow(overlapCRMs))

# 5.2. Selection of CRMs active in breast tissue ----

# Supplementary information of identified CRMs
CRMs_info_list <- vector("list", nrow(overlapCRMs))
CRMs_add_info_list <- vector("list", nrow(overlapCRMs))

for (i in seq_len(nrow(overlapCRMs))) {
  
  crm <- overlapCRMs$crm_name[i]
  
  CRMs_info_list[[i]] <- as.data.frame(
    lapply(getCRM_info(crm), paste, collapse = "; "))
  
  CRMs_add_info_list[[i]] <- as.data.frame(
    lapply(getCRM_add_info(crm), paste, collapse = "; "))
}

CRMs_info <- bind_rows(CRMs_info_list)
CRMs_add_info <- bind_rows(CRMs_add_info_list)

allCRMs <- cbind(overlapCRMs, CRMs_info[, -c(1, 2)], CRMs_add_info)

# CRMs active in breast tissue
breastCRMs <- allCRMs[grepl("UBERON_0000310", allCRMs$biological_samples),]

cat("CRMs active in breast tissue:", nrow(breastCRMs))

# 5.3. Verification of CMRs with BRCA-associated phenotype ----

breastCRMs$phenotype <- NA_character_

for (i in seq_len(nrow(breastCRMs))) {
  
  phen <- as.data.frame(crm2phen(breastCRMs$crm_name[i]))$phenotype
  
  if (length(phen) > 0) {
    breastCRMs$phenotype[i] <- paste(phen, collapse = "; ")
  }
}

print(breastCRMs$phenotype)

saveRDS(breastCRMs, file = "data/raw/breastCRMs.rds")

# 5.4. Identification of target genes of identified CRMs ----

breastCRMs$target_genes <- NA_character_

for (i in seq_len(nrow(breastCRMs))) {
  
  genes <- as.data.frame(crm2gene(breastCRMs$crm_name[i]))$gene_name
  
  if (length(genes) > 0) {breastCRMs$target_genes[i] <- paste(genes, collapse = "; ")}
}

cat("CRMs with identified target genes:", sum(!is.na(breastCRMs$target_genes)))
cat("CRMs without identified target genes:", sum(is.na(breastCRMs$target_genes)))

filtCRMs <- breastCRMs %>% dplyr::filter(!is.na(target_genes))

# 5.5. Identification of located genes of identified CRMs ----

# Ensembl chromosome format
filtCRMs$chr <- sub("^chr-?", "", filtCRMs$chr)

# Convert CRMs to GRanges
crm_gr <- makeGRangesFromDataFrame(
  filtCRMs,
  seqnames.field = "chr",
  start.field = "start",
  end.field = "end",
)

# Get annotated genes
gene_gr <- genes(EnsDb.Hsapiens.v86, columns = c("gene_id", "gene_name"))

# Find CRMs-genes overlaps
gene_hits <- findOverlaps(crm_gr, gene_gr, ignore.strand = TRUE)

# Use gene_name and if missing, use gene_id
gene_labels <- mcols(gene_gr)$gene_name
gene_labels[is.na(gene_labels) | gene_labels == ""] <- 
  mcols(gene_gr)$gene_id[is.na(gene_labels) | gene_labels == ""]

# Collapse genes per CRM
located_gene <- rep(NA_character_, length(crm_gr))

for (i in seq_len(length(crm_gr))) {
  
  gene_index <- subjectHits(gene_hits)[queryHits(gene_hits) == i]
  
  if (length(gene_index) > 0) {
    
    genes <- gene_labels[gene_index]
    
    located_gene[i] <- paste(sort(genes), collapse = "; ")
  }
}

filtCRMs$located_gene <- located_gene

# 5.6. Identification of self-targeting CRMs in mutated DEGs ----

filtCRMs$common_genes <- NA_character_

for (i in seq_len(nrow(filtCRMs))) {
  
  mut_DEGs <- trimws(unlist(strsplit(filtCRMs$mut_DEGs[i], ";")))
  target_genes <- trimws(unlist(strsplit(filtCRMs$target_genes[i], ";")))
  located_genes <- trimws(unlist(strsplit(filtCRMs$located_gene[i], ";")))
  
  common_genes <- Reduce(intersect, list(mut_DEGs, target_genes, located_genes))
  
  if (length(common_genes) > 0) {
    filtCRMs$common_genes[i] <- paste(common_genes, collapse = "; ")
  }
}

# Keep only CRMs where target gene, located gene and mutated gene match
selfTargetingCRMs <- filtCRMs[!is.na(filtCRMs$common_genes), ]

cat("CRMs located in and targeting the same mutated DEG:", nrow(selfTargetingCRMs))
cat("Mutated DEGs potentially regulated by these CRMs:", 
    paste(unique(selfTargetingCRMs$common_genes), collapse = "; "))

write.csv(selfTargetingCRMs, "results/tables/selfTargetingCRMs.csv", row.names = FALSE)
saveRDS(selfTargetingCRMs, "data/processed/selfTargetingCRMs.rds")

# 5.7. TFs associated with identified CRMs ----

CRMs_tfac_list <- vector("list", nrow(selfTargetingCRMs))

for (i in seq_len(nrow(selfTargetingCRMs))) {
  
  crm <- selfTargetingCRMs$crm_name[i]
  
  tfac <- crm2tfac(crm)
  
  # Skip CRMs without TFs
  if (is.null(tfac) || length(tfac) == 0) { next }
  
  tfac <- as.data.frame(tfac)
  
  # Keep only TFs associated with breast tissue
  if (!"biological_samples" %in% colnames(tfac)) { next }
  
  tfac <- tfac[grepl("UBERON_0000310", tfac$biological_samples), ]
  
  tfac$crm_name <- crm
  
  CRMs_tfac_list[[i]] <- tfac
}

CRMs_tf <- bind_rows(CRMs_tfac_list)

CRMs_tf <- CRMs_tf[, c("crm_name", setdiff(names(CRMs_tf), "crm_name"))]

cat("CRM-TF associations in breast tissue:", nrow(CRMs_tf))
cat("CRMs with breast-associated TFs:", length(unique(CRMs_tf$crm_name)))

write.csv(CRMs_tf, "results/tables/CRMs_tf.csv", row.names = FALSE)
saveRDS(CRMs_tf, "data/processed/CRMs_tf.rds")
