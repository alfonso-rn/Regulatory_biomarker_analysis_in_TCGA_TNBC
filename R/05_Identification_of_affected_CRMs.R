# 5. Identification of CRMs affected by mutations in DEGs ----

library(RBioGateway)
library(dplyr)
library(GenomicRanges)
library(ensembldb)
library(EnsDb.Hsapiens.v86)

maf_intronic_TNBC <- readRDS("data/processed/maf_intronic_TNBC.rds")

# 5.1. Identification of CRMs overlapping with mutation coordintes ----

maf_intronic_TNBC$Chromosome <- sub("^chr", "chr-", maf_intronic_TNBC$Chromosome)

overlapCRMs_list <- vector("list", nrow(maf_intronic_TNBC))

for (i in seq_len(nrow(maf_intronic_TNBC))) {
 
   overlapCRMs_list[[i]] <- as.data.frame(
     
     getCRMs_by_overlap(
      maf_intronic_TNBC$Chromosome[i],
      maf_intronic_TNBC$Start_Position[i],
      maf_intronic_TNBC$End_Position[i] 
    )
  )
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

write.csv(breastCRMs, file = "results/tables/breastCRMs.csv", row.names = FALSE)
saveRDS(breastCRMs, file = "data/raw/breastCRMs.rds")

# 5.3. Verification of CMRs with BRCA-associated phenotype ----

breastCRMs$phenotype <- NA_character_

for (i in seq_len(nrow(breastCRMs))) {
  
  phen <- as.data.frame(crm2phen(breastCRMs$crm_name[i]))$phenotype
  
  if (length(phen) > 0) {
    breastCRMs$phenotype[i] <- paste(phen, collapse = "; ")
  }
}

print(breastCRMs$phenotype)

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

# 5.6. Identification of CRMs active on the host gene ----
filtCRMs$target_located_match <- FALSE

for (i in seq_len(nrow(filtCRMs))) {
  
  target_genes <- trimws(unlist(strsplit(filtCRMs$target_genes[i], ";")))
  located_genes <- trimws(unlist(strsplit(filtCRMs$located_gene[i], ";")))
  
  if (any(target_genes %in% located_genes)) {
    filtCRMs$target_located_match[i] <- TRUE
  }
}

# Keep only CRMs where target_genes and located_gene match
hostgeneCRMs <- filtCRMs[filtCRMs$target_located_match == TRUE, ]

cat("CRMs with activity on the host gene:", nrow(hostgeneCRMs))

# Save results
write.csv(hostgeneCRMs, "results/tables/hostgeneCRMs.csv", row.names = FALSE)

# 5.7. TFs associated with identified CRMs ----

CRMs_tfac_list <- vector("list", nrow(hostgeneCRMs))

for (i in seq_len(nrow(hostgeneCRMs))) {
  
  crm <- hostgeneCRMs$crm_name[i]
  
  tfac <- crm2tfac(crm)
  
  if (is.list(tfac) && length(tfac) > 0) {
    CRMs_tfac_list[[i]] <- as.data.frame(
      lapply(tfac, paste, collapse = "; ")
    )
    
    CRMs_tfac_list[[i]]$crm_name <- crm
    
  } else {
    
    CRMs_tfac_list[[i]] <- data.frame(
      crm_name = crm,
      tfac_name = NA_character_
    )
  }
}

CRMs_tfac <- bind_rows(CRMs_tfac_list)

cat("CRMs with associated TFs:", sum(!is.na(CRMs_tfac$tfac_name)))
