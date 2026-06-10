# 5. Identification of CRMs affected by mutations in DEGs ----

library(RBioGateway)
library(dplyr)

maf_intronic_TNBC <- readRDS("data/processed/maf_intronic_TNBC.rds")

# 5.1. Identification of CRMs overlapping with the mutation position ----

maf_intronic_TNBC$Chromosome <- sub("^chr", "chr-", maf_intronic_TNBC$Chromosome)

maf_intronic_TNBC$CRM_overlap <- mapply(
  getCRMs_by_overlap,
  maf_intronic_TNBC$Chromosome,
  maf_intronic_TNBC$Start_Position,
  maf_intronic_TNBC$End_Position,
  SIMPLIFY = FALSE
)

overlapCRMs <- bind_rows(
  lapply(maf_intronic_TNBC$CRM_overlap, as.data.frame)) %>%
  filter(!is.na(crm_name))

cat("CRMs overlapping intronic mutations:", nrow(overlapCRMs))

# 5.2. Verification of CMRs with BRCA-associated phenotype ----

overlapCRMs$phenotype <- mapply(
  crm2phen,
  overlapCRMs$crm_name,
  SIMPLIFY = FALSE
)

print(overlapCRMs$phenotype)

# 5.3. Supplementary information of identified CRMs ----

CRMs_info <- bind_rows(
  lapply(overlapCRMs$crm_name, function(crm) {
    as.data.frame(lapply(
      getCRM_info(crm), paste, collapse = "; "))}
  ))

CRMs_add_info <- bind_rows(
  lapply(overlapCRMs$crm_name, function(crm) {
    as.data.frame(lapply(
      getCRM_add_info(crm), paste, collapse = "; "))}
    ))

allCRMs <- cbind(overlapCRMs, CRMs_info[, -c(1, 2)], CRMs_add_info)

write.csv(allCRMs, file = "results/tables/allCRMs.csv", row.names = FALSE)

# 5.4. Prioritize active CRMs in breast tissue ----

breastCRMs <- allCRMs[grepl("UBERON_0000310", allCRMs$biological_samples),]

cat("CRMs active in breast tissue:", nrow(breastCRMs))

# 5.5. Classification of CRMs based on genomic location: exon vs intron ----

download.file(
  "https://ftp.ensembl.org/pub/release-116/gtf/homo_sapiens/Homo_sapiens.GRCh38.116.chr.gtf.gz",
  "Homo_sapiens.GRCh38.116.chr.gtf.gz")

file.rename("Homo_sapiens.GRCh38.116.chr.gtf.gz", 
            "data/raw/Homo_sapiens.GRCh38.116.chr.gtf.gz")

gtf <- read.delim(
  gzfile("data/raw/Homo_sapiens.GRCh38.116.chr.gtf.gz"),
  sep = "\t",
  header = FALSE,
  comment.char = "#",
  quote = "",
  col.names = c("chr", "source", "feature", "start", "end", "score", "strand", 
                "frame", "attribute")
)
