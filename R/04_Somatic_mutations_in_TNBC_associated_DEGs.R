# 4. Somatic mutations in differentially expressed genes ----
#
# Masked somatic mutation data from TCGA-BRCA are downloaded and prepared from
# the GDC Simple Nucleotide Variation data category.
#
# Differentially expressed genes previously identified in the TNBC vs non-TNBC
# analysis are used to select mutation records affecting DEG-associated genes,
# based on Ensembl gene identifiers in ENSG format.
#
# The MAF table is then filtered to retain intronic somatic mutations located
# in DEGs.

library(TCGAbiolinks)
library(dplyr)
library(maftools)

samplesTNBC <- readRDS("data/processed/samplesTNBC.rds")
samplesNonTNBC <- readRDS("data/processed/samplesNonTNBC.rds")
dataDEGsLevel <- readRDS("data/processed/dataDEGsLevel.rds")

# 4.1. Query, download and prepare TCGA-BRCA somatic mutation data ----

query_mut <- GDCquery(
  project = "TCGA-BRCA",
  data.category = "Simple Nucleotide Variation",
  data.type = "Masked Somatic Mutation",
  workflow.type = "Aliquot Ensemble Somatic Variant Merging and Masking",
  access = "open"
)

GDCdownload(
  query_mut,
  method = "api",
  files.per.chunk = 100
)

TCGA_BRCA_mut <- GDCprepare(query_mut)

saveRDS(TCGA_BRCA_mut, "data/raw/TCGA_BRCA_mut.rds")

cat("Total somatic mutations downloaded from TCGA-BRCA:", nrow(TCGA_BRCA_mut))

# 4.2. Comparison of TNBC vs Non-TNBC mutational profile ----

# Extract TCGA sample-level barcodes for TNBC samples
tnbc_ids <- substr(samplesTNBC, 1, 16)

maf_TNBC <- TCGA_BRCA_mut %>%
  filter(substr(Tumor_Sample_Barcode, 1, 16) %in% tnbc_ids)

cat("Somatic mutations identified in TNBC samples:", nrow(maf_TNBC))

saveRDS(maf_TNBC, "data/processed/maf_TNBC.rds")

maf_TNBC_ids <- substr(maf_TNBC$Tumor_Sample_Barcode, 1, 16)

# Extract TCGA sample-level barcodes for Non-TNBC samples
nontnbc_ids <- substr(samplesNonTNBC, 1, 16)

maf_NonTNBC <- TCGA_BRCA_mut %>%
  filter(substr(Tumor_Sample_Barcode, 1, 16) %in% nontnbc_ids)

cat("Somatic mutations identified in Non-TNBC samples:", nrow(maf_NonTNBC))

saveRDS(maf_NonTNBC, "data/processed/maf_NonTNBC.rds")

maf_NonTNBC_ids <- substr(maf_NonTNBC$Tumor_Sample_Barcode, 1, 16)

# Summarize the number of samples with and without mutation data in the MAF
maf_sampleTNBCvsNonTNBC_coverage <- data.frame(
  Group = c("TNBC", "Non-TNBC"),
  
  Total_samples = c(length(tnbc_ids), length(nontnbc_ids)),
  
  Samples_with_mutations_in_MAF = c(
    sum(tnbc_ids %in% maf_TNBC_ids),
    sum(nontnbc_ids %in% maf_NonTNBC_ids)),
  
  Samples_absent_from_MAF = c(
    sum(!(tnbc_ids %in% maf_TNBC_ids)),
    sum(!(nontnbc_ids %in% maf_NonTNBC_ids)))
)

write.csv(maf_sampleTNBCvsNonTNBC_coverage, 
          "results/tables/maf_sampleTNBCvsNonTNBC_coverage.csv", 
          row.names = FALSE)

# 4.3. Visual representation of MAFs TNBC vs Non-TNBC distribution ----

MAF_TNBC <- read.maf(maf_TNBC, vc_nonSyn = unique(TCGA_BRCA_mut$Variant_Classification))
MAF_NonTNBC <- read.maf(maf_NonTNBC, vc_nonSyn = unique(TCGA_BRCA_mut$Variant_Classification))

pdf(file = "results/figures/MAFoncoplots.pdf", width = 12, height = 8)

oncoplot(
  maf = MAF_TNBC, 
  titleText = "TNBC",
  top = 20,
  logColBar = TRUE,
  removeNonMutated = TRUE,
  draw_titv = TRUE, 
  titv_col = c("C>T" = "#7F0000", "C>G" = "#B30000", "C>A" = "#E34A33", 
               "T>A" = "#FC8D59", "T>C" = "#FDBB84", "T>G" = "#FDD49E")
)

oncoplot(
  maf = MAF_NonTNBC, 
  titleText = "NonTNBC",
  top = 20,
  logColBar = TRUE,
  removeNonMutated = TRUE,
  draw_titv = TRUE, 
  titv_col = c("C>T" = "#7F0000", "C>G" = "#B30000", "C>A" = "#E34A33", 
               "T>A" = "#FC8D59", "T>C" = "#FDBB84", "T>G" = "#FDD49E")
)

dev.off()

# 4.4. Somatic mutations in TNBC-associated DEGs ----

DEGs_ID <- rownames(dataDEGsLevel)

maf_TNBC_DEGs <- maf_TNBC[maf_TNBC$Gene %in% DEGs_ID, ]

cat("Mutations in DEGs from TNBC patient:", nrow(maf_TNBC_DEGs))

nonMutDEGs <- dataDEGsLevel[!DEGs_ID %in% maf_TNBC_DEGs$Gene, ]
saveRDS(nonMutDEGs, "data/processed/nonMutDEGs.rds")

MutDEGs <- dataDEGsLevel[DEGs_ID %in% maf_TNBC_DEGs$Gene, ]

cat("DEGs from TNBC samples with any mutations:", nrow(MutDEGs), "\n",
    "DEGs from TNBC samples without any mutations:", nrow(nonMutDEGs))

# 4.5. Selection of non-truncating somatic mutations in TNBC-associated DEGs ----

maf_non_trunc <- maf_TNBC_DEGs %>% 
  filter(
    # Remove truncating variants
    !(Variant_Classification %in% c(
      "Frame_Shift_Del", 
      "Frame_Shift_Ins", 
      "Nonsense_Mutation",
      "Splice_Site", 
      "Nonstop_Mutation")
    
    # Remove truncating consequences
    | grepl(
      "frameshift_variant|stop_gained|splice_acceptor_variant|splice_donor_variant",
      Consequence)))

cat("Non-truncating mutations in DEGs from TNBC patient:", nrow(maf_non_trunc))

# 4.6. Categorizing exonic vs intronic mutations in TNBC-associated DEGs ----

# Intronic mutations
maf_intronic_TNBC <- maf_non_trunc[!is.na(maf_non_trunc$INTRON), ]

cat("Identification of intronic somatic mutations in TNBC-associated DEGs:", 
    nrow(maf_intronic_TNBC))

# Add expression TNBC status level column from DEA
maf_intronic_TNBC$expression_status <- dataDEGsLevel$expression_status[
  match(maf_intronic_TNBC$Gene, rownames(dataDEGsLevel))]

saveRDS(maf_intronic_TNBC, "data/processed/maf_intronic_TNBC.rds")

# Exonic mutations
maf_exonic_TNBC <- maf_non_trunc[!is.na(maf_non_trunc$EXON), ]

cat("Identification of exonic somatic mutations in TNBC-associated DEGs:", 
    nrow(maf_exonic_TNBC))

# Add expression TNBC status level column from DEA
maf_exonic_TNBC$expression_status <- dataDEGsLevel$expression_status[
  match(maf_exonic_TNBC$Gene, rownames(dataDEGsLevel))]

saveRDS(maf_exonic_TNBC, "data/processed/maf_exonic_TNBC.rds")
