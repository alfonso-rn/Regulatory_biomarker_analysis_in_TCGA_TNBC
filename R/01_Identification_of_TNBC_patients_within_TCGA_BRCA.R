# 1. Identification of TNBC patients within TCGA-BRCA ----
#
# TCGA-BRCA clinical data are downloaded and used to identify patients with
# clinically defined triple-negative breast cancer (TNBC).
#
# TNBC patients are selected using clinical immunohistochemistry (IHC) and
# HER2-FISH information available in the TCGA-BRCA BCR Biotab clinical tables.
#
# Patients are first required to have negative estrogen receptor (ER) and
# progesterone receptor (PR) status by IHC.
#
# HER2-positive cases by FISH are then excluded. The remaining cases are kept
# as clinically defined TNBC when HER2 evidence is compatible with a negative
# or non-positive HER2 status, based on HER2-FISH status, HER2 status by IHC
# and HER2 IHC score.
#
# The final output consists of the clinical table for the selected TNBC patients
# and their corresponding TCGA patient barcodes.

library(TCGAbiolinks)
library(dplyr)

# 1.1. Query, download and prepare TCGA-BRCA clinical data ----

query_clin <- GDCquery(
  project = "TCGA-BRCA",
  data.category = "Clinical",
  data.type = "Clinical Supplement",
  data.format = "BCR Biotab"
)

GDCdownload(query_clin)

TCGA_BRCA_clinic <- GDCprepare(query_clin)

saveRDS(TCGA_BRCA_clinic, "data/raw/TCGA_BRCA_clinic.rds")

# 1.2. Extract patient-level clinical information ----

# Display the available clinical tables contained in the downloaded object.
names(TCGA_BRCA_clinic)

# Extract the patient-level clinical table.
clinical_patient <- TCGA_BRCA_clinic[["clinical_patient_brca"]][-c(1, 2), ]

cat("Total number of clinical BRCA patient:", nrow(clinical_patient))

# Display the available clinical variables.
names(clinical_patient)

# Results of the columns used to identify clinically TNBC patients.
sapply(clinical_patient[,c( 
  "er_status_by_ihc",
  "pr_status_by_ihc",
  "her2_status_by_ihc",
  "her2_fish_status",
  "her2_ihc_score")],
table,
useNA = "ifany"
)

# 1.3. Identify clinically defined TNBC patients ----

# Keep ER & PR negative cases
tnbc <- clinical_patient %>%
  filter(
    er_status_by_ihc == "Negative",
    pr_status_by_ihc == "Negative",
  )

cat("N° patients at the first level of selection:", nrow(tnbc))

# Exclude HER2-FISH positive cases.
tnbc <- tnbc %>%
  filter(
    her2_fish_status %in% c(
      "Negative", 
      "Indeterminate", 
      "Equivocal", 
      "[Not Evaluated]", 
      "[Not Available]"
    ) 
    | is.na(her2_fish_status),
  )

cat("N° patients at the second level of selection:", nrow(tnbc))

# Keep remaining cases if one of the two criteria is met:
tnbc <- tnbc %>%
  filter(
    # Criterion 1: HER2-FISH negative & HER2-IHC not positive.
    ( her2_fish_status == "Negative" 
      &
        ( her2_status_by_ihc %in% c(
          "Negative", 
          "Indeterminate", 
          "Equivocal", 
          "[Not Evaluated]", 
          "[Not Available]"
        ) 
        | is.na(her2_status_by_ihc))) 
    |
      # Criterion 2: HER2-IHC negative & IHC score 0, 1+ or not available.
      (her2_status_by_ihc == "Negative" 
       &
         (her2_ihc_score %in% c("[Not Available]", "0", "1+", "") 
          | is.na(her2_ihc_score))) 
  )

cat("Number of clinically defined TNBC patients identified:", nrow(tnbc))

saveRDS(tnbc, "data/processed/tnbc_clinical_table.rds")

# 1.4. Extract TCGA-BRCA barcodes for the identified TNBC patients.
tnbc_barcodes <- tnbc$bcr_patient_barcode

saveRDS(tnbc_barcodes, "data/processed/tnbc_barcodes.rds")
