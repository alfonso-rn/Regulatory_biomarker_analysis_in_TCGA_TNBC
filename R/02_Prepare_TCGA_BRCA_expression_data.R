# 2. Preparation of TCGA-BRCA expression data for TNBC vs non-TNBC analysis ----
#
# TCGA-BRCA RNA-seq gene expression data are downloaded from primary tumor
# samples using open-access STAR-Counts data.
#
# The expression matrix is preprocessed to remove poorly correlated samples,
# normalized for GC-content bias, and filtered to exclude lowly expressed genes.
#
# Tumor samples are then classified into TNBC and non-TNBC groups according to
# the previously identified TNBC clinical barcodes.

library(TCGAbiolinks)
library(SummarizedExperiment)

# 2.1. Query, Download & Prepare TCGA-BRCA RNA-seq expression data ----

query_RNAseq <- GDCquery(
  project = "TCGA-BRCA",
  data.category = "Transcriptome Profiling",
  data.type = "Gene Expression Quantification",
  workflow.type= "STAR - Counts",
  experimental.strategy = "RNA-Seq",
  access = "open",
  sample.type = "Primary Tumor"
)

GDCdownload(
  query_RNAseq, 
  method = "api",
  files.per.chunk = 100, 
)

TCGA_BRCA_primarytumors <- GDCprepare(query_RNAseq)

saveRDS(TCGA_BRCA_primarytumors, "data/raw/TCGA_BRCA_primarytumors.rds")

# 2.2. Preprocess RNA-seq expression data ----

dataPrep <- TCGAanalyze_Preprocessing(
  object = TCGA_BRCA_primarytumors,
  cor.cut = 0.6
)

saveRDS(dataPrep, "data/processed/dataPrep.rds")

# 2.3. Normalize gene expression data ----

dataNorm <- TCGAanalyze_Normalization(
  tabDF = dataPrep,
  geneInfo = geneInfoHT,
  method = "gcContent"
)

saveRDS(dataNorm, "data/processed/dataNorm.rds")

# 2.4. Quantil filter of expressed genes ----

dataFilt <- TCGAanalyze_Filtering(
  tabDF = dataNorm,
  method = "quantile",
  qnt.cut = 0.25
)

saveRDS(dataFilt, "data/processed/dataFilt.rds")

print(data.frame(
  Matrix = c("Preprocessed", "Normalized", "Filtered"),
  Samples = c(ncol(dataPrep), ncol(dataNorm), ncol(dataFilt)),
  Genes = c(nrow(dataPrep), nrow(dataNorm), nrow(dataFilt))))

# 2.5. Selection of sample groups: TNBC vs non-TNBC ----

patient_barcodes <- substr(colnames(dataFilt), 1, 12)

tnbc_barcodes <- readRDS("data/processed/tnbc_barcodes.rds")

samplesTNBC <- colnames(dataFilt)[patient_barcodes %in% tnbc_barcodes]
samplesNonTNBC <- colnames(dataFilt)[!(patient_barcodes %in% tnbc_barcodes)]

cat("TNBC samples:", length(samplesTNBC))
cat("Non-TNBC samples:", length(samplesNonTNBC))

saveRDS(samplesTNBC, "data/processed/samplesTNBC.rds")
saveRDS(samplesNonTNBC, "data/processed/samplesNonTNBC.rds")
