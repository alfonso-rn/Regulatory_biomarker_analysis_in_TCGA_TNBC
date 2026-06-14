# 6. Identification of TFBS factors affected by mutations in DEGs ----

library(dplyr)
library(jsonlite)

maf_intronic_TNBC <- readRDS("data/processed/maf_intronic_TNBC.rds")
CRMs_tf <- readRDS("data/processed/CRMs_tf.rds")

# 6.1. Identification of TFBS overlapping with mutation coordinates ----

mut <- maf_intronic_TNBC %>%
  transmute( 
    mutation_id = paste0(Hugo_Symbol, "_", Chromosome, ":", Start_Position, "-", End_Position, "_", 
                         Reference_Allele, ">", Tumor_Seq_Allele2, "_", Tumor_Sample_Barcode),
    Tumor_Sample_Barcode,
    gene = Hugo_Symbol,
    chrom = Chromosome,
    start0 = Start_Position - 1, # convert start:MAF 1-based to UCSC 0-based
    start = Start_Position,
    end = End_Position,
    ref = Reference_Allele,
    alt = Tumor_Seq_Allele2,
    variant_type = Variant_Type,
    expression_status
  ) %>%
  mutate(ref = ifelse(ref == "-", "", ref), alt = ifelse(alt == "-", "", alt))

query_jaspar <- function(x) {
  
  url <- paste0(
    "https://api.genome.ucsc.edu/getData/track?genome=hg38;track=jaspar2026",
    ";chrom=", x$chrom,
    ";start=", x$start0,
    ";end=", x$end
  )
  
  tfbs <- fromJSON(url)$jaspar2026
  
  if (is.null(tfbs) || length(tfbs) == 0) return(NULL)
  
  as_tibble(tfbs) %>%
    mutate(
      mutation_id = x$mutation_id,
      Tumor_Sample_Barcode = x$Tumor_Sample_Barcode,
      gene = x$gene,
      mut_chrom = x$chrom,
      start0 = x$start0,
      mut_start = x$start,
      mut_end = x$end,
      ref = x$ref,
      alt = x$alt,
      variant_type = x$variant_type,
      expression_status = x$expression_status,
      .before = 1
    )
}

mut_overlapTFBS <- list()

for (i in seq_len(nrow(mut))) {
  mut_overlapTFBS[[i]] <- query_jaspar(mut[i, ])
}

mut_overlapTFBS <- bind_rows(mut_overlapTFBS) 

cat("Total of TFBS overlapping mutations:", nrow(mut_overlapTFBS))

overlapTFBS <- mut_overlapTFBS %>% dplyr::filter(mut_overlapTFBS$score >= 300)

cat("TFBS overlapping mutations with score ≥ 300:", nrow(overlapTFBS), "\n",
    "Mutations with at least one overlapping TFBS:",n_distinct(overlapTFBS$mutation_id), "\n",
    "Unique TFs identified:", n_distinct(overlapTFBS$TFName))

head(overlapTFBS)
