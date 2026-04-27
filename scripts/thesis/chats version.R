# ===
# IMPORTING DATA
# ===
library(readr)
library(readxl)
library(tidyr)
library(dplyr)
library(phyloseq)
library(stringr)

ASV_raw <- read_tsv('/Users/chloecamacho/Desktop/thesis/Bioinformatics/raw data/library_g3_trim_prok_asvs.df.tsv')
TAXA <- read_tsv('/Users/chloecamacho/Desktop/thesis/Bioinformatics/raw data/library_g3_trim_prok_asvs.taxa.tsv')
meta <- read_excel("~/Desktop/thesis/Bioinformatics/meta?.xlsx")

# ===
# CLEANING DATA - CORRECTED ASV FORMAT
# ===
# CRITICAL: ASV_matrix must have ASVs as ROWS, samples as COLUMNS for phyloseq
ASV <- ASV_raw %>% 
  pivot_wider(
    id_cols    = asv_id,          # ROWS = ASVs ✓
    names_from = sample_name,     # COLUMNS = samples ✓
    values_from = counts,
    values_fill = 0
  ) %>%
  column_to_rownames("asv_id") %>%
  mutate(across(everything(), as.numeric)) %>%
  as.matrix()
ASV_matrix <- ASV
class(ASV_matrix)  # "matrix"

# Taxonomy matrix (ASVs as rownames)
TAXA_MAT <- TAXA %>%
  column_to_rownames("asv_id") %>%
  as.matrix()
class(TAXA_MAT)

# ===
# METADATA PROCESSING
# ===
meta <- meta %>% 
  distinct(CamachoName, .keep_all = TRUE) %>% 
  rename(SamplingDate = `Sampling Date`) %>%
  mutate(CamachoName = str_replace(CamachoName, "PSV", "PSU"))

# Make rownames from CamachoName
meta_df <- as.data.frame(meta)
rownames(meta_df) <- meta_df$CamachoName
meta_df$CamachoName <- NULL

# ===
# ALIGN SAMPLES - CRITICAL STEP
# ===
# Find common samples and align ORDER
common_samples <- intersect(colnames(ASV_matrix), rownames(meta_df))

# Subset both to matching samples, SAME ORDER
ASV_matrix <- ASV_matrix[, common_samples, drop = FALSE]
meta_df    <- meta_df[common_samples, , drop = FALSE]

# VERIFY perfect alignment
identical(colnames(ASV_matrix), rownames(meta_df))  # Must be TRUE

# ===
# CREATE PHYLOSEQ OBJECT
# ===
OTU  <- otu_table(ASV_matrix, taxa_are_rows = TRUE)
TAX  <- tax_table(TAXA_MAT)
META <- sample_data(meta_df)

physeq <- phyloseq(OTU, TAX, META)

# ===
# CHECK RESULTS
# ===
physeq
sample_names(physeq)    # Should match your CamachoNames
taxa_names(physeq)     # Should match ASV IDs
sample_variables(physeq)  # Your metadata columns
