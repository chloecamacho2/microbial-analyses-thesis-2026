library(tidyverse)
library(phyloseq)
library(readxl)

# 1. READ DATA
ASV_coral_raw <- readRDS("~/Desktop/thesis/Bioinformatics/raw data/camacho_g3_h3_prok_asvs.df.rds")
TAXA_coral <- read_tsv("~/Desktop/thesis/Bioinformatics/raw data/camacho_g3_h3_prok_asvs.taxa.df.tsv")
meta_coral <- read_excel("~/Desktop/thesis/Bioinformatics/meta?.xlsx")
contam_ASV_coral <- read.delim("~/Desktop/thesis/Bioinformatics/raw data/potential_contam_asvs.camacho_g3_h3.taxa.tsv")

# 2. CONTAMINANT FILTERING (documentation)
tax_cols <- c("Domain", "Phylum", "Class", "Order", "Family", "Genus", "Species")  
bad_lineage <- contam_ASV_coral %>% select(all_of(tax_cols)) %>% distinct()
matching_taxa <- TAXA_coral %>%
  inner_join(bad_lineage, by = tax_cols, suffix = c("_main", "_bad")) %>%
  pull(asv_id)
cat("Taxonomy matches to contaminants:", length(matching_taxa), "\n")

# 3. BUILD CLEAN ASV MATRIX
ASV_coral_wide <- ASV_coral_raw %>% 
  select(asv_id, relabeled, counts) %>%
  pivot_wider(id_cols = asv_id, names_from = relabeled, values_from = counts, values_fill = 0) %>%
  column_to_rownames("asv_id") %>%
  mutate(across(everything(), as.numeric)) %>%
  as.matrix()

# Remove contaminants 
ASV_coral_final <- ASV_coral_wide[!rownames(ASV_coral_wide) %in% contam_ASV_coral$asv_id, ]

#prevalence filter
coral_cols <- colnames(ASV_coral_final)
prevalence_threshold <- 0.1 * length(coral_cols)
keep_asvs <- rowSums(ASV_coral_final > 0) >= prevalence_threshold
cat("High-prevalence ASVs kept:", sum(keep_asvs), "out of", nrow(ASV_coral_final), "\n")
ASV_coral_final <- ASV_coral_final[keep_asvs, ]

cat("Final ASVs:", nrow(ASV_coral_final), "x", ncol(ASV_coral_final), "samples\n") #2105 x 233 samples

# 4. CLEAN & MATCH METADATA
sample_names <- colnames(ASV_coral_final)
meta_coral_clean <- meta_coral %>% 
  distinct(CamachoName, .keep_all = TRUE) %>%
  rename(SamplingDate = `Sampling Date`) %>%
  mutate(CamachoName = str_replace(CamachoName, "PSV", "PSU")) %>%
  mutate(SamplingDate_fixed = case_when(
    SamplingDate == "Aug-14-2025" | SamplingDate == "August-14-2025" ~ "2025-08-14",
    SamplingDate == "July-14-2025" ~ "2025-07-14",
    SamplingDate == "July-15-2025" ~ "2025-07-15", 
    SamplingDate == "July-16-2025" ~ "2025-07-16",
    TRUE ~ as.character(SamplingDate)
  )) %>%
  mutate(Timepoint = case_when(
    SamplingDate_fixed %in% c("2025-07-14", "2025-07-15", "2025-07-16") ~ "T1 (July 2025)",
    SamplingDate_fixed %in% c("2025-08-14") ~ "T2 (August 2025)",
    TRUE ~ NA_character_
  ))

# Match metadata to samples
meta_filtered <- meta_coral_clean %>% filter(CamachoName %in% sample_names)
meta_ordered <- meta_filtered[match(sample_names, meta_filtered$CamachoName), ]
meta_coral_df <- as.data.frame(meta_ordered)
rownames(meta_coral_df) <- sample_names
meta_coral_df$CamachoName <- NULL

# 5. BUILD TAXONOMY & PHYLOSEQ
TAXA_matrix_coral <- TAXA_coral %>%
  filter(asv_id %in% rownames(ASV_coral_final)) %>%
  column_to_rownames("asv_id") %>%
  select(all_of(tax_cols)) %>%
  as.matrix()

ps_coral <- phyloseq(
  otu_table(ASV_coral_final, taxa_are_rows = TRUE),
  sample_data(meta_coral_df),
  tax_table(TAXA_matrix_coral)
)

# 6. FILTER TO CORAL SAMPLES ONLY (remove controls + Block samples)
# Keep APAL_PSV1/3 but rename to PSU, reassign to T1 UVI Nursery
# Remove ALL other controls + ANY "Block" samples
all_samples <- sample_names(ps_coral)

# Define samples to KEEP (corals + your 2 special APALs)
keep_apal <- c("APAL_PSV1", "APAL_PSV3")
remove_controls <- all_samples[grepl("Kit|Negative|Mock|BS_|WT_SP", all_samples)]
remove_blocks <- all_samples[grepl("Block", all_samples)]

# Combined removal list (everything except corals + your 2 APALs)
samples_to_drop <- c(remove_controls, remove_blocks)
samples_to_drop <- samples_to_drop[!samples_to_drop %in% keep_apal]

ps_coral_clean <- prune_samples(!sample_names(ps_coral) %in% samples_to_drop, ps_coral)

# 7. RENAME & REASSIGN APAL SAMPLES TO T1 CORALS
sample_names(ps_coral_clean) <- str_replace_all(
  sample_names(ps_coral_clean), 
  c("APAL_PSV1" = "APAL_PSU1", "APAL_PSV3" = "APAL_PSU3")
)

# Reassign APALs to T1 UVI Nursery CORALS
apal_samples <- c("APAL_PSU1", "APAL_PSU3")
sample_data(ps_coral_clean)$Timepoint[sample_names(ps_coral_clean) %in% apal_samples] <- "T1 (July 2025)"
sample_data(ps_coral_clean)$Site[sample_names(ps_coral_clean) %in% apal_samples] <- "UVI Nursery"

# 8. FINAL VALIDATION - CORALS ONLY
cat("\n=== FINAL CORAL-ONLY PHYLOSEQ ===\n")
print(ps_coral_clean)
cat("Samples:", nsamples(ps_coral_clean), "| ASVs:", ntaxa(ps_coral_clean), "\n")
print("T1/T2 CORALS by Site:")
print(table(sample_data(ps_coral_clean)$Timepoint, sample_data(ps_coral_clean)$Site))

# Check no blocks/controls remain
remaining_samples <- sample_names(ps_coral_clean)
cat("Block samples remaining:", sum(grepl("Block", remaining_samples)), "\n")
cat("Control patterns remaining:", sum(grepl("Kit|Negative|Mock|BS_|WT_SP", remaining_samples)), "\n")

#merging phyloseq objects 

physeq_filtered <- readRDS("~/Desktop/thesis/Bioinformatics/physeq_filtered.rds")
coral_physeq    <- readRDS("/Users/chloecamacho/Desktop/thesis/Bioinformatics/Microbial Analysis/coral_physeq.rds")

PHYLOSEQ <- merge_phyloseq(physeq_filtered, coral_physeq)

#double checking theyre both phyloseq objects
class(physeq_filtered)
class(coral_physeq)

#verifying counts in each object
nsamples(physeq_filtered) #168 total samples
# For physeq_filtered
md_filt <- sample_data(physeq_filtered)
ncoral_filt <- sum(md_filt$SampleType == "Coral", na.rm = TRUE)
cat("Coral samples in physeq_filtered:", ncoral_filt, "\n") #68 coral samples

nsamples(coral_physeq) #122 total samples
# For coral_physeq
md_coral <- sample_data(coral_physeq)
ncoral_coral <- sum(md_coral$SampleType == "Coral", na.rm = TRUE)
cat("Coral samples in coral_physeq:", ncoral_coral, "\n") #120 coral samples


#checking counts in combined object 
md_PHYLO <- sample_data(PHYLOSEQ)
ncoral <- sum(md_PHYLO$SampleType == "Coral", na.rm = TRUE)
cat("Coral samples in PHYLOSEQ:", ncoral, "\n") #120 coral samples

nsamples(PHYLOSEQ) #222 total samples
sample_names(PHYLOSEQ)

#saving combined phyloseq object
saveRDS(PHYLOSEQ, "PHYLOSEQ.rds")
