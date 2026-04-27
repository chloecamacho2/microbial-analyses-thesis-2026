library(vegan)
library(phyloseq)


# Load full phyloseq
physeq <- readRDS("~/Desktop/thesis/Bioinformatics/physeq_filtered.rds")

# Create distance matrix for FULL dataset first (for reference)
dist_mat_full <- phyloseq::distance(physeq, "bray")

# Subset to blocks from your 5 sites (YOUR CODE - PERFECT ✓)
physeq_blocks <- subset_samples(physeq, Site %in% c("Brewer's Bay", "Saba", "Perseverance", "Krum Bay", "Rupert's Rock")) %>% 
  subset_samples(., SampleType == "Block") %>%  
  filter_taxa(function(x) sum(x) > 0, TRUE)

# Fix dates 
sample_data(physeq_blocks)$SamplingDate_fixed <- case_when(
  sample_data(physeq_blocks)$SamplingDate == "Aug-14-2025"      ~ "2025-08-14",
  sample_data(physeq_blocks)$SamplingDate == "August-14-2025"   ~ "2025-08-14",
  sample_data(physeq_blocks)$SamplingDate == "July-14-2025"     ~ "2025-07-14",
  sample_data(physeq_blocks)$SamplingDate == "July-15-2025"     ~ "2025-07-15", 
  sample_data(physeq_blocks)$SamplingDate == "July-16-2025"     ~ "2025-07-16",
  TRUE ~ NA_character_
)

# Create timepoints
sample_data(physeq_blocks)$Timepoint <- case_when(
  sample_data(physeq_blocks)$SamplingDate_fixed %in% c("2025-07-14", "2025-07-15", "2025-07-16") ~ "T1 (July 2025)",
  sample_data(physeq_blocks)$SamplingDate_fixed %in% c("2025-08-14") ~ "T2 (August 2025)",
  TRUE ~ NA_character_
)

# verify
table(sample_data(physeq_blocks)$Site, sample_data(physeq_blocks)$Timepoint)

# Create factors using Timepoint column name
sample_data(physeq_blocks)$Site <- factor(sample_data(physeq_blocks)$Site)
sample_data(physeq_blocks)$Time <- factor(sample_data(physeq_blocks)$Timepoint)  

# Distance matrix for physeq_blocks ONLY
dist_mat <- phyloseq::distance(physeq_blocks, "bray")

# Extract metadata as data.frame
metadata <- as.data.frame(sample_data(physeq_blocks))

# permanova
permanova_res <- adonis2(dist_mat ~ metadata$Site * metadata$Time, 
                         permutations = 9999)

print(permanova_res)

#main effects only 

main_effects_perm <- adonis2(dist_mat ~ metadata$Site + metadata$Time, 
                             permutations = 9999)
print(main_effects_perm)

# Pairwise site comparisons (all timepoints combined)
library(pairwiseAdonis) 
pairwise_res <- pairwise.adonis(dist_mat, metadata$Site)

print(pairwise_res)

#T1

t1_samples <- subset_samples(physeq_blocks, Timepoint == "T1 (July 2025)")
dist_t1 <- distance(t1_samples, "bray")
pairwise.adonis(dist_t1, sample_data(t1_samples)$Site)

#T2

t2_samples <- subset_samples(physeq_blocks, Timepoint == "T2 (August 2025)")
dist_t2 <- distance(t2_samples, "bray")
pairwise.adonis(dist_t2, sample_data(t2_samples)$Site)

