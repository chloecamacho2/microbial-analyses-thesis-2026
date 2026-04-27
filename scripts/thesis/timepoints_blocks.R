#Stacked bar chart of blocks - initial vs final

library("randomcoloR")
library("pals")
library("cowplot")
library("dplyr")
library("tibble")
library("phyloseq")

physeq <-readRDS("scripts/thesis/my_physeq_filtered.rds")

# Transform counts to relative abundance to visualize bar charts
# Create an individual phyloseq object for each coral or the SW control samples.

#getting relative abundance for only class level
ps_ra <- transform_sample_counts(physeq, function(OTU) OTU/sum(OTU)) %>%
  filter_taxa(function(x) sum(x) > 0.01, TRUE) %>%  # Keep abundant taxa
  tax_glom(taxrank = "Family") %>% 
  subset_samples(Site != "UVI Nursery")

# Create condition column based on your sites
sample_data(ps_ra)$condition <- dplyr::case_when(
  sample_data(ps_ra)$Site %in% c("Saba", "Perseverance") ~ "Unaffected",
  sample_data(ps_ra)$Site %in% c("Krum Bay", "Rupert's Rock") ~ "Compromised", 
  sample_data(ps_ra)$Site == "Brewer's Bay" ~ "Target")

# Set factor levels for proper ordering
sample_data(ps_ra)$condition <- factor(sample_data(ps_ra)$condition, 
                                       levels = c("Unaffected", "Compromised", "Target"))

# Optional: verify it worked
table(sample_data(ps_ra)$Site, sample_data(ps_ra)$condition)

#clean up 

# STEP 1: Extract as REGULAR data.frame (FORCE it)
sdat_df <- data.frame(sample_data(physeq), stringsAsFactors = FALSE)

# fixing inconsistent dates
sdat_df$SamplingDate_fixed <- case_when(
  sdat_df$SamplingDate == "Aug-14-2025"      ~ "2025-08-14",
  sdat_df$SamplingDate == "August-14-2025"   ~ "2025-08-14",
  sdat_df$SamplingDate == "July-14-2025"     ~ "2025-07-14",
  sdat_df$SamplingDate == "July-15-2025"     ~ "2025-07-15", 
  sdat_df$SamplingDate == "July-16-2025"     ~ "2025-07-16",
  TRUE ~ NA_character_
)

#creating T1 and T2
sdat_df$Timepoint <- case_when(
    sdat_df$SamplingDate_fixed %in% c("2025-07-14", "2025-07-15", "2025-07-16") ~ "T1 (July 2025)",
    sdat_df$SamplingDate_fixed %in% c("2025-08-14") ~ "T2 (August 2025)",
    TRUE ~ NA_character_
  )
  
table(sdat_df$Timepoint)
#T1 (July 2025) = 71     
#T2 (August 2025) = 97 

sdat_df <- sdat_df %>% 
  filter(., SampleType == "Block")

sdat_df$Species <- NULL
sdat_df$Genotype <- NULL

# Add Sample column for matching (use rownames if no Sample column exists)
if(!"Sample" %in% colnames(sdat_df)) {
  sdat_df$Sample <- rownames(sdat_df)
}

# Get current ps_ra metadata and merge Timepoint
current_meta <- as(sample_data(ps_ra), "data.frame")
current_meta$Sample <- rownames(current_meta)

# Merge (only Block samples get Timepoint values)
current_meta <- current_meta %>%
  left_join(sdat_df %>% select(Sample, Timepoint), by = "Sample")

# Update ps_ra
rownames(current_meta) <- current_meta$Sample
sample_data(ps_ra) <- sample_data(current_meta)

# Subset ps_ra to only Block samples (now with Timepoint)
ps_ra_block <- subset_samples(ps_ra, SampleType == "Block")

# Perfect table now
table(sample_data(ps_ra_block)$Timepoint, sample_data(ps_ra_block)$condition)

length(get_taxa_unique(ps_ra_block, "Family")) 
n <- 165
palette <- distinctColorPalette(n)

# Each sample gets its own bar (default behavior)
t1 <- plot_bar(ps_ra_block, fill = "Family") +
  facet_wrap(~Timepoint + condition, scales = "free_x") +
  scale_fill_manual(values=palette) +
  labs(x = "Individual Samples", y = "Relative Abundance", fill = "Family") +
  theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5, size = 6))

print(t1)

# Save if needed
ggsave("3_conditions_blocks_T1andT2.png", width = 12, height = 10, dpi = 300)

#top 20 families
top20 <- names(sort(taxa_sums(ps_ra_block), TRUE))[1:20]
ps_top20 <- prune_taxa(top20, ps_ra_block)
ps_top20 <- merge_taxa(ps_top20, taxa_names(ps_top20)[-(1:20)], "Other")

n <- length(get_taxa_unique(ps_top20, "Family"))  # ~20
palette20 <- distinctColorPalette(n)

# Then your t1 plot code
t20 <- plot_bar(ps_top20, fill = "Family") +
  facet_wrap(~Timepoint + condition, scales = "free_x") +
  scale_fill_manual(values=palette20) +
  labs(x = "Individual Samples", y = "Relative Abundance", fill = "Family") +
  theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5, size = 6))

print(t20)

# Save if needed
ggsave("3_conditions_blocks_T1andT2_top20fam.png", width = 12, height = 10, dpi = 300)









#all sites + T1 and T2

# Each sample gets its own bar (default behavior)
t2 <- plot_bar(ps_ra_block, fill = "Class") +
  facet_wrap(Site ~ Timepoint, scales = "free_x") +
  scale_fill_manual(values=palette) +
  labs(x = "Individual Samples", y = "Relative Abundance", fill = "Class") +
  theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5, size = 6))

print(t2)

ggsave("blocks_all_sites_T1T2.png", width = 12, height = 10, dpi = 300)
# Plot 1: Unaffected sites (Brewer's Bay + Saba + Perseverance)
ps_unaffected <- subset_samples(ps_ra_block, 
                                Site %in% c("Brewer's Bay", "Saba", "Perseverance"))

p1 <- plot_bar(ps_unaffected, fill = "Class") +
  facet_wrap(.~ Site + Timepoint, scales = "free_x", ncol=2) +
  scale_fill_manual(values=palette) +
  labs(x = "Individual Blocks", y = "Relative Abundance", fill = "Class",
       title = "Communities by Class on Blocks from Unaffected Sites") +
  theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5, size = 6)
        )

print(p1)

# Plot 2: Compromised sites (Krum Bay + Rupert's Rock + Brewer's Bay)
ps_compromised <- subset_samples(ps_ra_block, 
                                 Site %in% c("Krum Bay", "Rupert's Rock", "Brewer's Bay"))

p2 <- plot_bar(ps_compromised, fill = "Class") +
  facet_wrap(.~Site + Timepoint, scales = "free_x", ncol=2) +
  scale_fill_manual(values=palette) +
  labs(x = "Individual Blocks", y = "Relative Abundance", fill = "Class",
       title = "Communities by Class on Blocks from Compromised Sites + Target") +
  theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5, size = 6),
        strip.text.x = element_text(size = 8, angle = 90))

print(p2)

