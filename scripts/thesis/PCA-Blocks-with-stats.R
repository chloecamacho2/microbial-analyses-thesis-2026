# BLOCK PCA

library("dplyr")
library("phyloseq")
library("vegan")
library("ggplot2")
library("devtools")
library(pairwiseAdonis)

physeq <- readRDS("~/Desktop/thesis/Bioinformatics/physeq_filtered.rds")

site_colors <- c("Brewer's Bay"= "turquoise", "Saba"="lightpink", "Perseverance"="navy", "Krum Bay"="darkorange", "Rupert's Rock"="brown")

physeq_blocks <- subset_samples(physeq, Site %in% c("Brewer's Bay", "Saba", "Perseverance", "Krum Bay", "Rupert's Rock")) %>% 
  subset_samples(., SampleType == "Block") %>%  # ← Changed from "Coral" to "Block"
  filter_taxa(function(x) sum(x) > 0, TRUE)

physeq_blocks_metadata <- as(sample_data(physeq_blocks), "data.frame")

physeq_rel_abund_blocks <- transform_sample_counts(physeq_blocks, function(x) x/sum(x))

# Input matrix for CCA
physeq_blocks_test <- otu_table(physeq_rel_abund_blocks) %>%  
  t() %>% 
  as.matrix(.)

# Run PCA
block_pca <- vegan::cca(physeq_blocks_test ~ Site, data = physeq_blocks_metadata)

# Plot points
pca_block_points <- as.data.frame(vegan::scores(block_pca, display = "sites")) %>% 
  tibble::rownames_to_column(., var = "sample_id") %>% 
  left_join(., (physeq_blocks_metadata %>% tibble::rownames_to_column(var = "sample_id")), 
            by = c("sample_id"))

# Axis labels
cca_var_block <- block_pca$CCA$eig / sum(block_pca$CCA$eig) * 100
axis_names_block <- paste0(c("CCA1", "CCA2"), " (", round(cca_var_block[1:2], 1), "%)")

# Plot 
pca_blocks <- ggplot(pca_block_points, aes(x = CCA1, y = CCA2, color = Site)) +
  geom_point(size = 4) +
  labs(x = axis_names_block[1], y = axis_names_block[2],
       color = "Site", shape = "Species",
       title = "PCA of Block Samples By Site") +
  theme_minimal() +  
  scale_color_manual(values = site_colors)

print(pca_blocks)

# Save if needed
ggsave("PCA_Blocks_all_sites.png", width = 12, height = 10, dpi = 300)

############################## PERMANOVA ##############################

#Does Site structure communities?
physeq_sites_blocks <- otu_table(physeq_rel_abund_blocks) %>%  
  t() %>% #transpose it if necsessary
  #vegan::vegdist(., method = "euclidean", binary = FALSE, na.rm = TRUE) %>%  
  as.matrix(.)

adonis_asv.res <- vegan::adonis2(data = physeq_blocks_metadata, method = "euclidean", permutations = 999, formula = physeq_sites_blocks ~ Site, parallel = FALSE, by = "terms")

adonis_asv.res

#trying with bray distance instead 
bray_blocks <- vegdist(physeq_sites_blocks, method = "bray")
adonis_bray_blocks <- adonis2(bray_blocks ~ Site, data = physeq_blocks_metadata, permutations = 999)

adonis_bray_blocks

# Df SumOfSqs      R2      F Pr(>F)    
# Model     4    3.682 0.11495 3.0848  0.001 ***
#   Residual 95   28.352 0.88505                  
# Total    99   32.034 1.00000  

#Pairwise comparison (to see which sites differed)
pairwise.adonis(bray_blocks, physeq_blocks_metadata$Site, p.adjust.m = "BH")
          #Top Comparison
#                           Df SumsOfSqs  F.Model   R2      p.value
# Brewer's Bay vs Krum Bay  1 1.3641490 4.829426 0.11275953   0.001

#Check dispersion homogeneity (PERMANOVA assumption)
dist <- vegdist(physeq_sites_blocks, "bray")
anova(betadisper(dist, physeq_blocks_metadata$Site))

# Plot dispersion by site (shows which sites vary most)
boxplot(betadisper(dist, physeq_blocks_metadata$Site)$distances ~ physeq_blocks_metadata$Site)

##############################################################################

#BLOCK PCA w/ Time

physeq_blocks_metadata <- as(sample_data(physeq_blocks), "data.frame")

# Build Timepoint on that SAME object (no need for a separate sdat_blocks_df)
physeq_blocks_metadata$SamplingDate_fixed <- case_when(
  physeq_blocks_metadata$SamplingDate == "Aug-14-2025"    ~ "2025-08-14",
  physeq_blocks_metadata$SamplingDate == "August-14-2025" ~ "2025-08-14",
  physeq_blocks_metadata$SamplingDate == "July-14-2025"   ~ "2025-07-14",
  physeq_blocks_metadata$SamplingDate == "July-15-2025"   ~ "2025-07-15",
  physeq_blocks_metadata$SamplingDate == "July-16-2025"   ~ "2025-07-16",
  TRUE ~ NA_character_
)

physeq_blocks_metadata$Timepoint <- case_when(
  physeq_blocks_metadata$SamplingDate_fixed %in% c("2025-07-14", "2025-07-15", "2025-07-16") ~ "T1 (July 2025)",
  physeq_blocks_metadata$SamplingDate_fixed %in% c("2025-08-14") ~ "T2 (August 2025)",
  TRUE ~ NA_character_
)

physeq_blocks_metadata$Timepoint <- factor(
  physeq_blocks_metadata$Timepoint,
  levels = c("T1 (July 2025)", "T2 (August 2025)")
)

# Check that Timepoint is really there
str(physeq_blocks_metadata$Timepoint)
table(physeq_blocks_metadata$Timepoint, useNA = "always")
physeq_blocks_metadata$Species <- NULL
physeq_blocks_metadata$Genotype <- NULL


# Create physeq_sites_blocks (OTUs as rows, samples as columns, relative abundance)
physeq_rel_abund_blocks <- transform_sample_counts(physeq_blocks, function(x) x/sum(x))

physeq_sites_blocks <- otu_table(physeq_rel_abund_blocks) %>%  
  t() %>%  # Transpose: OTUs = rows, Samples = columns  
  as.matrix(.)

#run PCA
block_pca_time <- vegan::cca(physeq_sites_blocks ~ Site + Timepoint, 
                             data = physeq_blocks_metadata)

# getting points
pca_block_time_points <- as.data.frame(vegan::scores(block_pca_time, display = "sites")) %>% 
  tibble::rownames_to_column(., var = "sample_id") %>% 
  left_join(., (physeq_blocks_metadata %>% tibble::rownames_to_column(var = "sample_id")), 
            by = c("sample_id"))
# Axis labels
cca_var_block_time <- block_pca_time$CCA$eig / sum(block_pca_time$CCA$eig) * 100
axis_names_block_time <- paste0(c("CCA1", "CCA2"), " (", round(cca_var_block_time[1:2], 1), "%)")

# Plot 
pca_blocks_time <- ggplot(pca_block_time_points, 
                          aes(x = CCA1, y = CCA2, color = Site)) +
  geom_point(size = 4, stroke = 1) +
  facet_wrap(~ Timepoint, ncol = 2) +  # T1 | T2 side-by-side
  labs(x = axis_names_block_time[1], y = axis_names_block_time[2],
       color = "Site", title = "Block Samples by Site") +
  theme_minimal() +  
  scale_color_manual(values = site_colors)

print(pca_blocks_time)

# Save if needed
ggsave("PCA_Blocks_site_time.png", width = 12, height = 10, dpi = 300)

####### QUANTIFYING DIFFERENCES #######

#comparing bray-curtis dissimilarities for Saba T1 vs T2

#Using physeq_blocks_test from earlier as the matrix 

# Saba subset (now works - samples as rows)
saba_data <- physeq_blocks_test[physeq_blocks_metadata$Site == "Saba", ]
saba_meta <- physeq_blocks_metadata[physeq_blocks_metadata$Site == "Saba", ]

# T1/T2 indices
t1_idx <- which(saba_meta$Timepoint == "T1 (July 2025)")
t2_idx <- which(saba_meta$Timepoint == "T2 (August 2025)")

# T1-T2 pairwise BC
saba_dist <- vegdist(saba_data, "bray")
saba_t1t2_bc <- as.matrix(saba_dist)[t1_idx, t2_idx]
mean(saba_t1t2_bc)  # Saba turnover metric = 0.8183283

#comparing bray-curtis dissimilarities for Brewer's Bay T1 vs T2

bb_data <- physeq_blocks_test[physeq_blocks_metadata$Site == "Brewer's Bay", ]
bb_meta <- physeq_blocks_metadata[physeq_blocks_metadata$Site == "Brewer's Bay", ]

# T1/T2 indices
bb_t1_idx <- which(bb_meta$Timepoint == "T1 (July 2025)")
bb_t2_idx <- which(bb_meta$Timepoint == "T2 (August 2025)")

# T1-T2 pairwise BC
bb_dist <- vegdist(bb_data, "bray")
bb_t1t2_bc <- as.matrix(bb_dist)[bb_t1_idx, bb_t2_idx]
mean(bb_t1t2_bc)  # Brewer's Bay turnover metric = 0.7230829

#comparing bray-curtis dissimilarities for Rupert's Rock T1 vs T2

rr_data <- physeq_blocks_test[physeq_blocks_metadata$Site == "Rupert's Rock", ]
rr_meta <- physeq_blocks_metadata[physeq_blocks_metadata$Site == "Rupert's Rock", ]

# T1/T2 indices
rr_t1_idx <- which(rr_meta$Timepoint == "T1 (July 2025)")
rr_t2_idx <- which(rr_meta$Timepoint == "T2 (August 2025)")

# T1-T2 pairwise BC
rr_dist <- vegdist(rr_data, "bray")
rr_t1t2_bc <- as.matrix(rr_dist)[rr_t1_idx, rr_t2_idx]
mean(rr_t1t2_bc)  # Rupert's Rock turnover metric = 0.8245728

#comparing bray-curtis dissimilarities for Krum Bay T1 vs T2

kb_data <- physeq_blocks_test[physeq_blocks_metadata$Site == "Krum Bay", ]
kb_meta <- physeq_blocks_metadata[physeq_blocks_metadata$Site == "Krum Bay", ]

# T1/T2 indices
kb_t1_idx <- which(kb_meta$Timepoint == "T1 (July 2025)")
kb_t2_idx <- which(kb_meta$Timepoint == "T2 (August 2025)")

# T1-T2 pairwise BC
kb_dist <- vegdist(kb_data, "bray")
kb_t1t2_bc <- as.matrix(kb_dist)[kb_t1_idx, kb_t2_idx]
mean(kb_t1t2_bc)  # Krum Bay turnover metric = 0.875306

#comparing bray-curtis dissimilarities for Perseverance T1 vs T2

pv_data <- physeq_blocks_test[physeq_blocks_metadata$Site == "Perseverance", ]
pv_meta <- physeq_blocks_metadata[physeq_blocks_metadata$Site == "Perseverance", ]

# T1/T2 indices
pv_t1_idx <- which(pv_meta$Timepoint == "T1 (July 2025)")
pv_t2_idx <- which(pv_meta$Timepoint == "T2 (August 2025)")

# T1-T2 pairwise BC
pv_dist <- vegdist(pv_data, "bray")
pv_t1t2_bc <- as.matrix(pv_dist)[pv_t1_idx, pv_t2_idx]
mean(pv_t1t2_bc)  # Perseverance turnover metric = 0.8075703

