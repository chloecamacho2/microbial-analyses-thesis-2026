library("dplyr")
library("phyloseq")
library("vegan")
library("ggplot2")

physeq <- readRDS("~/Desktop/thesis/Bioinformatics/physeq_filtered.rds")

site_colors <- c("Brewer's Bay"= "turquoise", "Saba"="lightpink", "Perseverance"="navy", "Krum Bay"="darkorange", "Rupert's Rock"="brown")

coral_colors <- c("Acropora palmata"="gold3", "Montastraea cavernosa"="darkorchid","Porites astreoides"="darkgreen")

sample_data(physeq)$Site <-factor(sample_data(physeq)$Site, levels = c("Brewer's Bay", "Saba", "Perseverance", "Krum Bay", "Rupert's Rock"))

physeq_sites <- subset_samples(physeq, Site %in% c("Brewer's Bay", "Saba", "Perseverance", "Krum Bay", "Rupert's Rock")) %>% subset_samples(., SampleType == "Coral") %>% 
  filter_taxa(function(x) sum(x) > 0, TRUE) #subsetted for site and just coral samples 

physeq_sites_metadata <- as(sample_data(physeq_sites), "data.frame")  #saving metadata from physeq as a data frame to pull Site and Species info; needs to be dataframe and not phyloseq object

physeq_rel_abund <- transform_sample_counts(physeq_sites, function(x) x/sum(x)) #getting relative abundances instead of count

#input for statistical test (can use for nMDS or others)
physeq_sites_test <- otu_table(physeq_rel_abund) %>%  
  t() %>% #transpose it if necsessary
  #vegan::vegdist(., method = "euclidean", binary = FALSE, na.rm = TRUE) %>%  
  as.matrix(.)

#calculates PCoA 
coral_pca <- vegan::cca(physeq_sites_test ~ Site, dist = "euclidean", data = physeq_sites_metadata)

coral_pca
summary(coral_pca)

# plot PCA

pca_points <- as.data.frame(vegan::scores(coral_pca,display = "sites")) %>% 
  tibble::rownames_to_column(.,var = "sample_id") %>% 
  left_join(., (physeq_sites_metadata %>% 
                  tibble::rownames_to_column(var= "sample_id")), by = c("sample_id"))
    

pca <- ggplot(pca_points, aes(x= CCA1, y= CCA2, color = Site, shape = Species)) +
         geom_point(size = 4) +
         labs(color = "Site", shape = "Species")+
         theme_minimal() +  
        scale_color_manual(values = site_colors)
pca

# PERMANOVA: Does Site structure communities?
physeq_sites_test <- otu_table(physeq_rel_abund) %>%  
  t() %>% #transpose it if necsessary
  #vegan::vegdist(., method = "euclidean", binary = FALSE, na.rm = TRUE) %>%  
  as.matrix(.)

adonis_asv.res <- vegan::adonis2(data = physeq_sites_metadata, method = "euclidean", permutations = 999, formula = physeq_sites_test ~ Site, parallel = FALSE, by = "terms")

adonis_asv.res

# Save if needed
ggsave("PCA_coral.png", width = 12, height = 10, dpi = 300)

#trying with bray distance instead 
bray <- vegdist(physeq_sites_test, method = "bray")
adonis2(bray ~ Site, data = physeq_sites_metadata, permutations = 999)

#testing effect of species
adonis2(bray ~ Site + Species, data = physeq_sites_metadata, permutations = 999)

#testing interactions

#Does Site effect differ BY species?
adonis2(bray ~ Site * Species, data = physeq_sites_metadata, permutations = 999)

# 1. Check dispersion assumption (should be non-significant)
disp <- betadisper(bray, interaction(physeq_sites_metadata$Site, physeq_sites_metadata$Species))
anova(disp)

# 2. Plot interaction (species colored, sites shaped)
plot_ordination(physeq_sites, ordinate(physeq_sites, "PCoA", "bray"), 
                color = "Species", shape = "Site") +
  scale_color_manual(values = coral_colors) +
  facet_wrap(~Species)  # See site separation BY species

# Test dispersion by Species only
disp_species <- betadisper(bray, physeq_sites_metadata$Species) #FAILED --> unequal dispersion 
anova(disp_species)
plot(disp_species)  # Visualizes group spreads
