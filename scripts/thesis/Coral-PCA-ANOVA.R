library("dplyr")
library("phyloseq")
library("vegan")
library("ggplot2")
library("tibble")
library("readxl")

physeq <- readRDS("/Users/chloecamacho/Desktop/thesis/Bioinformatics/Microbial Analysis/coral_physeq.rds")

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


#get % constrained variance explained by CCA axes
cca_var <- coral_pca$CCA$eig / sum(coral_pca$CCA$eig) * 100
axis_names <- paste0(c("CCA1", "CCA2"), 
                     " (",
                     round(cca_var[1:2], 1),
                     "%)")

pca <- ggplot(pca_points, aes(x= CCA1, y= CCA2, color = Site, shape = Species)) +
  geom_point(size = 4) +
  labs(x = axis_names[1], y = axis_names[2],  # ← Adds coefficients here!
       color = "Site", shape = "Species") +
  theme_minimal() +  
  scale_color_manual(values = site_colors)

pca

# Save if needed
ggsave("PCA_allcoral.png", width = 12, height = 10, dpi = 300)

############################## PERMANOVA ##############################

#Does Site structure communities?
physeq_sites_test <- otu_table(physeq_rel_abund) %>%  
  t() %>% #transpose it if necsessary
  #vegan::vegdist(., method = "euclidean", binary = FALSE, na.rm = TRUE) %>%  
  as.matrix(.)

adonis_euclidean <- adonis2(
  physeq_sites_test ~ Site, 
  data = physeq_sites_metadata, 
  permutations = 999
)

adonis_euclidean

#trying with bray distance instead 
bray <- vegdist(physeq_sites_test, method = "bray")

#Site only
adonis2(bray ~ Site, data = physeq_sites_metadata, permutations = 999)

#Site and Species
adonis2(bray ~ Site + Species, data = physeq_sites_metadata, permutations = 999)

# Interaction of Site and Species
adonis2(bray ~ Site * Species, data = physeq_sites_metadata, permutations = 999)

# 1. Check dispersion assumption (should be non-significant)
disp <- betadisper(bray, interaction(physeq_sites_metadata$Site, physeq_sites_metadata$Species))
anova(disp)

# 2. Plot interaction (species colored, sites shaped)

#NO JITTER
plot_ordination(physeq_sites, ordinate(physeq_sites, "PCoA", "bray"), 
                color = "Species", shape = "Site") +
  scale_color_manual(values = coral_colors) +
  facet_wrap(~Species)  # See site separation BY species

#WITH JITTER
plot_ordination(physeq_phylo, ordinate(physeq_phylo, "PCoA", "bray"), 
                color = "Species", shape = "Site") +
  geom_jitter(width = 0.02, height = 0.02, alpha = 0.7) +
  scale_color_manual(values = coral_colors) +
  facet_wrap(~Species)


# Test dispersion by Species only
disp_species <- betadisper(bray, physeq_sites_metadata$Species) #PASSED 
anova(disp_species)
plot(disp_species)  # Visualizes group spreads

plot_ordination(physeq_sites, ordinate(physeq_sites, "PCoA", "bray"), 
                color = "Species", shape = "Site") +
  scale_color_manual(values = coral_colors) +
  facet_wrap(~Species)



################  3 FACTOR ANOVA  ############

physeq1 <- readRDS("coral_physeq.rds")

md <- read_excel("~/Desktop/thesis/Bioinformatics/meta?.xlsx")
names(md)  # check it has CamachoName, Species, Site, Timepoint

physeq_phylo <- physeq1 %>% 
  subset_samples(Site %in% c("Brewer's Bay", "Saba", "Perseverance", "Krum Bay", "Rupert's Rock", "UVI Nursery")) %>%
  subset_samples(SampleType == "Coral") %>%
  filter_taxa(function(x) sum(x) > 0, TRUE)

class(physeq_phylo) #phyloseq package 
nrow(sample_data(physeq_phylo))  # 120

#Keep only samples that are in physeq_phylo
samples_in_physeq <- sample_names(physeq_phylo)
md <- md %>%
  filter(CamachoName %in% samples_in_physeq)  # CamachoName is your sample ID

# Rename CamachoName to match shannon_df
md <- md %>%
  rename(sample_id = CamachoName)

# 4. Ensure factors
md <- md %>%
  mutate(
    Species   = factor(Species),
    Site      = factor(Site),
    Timepoint = factor(Timepoint)
  )

#Compute Shannon and add to md

shannon_df <- estimate_richness(physeq_phylo, measures = "Shannon") %>%
  tibble::rownames_to_column("sample_id")

#adding md and shannon_df using sample_id column
md <- md %>%
  left_join(shannon_df, by = "sample_id")

# Check Timepoint and drop NAs
unique(md$Timepoint)   
md <- md %>%
  filter(!is.na(Shannon), !is.na(Timepoint))

# make dataframe for 3 factor anova 
df <- data.frame(
  y         = md$Shannon,
  Species   = md$Species,
  Site      = md$Site,
  Timepoint = md$Timepoint   
)

#all 6 comparisons
model <- aov(y ~ Site * Species * Timepoint, data = df)
summary(model)

#Df Sum Sq Mean Sq F value   Pr(>F)    
#Site           5  83.62  16.724   9.446 1.95e-07 ***
#Species        2   5.51   2.757   1.557    0.216    
#Site:Species   9  13.87   1.541   0.870    0.554    
#Residuals    103 182.37   1.771                     
#---
  #Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1

#just main effect comparisons
model2 <- aov(y ~ Site + Species + Timepoint, data = df)
summary(model2)

#Df Sum Sq Mean Sq F value   Pr(>F)    
#Site          5  83.62  16.724   9.545 1.32e-07 ***
#Species       2   5.51   2.757   1.574    0.212    
#Residuals   112 196.24   1.752                     
#---
#Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1

#Tukey Test to see which Sites differ
TukeyHSD(model, which = "Site")
plot(TukeyHSD(model, "Site"))


#JUST reef sites
df_reefs <- df %>%
  filter(Site != "UVI Nursery")

model_reefs <- aov(y ~ Site * Species, data = df_reefs)
summary(model_reefs)

#Df Sum Sq Mean Sq F value Pr(>F)
#Site          4   3.63  0.9064   0.462  0.763
#Species       2   1.76  0.8789   0.448  0.640
#Site:Species  7   8.65  1.2350   0.630  0.730
#Residuals    85 166.75  1.9617

TukeyHSD(model_reefs, which = "Site")


# SAMPLE SIZE BAR PLOT ##
library(ggplot2)

# Count samples per Species from the phyloseq object
ns <- as.data.frame(
  table(sample_data(physeq1)$Species)
)
names(ns) <- c("Species", "n")

# Sort so bars are ordered nicely (optional)
ns <- ns[order(ns$n, decreasing = TRUE), ]

# Draw the barplot
ggplot(ns, aes(x = Species, y = n, fill = Species)) +
  geom_col() +
  geom_text(aes(label = n), vjust = -0.5, size = 4) +
  labs(
    x = "Coral species",
    y = "Number of samples",
    title = "Sample size per species"
  ) +
  theme_minimal() +
  theme(legend.position = "none")

