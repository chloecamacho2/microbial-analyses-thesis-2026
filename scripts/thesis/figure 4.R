#Stacked bar chart of 3 coral species by site
library("randomcoloR")
library("pals")
library("cowplot")
library("dplyr")
library("phyloseq")
library("ggplot2")

physeq <- readRDS("/Users/chloecamacho/Desktop/thesis/Bioinformatics/Microbial Analysis/coral_physeq.rds")

# Transform counts to relative abundance to visualize bar charts
# Create an individual phyloseq object for each coral or the SW control samples.

#getting relative abundance for only class level
ps_ra <- transform_sample_counts(physeq, function(OTU) OTU/sum(OTU)) %>%
  #filter_taxa(function(x) sum(x) > 0.01, TRUE) %>%  # Keep abundant taxa
  tax_glom(taxrank = "Family") %>% 
  subset_samples(Site != "UVI Nursery")

ps_ra <- transform_sample_counts(physeq, function(OTU) OTU/sum(OTU)) %>%
  tax_glom(taxrank = "Family", NArm = FALSE) %>%  # Keep NAs as zero
  subset_samples(Site != "UVI Nursery")

# Create condition column based on your sites
sample_data(ps_ra)$condition <- dplyr::case_when(
  sample_data(ps_ra)$Site %in% c("Saba", "Perseverance") ~ "Unimpaired",
  sample_data(ps_ra)$Site %in% c("Krum Bay", "Rupert's Rock") ~ "Impaired", 
  sample_data(ps_ra)$Site == "Brewer's Bay" ~ "Target")

# Set factor levels for proper ordering
sample_data(ps_ra)$condition <- factor(sample_data(ps_ra)$condition, 
                                       levels = c("Unimpaired", "Target","Impaired" ))

# Optional: verify it worked
table(sample_data(ps_ra)$Site, sample_data(ps_ra)$condition)

ps_ra_mcav = subset_samples(ps_ra, Species == "Montastraea cavernosa")
ps_ra_past = subset_samples(ps_ra, Species == "Porites astreoides")
ps_ra_apal = subset_samples(ps_ra, Species == "Acropora palmata")

#figure out how many colors you need
length(get_taxa_unique(ps_ra, "Family")) #169
length(get_taxa_unique(ps_ra, "Class")) #54
length(get_taxa_unique(ps_ra, "Phylum")) #31

#you can make n any number of colors you want; with as much difference between the colors as possible (distinct colors)
n <- 169
palette <- distinctColorPalette(n)

# Create all coral graphs from the level of Family

p1 <- plot_bar(ps_ra_mcav, fill="Family") +
  geom_bar(aes(fill=Family), stat="identity",position="stack",na.rm =TRUE)  +
  theme(strip.text=element_text(face="bold")) +
  scale_fill_manual(values=palette) +
  facet_grid(. ~ condition, scales="free", space="free") +
  ggtitle("Montastraea cavernosa") +
  labs(y = "Relative Abundance") +
  theme(plot.title = element_text(face="italic")) +
  theme(legend.position = "none") +
  theme(axis.title.x = element_blank(), axis.text.x = element_blank())

p2 <- plot_bar(ps_ra_past, fill="Family") +
  geom_bar(aes(fill=Family), stat="identity",position="stack", na.rm = TRUE) +
  theme(strip.text=element_text(face="bold")) +
  scale_fill_manual(values=palette) +
  facet_grid(. ~ condition, scales="free", space="free") +
  ggtitle("Porites astreoides") +
  labs(y = "Relative Abundance") +
  theme(plot.title = element_text(face="italic")) +
  theme(legend.position = "none") +
  theme(axis.title.x = element_blank(), axis.text.x = element_blank())

p3 <- plot_bar(ps_ra_apal, fill="Family") +
  geom_bar(aes(fill=Family), stat="identity",position="stack", na.rm = TRUE) +
  theme(strip.text=element_text(face="bold")) +
  scale_fill_manual(values=palette) +
  facet_grid(. ~ condition, scales="free", space="free") +
  ggtitle("Acropora palmata") +
  labs(y = "Relative Abundance") +
  theme(plot.title = element_text(face="italic")) +
  theme(legend.position = "none") +
  theme(axis.title.x = element_blank(), axis.text.x = element_blank())


plot_grid(p3,p1,p2, labels = c("A","B","C"), ncol = 2, nrow = 2)

library(cowplot)

# Extract legend from p1 (make sure it has one first)
p_legend <- p1 + theme(legend.position = "right")  # Temp plot with legend
legend <- get_legend(p_legend)

print(legend)

# Create 2x2 grid: plots in A,B,C + legend in D position
plot_grid(p3, p1, p2, legend, 
          labels = c("A", "B", "C", ""),  # Empty label for legend
          ncol = 2, nrow = 2,
          rel_heights = c(1, 1), 
          rel_widths = c(1, 1))

# Save if needed
ggsave("3coral_plots_with_legend.png", width = 12, height = 10, dpi = 300)

#Save plot for higher res
ggsave(
  filename = "Familystackedbar.png",   # File name
  plot = p_site2,            # The plot to save
  width = 40,                      # Width in inches
  height = 20,                      # Height in inches
  dpi = 300,                       # Resolution
  units = "in"                     # Unit for width and height
)



# Extract legend with MAXIMUM text width
p_legend <- p1 + 
  theme_void() + 
  guides(fill = guide_legend(
    keywidth = 1.0,      # ← Wide legend keys  
    keyheight = 1.0,
    label.position = "right",
    title.position = "top",
    nrow = NULL,         # ← Auto vertical layout
    label.hjust = 0      # ← Left-align text (no cutoff)
  )) +
  theme(legend.position = "right",
        legend.text = element_text(size = 6, margin = margin(r = 5)),  # Tiny text + right margin
        legend.title = element_text(size = 8),
        legend.margin = margin(l = 10, r = 20),
        legend.box.margin = margin(20, 40, 20, 20))

legend <- get_legend(p_legend)

# Visualize first
grid.draw(legend)

#NEED TO FIRGURE OUT HOW TO SAVE JUST LEGEND WITH GOOD DIMENSIONS








############# AVERAGE PER CATEGORY ##########
library(phyloseq)
library(dplyr)
library(ggplot2)
library(cowplot)

# Load physeq
physeq <- readRDS("/Users/chloecamacho/Desktop/thesis/Bioinformatics/Microbial Analysis/coral_physeq.rds")

# Relative abundance at Family level (no UVI Nursery)
ps_ra <- transform_sample_counts(physeq, function(OTU) OTU/sum(OTU)) %>%
  tax_glom(taxrank = "Family", NArm = FALSE) %>%
  subset_samples(Site != "UVI Nursery")

# Define condition
sample_data(ps_ra)$condition <- dplyr::case_when(
  sample_data(ps_ra)$Site %in% c("Saba", "Perseverance") ~ "Unimpaired",
  sample_data(ps_ra)$Site %in% c("Krum Bay", "Rupert's Rock") ~ "Impaired",
  sample_data(ps_ra)$Site == "Brewer's Bay" ~ "Target"
)

# Set factor levels (order in plots)
sample_data(ps_ra)$condition <- factor(sample_data(ps_ra)$condition,
                                       levels = c("Unimpaired", "Target", "Impaired"))

# Optional: check
table(sample_data(ps_ra)$Site, sample_data(ps_ra)$condition)



# --- Helper: compute mean rel. abund by Family per condition ---
compute_mean_by_condition <- function(ps_species) {
  # Relative abundances
  ps_ra <- transform_sample_counts(ps_species, function(OTU) OTU/sum(OTU))
  # Melt to long (Family column already exists)
  ps_long <- psmelt(ps_ra) %>%
    select(Sample, condition, Family, Abundance)
  # Mean per Family per condition
  avg <- ps_long %>%
    group_by(condition, Family) %>%
    summarise(mean_abund = mean(Abundance, na.rm = TRUE), .groups = "drop")
  # Ensure all 3 conditions appear for every Family (fill 0)
  conditions <- c("Unimpaired", "Target", "Impaired")
  avg <- avg %>%
    complete(condition = conditions, Family, fill = list(mean_abund = 0))
  return(avg)
}

# --- Compute mean‑by‑condition tables for each species ---
ps_ra_mcav <- subset_samples(ps_ra, Species == "Montastraea cavernosa")
ps_ra_past <- subset_samples(ps_ra, Species == "Porites astreoides")
ps_ra_apal <- subset_samples(ps_ra, Species == "Acropora palmata")

# Use exactly this line (no extra string)
mcav_avg <- compute_mean_by_condition(ps_ra_mcav)
past_avg <- compute_mean_by_condition(ps_ra_past)
apal_avg <- compute_mean_by_condition(ps_ra_apal)
# --- Create a distinct color palette for as many families as exist in all 3 ---
families <- unique(c(mcav_avg$Family, past_avg$Family, apal_avg$Family))
n_fam <- length(families)

# Function to get n very distinct colors (if you don't have distinctColorPalette already)
set.seed(123)
palette <- grDevices::hcl.colors(n_fam, palette = "Viridis")

# Map family → color
family_palette <- setNames(palette, families)


# --- Make the averaged bar plots ---

# Montastraea cavernosa (M. cavernosa)
p1_avg <- ggplot(mcav_avg, aes(x = condition, y = mean_abund, fill = Family)) +
  geom_bar(stat = "identity", position = "stack") +
  scale_fill_manual(values = family_palette, na.value = "grey90") +
  ggtitle("Montastraea cavernosa (mean by condition)") +
  labs(y = "Mean Relative Abundance",
       x = "Condition") +
  theme(strip.text = element_text(face = "bold"),
        plot.title = element_text(face = "italic"),
        legend.position = "none",
        axis.title.x = element_blank(),
        axis.text.x = element_blank())

# Porites astreoides
p2_avg <- ggplot(past_avg, aes(x = condition, y = mean_abund, fill = Family)) +
  geom_bar(stat = "identity", position = "stack") +
  scale_fill_manual(values = family_palette, na.value = "grey90") +
  ggtitle("Porites astreoides (mean by condition)") +
  labs(y = "Mean Relative Abundance",
       x = "Condition") +
  theme(strip.text = element_text(face = "bold"),
        plot.title = element_text(face = "italic"),
        legend.position = "none",
        axis.title.x = element_blank(),
        axis.text.x = element_blank())

# Acropora palmata
p3_avg <- ggplot(apal_avg, aes(x = condition, y = mean_abund, fill = Family)) +
  geom_bar(stat = "identity", position = "stack") +
  scale_fill_manual(values = family_palette, na.value = "grey90") +
  ggtitle("Acropora palmata (mean by condition)") +
  labs(y = "Mean Relative Abundance",
       x = "Condition") +
  theme(strip.text = element_text(face = "bold"),
        plot.title = element_text(face = "italic"),
        legend.position = "none",
        axis.title.x = element_blank(),
        axis.text.x = element_blank())


# --- Arrange all three plots ---
plot_grid(p3_avg, p1_avg, p2_avg, labels = c("A", "B", "C"), ncol = 2, nrow = 2)