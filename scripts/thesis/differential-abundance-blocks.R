#differential abundance plot
library("corncob")
library("phyloseq")
library("dplyr")
library("ggplot2")
library("ggthemes")


ps<- readRDS("scripts/thesis/my_physeq_filtered.rds")

ps2 <- ps %>%
  phyloseq::subset_samples(SampleType == "Block")

# Create condition column based on your sites
sample_data(ps2)$condition <- dplyr::case_when(
  sample_data(ps2)$Site %in% c("Saba", "Perseverance") ~ "Unaffected",
  sample_data(ps2)$Site %in% c("Krum Bay", "Rupert's Rock") ~ "Compromised", 
  sample_data(ps2)$Site == "Brewer's Bay" ~ "Target")

#UNAFFECTED SITES VS BB

# Subset to just Saba/Perseverance + BB
ps_saba_bb <- ps2 %>% subset_samples(Site %in% c("Saba", "Perseverance", "Brewer's Bay"))

# Create binary factor (Saba/Perseverance as reference)
sample_data(ps_saba_bb)$group <- factor(ifelse(sample_data(ps_saba_bb)$Site %in% c("Saba", "Perseverance"), "S&P", "BB"))

# Check balance
table(sample_data(ps_saba_bb)$group)

set.seed(1)
saba_bb_test <- differentialTest(formula = ~ group, 
                                 phi.formula = ~ group,
                                 formula_null = ~ 1,
                                 phi.formula_null = ~ 1,  # Proper DV test
                                 test = "Wald", boot = FALSE,
                                 data = ps_saba_bb,
                                 fdr_cutoff = 0.05,
                                 verbose = TRUE)

# Results
saba_bb_test$significant_taxa  # 95 sig ASVs

#creating function to look at coefficients

sigtaxto_df <- function(daout, psobj) {
  if (length(daout$significant_models) == 0) return(data.frame())
  df.out <- do.call(rbind, lapply(daout$significant_models, function(m) m$coefficients[2, 1:4]))
  asv <- daout$significant_taxa
  tax_mat <- as(tax_table(psobj), "matrix")[asv, , drop = FALSE]
  df.bind <- data.frame(
    Estimate = as.numeric(df.out[,1]), StdError = as.numeric(df.out[,2]),
    tvalue = as.numeric(df.out[,3]), Pr = as.numeric(df.out[,4]),
    asv = asv,
    Class = ifelse(tax_mat[,"Class"] == "" | is.na(tax_mat[,"Class"]), "unclassified", tax_mat[,"Class"]),
    Family = ifelse(tax_mat[,"Family"] == "" | is.na(tax_mat[,"Family"]), "unclassified", tax_mat[,"Family"]),
    Genus = ifelse(tax_mat[,"Genus"] == "" | is.na(tax_mat[,"Genus"]), "unclassified", tax_mat[,"Genus"]),
    stringsAsFactors = FALSE
  )
  df.bind$asvtaxa <- paste0(df.bind$Genus, " (", df.bind$asv, ")")
  return(df.bind)
}

#using function
saba_bb_forgraph <- sigtaxto_df(saba_bb_test, ps_saba_bb)
head(saba_bb_forgraph)


# reorder the plot by decreasing coefficients

saba_bb_forgraph <- saba_bb_forgraph[order(-saba_bb_forgraph$Estimate), ]
saba_bb_forgraph$asvtaxa <- factor(saba_bb_forgraph$asvtaxa, levels = as.vector(saba_bb_forgraph$asvtaxa))

# Count your families and create exact palette
u_vs_bb_families <- length(unique(saba_bb_forgraph$Family)) #31
u_vs_bb_distinct_colors <- colorRampPalette(RColorBrewer::brewer.pal(8, "Dark2"))(u_vs_bb_families)

#graph

ggplot(saba_bb_forgraph, aes(x = asvtaxa, y = Estimate, fill = Family)) +
  geom_errorbar(aes(ymin = Estimate-StdError, ymax = Estimate+StdError), 
                color = "black", width = .3, position=position_dodge(.9)) +
  geom_point(size = 4, pch = 21, position = position_dodge(.9)) + 
  scale_x_discrete(limits = rev(levels(saba_bb_forgraph$asvtaxa))) +
  scale_alpha_manual(values = u_vs_bb_distinct_colors) + 
  coord_flip() +
  theme_bw() +
  labs(x = "Taxa", y = "Coefficient", title = "Blocks from Unaffected (Saba & Perseverance Bay) vs Target (Brewer's Bay) Locations") +
  geom_hline(yintercept = 0, linetype = "dashed", color = "darkgray")+
  theme(axis.text.y = element_text(size = 6),    # ← SMALLER TEXT
        axis.text.x = element_text(size = 10, angle = 0, hjust = 1))  # ← OPTIONAL: bottom rotation

#ggsave("blocks_diff_abund_unaffected_vs_target.png", width = 12, height = 10, dpi = 300)

#counting ASVs

bb_enriched <- sum(saba_bb_forgraph$Estimate > 0)    # BB = 8
sp_enriched   <- sum(saba_bb_forgraph$Estimate < 0)    # Saba/Perseverance = 87

c(SABA_PERS = sp_enriched, BB = bb_enriched)


####################### COMPROMISED SITES VS BB ###################### 

rr_krum_bb <- ps2 %>%
  subset_samples(Site %in% c("Krum Bay", "Rupert's Rock", "Brewer's Bay"))

sample_data(rr_krum_bb)$group <- factor(ifelse(sample_data(rr_krum_bb)$Site %in% c("Krum Bay", "Rupert's Rock"), "K&R", "BB"))

# Check balance
table(sample_data(rr_krum_bb)$group)


set.seed(1)
krum_bb_test <- differentialTest(formula = ~ group, 
                                 phi.formula = ~ group,
                                 formula_null = ~ 1,
                                 phi.formula_null = ~ 1,
                                 test = "Wald", boot = FALSE,
                                 data = rr_krum_bb,
                                 fdr_cutoff = 0.05,
                                 verbose = TRUE)

# Results
krum_bb_test$significant_taxa  #121

#using function
krum_bb_forgraph <- sigtaxto_df(krum_bb_test, rr_krum_bb)
head(krum_bb_forgraph)

# Count your families and create exact palette
c_vs_bb_families <- length(unique(krum_bb_forgraph$Family)) #37
c_vs_bb_distinct_colors <- colorRampPalette(RColorBrewer::brewer.pal(8, "Dark2"))(c_vs_bb_families)


#reorder the plot by decreasing coefficients

krum_bb_forgraph <- krum_bb_forgraph[order(-krum_bb_forgraph$Estimate), ]
krum_bb_forgraph$asvtaxa <- factor(krum_bb_forgraph$asvtaxa, levels = as.vector(krum_bb_forgraph$asvtaxa))

ggplot(krum_bb_forgraph, aes(x = asvtaxa, y = Estimate, fill = Family)) +
  geom_errorbar(aes(ymin = Estimate-StdError, ymax = Estimate+StdError), 
                color = "black", width = .3, position=position_dodge(.9)) +
  geom_point(size = 4, pch = 21, position = position_dodge(.9)) + 
  scale_x_discrete(limits = rev(levels(krum_bb_forgraph$asvtaxa))) +
  scale_fill_manual(values = c_vs_bb_distinct_colors) +  
  coord_flip() +
  theme_bw() +
  labs(x = "Taxa", y = "Coefficient", title = "Blocks from Compromised (Krum Bay & Rupert's Rock) vs Target (Brewer's Bay) Locations") +
  geom_hline(yintercept = 0, linetype = "dashed", color = "darkgray")+
  theme(axis.text.y = element_text(size = 6),   
        axis.text.x = element_text(size = 10, angle = 0, hjust = 1))  

ggsave("blocks_diff_abund_compromised_vs_target.png", width = 12, height = 10, dpi = 300)

bb_enriched2 <- sum(krum_bb_forgraph$Estimate > 0)    # BB = 13
krr_enriched   <- sum(krum_bb_forgraph$Estimate < 0)    # KRUM / RR = 108

c(KRUM_RR = krr_enriched, BB = bb_enriched2)


###################### COMPROMISED VS UNAFFECTED ###################### 

# Subset to Compromised + Unaffected sites only (exclude BB)
ps_comp_unaffected <- ps2 %>%
  subset_samples(Site %in% c("Krum Bay", "Rupert's Rock", "Saba", "Perseverance"))

# Binary factor: Compromised vs Unaffected  
sample_data(ps_comp_unaffected)$group <- factor(
  ifelse(sample_data(ps_comp_unaffected)$Site %in% c("Krum Bay", "Rupert's Rock"), 
         "Compromised", "Unaffected")
)

# Check balance
table(sample_data(ps_comp_unaffected)$group)

set.seed(1)
comp_vs_unaffected_test <- differentialTest(
                          formula = ~ group, 
                          phi.formula = ~ group,
                          formula_null = ~ 1,
                          phi.formula_null = ~ 1,
                          test = "Wald", boot = FALSE,
                          data = ps_comp_unaffected,
                          fdr_cutoff = 0.05,
                          verbose = TRUE
)

# Extract results
comp_vs_unaffected_test$significant_taxa  # 59
comp_vs_unaffected_forgraph <- sigtaxto_df(comp_vs_unaffected_test, ps_comp_unaffected)

#reorder the plot by decreasing coefficients

comp_vs_unaffected_forgraph <- comp_vs_unaffected_forgraph[order(-comp_vs_unaffected_forgraph$Estimate), ]
comp_vs_unaffected_forgraph$asvtaxa <- factor(comp_vs_unaffected_forgraph$asvtaxa, levels = as.vector(comp_vs_unaffected_forgraph$asvtaxa))

# Count your families and create exact palette
c_vs_u_families <- length(unique(comp_vs_unaffected_forgraph$Family))
c_vs_u_distinct_colors <- colorRampPalette(RColorBrewer::brewer.pal(8, "Dark2"))(c_vs_u_families)

# plot
ggplot(comp_vs_unaffected_forgraph, aes(x = asvtaxa, y = Estimate, fill = Family)) +
  geom_errorbar(aes(ymin = Estimate-StdError, ymax = Estimate+StdError), 
                color = "black", width = .3, position=position_dodge(.9)) +
  geom_point(size = 4, pch = 21, position = position_dodge(.9)) + 
  scale_x_discrete(limits = rev(levels(comp_vs_unaffected_forgraph$asvtaxa))) +
  scale_fill_manual(values = c_vs_u_distinct_colors) +  
  coord_flip() +
  theme_bw() +
  labs(x = "Taxa", y = "Logit Coefficient", title = "Blocks from Compromised (Krum Bay & Rupert's Rock) vs Unaffected (Saba & Perseverance) Locations") +
  geom_hline(yintercept = 0, linetype = "dashed", color = "darkgray") +
  theme(axis.text.y = element_text(size = 6),
        legend.position = "right")

ggsave("blocks_diff_abund_compromised_vs_unaffected.png", width = 12, height = 10, dpi = 300)

#counting ASVs

comp_enriched <- sum(comp_vs_unaffected_forgraph$Estimate > 0) #Compromised = 
unaffected_enriched   <- sum(comp_vs_unaffected_forgraph$Estimate < 0) #Unaffected = 

c(T2_Saba = t2_saba_enriched, T2_BB = t2_bb_enriched)


####################### SABA VS BB - T1 vs T2 ###################### 

#fixing dates
sample_data(ps2)$SamplingDate_fixed <- case_when(
  sample_data(ps2)$SamplingDate == "Aug-14-2025"      ~ "2025-08-14",
  sample_data(ps2)$SamplingDate == "August-14-2025"   ~ "2025-08-14",
  sample_data(ps2)$SamplingDate == "July-14-2025"     ~ "2025-07-14",
  sample_data(ps2)$SamplingDate == "July-15-2025"     ~ "2025-07-15", 
  sample_data(ps2)$SamplingDate == "July-16-2025"     ~ "2025-07-16",
  TRUE ~ NA_character_
)

#creating timepoints
sample_data(ps2)$Timepoint <- case_when(
  sample_data(ps2)$SamplingDate_fixed %in% c("2025-07-14", "2025-07-15", "2025-07-16") ~ "T1 (July 2025)",
  sample_data(ps2)$SamplingDate_fixed %in% c("2025-08-14") ~ "T2 (August 2025)",
  TRUE ~ NA_character_
)

# Verify it worked
table(sample_data(ps2)$Site, sample_data(ps2)$Timepoint)


# T1: Saba vs BB (July 2025)
saba_vs_bb_t1 <- ps2 %>% 
  subset_samples(Timepoint == "T1 (July 2025)" & 
                   Site %in% c("Saba", "Brewer's Bay"))
sample_data(saba_vs_bb_t1)$group <- factor(
  ifelse(sample_data(saba_vs_bb_t1)$Site %in% c("Saba"), "Saba", "BB")
)

set.seed(1)
saba_vs_bb_t1_test <- differentialTest(formula = ~ group, 
                            phi.formula = ~ group,
                            formula_null = ~ 1, 
                            phi.formula_null = ~ 1,
                            test = "Wald", 
                            boot = FALSE, 
                            data = saba_vs_bb_t1,
                            fdr_cutoff = 0.05,
                            verbose = TRUE)

# T2: Saba vs BB (August 2025)  
saba_vs_bb_t2 <- ps2 %>% 
  subset_samples(Timepoint == "T2 (August 2025)" & 
                   Site %in% c("Saba", "Brewer's Bay"))
sample_data(saba_vs_bb_t2)$group <- factor(
  ifelse(sample_data(saba_vs_bb_t2)$Site %in% c("Saba"), "Saba", "BB")
)

set.seed(1)
saba_vs_bb_t2_test <- differentialTest(formula = ~ group, 
                                      phi.formula = ~ group,
                                      formula_null = ~ 1, 
                                      phi.formula_null = ~ 1,
                                      test = "Wald", boot = FALSE, 
                                      data = saba_vs_bb_t2,
                                      fdr_cutoff = 0.05,
                                      verbose = TRUE)


saba_vs_bb_t1_test$significant_taxa #79
saba_vs_bb_t1_forgraph <- sigtaxto_df(saba_vs_bb_t1_test, saba_vs_bb_t1)

saba_vs_bb_t1_forgraph <- saba_vs_bb_t1_forgraph[order(-saba_vs_bb_t1_forgraph$Estimate), ]
saba_vs_bb_t1_forgraph$asvtaxa <- factor(saba_vs_bb_t1_forgraph$asvtaxa, levels = as.vector(saba_vs_bb_t1_forgraph$asvtaxa))

# Count your families and create exact palette
s_vs_b_t1_families <- length(unique(saba_vs_bb_t1_forgraph$Family)) #33
s_vs_b_t1_distinct_colors <- colorRampPalette(RColorBrewer::brewer.pal(8, "Dark2"))(s_vs_b_t1_families)

# plot
ggplot(saba_vs_bb_t1_forgraph, aes(x = asvtaxa, y = Estimate, fill = Family)) +
  geom_errorbar(aes(ymin = Estimate-StdError, ymax = Estimate+StdError), 
                color = "black", width = .3, position=position_dodge(.9)) +
  geom_point(size = 4, pch = 21, position = position_dodge(.9)) + 
  scale_x_discrete(limits = rev(levels(saba_vs_bb_t1_forgraph$asvtaxa))) +
  scale_fill_manual(values = s_vs_b_t1_distinct_colors) +  
  coord_flip() +
  theme_bw() +
  labs(x = "Taxa", y = "Logit Coefficient", title = "Blocks from Saba (Unimpaired) vs. Brewer's Bay (Target) at T1") +
  geom_hline(yintercept = 0, linetype = "dashed", color = "darkgray") +
  theme(axis.text.y = element_text(size = 6),
        legend.position = "right")

ggsave("blocks_diff_abund_saba_vs_bbT1.png", width = 12, height = 10, dpi = 300)

saba_vs_bb_t2_test$significant_taxa #25
saba_vs_bb_t2_forgraph <- sigtaxto_df(saba_vs_bb_t2_test, saba_vs_bb_t2)

saba_vs_bb_t2_forgraph <- saba_vs_bb_t2_forgraph[order(-saba_vs_bb_t2_forgraph$Estimate), ]
saba_vs_bb_t2_forgraph$asvtaxa <- factor(saba_vs_bb_t2_forgraph$asvtaxa, levels = as.vector(saba_vs_bb_t2_forgraph$asvtaxa))

# Count your families and create exact palette
s_vs_b_t2_families <- length(unique(saba_vs_bb_t2_forgraph$Family)) #11
s_vs_b_t2_distinct_colors <- colorRampPalette(RColorBrewer::brewer.pal(8, "Dark2"))(s_vs_b_t2_families)

# plot
ggplot(saba_vs_bb_t2_forgraph, aes(x = asvtaxa, y = Estimate, fill = Family)) +
  geom_errorbar(aes(ymin = Estimate-StdError, ymax = Estimate+StdError), 
                color = "black", width = .3, position=position_dodge(.9)) +
  geom_point(size = 4, pch = 21, position = position_dodge(.9)) + 
  scale_x_discrete(limits = rev(levels(saba_vs_bb_t2_forgraph$asvtaxa))) +
  scale_fill_manual(values = s_vs_b_t2_distinct_colors) +  
  coord_flip() +
  theme_bw() +
  labs(x = "Taxa", y = "Logit Coefficient", title = "Blocks from Saba (Unimpaired) vs. Brewer's Bay (Target) at T2") +
  geom_hline(yintercept = 0, linetype = "dashed", color = "darkgray") +
  theme(axis.text.y = element_text(size = 6),
        legend.position = "right")

ggsave("blocks_diff_abund_saba_vs_bbT2.png", width = 12, height = 10, dpi = 300)


#counting ASVs

t1_saba_enriched <- sum(saba_vs_bb_t1_forgraph$Estimate < 0)    # Saba = 54
t1_bb_enriched_svsbb   <- sum(saba_vs_bb_t1_forgraph$Estimate > 0)    # BB = 25

c(T1_Saba = t1_saba_enriched, T1_BB = t1_bb_enriched_svsbb)


t2_saba_enriched <- sum(saba_vs_bb_t2_forgraph$Estimate < 0) #Saba = 15
t2_bb_enriched_svsbb   <- sum(saba_vs_bb_t2_forgraph$Estimate > 0) #BB = 10

c(T2_Saba = t2_saba_enriched, T2_BB = t2_bb_enriched_svsbb)

####################### PERSEVERANCE VS BB - T1 vs T2 ###################### 

#fixing dates
sample_data(ps2)$SamplingDate_fixed <- case_when(
  sample_data(ps2)$SamplingDate == "Aug-14-2025"      ~ "2025-08-14",
  sample_data(ps2)$SamplingDate == "August-14-2025"   ~ "2025-08-14",
  sample_data(ps2)$SamplingDate == "July-14-2025"     ~ "2025-07-14",
  sample_data(ps2)$SamplingDate == "July-15-2025"     ~ "2025-07-15", 
  sample_data(ps2)$SamplingDate == "July-16-2025"     ~ "2025-07-16",
  TRUE ~ NA_character_
)

#creating timepoints
sample_data(ps2)$Timepoint <- case_when(
  sample_data(ps2)$SamplingDate_fixed %in% c("2025-07-14", "2025-07-15", "2025-07-16") ~ "T1 (July 2025)",
  sample_data(ps2)$SamplingDate_fixed %in% c("2025-08-14") ~ "T2 (August 2025)",
  TRUE ~ NA_character_
)

# Verify it worked
table(sample_data(ps2)$Site, sample_data(ps2)$Timepoint)


# T1: Perseverance vs BB (July 2025)
pers_vs_bb_t1 <- ps2 %>% 
  subset_samples(Timepoint == "T1 (July 2025)" & 
                   Site %in% c("Perseverance", "Brewer's Bay"))
sample_data(pers_vs_bb_t1)$group <- factor(
  ifelse(sample_data(pers_vs_bb_t1)$Site %in% c("Perseverance"), "Perseverance", "BB")
)

set.seed(1)
pers_vs_bb_t1_test <- differentialTest(formula = ~ group, 
                                       phi.formula = ~ group,
                                       formula_null = ~ 1, 
                                       phi.formula_null = ~ 1,
                                       test = "Wald", 
                                       boot = FALSE, 
                                       data = pers_vs_bb_t1,
                                       fdr_cutoff = 0.05,
                                       verbose = TRUE)

# T2: Perseverance vs BB (August 2025)  
pers_vs_bb_t2 <- ps2 %>% 
  subset_samples(Timepoint == "T2 (August 2025)" & 
                   Site %in% c("Perseverance", "Brewer's Bay"))
sample_data(pers_vs_bb_t2)$group <- factor(
  ifelse(sample_data(pers_vs_bb_t2)$Site %in% c("Perseverance"), "Perseverance", "BB")
)

set.seed(1)
pers_vs_bb_t2_test <- differentialTest(formula = ~ group, 
                                       phi.formula = ~ group,
                                       formula_null = ~ 1, 
                                       phi.formula_null = ~ 1,
                                       test = "Wald", 
                                       boot = FALSE, 
                                       data = pers_vs_bb_t2,
                                       fdr_cutoff = 0.05,
                                       verbose = TRUE)

#counting taxa and creating forgraph
pers_vs_bb_t1_test$significant_taxa #58
pers_vs_bb_t1_forgraph <- sigtaxto_df(pers_vs_bb_t1_test, pers_vs_bb_t1)

pers_vs_bb_t2_test$significant_taxa #19
pers_vs_bb_t2_forgraph <- sigtaxto_df(pers_vs_bb_t2_test, pers_vs_bb_t2)

#reorder
pers_vs_bb_t1_forgraph <- pers_vs_bb_t1_forgraph[order(-pers_vs_bb_t1_forgraph$Estimate), ]
pers_vs_bb_t1_forgraph$asvtaxa <- factor(pers_vs_bb_t1_forgraph$asvtaxa, levels = as.vector(pers_vs_bb_t1_forgraph$asvtaxa))

pers_vs_bb_t2_forgraph <- pers_vs_bb_t2_forgraph[order(-pers_vs_bb_t2_forgraph$Estimate), ]
pers_vs_bb_t2_forgraph$asvtaxa <- factor(pers_vs_bb_t2_forgraph$asvtaxa, levels = as.vector(pers_vs_bb_t2_forgraph$asvtaxa))


# Count your families and create exact palette
p_vs_b_t1_families <- length(unique(pers_vs_bb_t1_forgraph$Family)) #26
p_vs_b_t1_distinct_colors <- colorRampPalette(RColorBrewer::brewer.pal(8, "Dark2"))(p_vs_b_t1_families)

p_vs_b_t2_families <- length(unique(pers_vs_bb_t2_forgraph$Family)) #8
p_vs_b_t2_distinct_colors <- colorRampPalette(RColorBrewer::brewer.pal(8, "Dark2"))(p_vs_b_t2_families)

# plot pers vs bb at T1
ggplot(pers_vs_bb_t1_forgraph, aes(x = asvtaxa, y = Estimate, fill = Family)) +
  geom_errorbar(aes(ymin = Estimate-StdError, ymax = Estimate+StdError), 
                color = "black", width = .3, position=position_dodge(.9)) +
  geom_point(size = 4, pch = 21, position = position_dodge(.9)) + 
  scale_x_discrete(limits = rev(levels(pers_vs_bb_t1_forgraph$asvtaxa))) +
  scale_fill_manual(values = p_vs_b_t1_distinct_colors) +  
  coord_flip() +
  theme_bw() +
  labs(x = "Taxa", y = "Logit Coefficient", title = "Blocks from Perseverance (Unimpaired) vs. Brewer's Bay (Target) at T1") +
  geom_hline(yintercept = 0, linetype = "dashed", color = "darkgray") +
  theme(axis.text.y = element_text(size = 6),
        legend.position = "right")

ggsave("blocks_diff_abund_pers_vs_bbT1.png", width = 12, height = 10, dpi = 300)

# plot pers at T2
ggplot(pers_vs_bb_t2_forgraph, aes(x = asvtaxa, y = Estimate, fill = Family)) +
  geom_errorbar(aes(ymin = Estimate-StdError, ymax = Estimate+StdError), 
                color = "black", width = .3, position=position_dodge(.9)) +
  geom_point(size = 4, pch = 21, position = position_dodge(.9)) + 
  scale_x_discrete(limits = rev(levels(pers_vs_bb_t2_forgraph$asvtaxa))) +
  scale_fill_manual(values = p_vs_b_t2_distinct_colors) +  
  coord_flip() +
  theme_bw() +
  labs(x = "Taxa", y = "Logit Coefficient", title = "Blocks from Perseverance (Unimpaired) vs. Brewer's Bay (Target) at T2") +
  geom_hline(yintercept = 0, linetype = "dashed", color = "darkgray") +
  theme(axis.text.y = element_text(size = 6),
        legend.position = "right")

ggsave("blocks_diff_abund_pers_vs_bbT2.png", width = 12, height = 10, dpi = 300)


#counting ASVs

t1_pers_enriched <- sum(pers_vs_bb_t1_forgraph$Estimate < 0)    # Perseverance = 36
t1_bb_enriched_pvsbb   <- sum(pers_vs_bb_t1_forgraph$Estimate > 0)    # BB = 22

c(T1_Pers = t1_pers_enriched, T1_BB = t1_bb_enriched_pvsbb)

t2_pers_enriched <- sum(pers_vs_bb_t2_forgraph$Estimate < 0) #Perseverance = 14
t2_bb_enriched_pvsbb   <- sum(pers_vs_bb_t2_forgraph$Estimate > 0) #BB = 5

c(T2_Pers = t2_pers_enriched, T2_BB = t2_bb_enriched_pvsbb)


######################## KRUM BAY VS BREWERS T1 and T2 #########################

# T1: Krum Bay vs BB (July 2025)


krum_vs_bb_t1 <- ps2 %>% 
  subset_samples(Timepoint == "T1 (July 2025)" & 
                   Site %in% c("Krum Bay", "Brewer's Bay"))
sample_data(krum_vs_bb_t1)$group <- factor(
  ifelse(sample_data(krum_vs_bb_t1)$Site %in% c("Krum Bay"), "Krum Bay", "BB")
)

set.seed(1)
krum_vs_bb_t1_test <- differentialTest(formula = ~ group, 
                                       phi.formula = ~ group,
                                       formula_null = ~ 1, 
                                       phi.formula_null = ~ 1,
                                       test = "Wald", 
                                       boot = FALSE, 
                                       data = krum_vs_bb_t1,
                                       fdr_cutoff = 0.05,
                                       verbose = TRUE)

# T2: Krum Bay vs BB (August 2025)  
krum_vs_bb_t2 <- ps2 %>% 
  subset_samples(Timepoint == "T2 (August 2025)" & 
                   Site %in% c("Krum Bay", "Brewer's Bay"))
sample_data(krum_vs_bb_t2)$group <- factor(
  ifelse(sample_data(krum_vs_bb_t2)$Site %in% c("Krum Bay"), "Krum Bay", "BB")
)

set.seed(1)
krum_vs_bb_t2_test <- differentialTest(formula = ~ group, 
                                       phi.formula = ~ group,
                                       formula_null = ~ 1, 
                                       phi.formula_null = ~ 1,
                                       test = "Wald", 
                                       boot = FALSE, 
                                       data = krum_vs_bb_t2,
                                       fdr_cutoff = 0.05,
                                       verbose = TRUE)

#counting taxa and creating forgraph for both timepoints

krum_vs_bb_t1_test$significant_taxa 
krum_vs_bb_t1_forgraph <- sigtaxto_df(krum_vs_bb_t1_test, krum_vs_bb_t1)

krum_vs_bb_t2_test$significant_taxa #19
krum_vs_bb_t2_forgraph <- sigtaxto_df(krum_vs_bb_t2_test, krum_vs_bb_t2)

#reorder
krum_vs_bb_t1_forgraph <- krum_vs_bb_t1_forgraph[order(-krum_vs_bb_t1_forgraph$Estimate), ]
krum_vs_bb_t1_forgraph$asvtaxa <- factor(krum_vs_bb_t1_forgraph$asvtaxa, levels = as.vector(krum_vs_bb_t1_forgraph$asvtaxa))

krum_vs_bb_t2_forgraph <- krum_vs_bb_t2_forgraph[order(-krum_vs_bb_t2_forgraph$Estimate), ]
krum_vs_bb_t2_forgraph$asvtaxa <- factor(krum_vs_bb_t2_forgraph$asvtaxa, levels = as.vector(krum_vs_bb_t2_forgraph$asvtaxa))


# Count your families and create exact palette
k_vs_b_t1_families <- length(unique(krum_vs_bb_t1_forgraph$Family)) #42
k_vs_b_t1_distinct_colors <- colorRampPalette(RColorBrewer::brewer.pal(8, "Dark2"))(k_vs_b_t1_families)

k_vs_b_t2_families <- length(unique(krum_vs_bb_t2_forgraph$Family)) #9
k_vs_b_t2_distinct_colors <- colorRampPalette(RColorBrewer::brewer.pal(8, "Dark2"))(k_vs_b_t2_families)

# plot krum vs bb at T1
ggplot(krum_vs_bb_t1_forgraph, aes(x = asvtaxa, y = Estimate, fill = Family)) +
  geom_errorbar(aes(ymin = Estimate-StdError, ymax = Estimate+StdError), 
                color = "black", width = .3, position=position_dodge(.9)) +
  geom_point(size = 4, pch = 21, position = position_dodge(.9)) + 
  scale_x_discrete(limits = rev(levels(krum_vs_bb_t1_forgraph$asvtaxa))) +
  scale_fill_manual(values = k_vs_b_t1_distinct_colors) +  
  coord_flip() +
  theme_bw() +
  labs(x = "Taxa", y = "Logit Coefficient", title = "Blocks from Krum Bay (Impaired) vs. Brewer's Bay (Target) at T1") +
  geom_hline(yintercept = 0, linetype = "dashed", color = "darkgray") +
  theme(axis.text.y = element_text(size = 6),
        legend.position = "right")

ggsave("blocks_diff_abund_krum_vs_bbT1.png", width = 12, height = 10, dpi = 300)

# plot krum vs bb at T2
ggplot(krum_vs_bb_t2_forgraph, aes(x = asvtaxa, y = Estimate, fill = Family)) +
  geom_errorbar(aes(ymin = Estimate-StdError, ymax = Estimate+StdError), 
                color = "black", width = .3, position=position_dodge(.9)) +
  geom_point(size = 4, pch = 21, position = position_dodge(.9)) + 
  scale_x_discrete(limits = rev(levels(krum_vs_bb_t2_forgraph$asvtaxa))) +
  scale_fill_manual(values = k_vs_b_t2_distinct_colors) +  
  coord_flip() +
  theme_bw() +
  labs(x = "Taxa", y = "Logit Coefficient", title = "Blocks from Krum Bay (Impaired) vs. Brewer's Bay (Target) at T2") +
  geom_hline(yintercept = 0, linetype = "dashed", color = "darkgray") +
  theme(axis.text.y = element_text(size = 6),
        legend.position = "right")

ggsave("blocks_diff_abund_krum_vs_bbT2.png", width = 12, height = 10, dpi = 300)


#counting ASVs

t1_krum_enriched <- sum(krum_vs_bb_t1_forgraph$Estimate < 0)    # krum = 58
bb_enriched3   <- sum(krum_vs_bb_t1_forgraph$Estimate > 0)    # BB = 68

c(T1_KRUM = t1_krum_enriched, T1_BB = bb_enriched3)

t2_krum_enriched <- sum(krum_vs_bb_t2_forgraph$Estimate < 0) #Krum = 16
t2_bb_enriched3   <- sum(krum_vs_bb_t2_forgraph$Estimate > 0) #BB = 4

c(T2_KRUM = t2_krum_enriched, T2_BB = t2_bb_enriched3)

######################## RR VS BREWERS T1 and T2 #########################

# T1: Krum Bay vs BB (July 2025)

rr_vs_bb_t1 <- ps2 %>% 
  subset_samples(Timepoint == "T1 (July 2025)" & 
                   Site %in% c("Rupert's Rock", "Brewer's Bay"))
sample_data(rr_vs_bb_t1)$group <- factor(
  ifelse(sample_data(rr_vs_bb_t1)$Site %in% c("Rupert's Rock"), "RR", "BB")
)

set.seed(1)
rr_vs_bb_t1_test <- differentialTest(formula = ~ group, 
                                       phi.formula = ~ group,
                                       formula_null = ~ 1, 
                                       phi.formula_null = ~ 1,
                                       test = "Wald", 
                                       boot = FALSE, 
                                       data = rr_vs_bb_t1,
                                       fdr_cutoff = 0.05,
                                       verbose = TRUE)

# T2: RR vs BB (August 2025)  
rr_vs_bb_t2 <- ps2 %>% 
  subset_samples(Timepoint == "T2 (August 2025)" & 
                   Site %in% c("Rupert's Rock", "Brewer's Bay"))
sample_data(rr_vs_bb_t2)$group <- factor(
  ifelse(sample_data(rr_vs_bb_t2)$Site %in% c("Rupert's Rock"), "RR", "BB")
)

set.seed(1)
rr_vs_bb_t2_test <- differentialTest(formula = ~ group, 
                                       phi.formula = ~ group,
                                       formula_null = ~ 1, 
                                       phi.formula_null = ~ 1,
                                       test = "Wald", 
                                       boot = FALSE, 
                                       data = rr_vs_bb_t2,
                                       fdr_cutoff = 0.05,
                                       verbose = TRUE)

#counting taxa and creating forgraph for both timepoints

rr_vs_bb_t1_test$significant_taxa #67
rr_vs_bb_t1_forgraph <- sigtaxto_df(rr_vs_bb_t1_test, rr_vs_bb_t1)

rr_vs_bb_t2_test$significant_taxa #21
rr_vs_bb_t2_forgraph <- sigtaxto_df(rr_vs_bb_t2_test, rr_vs_bb_t2)

#reorder
rr_vs_bb_t1_forgraph <- rr_vs_bb_t1_forgraph[order(-rr_vs_bb_t1_forgraph$Estimate), ]
rr_vs_bb_t1_forgraph$asvtaxa <- factor(rr_vs_bb_t1_forgraph$asvtaxa, levels = as.vector(rr_vs_bb_t1_forgraph$asvtaxa))

rr_vs_bb_t2_forgraph <- rr_vs_bb_t2_forgraph[order(-rr_vs_bb_t2_forgraph$Estimate), ]
rr_vs_bb_t2_forgraph$asvtaxa <- factor(rr_vs_bb_t2_forgraph$asvtaxa, levels = as.vector(rr_vs_bb_t2_forgraph$asvtaxa))


# Count your families and create exact palette
rr_vs_b_t1_families <- length(unique(rr_vs_bb_t1_forgraph$Family)) #26
rr_vs_b_t1_distinct_colors <- colorRampPalette(RColorBrewer::brewer.pal(8, "Dark2"))(rr_vs_b_t1_families)

rr_vs_b_t2_families <- length(unique(rr_vs_bb_t2_forgraph$Family)) #10
rr_vs_b_t2_distinct_colors <- colorRampPalette(RColorBrewer::brewer.pal(8, "Dark2"))(rr_vs_b_t2_families)

# plot rr vs bb at T1
ggplot(rr_vs_bb_t1_forgraph, aes(x = asvtaxa, y = Estimate, fill = Family)) +
  geom_errorbar(aes(ymin = Estimate-StdError, ymax = Estimate+StdError), 
                color = "black", width = .3, position=position_dodge(.9)) +
  geom_point(size = 4, pch = 21, position = position_dodge(.9)) + 
  scale_x_discrete(limits = rev(levels(rr_vs_bb_t1_forgraph$asvtaxa))) +
  scale_fill_manual(values = rr_vs_b_t1_distinct_colors) +  
  coord_flip() +
  theme_bw() +
  labs(x = "Taxa", y = "Logit Coefficient", title = "Blocks from Rupert's Rock (Impaired) vs. Brewer's Bay (Target) at T1") +
  geom_hline(yintercept = 0, linetype = "dashed", color = "darkgray") +
  theme(axis.text.y = element_text(size = 6),
        legend.position = "right")

ggsave("blocks_diff_abund_rr_vs_bbT1.png", width = 12, height = 10, dpi = 300)

# plot rr vs bb at T2
ggplot(rr_vs_bb_t2_forgraph, aes(x = asvtaxa, y = Estimate, fill = Family)) +
  geom_errorbar(aes(ymin = Estimate-StdError, ymax = Estimate+StdError), 
                color = "black", width = .3, position=position_dodge(.9)) +
  geom_point(size = 4, pch = 21, position = position_dodge(.9)) + 
  scale_x_discrete(limits = rev(levels(rr_vs_bb_t2_forgraph$asvtaxa))) +
  scale_fill_manual(values = rr_vs_b_t2_distinct_colors) +  
  coord_flip() +
  theme_bw() +
  labs(x = "Taxa", y = "Logit Coefficient", title = "Blocks from Rupert's Rock (Impaired) vs. Brewer's Bay (Target) at T2") +
  geom_hline(yintercept = 0, linetype = "dashed", color = "darkgray") +
  theme(axis.text.y = element_text(size = 6),
        legend.position = "right")

ggsave("blocks_diff_abund_rr_vs_bbT2.png", width = 12, height = 10, dpi = 300)


#counting ASVs

t1_rr_enriched <- sum(rr_vs_bb_t1_forgraph$Estimate < 0)    # rr = 44
bb_enriched4   <- sum(rr_vs_bb_t1_forgraph$Estimate > 0)    # BB = 23

c(T1_RR = t1_rr_enriched, T1_BB = bb_enriched4)

t2_rr_enriched <- sum(rr_vs_bb_t2_forgraph$Estimate < 0) # rr = 17
t2_bb_enriched4  <- sum(rr_vs_bb_t2_forgraph$Estimate > 0) #BB = 4

c(T2_RR = t2_rr_enriched, T2_BB = t2_bb_enriched4)


######################### CREATING ASV TABLE ####################################


# Full table with TOP ASV for BOTH sites in each comparison
site_vs_bb_table <- data.frame(
  Comparison = c("Saba vs BB (T1)", "Saba vs BB (T2)", 
                 "Perseverance vs BB (T1)", "Perseverance vs BB (T2)",
                 "Krum Bay vs BB (T1)", "Krum Bay vs BB (T2)", 
                 "Rupert's Rock vs BB (T1)", "Rupert's Rock vs BB (T2)"),
  
  Total_Sig_ASVs = c(
    nrow(saba_vs_bb_t1_forgraph), nrow(saba_vs_bb_t2_forgraph),
    nrow(pers_vs_bb_t1_forgraph), nrow(pers_vs_bb_t2_forgraph),
    nrow(krum_vs_bb_t1_forgraph), nrow(krum_vs_bb_t2_forgraph),
    nrow(rr_vs_bb_t1_forgraph), nrow(rr_vs_bb_t2_forgraph)
  ),
  
  BB_Enriched = c(
    sum(saba_vs_bb_t1_forgraph$Estimate > 0), sum(saba_vs_bb_t2_forgraph$Estimate > 0),
    sum(pers_vs_bb_t1_forgraph$Estimate > 0), sum(pers_vs_bb_t2_forgraph$Estimate > 0),
    sum(krum_vs_bb_t1_forgraph$Estimate > 0), sum(krum_vs_bb_t2_forgraph$Estimate > 0),
    sum(rr_vs_bb_t1_forgraph$Estimate > 0), sum(rr_vs_bb_t2_forgraph$Estimate > 0)
  ),
  
  Site_Enriched = c(
    sum(saba_vs_bb_t1_forgraph$Estimate < 0), sum(saba_vs_bb_t2_forgraph$Estimate < 0),
    sum(pers_vs_bb_t1_forgraph$Estimate < 0), sum(pers_vs_bb_t2_forgraph$Estimate < 0),
    sum(krum_vs_bb_t1_forgraph$Estimate < 0), sum(krum_vs_bb_t2_forgraph$Estimate < 0),
    sum(rr_vs_bb_t1_forgraph$Estimate < 0), sum(rr_vs_bb_t2_forgraph$Estimate < 0)
  ),
  
  Top_BB_ASV = c(
    ifelse(sum(saba_vs_bb_t1_forgraph$Estimate > 0) > 0,
           as.character(saba_vs_bb_t1_forgraph$asvtaxa[which.max(saba_vs_bb_t1_forgraph$Estimate)]),
           "None"),
    ifelse(sum(saba_vs_bb_t2_forgraph$Estimate > 0) > 0,
           as.character(saba_vs_bb_t2_forgraph$asvtaxa[which.max(saba_vs_bb_t2_forgraph$Estimate)]),
           "None"),
    ifelse(sum(pers_vs_bb_t1_forgraph$Estimate > 0) > 0,
           as.character(pers_vs_bb_t1_forgraph$asvtaxa[which.max(pers_vs_bb_t1_forgraph$Estimate)]),
           "None"),
    ifelse(sum(pers_vs_bb_t2_forgraph$Estimate > 0) > 0,
           as.character(pers_vs_bb_t2_forgraph$asvtaxa[which.max(pers_vs_bb_t2_forgraph$Estimate)]),
           "None"),
    ifelse(sum(krum_vs_bb_t1_forgraph$Estimate > 0) > 0,
           as.character(krum_vs_bb_t1_forgraph$asvtaxa[which.max(krum_vs_bb_t1_forgraph$Estimate)]),
           "None"),
    ifelse(sum(krum_vs_bb_t2_forgraph$Estimate > 0) > 0,
           as.character(krum_vs_bb_t2_forgraph$asvtaxa[which.max(krum_vs_bb_t2_forgraph$Estimate)]),
           "None"),
    ifelse(sum(rr_vs_bb_t1_forgraph$Estimate > 0) > 0,
           as.character(rr_vs_bb_t1_forgraph$asvtaxa[which.max(rr_vs_bb_t1_forgraph$Estimate)]),
           "None"),
    ifelse(sum(rr_vs_bb_t2_forgraph$Estimate > 0) > 0,
           as.character(rr_vs_bb_t2_forgraph$asvtaxa[which.max(rr_vs_bb_t2_forgraph$Estimate)]),
           "None")
  ),
  
  Top_Site_ASV = c(
    ifelse(sum(saba_vs_bb_t1_forgraph$Estimate < 0) > 0,
           as.character(saba_vs_bb_t1_forgraph$asvtaxa[which.min(saba_vs_bb_t1_forgraph$Estimate)]),
           "None"),
    ifelse(sum(saba_vs_bb_t2_forgraph$Estimate < 0) > 0,
           as.character(saba_vs_bb_t2_forgraph$asvtaxa[which.min(saba_vs_bb_t2_forgraph$Estimate)]),
           "None"),
    ifelse(sum(pers_vs_bb_t1_forgraph$Estimate < 0) > 0,
           as.character(pers_vs_bb_t1_forgraph$asvtaxa[which.min(pers_vs_bb_t1_forgraph$Estimate)]),
           "None"),
    ifelse(sum(pers_vs_bb_t2_forgraph$Estimate < 0) > 0,
           as.character(pers_vs_bb_t2_forgraph$asvtaxa[which.min(pers_vs_bb_t2_forgraph$Estimate)]),
           "None"),
    ifelse(sum(krum_vs_bb_t1_forgraph$Estimate < 0) > 0,
           as.character(krum_vs_bb_t1_forgraph$asvtaxa[which.min(krum_vs_bb_t1_forgraph$Estimate)]),
           "None"),
    ifelse(sum(krum_vs_bb_t2_forgraph$Estimate < 0) > 0,
           as.character(krum_vs_bb_t2_forgraph$asvtaxa[which.min(krum_vs_bb_t2_forgraph$Estimate)]),
           "None"),
    ifelse(sum(rr_vs_bb_t1_forgraph$Estimate < 0) > 0,
           as.character(rr_vs_bb_t1_forgraph$asvtaxa[which.min(rr_vs_bb_t1_forgraph$Estimate)]),
           "None"),
    ifelse(sum(rr_vs_bb_t2_forgraph$Estimate < 0) > 0,
           as.character(rr_vs_bb_t2_forgraph$asvtaxa[which.min(rr_vs_bb_t2_forgraph$Estimate)]),
           "None")
  )
)

# Display and save
print(site_vs_bb_table)
write_csv(site_vs_bb_table, "site_vs_bb_both_top_ASVs.csv")

# Top 3 ASVs for BOTH Site AND Brewers Bay per comparison
top3_table <- bind_rows(
  # Saba vs BB (T1)
  saba_vs_bb_t1_forgraph %>% 
    filter(Estimate > 0) %>% arrange(desc(Estimate)) %>% slice_head(n=3) %>% 
    mutate(Comparison = "Saba vs BB (T1)", Direction = "Brewers Bay"),
  saba_vs_bb_t1_forgraph %>% 
    filter(Estimate < 0) %>% arrange(Estimate) %>% slice_head(n=3) %>% 
    mutate(Comparison = "Saba vs BB (T1)", Direction = "Saba"),
  
  # Saba vs BB (T2)  
  saba_vs_bb_t2_forgraph %>% 
    filter(Estimate > 0) %>% arrange(desc(Estimate)) %>% slice_head(n=3) %>% 
    mutate(Comparison = "Saba vs BB (T2)", Direction = "Brewers Bay"),
  saba_vs_bb_t2_forgraph %>% 
    filter(Estimate < 0) %>% arrange(Estimate) %>% slice_head(n=3) %>% 
    mutate(Comparison = "Saba vs BB (T2)", Direction = "Saba"),
  
  # Perseverance vs BB (T1)
  pers_vs_bb_t1_forgraph %>% 
    filter(Estimate > 0) %>% arrange(desc(Estimate)) %>% slice_head(n=3) %>% 
    mutate(Comparison = "Perseverance vs BB (T1)", Direction = "Brewers Bay"),
  pers_vs_bb_t1_forgraph %>% 
    filter(Estimate < 0) %>% arrange(Estimate) %>% slice_head(n=3) %>% 
    mutate(Comparison = "Perseverance vs BB (T1)", Direction = "Perseverance"),
  
  # Perseverance vs BB (T2)
  pers_vs_bb_t2_forgraph %>% 
    filter(Estimate > 0) %>% arrange(desc(Estimate)) %>% slice_head(n=3) %>% 
    mutate(Comparison = "Perseverance vs BB (T2)", Direction = "Brewers Bay"),
  pers_vs_bb_t2_forgraph %>% 
    filter(Estimate < 0) %>% arrange(Estimate) %>% slice_head(n=3) %>% 
    mutate(Comparison = "Perseverance vs BB (T2)", Direction = "Perseverance"),
  
  # Krum Bay vs BB (T1)
  krum_vs_bb_t1_forgraph %>% 
    filter(Estimate > 0) %>% arrange(desc(Estimate)) %>% slice_head(n=3) %>% 
    mutate(Comparison = "Krum Bay vs BB (T1)", Direction = "Brewers Bay"),
  krum_vs_bb_t1_forgraph %>% 
    filter(Estimate < 0) %>% arrange(Estimate) %>% slice_head(n=3) %>% 
    mutate(Comparison = "Krum Bay vs BB (T1)", Direction = "Krum Bay"),
  
  # Krum Bay vs BB (T2)
  krum_vs_bb_t2_forgraph %>% 
    filter(Estimate > 0) %>% arrange(desc(Estimate)) %>% slice_head(n=3) %>% 
    mutate(Comparison = "Krum Bay vs BB (T2)", Direction = "Brewers Bay"),
  krum_vs_bb_t2_forgraph %>% 
    filter(Estimate < 0) %>% arrange(Estimate) %>% slice_head(n=3) %>% 
    mutate(Comparison = "Krum Bay vs BB (T2)", Direction = "Krum Bay"),
  
  # Rupert's Rock vs BB (T1)
  rr_vs_bb_t1_forgraph %>% 
    filter(Estimate > 0) %>% arrange(desc(Estimate)) %>% slice_head(n=3) %>% 
    mutate(Comparison = "Rupert's Rock vs BB (T1)", Direction = "Brewers Bay"),
  rr_vs_bb_t1_forgraph %>% 
    filter(Estimate < 0) %>% arrange(Estimate) %>% slice_head(n=3) %>% 
    mutate(Comparison = "Rupert's Rock vs BB (T1)", Direction = "Rupert's Rock"),
  
  # Rupert's Rock vs BB (T2)
  rr_vs_bb_t2_forgraph %>% 
    filter(Estimate > 0) %>% arrange(desc(Estimate)) %>% slice_head(n=3) %>% 
    mutate(Comparison = "Rupert's Rock vs BB (T2)", Direction = "Brewers Bay"),
  rr_vs_bb_t2_forgraph %>% 
    filter(Estimate < 0) %>% arrange(Estimate) %>% slice_head(n=3) %>% 
    mutate(Comparison = "Rupert's Rock vs BB (T2)", Direction = "Rupert's Rock")
) %>%
  select(Comparison, Direction, asvtaxa, Estimate, Pr, Family) %>%
  arrange(Comparison, Direction, desc(Estimate))

# Display and save
print(top3_table)
write_csv(top3_table, "top3_ASVs_by_direction_site_vs_bb.csv")