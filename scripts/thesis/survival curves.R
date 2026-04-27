#SURVIVAL CURVES

library(readxl)
library(dplyr)
library(survival)
library(survminer)

# List sheets to confirm names
excel_sheets("~/Desktop/thesis/DATA/Outplant Data (Months 1-3).xlsx")

# Read both sheets and combine
survival_data <- read_excel("~/Desktop/thesis/DATA/Outplant Data (Months 1-3).xlsx", 
                            sheet = "survival 1-3 ") %>%
  bind_rows(
    read_excel("~/Desktop/thesis/DATA/Outplant Data (Months 4-6).xlsx", 
               sheet = "survival 4-6")
  ) %>%
  mutate(
    site = case_when(
      str_detect(`Block ID`, "^RR") ~ "Rupert's Rock",
      str_detect(`Block ID`, "^PV") ~ "Perseverance", 
      str_detect(`Block ID`, "^SB") ~ "Saba",
      str_detect(`Block ID`, "^BB") ~ "Brewer's Bay",
      str_detect(`Block ID`, "^KB") ~ "Krum Bay",
      TRUE ~ "Other"
    ),
    status = case_when(
      Survival == "Dead" ~ 1,
      Survival == "Alive" ~ 0,
      Survival == "Missing" ~ NA,
      TRUE ~ NA
    ),
    Timepoint = as.numeric(Timepoint)
  ) %>%
  filter(!is.na(status))

# Rest of your plotting loop stays identical
species_list <- unique(survival_data_123$Species)
for(sp in species_list) {
  cat("\n=== ", sp, " Survival by Site ===\n")
  sp_data <- filter(survival_data, Species == sp)
  fit <- survfit(Surv(Timepoint, status) ~ site, data = sp_data)
  
  p <- ggsurvplot(fit, data = sp_data, pval = TRUE, risk.table = TRUE,
                  risk.table.height = 0.25, break.time.by = 1,
                  xlab = "Timepoint (Months)", ylab = "Survival Probability",
                  title = paste(sp, "Survival by Site (Months 1-6)"),
                  palette = "Set1", legend.title = "Site",
                  legend.labs = c("BB", "KB", "PV", "RR", "SB"),
                  ggtheme = theme_bw())
  print(p)
  ggsave(paste0(sp, "_survival_by_site_full.png"), width = 10, height = 8)
}


########## COX REGRESSION MODEL #########

library(survival)

# Per species Cox regression (site as factor, first level = reference)
for(sp in species_list) {
  cat("\n=== ", sp, " CoxPH by Site ===\n")
  
  sp_data <- filter(survival_data, Species == sp)
  
  # Fit model: site reference = first alphabetical ("BB" typically)
  cox_model <- coxph(Surv(Timepoint, status) ~ site, data = sp_data)
  
  # Summary: hazard ratios, p-values, CIs
  print(summary(cox_model))
  
  # Save results table
  cox_table <- summary(cox_model)$coefficients %>%
    as.data.frame() %>%
    mutate(
      HR = exp(`coef`),  # Hazard ratio
      CI_low = exp(`coef` - 1.96*`se(coef)`),
      CI_high = exp(`coef` + 1.96*`se(coef)`),
      Site = rownames(.)
    )
  write.csv(cox_table, paste0(sp, "_coxph_results.csv"), row.names = FALSE)
}
