################### 5. Calculation of microbial risk scores ###################
###5-1. Abundance-based microbial risk score (aMRS)
cvd_log2abundance_data <- read.delim("taxa_log2transformed_feature_table.txt", sep = "\t", row.names = 1, check.names = F)  ##row = SampleID, col = taxa
taxa_sub <- cvd_abundance_data[, ##significant abundance bacterial species list]
coef_sub <- data.frame(feature = ##significant abundance bacterial species list, coef = ##significant abundance bacterial species coefficient)

plot_long <- taxa_sub %>%
  tibble::rownames_to_column("SampleID") %>%
  pivot_longer(
    -SampleID,
    names_to = "feature",
    values_to = "abundance"
  )

merged <- plot_long %>%
  left_join(coef_sub, by = "feature") %>%
  mutate(
    weighted_value = abundance * coef
  )

aMRS <- merged %>%
  group_by(SampleID) %>%
  summarise(
    aMRS = sum(weighted_value, na.rm = TRUE),
    .groups = "drop"
  )


###5-2. Presence/absence–based microbial risk score (paMRS)
cvd_abundance_data <- read.delim("taxa_feature_table.txt", sep = "\t", row.names = 1, check.names = F)
taxa_sub <- cvd_abundance_data[, ##significant presence/absence bacterial species list]
taxa_binary <- taxa_sub %>% mutate(across(everything(), ~ ifelse(is.na(.), 0, 1)))
coef_sub <- data.frame(feature = ##significant presence/absence bacterial species list, coef = ##significant presence/absence bacterial species coefficient)

plot_long <- taxa_binary %>%
  tibble::rownames_to_column("SampleID") %>%
  pivot_longer(
    -SampleID,
    names_to = "feature",
    values_to = "prevalence"
  )

merged <- plot_long %>%
  left_join(coef_sub, by = "feature") %>%
  mutate(
    weighted_value = prevalence * coef
  )

paMRS <- merged %>%
  group_by(SampleID) %>%
  summarise(
    paMRS = sum(weighted_value, na.rm = TRUE),
    .groups = "drop"
  )
