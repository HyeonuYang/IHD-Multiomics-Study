################### 4. Statistical analysis ###################
###4-1. Alpha diversity analysis
library(vegan)

cvd_metadata <- read.delim("metadata.txt", sep = "\t", row.names = 1, check.names = F)
cvd_abundance_data <- read.delim("taxa_feature_table.txt", sep = "\t", row.names = 1, check.names = F)

Shannon_alpha <- diversity(t(cvd_abundance_data), index = "shannon")
Richness_alpha <- specnumber(t(cvd_abundance_data))

Alpha_diversity_df <- data.frame(SampleID = names(Shannon_alpha), Shannon = Shannon_alpha, Richness = Richness_alpha, Group = cvd_metadata$Group, age = cvd_metadata$age, sex = cvd_metadata$sex, BMI = cvd_metadata$BMI, PPI = cvd_metadata$PPI)

summary(aov(lm(Shannon ~ age + sex + BMI + PPI + Group, data = Alpha_diversity_df)))
summary(aov(lm(Richness ~ age + sex + BMI + PPI + Group, data = Alpha_diversity_df)))


###4-2. Beta diversity analysis
library(vegan)

cvd_metadata <- read.delim("metadata.txt", sep = "\t", row.names = 1, check.names = F)
cvd_abundance_data <- read.delim("taxa_feature_table.txt", sep = "\t", row.names = 1, check.names = F)

## aitchison
pseudocount_half <- min(cvd_abundance_data[cvd_abundance_data > 0]) / 2

aitchison_dist <- vegdist(t(cvd_abundance_data), method = "aitchison", pseudocount = pseudocount_half)
aitchison_view <- as.matrix(aitchison_dist)

aitchison_per <- read.delim("aitchison_distance.txt")
meta_aitchison <- inner_join(cvd_metadata, aitchison_per, by="SampleID")
all_dist_aitchison <- meta_aitchison %>%
  select(all_of(.[["SampleID"]])) %>%
  as.dist()

## PERMANOVA aitchison
set.seed(12345)
permanova_aitchison <- adonis2(all_dist_aitchison ~ age + sex + BMI + PPI + Group, 
                                  data = meta_aitchison, 
                                  permutations = 999, by = "term")


## bray
bray_dist <- vegdist(t(cvd_abundance_data), method = "bray")
bray_view <- as.matrix(bray_dist)

bray_per <- read.delim("bray_distance.txt")
meta_bray <- inner_join(cvd_metadata, bray_per, by="SampleID")
all_dist_bray <- meta_bray %>%
  select(all_of(.[["SampleID"]])) %>%
  as.dist()

## PERMANOVA bray
set.seed(12345)
permanova_bray <- adonis2(all_dist_bray ~ age + sex + BMI + PPI + Group, 
                             data = meta_bray, 
                             permutations = 999, by = "term")


###4-3. Differential analysis for bacterial species, MetaCyc metabolic pathways, GMMs using MaAsLin3
##4-3-1. Ordinal differential analysis across CVD severity
library(maaslin3)
cvd_metadata <- read.delim("metadata.txt", sep = "\t", row.names = 1, check.names = F)
cvd_abundance_data <- read.delim("taxa_feature_table.txt", sep = "\t", row.names = 1, check.names = F)

cvd_metadata$Group_ord <- factor(
  cvd_metadata$Group,
  levels = c("SA", "ACS", "HFrEF"),
  ordered = TRUE
)

fit_data <- maaslin3(
  input_data      = as.data.frame(t(cvd_abundance_data)),
  input_metadata  = cvd_metadata,
  output          = "/path/to/output",
  formula         = "~ age + sex + BMI + PPI + Group_ord",    ### add statin or/and drinking for sensitivity analysis
  min_abundance   = ###, 		## 0.01 for bacterial species and 0 for functional pathways
  min_prevalence  = 0.20,
  normalization   = "NONE",
  transform       = "LOG",
  standardize     = TRUE,
  correction      = "BH",
  cores           = 7
)

##4-3-2. Pairwise differential analysis for ACS-specific bacterial species
library(maaslin3)
cvd_metadata <- read.delim("metadata.txt", sep = "\t", row.names = 1, check.names = F)
cvd_abundance_data <- read.delim("taxa_feature_table.txt", sep = "\t", row.names = 1, check.names = F)


# 1️⃣ SA reference
cvd_metadata$Group_refSA <- factor(
  cvd_metadata$Group,
  levels = c("SA", "ACS", "HFrEF")
)

fit_SARef <- maaslin3(
  input_data      = as.data.frame(t(cvd_abundance_data)),
  input_metadata  = cvd_metadata,
  output          = "/path/to/output",
  formula         = "~ age + sex + BMI + PPI + Group_refSA",
  min_abundance   = 0.01,
  min_prevalence  = 0.20,
  normalization   = "NONE",
  transform       = "LOG",
  standardize     = TRUE,
  correction      = "BH",
  cores           = 7
)

refSA_df <- read.delim("all_results_refSA.tsv", header = T, sep = "\t")
refSA_df_filtered <- refSA_df %>% filter(metadata == "Group_refSA") %>% filter(model == "abundance") %>% arrange(pval_individual) %>% select(feature, value, coef, pval_individual, model)

refSA_df_filtered <- refSA_df_filtered %>%
  mutate(
    value = case_when(
      value == "ACS" ~ "SA vs ACS",
      value == "HFrEF" ~ "SA vs HFrEF",
      TRUE ~ value
    )
  )

refSA_df_filtered_p <- refSA_df_filtered %>% filter(pval_individual < 0.05)


# 2️⃣ ACS reference
cvd_metadata$Group_refACS <- factor(
  cvd_metadata$Group,
  levels = c("ACS", "SA", "HFrEF")
)

fit_ACSRef <- maaslin3(
  input_data      = as.data.frame(t(cvd_abundance_data)),
  input_metadata  = cvd_metadata,
  output          = "/path/to/output",
  formula         = "~ age + sex + BMI + PPI + Group_refACS",
  min_abundance   = 0.01,
  min_prevalence  = 0.20,
  normalization   = "NONE",
  transform       = "LOG",
  standardize     = TRUE,
  correction      = "BH",
  cores           = 7
)

refACS_df <- read.delim("all_results_refACS.tsv", header = T, sep = "\t")
refACS_df_filtered <- refACS_df %>% filter(metadata == "Group_refACS") %>% filter(value == "HFrEF") %>% filter(model == "abundance") %>% arrange(pval_individual) %>% select(feature, value, coef, pval_individual, model)

refACS_df_filtered <- refACS_df_filtered %>%
  mutate(
    value = case_when(
      value == "HFrEF" ~ "ACS vs HFrEF",
      TRUE ~ value
    )
  )

refACS_df_filtered_p <- refACS_df_filtered %>% filter(pval_individual < 0.05)

feature_list <- union(refSA_df_filtered_p$feature, refACS_df_filtered_p$feature)

summary_df <- data.frame(
  feature = feature_list,
  p1 = NA_real_,
  p2 = NA_real_,
  p3 = NA_real_,
  stringsAsFactors = FALSE
)

tmp_SAACS <- refSA_df_filtered %>%
  filter(value == "SA vs ACS") %>%
  select(feature, coef, pval_individual)

summary_df <- summary_df %>%
  left_join(tmp_SAACS, by = "feature") %>%
  mutate(
    "SA vs ACS" = coef,
    p1 = pval_individual
  ) %>%
  select(-coef, -pval_individual)


tmp_SAHF <- refSA_df_filtered %>%
  filter(value == "SA vs HFrEF") %>%
  select(feature, coef, pval_individual)

summary_df <- summary_df %>%
  left_join(tmp_SAHF, by = "feature") %>%
  mutate(
    "SA vs HFrEF" = coef,
    p2 = pval_individual
  ) %>%
  select(-coef, -pval_individual)


tmp_ACSHF <- refACS_df_filtered %>%
  select(feature, coef, pval_individual)

summary_df <- summary_df %>%
  left_join(tmp_ACSHF, by = "feature") %>%
  mutate(
    "ACS vs HFrEF" = coef,
    p3 = pval_individual
  ) %>%
  select(-coef, -pval_individual)


summary_df <- summary_df %>% arrange(feature) %>% select(feature, "SA vs ACS", p1, "SA vs HFrEF", p2, "ACS vs HFrEF", p3)

ACS_specific_species <- summary_df %>% filter(!is.na(p1), !is.na(p3), p1 < 0.05, p3 < 0.05, ("SA vs ACS" * "ACS vs HFrEF") < 0)


###4-4. Functional bacterial group analysis
##4-4-1. Literature-defined functional groups
cvd_metadata <- read.delim("metadata.txt", sep = "\t", row.names = 1, check.names = F)
cvd_abundance_data <- read.delim("taxa_feature_table.txt", sep = "\t", row.names = 1, check.names = F)
taxa_functional_group <- read.delim("taxa_functional_group.txt", sep = "\t", check.names = F)
group_taxa <- readLines("taxa_functional_group.txt") %>% trimws()
group_taxa <- intersect(group_taxa, rownames(cvd_abundance_data))
taxa_sub <- cvd_abundance_data[group_taxa, , drop = FALSE]

group_sum <- colSums(taxa_sub, na.rm = TRUE)
group_df <- tibble(
  SampleID = names(group_sum),
  abundance = as.numeric(group_sum)
) %>%
  left_join(
    cvd_metadata %>% select(SampleID, Group),
    by = "SampleID"
  ) %>%
  mutate(
    Group = factor(Group, levels = c("SA", "ACS", "HFrEF")),
    log_abundance = log10(abundance)
  )

kw_res <- group_df %>% kruskal_test(log_abundance ~ Group)
dunn_res <- group_df %>% dunn_test(log_abundance ~ Group)

##4-4-2. Prevotella-to-Bacteroides ratio
cvd_metadata <- read.delim("metadata.txt", sep = "\t", row.names = 1, check.names = F)
cvd_genus_data <- read.delim("genus_feature_table.txt", sep = "\t", row.names = 1, check.names = F)
pb_genus <- cvd_genus_data[c("g__Prevotella", "g__Bacteroides"), , drop = FALSE]
pb_df <- as.data.frame(t(pb_genus))

min_positive <- min(pb_df[pb_df > 0], na.rm = TRUE)
pseudocount <- min_positive / 2

pb_df <- pb_df %>%
  mutate(
    SampleID = rownames(.),
    PB_ratio = (g__Prevotella + pseudocount) / (g__Bacteroides + pseudocount),
    log10_PB_ratio = log10(PB_ratio)
  ) %>%
  left_join(
    cvd_metadata %>%
      mutate(SampleID = rownames(cvd_metadata)) %>%
      select(SampleID, Group),
    by = "SampleID"
  )

kw_res <- pb_df %>% kruskal_test(log10_PB_ratio ~ Group)
dunn_res <- pb_df %>% dunn_test(log10_PB_ratio ~ Group)


###4-5. Spearman correlation of bacterial species and functional pathways
cvd_abundance_data <- read.delim("taxa_feature_table.txt", sep = "\t", row.names = 1, check.names = F)
taxa_sub <- cvd_abundance_data[##significant abundance bacterial species list,]

cvd_metacyc_data <- read.delim("MetaCyc_abundance.txt", sep = "\t", row.names = 1, check.names = F)
metacyc_sub <- cvd_abundance_data[##significant MetaCyc list,]

cvd_gmm_data <- read.delim("GMM_abundance.txt", sep = "\t", row.names = 1, check.names = F)
gmm_sub <- cvd_abundance_data[##significant GMM list,]

taxa_mat    <- t(taxa_sub)
pathway_mat <- cbind(t(metacyc_sub), t(gmm_sub))

cor_mat <- matrix(NA, nrow = ncol(taxa_mat), ncol = ncol(pathway_mat), dimnames = list(colnames(taxa_mat), colnames(pathway_mat)))

p_mat <- cor_mat

for (i in seq_len(ncol(taxa_mat))) {
  for (j in seq_len(ncol(pathway_mat))) {
    test <- cor.test(taxa_mat[, i], pathway_mat[, j], method = "spearman")
    cor_mat[i, j] <- test$estimate
    p_mat[i, j]   <- test$p.value
  }
}


###4-6. Gene-level analysis of metabolite-producing pathways
cvd_metadata <- read.delim("metadata.txt", sep = "\t", row.names = 1, check.names = F)
cvd_EC_data <- read.delim("EC_feature_table.txt", sep = "\t", row.names = 1)
cvd_EC_data <- as.data.frame(t(cvd_EC_data))

min_positive <- min(cvd_EC_data[cvd_EC_data > 0], na.rm = TRUE)
pseudocount <- min_positive / 2
cvd_EC_data_log10 <- log10(cvd_EC_data + pseudocount)

EC_sub <- cvd_EC_data_log10[, ###enzyme_encoding_genes, drop = FALSE]
EC_sub$SampleID <- rownames(EC_sub)

plot_df <- EC_sub %>%
  left_join(
    cvd_metadata %>%
      mutate(SampleID = rownames(cvd_metadata)) %>%
      select(SampleID, Group),
    by = "SampleID"
  )

plot_long <- plot_df %>%
  pivot_longer(
    cols = all_of(sig_features),
    names_to = "feature",
    values_to = "abundance"
  ) %>%
  filter(!is.na(abundance)) %>%
  mutate(
    feature = factor(
      feature,
      levels = sig_features,
      labels = ec_labels[sig_features]
    ),
    Group = factor(Group, levels = c("SA", "ACS", "HFrEF"))
  )

# pairwise Wilcoxon rank-sum tests
pairwise_res <- plot_long %>%
  group_by(feature) %>%
  pairwise_wilcox_test(
    abundance ~ Group,
    p.adjust.method = "none"
  )