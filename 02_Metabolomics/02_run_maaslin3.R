library(maaslin3)
library(dplyr)

dir.create("MaAsLin3_results", showWarnings = FALSE)

raw_df <- read.csv(
    "/path/to/input/data.csv",
    check.names = FALSE
)

rownames(raw_df) <- raw_df$SampleID
raw_df$SampleID  <- NULL

metabolites_df <- raw_df %>% select(-Group)
metadata_df    <- raw_df %>% select(Group)

orig_names  <- colnames(metabolites_df)
clean_names <- make.names(orig_names, unique = TRUE)
colnames(metabolites_df) <- clean_names

feat_map <- data.frame(
    feature       = clean_names,
    feature_names = orig_names,
    stringsAsFactors = FALSE
)

write.table(
    feat_map,
    file      = "MaAsLin3_results/Metabolomics_feature_name_map.tsv",
    sep       = "\t",
    row.names = FALSE,
    quote     = FALSE
)

add_feature_names <- function(result_dir, map_file) {
    res_file <- file.path(result_dir, "all_results.tsv")
    if (!file.exists(res_file)) return(invisible(NULL))

    res <- read.delim(res_file, check.names = FALSE)
    feat_map <- read.delim(map_file, check.names = FALSE)

    res2 <- dplyr::left_join(res, feat_map, by = "feature")
    cols_order <- c("feature", "feature_names",
                    setdiff(colnames(res2), c("feature", "feature_names")))
    res2 <- res2[, cols_order]

    out_file <- file.path(result_dir, "all_results_with_names.tsv")
    write.table(
        res2,
        file      = out_file,
        sep       = "\t",
        quote     = FALSE,
        row.names = FALSE
    )

    invisible(NULL)
}

map_file <- "MaAsLin3_results/Metabolomics_feature_name_map.tsv"

meta2 <- read.delim(
    "/path/to/input/metadata2.tsv",
    row.names   = 1,
    check.names = TRUE
)

metadata <- cbind(meta2, metadata_df[rownames(meta2), , drop = FALSE])

metadata$severity_num <- ifelse(metadata$Group == "Angina", 0,
                         ifelse(metadata$Group == "MI", 1, 2))
metadata$Group <- factor(metadata$Group, levels = c("Angina", "MI", "HF"))


# Specify ordered = TRUE in the Group factor
metadata$Group_ord <- factor(
  metadata$Group,
  levels = c("Angina", "MI", "HF"),
  ordered = TRUE
)

common_samples <- intersect(rownames(metadata), rownames(metabolites_df))
metadata       <- metadata[common_samples, ]
metabolites_df <- metabolites_df[common_samples, ]

dir.create("MaAsLin3_results/Metabolomics_ordinal_MaAsLin3", showWarnings = FALSE)

# Define the formula: "~ Group_ord + age + sex + BMI + antacid_ppi"
maaslin3(
    input_data     = metabolites_df,
    input_metadata = metadata,
    output         = "MaAsLin3_results/Metabolomics_ordinal_MaAsLin3",

    formula  = "~ Group_ord + age + sex + BMI + antacid_ppi",

    normalization  = "NONE",
    transform      = "LOG",
    standardize    = FALSE,

    min_prevalence = 0,
    min_abundance  = 0,

    correction     = "BH",
    cores          = 7
)

add_feature_names(
    result_dir = "MaAsLin3_results/Metabolomics_ordinal_MaAsLin3",
    map_file   = map_file
)

pair_list <- list(
    Angina_vs_MI = c("Angina", "MI"),
    Angina_vs_HF = c("Angina", "HF"),
    MI_vs_HF     = c("MI", "HF")
)

for (nm in names(pair_list)) {

    gr <- pair_list[[nm]]
    samples_sub <- rownames(metadata)[metadata$Group %in% gr]

    meta_sub  <- metadata[samples_sub, , drop = FALSE]
    metab_sub <- metabolites_df[samples_sub, , drop = FALSE]

    meta_sub$Group <- trimws(as.character(meta_sub$Group))

    meta_sub$Group_pair <- factor(meta_sub$Group, levels = gr)

    out_dir <- paste0("MaAsLin3_results/Metabolomics_pairwise_", nm, "_MaAsLin3")
    dir.create(out_dir, showWarnings = FALSE)
    maaslin3(
        input_data     = metab_sub,
        input_metadata = meta_sub,
        output         = out_dir,

        fixed_effects  = c("Group_pair", "age", "sex", "BMI", "antacid_ppi"),

        normalization  = "NONE",
        transform      = "LOG",
        standardize    = FALSE,

        min_prevalence = 0,
        min_abundance  = 0,

        correction     = "BH",
        cores          = 7
    )

    add_feature_names(result_dir = out_dir, map_file = map_file)
}

meta_A_vs_MH <- metadata
meta_A_vs_MH$Group_bin <- ifelse(meta_A_vs_MH$Group == "Angina",
                                 "Angina", "MI_HF")
meta_A_vs_MH$Group_bin <- factor(meta_A_vs_MH$Group_bin,
                                 levels = c("Angina", "MI_HF"))

dir.create("MaAsLin3_results/Metabolomics_Angina_vs_MIHF_MaAsLin3", showWarnings = FALSE)

maaslin3(
    input_data     = metabolites_df,
    input_metadata = meta_A_vs_MH,
    output         = "MaAsLin3_results/Metabolomics_Angina_vs_MIHF_MaAsLin3",

    fixed_effects  = c("Group_bin", "age", "sex", "BMI", "antacid_ppi"),

    normalization  = "NONE",
    transform      = "LOG",
    standardize    = FALSE,

    min_prevalence = 0,
    min_abundance  = 0,

    correction     = "BH",
    cores          = 7
)

add_feature_names(
    result_dir = "MaAsLin3_results/Metabolomics_Angina_vs_MIHF_MaAsLin3",
    map_file   = map_file
)

meta_AM_vs_HF <- metadata
meta_AM_vs_HF$Group_bin <- ifelse(meta_AM_vs_HF$Group %in% c("Angina", "MI"),
                                  "Angina_MI", "HF")
meta_AM_vs_HF$Group_bin <- factor(meta_AM_vs_HF$Group_bin,
                                  levels = c("Angina_MI", "HF"))

dir.create("MaAsLin3_results/Metabolomics_AnginaMI_vs_HF_MaAsLin3", showWarnings = FALSE)

maaslin3(
    input_data     = metabolites_df,
    input_metadata = meta_AM_vs_HF,
    output         = "MaAsLin3_results/Metabolomics_AnginaMI_vs_HF_MaAsLin3",

    fixed_effects  = c("Group_bin", "age", "sex", "BMI", "antacid_ppi"),

    normalization  = "NONE",
    transform      = "LOG",
    standardize    = FALSE,

    min_prevalence = 0,
    min_abundance  = 0,

    correction     = "BH",
    cores          = 7
)

add_feature_names(
    result_dir = "MaAsLin3_results/Metabolomics_AnginaMI_vs_HF_MaAsLin3",
    map_file   = map_file
)