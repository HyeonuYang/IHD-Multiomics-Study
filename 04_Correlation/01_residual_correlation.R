suppressPackageStartupMessages({
  library(dplyr)
  library(readr)
  library(stringr)
})

# =========================
# 0) PATHS
# =========================
taxa_all_results_fp <- '/path/to/taxa_main012_maaslin3_noTSS/all_results.tsv'
met_all_results_fp  <- '/path/to/MaAsLin3_results/Metabolomics_ordinal_MaAsLin3/all_results.tsv'

taxa_matrix_fp <- '/path/to/taxa_main012_maaslin3_noTSS/features/data_transformed.tsv'
met_matrix_fp  <- '/path/to/MaAsLin3_results/Metabolomics_ordinal_MaAsLin3/features/data_transformed.tsv'

meta_fp <- '/path/to/input/metadata_filt.csv'

out_dir <- 'results/taxa_metabolite_corr'
dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)

# Preferred column name for disease stage
stage_meta_preferred <- 'Group_ord'

# =========================
# 0.5) CUTS
# =========================
# 1) Retain only abundance models for both taxa and metabolites
model_keep <- 'abundance'

# 3) Significance thresholds: p < 0.05 for taxa, q < 0.05 for metabolites
taxa_stage_p_cut <- 0.05
met_stage_q_cut  <- 0.05

# 2) Output thresholds for correlation analysis
corr_cuts <- tibble::tibble(
  stat = c('p', 'q', 'q'),
  cut  = c(0.05, 0.05, 0.1)
)

fmt_cut <- function(x) {
  s <- formatC(x, format = 'fg', digits = 3)
  sub('\\.?0+$', '', s)
}

# =========================
# 1) Load Metadata
# =========================
meta <- readr::read_csv(meta_fp, show_col_types = FALSE) %>%
  mutate(SampleID = as.character(SampleID))

if (!'Group_ord' %in% names(meta)) {
  if (!'Group' %in% names(meta)) stop('Metadata lacks both Group_ord and Group columns. Cannot proceed with stage definition.')
  meta <- meta %>% mutate(Group_ord = Group)
}

meta <- meta %>%
  mutate(
    Group_ord = factor(
      str_trim(as.character(Group_ord)),
      levels = c('Angina', 'MI', 'HF'),
      ordered = TRUE
    )
  ) %>%
  filter(!is.na(Group_ord))

stopifnot(all(c('SampleID', 'Group_ord') %in% names(meta)))

# =========================
# 2) Load Data Matrices
# =========================
read_data_transformed <- function(fp) {
  df <- readr::read_tsv(fp, show_col_types = FALSE)

  if ('SampleID' %in% names(df)) {
    return(df %>% mutate(SampleID = as.character(SampleID)))
  }

  if ('feature' %in% names(df)) {
    return(df %>%
      rename(SampleID = feature) %>%
      mutate(SampleID = as.character(SampleID))
    )
  }

  if (names(df)[1] %in% c('...1', 'X1', 'V1')) {
    return(df %>%
      rename(SampleID = 1) %>%
      mutate(SampleID = as.character(SampleID))
    )
  }

  stop(sprintf('Cannot locate SampleID column in: %s', fp))
}

taxa_mat <- read_data_transformed(taxa_matrix_fp)
met_mat  <- read_data_transformed(met_matrix_fp)

stopifnot('SampleID' %in% names(taxa_mat), 'SampleID' %in% names(met_mat))

common_ids <- Reduce(intersect, list(meta$SampleID, taxa_mat$SampleID, met_mat$SampleID))
if (length(common_ids) == 0) {
  stop(
    'No overlapping SampleIDs found across meta, taxa, and metabolite matrices.\n',
    'Example meta IDs: ', paste(head(meta$SampleID, 5), collapse = ', '), '\n',
    'Example taxa_mat IDs: ', paste(head(taxa_mat$SampleID, 5), collapse = ', '), '\n',
    'Example met_mat IDs: ', paste(head(met_mat$SampleID, 5), collapse = ', ')
  )
}

meta <- meta %>% filter(SampleID %in% common_ids) %>% arrange(SampleID)
taxa_mat <- taxa_mat %>% filter(SampleID %in% common_ids) %>% arrange(SampleID)
met_mat  <- met_mat  %>% filter(SampleID %in% common_ids) %>% arrange(SampleID)

cat(sprintf('✅ Common samples identified: %d\n', length(common_ids)))

# =========================
# 3) Load MaAsLin3 Results
# =========================
taxa_res <- readr::read_tsv(taxa_all_results_fp, show_col_types = FALSE)
met_res  <- readr::read_tsv(met_all_results_fp,  show_col_types = FALSE)

pick_first_col <- function(df, candidates) {
  hit <- candidates[candidates %in% names(df)]
  if (length(hit) == 0) return(NA_character_)
  hit[1]
}

taxa_p_col <- pick_first_col(taxa_res, c('pval', 'pval_individual', 'pval_individuals'))
met_q_col  <- pick_first_col(met_res,  c('qval_individual', 'qval', 'qval_individuals'))

if (is.na(taxa_p_col)) stop('Failed to locate p-value column in taxa results.')
if (is.na(met_q_col))  stop('Failed to locate q-value column in metabolite results.')

# =========================
# 4) Filter by Stage Term
# =========================
filter_stage_GroupL <- function(df, stage_meta_name, model_keep, stat_col, stat_cut) {
  req <- c('feature', 'metadata', 'name', 'coef', 'model', stat_col)
  miss <- setdiff(req, names(df))
  if (length(miss) > 0) stop('Missing required columns: ', paste(miss, collapse = ', '))

  target_name <- paste0(stage_meta_name, '.L')

  out <- df %>%
    mutate(model2 = tolower(as.character(model))) %>%
    filter(
      model2 == model_keep,
      metadata == stage_meta_name,
      name == target_name,
      !is.na(.data[[stat_col]]),
      .data[[stat_col]] < stat_cut,
      !is.na(coef)
    ) %>%
    transmute(
      feat = feature,
      coef = coef,
      stat_value = .data[[stat_col]],
      term = name,
      model = model2,
      stage_meta = stage_meta_name
    )

  cat(sprintf('🔎 Filtered to %d rows (name==%s, %s<%.3g, model==%s)\n',
              nrow(out), target_name, stat_col, stat_cut, model_keep))
  out
}

taxa_sig <- filter_stage_GroupL(taxa_res, stage_meta_preferred, model_keep, taxa_p_col, taxa_stage_p_cut) %>%
  rename(pv = stat_value)

met_sig  <- filter_stage_GroupL(met_res,  stage_meta_preferred, model_keep, met_q_col,  met_stage_q_cut) %>%
  rename(qv = stat_value)

readr::write_tsv(taxa_sig, file.path(out_dir, 'taxa_sig_stage_p0.05.tsv'))
readr::write_tsv(met_sig,  file.path(out_dir, 'met_sig_stage_q0.05.tsv'))

cat(sprintf('✅ Significant features retained - Taxa: %d, Metabolites: %d\n', nrow(taxa_sig), nrow(met_sig)))

if (nrow(taxa_sig) == 0) {
  stop(
    'No significant taxa features passed the filter criteria.\n',
    'name==', paste0(stage_meta_preferred, '.L'),
    ', model==', model_keep,
    ', ', taxa_p_col, '<', taxa_stage_p_cut, '\n'
  )
}
if (nrow(met_sig) == 0) {
  stop(
    'No significant metabolite features passed the filter criteria.\n',
    'name==', paste0(stage_meta_preferred, '.L'),
    ', model==', model_keep,
    ', ', met_q_col, '<', met_stage_q_cut, '\n'
  )
}

# =========================
# 5) Extract Features from Matrices
# =========================
taxa_cols_all <- setdiff(names(taxa_mat), 'SampleID')
met_cols_all  <- setdiff(names(met_mat),  'SampleID')

taxa_feats <- intersect(taxa_sig$feat, taxa_cols_all)
met_feats  <- intersect(met_sig$feat,  met_cols_all)

if (length(taxa_feats) == 0) stop('No overlapping significant taxa features found in matrix columns.')
if (length(met_feats) == 0)  stop('No overlapping significant metabolite features found in matrix columns.')

as_numeric_matrix <- function(df, cols) {
  X <- as.matrix(df[, cols, drop = FALSE])
  X <- apply(X, 2, function(v) as.numeric(v))
  X <- as.matrix(X)
  colnames(X) <- cols
  rownames(X) <- df$SampleID
  X
}

X_taxa <- as_numeric_matrix(taxa_mat, taxa_feats)
X_met  <- as_numeric_matrix(met_mat,  met_feats)

cat(sprintf('✅ Matrix dimensions ready - Taxa: %d features, Metabolites: %d features\n', ncol(X_taxa), ncol(X_met)))

# =========================
# 6) Residualize by Stage
# =========================
stage_num <- as.integer(meta$Group_ord)
names(stage_num) <- meta$SampleID

residualize_on_stage <- function(M, stage_num_named) {
  if (!all(rownames(M) %in% names(stage_num_named))) {
    stop('Stage vector names do not adequately cover matrix rownames.')
  }
  st <- stage_num_named[rownames(M)]

  out <- M
  for (j in seq_len(ncol(M))) {
    y <- M[, j]
    ok <- is.finite(y) & !is.na(y) & is.finite(st)
    if (sum(ok) < 5) {
      out[, j] <- NA_real_
      next
    }
    fit <- lm(y[ok] ~ st[ok])
    r <- rep(NA_real_, length(y))
    r[ok] <- residuals(fit)
    out[, j] <- r
  }
  out
}

X_taxa_resid <- residualize_on_stage(X_taxa, stage_num)
X_met_resid  <- residualize_on_stage(X_met,  stage_num)

# =========================
# 7) Compute Spearman Correlations
# =========================
pairwise_spearman <- function(A, B) {
  taxa <- colnames(A)
  mets <- colnames(B)

  r_mat <- matrix(NA_real_, nrow = length(taxa), ncol = length(mets),
                  dimnames = list(taxa, mets))
  p_mat <- matrix(NA_real_, nrow = length(taxa), ncol = length(mets),
                  dimnames = list(taxa, mets))

  for (i in seq_along(taxa)) {
    x <- A[, i]
    for (j in seq_along(mets)) {
      y <- B[, j]
      ok <- is.finite(x) & is.finite(y) & !is.na(x) & !is.na(y)
      if (sum(ok) < 5) next
      ct <- suppressWarnings(cor.test(x[ok], y[ok], method = 'spearman', exact = FALSE))
      r_mat[i, j] <- unname(ct$estimate)
      p_mat[i, j] <- ct$p.value
    }
  }

  list(r = r_mat, p = p_mat)
}

cat('⏳ Computing raw correlations...\n')
raw <- pairwise_spearman(X_taxa, X_met)

cat('⏳ Computing stage-residualized correlations...\n')
resid <- pairwise_spearman(X_taxa_resid, X_met_resid)

to_long <- function(r_mat, p_mat, tag) {
  df <- expand.grid(
    taxa = rownames(r_mat),
    metabolite = colnames(r_mat),
    stringsAsFactors = FALSE
  )
  df$r <- as.vector(r_mat)
  df$p <- as.vector(p_mat)
  df$analysis <- tag
  df
}

raw_long <- to_long(raw$r, raw$p, 'raw')
res_long <- to_long(resid$r, resid$p, 'residual_stage_removed')

all_long <- bind_rows(raw_long, res_long) %>%
  filter(is.finite(r), is.finite(p), !is.na(r), !is.na(p)) %>%
  group_by(analysis) %>%
  mutate(q = p.adjust(p, method = 'BH')) %>%
  ungroup()

taxa_coef <- taxa_sig %>% select(taxa = feat, taxa_coef = coef, taxa_pv = pv)
met_coef  <- met_sig  %>% select(metabolite = feat, met_coef = coef, met_qv = qv)

all_long2 <- all_long %>%
  left_join(taxa_coef, by = 'taxa') %>%
  left_join(met_coef,  by = 'metabolite') %>%
  mutate(concordant_trend = ifelse(is.na(taxa_coef) | is.na(met_coef), NA, sign(taxa_coef) == sign(met_coef)))

# =========================
# 8) Save Results
# =========================
save_full <- function(df, analysis_tag, out_dir) {
  out_fp <- file.path(out_dir, paste0('taxa_metabolite_corr_', analysis_tag, '_ALL.tsv'))
  out <- df %>%
    filter(analysis == analysis_tag) %>%
    arrange(taxa, metabolite)
  readr::write_tsv(out, out_fp)
  cat(sprintf('- %s (%s, unfiltered): %d rows\n', out_fp, analysis_tag, nrow(out)))
  invisible(out_fp)
}

save_filtered <- function(df, analysis_tag, stat, cut, out_dir) {
  stat_tag <- paste0(stat, fmt_cut(cut))
  out_fp <- file.path(out_dir, paste0('taxa_metabolite_corr_', analysis_tag, '_', stat_tag, '.tsv'))

  if (stat == 'p') {
    out <- df %>%
      filter(analysis == analysis_tag, p < cut) %>%
      arrange(p, desc(abs(r)))
  } else if (stat == 'q') {
    out <- df %>%
      filter(analysis == analysis_tag, q < cut) %>%
      arrange(q, desc(abs(r)))
  } else {
    stop('Invalid statistic specified. Must be "p" or "q".')
  }

  readr::write_tsv(out, out_fp)
  cat(sprintf('- %s (%s, %s < %s): %d rows\n', out_fp, analysis_tag, stat, fmt_cut(cut), nrow(out)))
  invisible(out_fp)
}

cat('Output files saved (Full and Filtered sets):\n')
analyses <- c('raw', 'residual_stage_removed')

# Export full result sets
for (a in analyses) {
  save_full(all_long2, a, out_dir)
}

# Export filtered result sets
for (a in analyses) {
  for (k in seq_len(nrow(corr_cuts))) {
    save_filtered(all_long2, a, corr_cuts$stat[k], corr_cuts$cut[k], out_dir)
  }
}