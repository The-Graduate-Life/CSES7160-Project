# ==============================================================================
# Prediction Diagnostics — Extended Accuracy Metrics (Independent Script)
# Genomic Dissection of Oleic Acid Concentration in Peanut
#
# Author  : Fritzner Pierre
# Course  : CSES 7160 — Genetic Data Analysis
#
# Description:
#   Computes extended prediction accuracy diagnostics for the additive GBLUP
#   model fitted in Scripts 01 and 02. Goes beyond the standard Pearson r
#   reported during cross-validation by evaluating multiple complementary
#   metrics that capture different aspects of prediction quality relevant
#   to a breeding program context.
#
# Metrics computed:
#   1. 95% CI (empirical)      — range of the middle 95% of fold-level CV r
#                                values; reflects true uncertainty in PA
#                                estimate across 5-fold × 10 rep CV
#   2. Spearman ρ              — rank-based accuracy; assesses whether the
#                                model correctly ranks accessions by oleic
#                                acid content, independent of scale
#   3. MAE (oleic acid units)  — mean absolute prediction error in % oleic;
#                                directly interpretable by breeders
#   4. Classification accuracy — proportion of accessions correctly classified
#                                as high- or low-oleic relative to the grand
#                                mean; answers the practical breeding question
#   5. SD ratio (pred/obs)     — ratio of predicted to observed standard
#                                deviations; quantifies GBLUP shrinkage
#                                toward the mean (expected < 1.0)
#   6. CV bias slope           — slope of obs ~ pred regression; values > 1.0
#                                indicate the model under-disperses predictions
#                                relative to observations
#   7. High-oleic group r      — Pearson r within the top tertile of observed
#                                values; critical for identifying elite lines
#   8. Low-oleic group r       — Pearson r within the bottom tertile; relevant
#                                for eliminating poor performers
#
# Requirements (produced by earlier pipeline scripts):
#   results/gp/gblup_cv_results.csv              — Script 01: fold-level CV
#   results/gp/gblup_full_predictions.csv        — Script 01: full model obs/pred
#   results/gp/gblup_deregressed_cv.csv          — Script 02: deBLUP fold CV
#   results/gp/gblup_deregressed_cv_pairs.csv    — Script 02: deBLUP obs/pred pairs
#
# Outputs:
#   results/gp/prediction_diagnostics_raw.csv
#   results/gp/prediction_diagnostics_deregressed.csv
#   results/gp/prediction_diagnostics_combined.csv
# 
#
# Run from the project root (Oleic_acid/):
#   source("scripts/prediction_diagnostics.R")
# ==============================================================================

suppressPackageStartupMessages(library(here))
PROJECT_DIR <- here::here()
setwd(PROJECT_DIR)

cat("======================================================================\n")
cat(" Extended Prediction Diagnostics — Oleic Acid GBLUP\n")
cat("======================================================================\n\n")


# ==============================================================================
# 1 — LOAD DATA
# ==============================================================================

# CV fold results (raw BLUPs)
cv_raw <- read.csv("results/gp/gblup_cv_results.csv",
                   stringsAsFactors = FALSE)

# Full model predictions (raw BLUPs)
full_pred <- read.csv("results/gp/gblup_full_predictions.csv",
                      stringsAsFactors = FALSE)

# Deregressed BLUPs and CV (if available)
deblup_file  <- "results/gp/deregressed_blups.csv"
dercv_file   <- "results/gp/gblup_deregressed_cv.csv"

has_deblup <- file.exists(deblup_file) && file.exists(dercv_file)
if (has_deblup) {
  deblup_df <- read.csv(deblup_file,  stringsAsFactors = FALSE)
  cv_der    <- read.csv(dercv_file,   stringsAsFactors = FALSE)
  cat("Deregressed BLUP files found — will compute diagnostics for both.\n\n")
} else {
  cat("Deregressed BLUP files not found — computing diagnostics for raw BLUPs only.\n\n")
}

# Column name check for full_pred
cat("Columns in gblup_full_predictions.csv:", paste(names(full_pred), collapse=", "), "\n\n")


# ==============================================================================
# 2 — HELPER FUNCTIONS
# ==============================================================================

compute_diagnostics <- function(obs, pred, cv_r_vec, label = "") {
  
  cat(sprintf("--- %s ---\n", label))
  
  # ── 1. 95% CI (empirical from CV fold r values) ───────────────────────────
  ci <- quantile(cv_r_vec, c(0.025, 0.975), na.rm = TRUE)
  cat(sprintf("Mean CV r               : %.4f\n", mean(cv_r_vec, na.rm = TRUE)))
  cat(sprintf("95%% CI (empirical)      : [%.4f, %.4f]\n", ci[1], ci[2]))
  
  # ── 2. Spearman ρ ─────────────────────────────────────────────────────────
  spearman <- cor(obs, pred, method = "spearman", use = "complete.obs")
  cat(sprintf("Spearman rho            : %.4f\n", spearman))
  
  # ── 3. MAE ────────────────────────────────────────────────────────────────
  mae <- mean(abs(pred - obs), na.rm = TRUE)
  cat(sprintf("MAE (oleic acid units)  : %.4f %%\n", mae))
  
  # ── 4. Classification accuracy ────────────────────────────────────────────
  threshold  <- mean(obs, na.rm = TRUE)   # grand mean as cut-off
  obs_class  <- ifelse(obs  >= threshold, "high", "low")
  pred_class <- ifelse(pred >= threshold, "high", "low")
  class_acc  <- mean(obs_class == pred_class, na.rm = TRUE)
  cat(sprintf("Classification accuracy : %.4f  (threshold = %.3f%%)\n",
              class_acc, threshold))
  
  # ── 5. SD ratio pred/obs ─────────────────────────────────────────────────
  sd_ratio <- sd(pred, na.rm = TRUE) / sd(obs, na.rm = TRUE)
  cat(sprintf("SD ratio (pred/obs)     : %.4f  (1.0 = no shrinkage)\n", sd_ratio))
  
  # ── 6. CV bias slope ─────────────────────────────────────────────────────
  bias_fit   <- lm(obs ~ pred)
  bias_slope <- coef(bias_fit)[2]
  cat(sprintf("CV bias slope           : %.4f  (1.0 = unbiased)\n", bias_slope))
  
  # ── 7 & 8. Group-stratified r ────────────────────────────────────────────
  tert       <- quantile(obs, c(1/3, 2/3), na.rm = TRUE)
  high_idx   <- obs >= tert[2]
  low_idx    <- obs <= tert[1]
  
  r_high <- if (sum(high_idx, na.rm = TRUE) > 3)
    cor(obs[high_idx], pred[high_idx], use = "complete.obs") else NA_real_
  r_low  <- if (sum(low_idx,  na.rm = TRUE) > 3)
    cor(obs[low_idx],  pred[low_idx],  use = "complete.obs") else NA_real_
  
  cat(sprintf("High-oleic group r      : %.4f  (n=%d, obs >= %.2f%%)\n",
              r_high, sum(high_idx, na.rm = TRUE), tert[2]))
  cat(sprintf("Low-oleic  group r      : %.4f  (n=%d, obs <= %.2f%%)\n",
              r_low,  sum(low_idx,  na.rm = TRUE), tert[1]))
  cat("\n")
  
  # Return as data frame
  data.frame(
    Metric = c(
      "Mean CV r",
      "95% CI lower (empirical)",
      "95% CI upper (empirical)",
      "Spearman rho",
      "MAE (oleic acid %)",
      "Classification accuracy",
      "SD ratio (pred/obs)",
      "CV bias slope",
      "High-oleic group r",
      "Low-oleic group r"
    ),
    Value = round(c(
      mean(cv_r_vec, na.rm = TRUE),
      ci[1], ci[2],
      spearman, mae, class_acc,
      sd_ratio, bias_slope,
      r_high, r_low
    ), 4),
    Note = c(
      "Pearson r across all CV folds",
      "2.5th percentile of fold-level r",
      "97.5th percentile of fold-level r",
      "Rank-based accuracy (robust to outliers)",
      "Mean |predicted - observed| in % oleic",
      paste0("Threshold = grand mean (", round(threshold, 2), "%)"),
      "< 1.0 = GBLUP shrinkage toward mean",
      "Slope of obs ~ pred regression (ideal = 1.0)",
      paste0("Top tertile: obs >= ", round(tert[2], 2), "%"),
      paste0("Bottom tertile: obs <= ", round(tert[1], 2), "%")
    ),
    stringsAsFactors = FALSE
  )
}


# ==============================================================================
# 3 — RAW BLUP DIAGNOSTICS
# ==============================================================================

cat("======================================================================\n")
cat("RAW BLUP MODEL\n")
cat("======================================================================\n")

# Identify observed and predicted columns
obs_col  <- if ("observed"  %in% names(full_pred)) {
  "observed"
} else if ("BLUP" %in% names(full_pred)) {
  "BLUP"
} else {
  names(full_pred)[2]
}
pred_col <- if ("predicted" %in% names(full_pred)) {
  "predicted"
} else {
  names(full_pred)[3]
}

cat(sprintf("Using columns: obs='%s', pred='%s'\n\n", obs_col, pred_col))

obs_raw  <- full_pred[[obs_col]]
pred_raw <- full_pred[[pred_col]]
cv_r_raw <- cv_raw$r

diag_raw <- compute_diagnostics(obs_raw, pred_raw, cv_r_raw,
                                label = "Raw BLUP GBLUP")

write.csv(diag_raw, "results/gp/prediction_diagnostics_raw.csv",
          row.names = FALSE)
cat("Saved: results/gp/prediction_diagnostics_raw.csv\n\n")


# ==============================================================================
# 4 — DEREGRESSED BLUP DIAGNOSTICS (if available)
# ==============================================================================

if (has_deblup) {
  
  cat("======================================================================\n")
  cat("DEREGRESSED BLUP MODEL\n")
  cat("======================================================================\n")
  
  # For deregressed BLUPs, obs = deBLUP, pred = weighted GBLUP prediction
  # We need to reconstruct full-model predictions using the saved deBLUPs
  # and the heritability file
  h2_file <- "results/tables/heritability.csv"
  if (file.exists(h2_file)) {
    h2_tbl <- read.csv(h2_file, stringsAsFactors = FALSE)
    h2g    <- h2_tbl$Estimate[h2_tbl$Parameter == "h2g"]
    Vu     <- h2_tbl$Estimate[h2_tbl$Parameter == "sigma2_g"]
    Ve     <- h2_tbl$Estimate[h2_tbl$Parameter == "sigma2_e"]
  } else {
    cat("heritability.csv not found — skipping deregressed full-model predictions.\n")
    has_deblup <- FALSE
  }
  
  cv_pairs_file <- "results/gp/gblup_deregressed_cv_pairs.csv"
  
  if (!file.exists(cv_pairs_file)) {
    cat("gblup_deregressed_cv_pairs.csv not found.\n")
    cat("Re-run 02_GP_Deregressed_BLUPs.R to generate it, then re-run this script.\n\n")
    has_deblup <- FALSE
  }
  
  if (has_deblup) {
    cv_pairs <- read.csv(cv_pairs_file, stringsAsFactors = FALSE)
    cat(sprintf("CV pairs loaded: %d observations across all folds\n\n",
                nrow(cv_pairs)))
    
    # Use all held-out obs/pred pairs — this is the honest out-of-sample estimate
    obs_der  <- cv_pairs$observed
    pred_der <- cv_pairs$predicted
    
    diag_der <- compute_diagnostics(obs_der, pred_der, cv_der$r,
                                    label = "Deregressed BLUP GBLUP (CV pairs)")
    
    write.csv(diag_der, "results/gp/prediction_diagnostics_deregressed.csv",
              row.names = FALSE)
    cat("Saved: results/gp/prediction_diagnostics_deregressed.csv\n\n")
  }
}


# ==============================================================================
# 5 — COMBINED SUMMARY TABLE
# ==============================================================================

cat("======================================================================\n")
cat("COMBINED SUMMARY\n")
cat("======================================================================\n\n")

if (has_deblup && exists("diag_der")) {
  combined <- merge(diag_raw[, c("Metric","Value")],
                    diag_der[, c("Metric","Value")],
                    by = "Metric", suffixes = c("_Raw", "_Deregressed"),
                    all = TRUE)
  print(combined, row.names = FALSE)
  write.csv(combined, "results/gp/prediction_diagnostics_combined.csv",
            row.names = FALSE)
  cat("\nSaved: results/gp/prediction_diagnostics_combined.csv\n")
} else {
  print(diag_raw[, c("Metric","Value","Note")], row.names = FALSE)
}

cat("\n======================================================================\n")
cat("Done.\n")