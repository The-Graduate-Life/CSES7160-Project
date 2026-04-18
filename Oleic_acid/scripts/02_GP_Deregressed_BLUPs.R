# ==============================================================================
# Genomic Prediction with Deregressed BLUPs (Independent Script)
# Genomic Dissection of Oleic Acid Concentration in Peanut
#
# Author  : Fritzner Pierre
# Course  : CSES 7160 — Genetic Data Analysis
#
# Revision applied:
#   The main script (CSES7160_Fritzner_Analysis.R Step 10) used raw lmer BLUPs
#   as the response in GBLUP cross-validation. This causes DOUBLE-SHRINKAGE:
#   BLUPs are already shrunk by lmer; passing them into a second mixed model
#   (rrBLUP::mixed.solve) shrinks them a second time, biasing variance
#   components downward and deflating prediction accuracy estimates.
#
#   Fix (Garrick et al. 2009, J. Dairy Sci.):
#     deBLUP_i = BLUP_i / r²_i
#     where r²_i = 1 - PEV_i / sigma²_g  (reliability of individual i's BLUP)
#     Garrick weights w_i = (1 - h²) / [(c + (1 - r²_i)/r²_i) * h²]
#     c = 0.1 (default Garrick constant; proportion of Vg not tagged by SNPs)
#
#   The deregressed BLUPs (and their weights) are then used as the response
#   in GBLUP cross-validation, giving unbiased prediction accuracy estimates.
#
# Requirements (produced by CSES7160_Fritzner_Analysis.R):
#   results/geno_num.rds      — post-QC genotype matrix {0,1,2}, taxa × SNPs
#   results/geno_map.rds      — SNP position map (SNP, Chr, Pos)
#   results/tables/BLUPs_OleicAcid.csv  — lmer BLUPs and GEBVs per accession
#
# Outputs:
#   results/gp/deregressed_blups.csv          — deBLUP, weight, reliability per accession
#   results/gp/gblup_deregressed_cv.csv       — fold-level CV results
#   results/plots/12_gblup_deregressed_cv.pdf — CV distribution + scatter
# ==============================================================================


# ==============================================================================
# 0 — WORKING DIRECTORY
# ==============================================================================
# ---- Package check ----
if (!requireNamespace("here", quietly = TRUE)) {
  stop("Package 'here' is required. Install with install.packages('here')")
}

suppressPackageStartupMessages(library(here))

# ---- Project root ----
PROJECT_DIR <- here::here()
setwd(PROJECT_DIR)  # ensure paths resolve from project root

cat("Project root:\n ", PROJECT_DIR, "\n\n")

# ---- Verify directory exists ----
if (!dir.exists(PROJECT_DIR)) {
  stop("Project directory not found:\n ", PROJECT_DIR)
}

# ==============================================================================
# 1 — PACKAGES
# ==============================================================================
# install.packages(c("rrBLUP", "ggplot2", "dplyr", "patchwork"))

pkgs <- c("rrBLUP","ggplot2","dplyr","patchwork")

missing <- pkgs[!sapply(pkgs, requireNamespace, quietly = TRUE)]

if (length(missing) > 0) {
  stop("Missing packages: ", paste(missing, collapse = ", "))
}

suppressPackageStartupMessages(
  lapply(pkgs, library, character.only = TRUE)
)

set.seed(2026)

for (d in c("results/gp", "results/plots"))
  dir.create(d, showWarnings = FALSE, recursive = TRUE)

cat("======================================================================\n")
cat(" Deregressed BLUP Genomic Prediction \n")
cat("======================================================================\n")
cat("R version :", R.version$version.string, "\n")
cat("Date/Time :", format(Sys.time(), "%Y-%m-%d %H:%M:%S"), "\n\n")


# ==============================================================================
# 2 — LOAD DATA
# ==============================================================================
cat("----------------------------------------------------------------------\n")
cat("STEP 2 | Load Data\n")
cat("----------------------------------------------------------------------\n")

# Genotype matrix {0,1,2}
if (!file.exists("results/geno_num.rds"))
  stop("oops! results/geno_num.rds not found.\nRun CSES7160_Fritzner_Analysis.R first.")
geno_num <- readRDS("results/geno_num.rds")
cat(" Loaded: results/geno_num.rds —", nrow(geno_num), "taxa x",
    ncol(geno_num), "SNPs\n")

# lmer BLUPs (per accession)
blup_file <- "results/tables/BLUPs_OleicAcid.csv"
if (!file.exists(blup_file))
  stop(blup_file, " Oops! not found.\nRun CSES7160_Fritzner_Analysis.R first.")
blup_df <- read.csv(blup_file, stringsAsFactors = FALSE)
cat(" Loaded:", blup_file, "—", nrow(blup_df), "accessions\n\n")

# Align taxa present in both objects
common_taxa  <- intersect(blup_df$Taxa, rownames(geno_num))
blup_aligned <- blup_df[match(common_taxa, blup_df$Taxa), ]
geno_num     <- geno_num[common_taxa, , drop = FALSE]
y            <- blup_aligned$BLUP
cat("Common taxa (pheno ∩ geno):", length(common_taxa), "\n\n")

# Center genotype matrix for rrBLUP: {0,1,2} -> {-1,0,1}
geno_rrblup <- geno_num - 1


# ==============================================================================
# 3 — GENOMIC RELATIONSHIP MATRIX (VanRaden 2008)
# ==============================================================================
cat("----------------------------------------------------------------------\n")
cat("STEP 3 | Build GRM\n")
cat("----------------------------------------------------------------------\n")

K_mat <- A.mat(geno_rrblup, min.MAF = 0)
cat("GRM:", nrow(K_mat), "x", ncol(K_mat), "\n\n")


# ==============================================================================
# 4 — FIT FULL ADDITIVE GBLUP (to extract PEV)
# ==============================================================================
# mixed.solve() returns:
#   $Vu  = sigma²_g (genomic variance)
#   $Ve  = sigma²_e (residual variance)
#   $u   = GEBVs (n×1 vector)
#   $LL  = log-likelihood
#
# Prediction error variance (PEV):
#   PEV_i = Var(u_i - u_hat_i) = sigma²_g * (1 - H_ii)
#   where H = K * (K + lambda*I)^{-1},  lambda = sigma²_e / sigma²_g
#   rrBLUP does NOT directly expose PEV, so we recompute it from H.

cat("----------------------------------------------------------------------\n")
cat("STEP 4 | Fit Full GBLUP and Compute PEV\n")
cat("----------------------------------------------------------------------\n")

lmm_full <- mixed.solve(y = y, K = K_mat, SE = FALSE, return.Hinv = TRUE)

Vu   <- lmm_full$Vu   # genomic variance
Ve   <- lmm_full$Ve   # residual variance
h2g  <- Vu / (Vu + Ve)
lam  <- Ve / Vu       # lambda = sigma²_e / sigma²_g

cat(sprintf("sigma²_g  = %.6f\n", Vu))
cat(sprintf("sigma²_e  = %.6f\n", Ve))
cat(sprintf("h²_g      = %.4f\n", h2g))
cat(sprintf("lambda    = %.4f\n\n", lam))

# Compute H = K (K + lambda I)^{-1}   using the returned H^{-1}
# mixed.solve returns Hinv = (K + lambda*I)^{-1} when return.Hinv = TRUE
Hinv <- lmm_full$Hinv                          # (K + lam*I)^{-1}
H    <- K_mat %*% Hinv                         # K * (K + lam*I)^{-1}

# PEV_i = sigma²_g * (1 - H_ii)
pev  <- Vu * (1 - diag(H))
pev  <- pmax(pev, 0)       # numerical floor

# Reliability r²_i = 1 - PEV_i / sigma²_g  (capped to [0.05, 0.99])
rel  <- 1 - pev / Vu
rel  <- pmin(pmax(rel, 0.05), 0.99)

cat("Reliability (r²) summary:\n")
print(round(summary(rel), 4))
cat("\n")


# ==============================================================================
# 5 — DEREGRESS BLUPs (Garrick et al. 2009)
# ==============================================================================
# deBLUP_i = BLUP_i / r²_i
# Garrick weight:
#   w_i = (1 - h²) / [(c + (1 - r²_i) / r²_i) * h²]
# c = 0.1 (default: ~10% of genetic variance not captured by SNPs)
# Larger w_i → higher confidence observation; passed as diagonal of
# residual covariance R = diag(1/w_i) * sigma²_e in the prediction model.

cat("----------------------------------------------------------------------\n")
cat("STEP 5 | Deregress BLUPs (Garrick et al. 2009)\n")
cat("----------------------------------------------------------------------\n")

c_const <- 0.1

gebv    <- as.numeric(lmm_full$u)             # GEBVs from full GBLUP
deblup  <- gebv / rel                         # deregressed BLUPs
w_i     <- (1 - h2g) / ((c_const + (1 - rel) / rel) * h2g)

deblup_df <- data.frame(
  Taxa   = common_taxa,
  BLUP   = y,
  GEBV   = gebv,
  PEV    = pev,
  rel    = rel,
  deBLUP = deblup,
  weight = w_i,
  stringsAsFactors = FALSE
)

cat("Deregressed BLUP summary:\n")
print(round(summary(deblup_df$deBLUP), 3))
cat("\nWeight summary:\n")
print(round(summary(deblup_df$weight), 4))

write.csv(deblup_df, "results/gp/deregressed_blups.csv", row.names = FALSE)
cat("\nTable saved: results/gp/deregressed_blups.csv\n\n")


# ==============================================================================
# 6 — 5-FOLD CV × 10 REPS  (Deregressed BLUPs as response)
# ==============================================================================
# Prediction model per fold:
#   deBLUP ~ 1 + g + e*
#   g  ~ N(0, K * sigma²_g)
#   e* ~ N(0, diag(1/w_i) * sigma²_e)   <- heterogeneous residuals
#
# rrBLUP::mixed.solve() does not natively accept a diagonal R matrix.
# We approximate it by scaling the phenotype and K:
#   Approach: fit  y_scaled = diag(sqrt(w)) * deBLUP
#             with K_scaled = diag(sqrt(w)) * K * diag(sqrt(w))
#   This is equivalent to WLS in the LMM context.
#
# Prediction accuracy (PA) = Pearson r(observed deBLUP, predicted deBLUP)
# in the test fold.

cat("----------------------------------------------------------------------\n")
cat("STEP 6 | 5-Fold CV × 10 Reps (Deregressed BLUPs)\n")
cat("----------------------------------------------------------------------\n")

N_FOLDS <- 5
N_REPS  <- 10
n_total <- length(common_taxa)
W_sqrt  <- sqrt(w_i)           # per-individual square root weights

cv_store <- data.frame(
  rep  = integer(),
  fold = integer(),
  r    = numeric(),
  RMSE = numeric(),
  bias = numeric()
)

for (rep in seq_len(N_REPS)) {
  fid <- sample(rep(seq_len(N_FOLDS), length.out = n_total))

  for (fold in seq_len(N_FOLDS)) {
    test_idx  <- which(fid == fold)
    train_idx <- which(fid != fold)

    # Scale response and K by sqrt(w) for WLS
    y_cv         <- deblup_df$deBLUP * W_sqrt
    y_cv[test_idx] <- NA                       # mask test individuals

    K_sc <- diag(W_sqrt) %*% K_mat %*% diag(W_sqrt)   # weighted K

    fv <- tryCatch(
      mixed.solve(y = y_cv, K = K_sc),
      error = function(e) NULL
    )

    if (!is.null(fv)) {
      # Predicted deBLUP: unscale by dividing by w_sqrt
      yp_scaled  <- as.numeric(fv$u) + as.numeric(fv$beta)
      yp         <- yp_scaled[test_idx] / W_sqrt[test_idx]   # unscale
      yo         <- deblup_df$deBLUP[test_idx]

      cv_store <- rbind(cv_store, data.frame(
        rep  = rep,
        fold = fold,
        r    = cor(yo, yp, use = "complete.obs"),
        RMSE = sqrt(mean((yo - yp)^2, na.rm = TRUE)),
        bias = coef(lm(yo ~ yp))[2]
      ))
    }
  }

  if (rep %% 2 == 0)
    cat(sprintf("  Rep %2d / %d  |  mean r = %.4f\n",
                rep, N_REPS, mean(cv_store$r, na.rm = TRUE)))
}

cat("\n--- Cross-Validation Summary (Deregressed BLUPs) ---\n")
cat("Mean r  :", round(mean(cv_store$r,    na.rm = TRUE), 4), "\n")
cat("SD r    :", round(sd(cv_store$r,      na.rm = TRUE), 4), "\n")
cat("95% CI  : [",
    round(quantile(cv_store$r, 0.025, na.rm = TRUE), 4), ",",
    round(quantile(cv_store$r, 0.975, na.rm = TRUE), 4), "]\n")
cat("Mean RMSE:", round(mean(cv_store$RMSE, na.rm = TRUE), 4), "\n")
cat("Bias (slope):", round(mean(cv_store$bias, na.rm = TRUE), 4),
    "(1.0 = unbiased)\n\n")

write.csv(cv_store, "results/gp/gblup_deregressed_cv.csv", row.names = FALSE)
cat("Table saved: results/gp/gblup_deregressed_cv.csv\n\n")


# ==============================================================================
# 7 — FULL MODEL PREDICTIONS (for scatter plot)
# ==============================================================================
cat("----------------------------------------------------------------------\n")
cat("STEP 7 | Full Model Predictions\n")
cat("----------------------------------------------------------------------\n")

y_scaled_full <- deblup_df$deBLUP * W_sqrt
K_sc_full     <- diag(W_sqrt) %*% K_mat %*% diag(W_sqrt)

fit_full <- mixed.solve(y = y_scaled_full, K = K_sc_full)

pred_scaled <- as.numeric(fit_full$u) + as.numeric(fit_full$beta)
pred_full   <- pred_scaled / W_sqrt      # unscale back to deBLUP scale

scatter_df <- data.frame(
  Taxa      = common_taxa,
  deBLUP    = deblup_df$deBLUP,
  predicted = pred_full
)
fr <- cor(scatter_df$deBLUP, scatter_df$predicted)
cat(sprintf("Full model r (deBLUP ~ predicted): %.4f\n\n", fr))


# ==============================================================================
# 8 — PLOTS
# ==============================================================================
cat("----------------------------------------------------------------------\n")
cat("STEP 8 | Plots\n")
cat("----------------------------------------------------------------------\n")

# -- CV distribution ----------------------------------------------------------
p_cv <- ggplot(cv_store, aes(x = r)) +
  geom_histogram(aes(y = after_stat(density)), binwidth = 0.03,
                 fill = "#2E6B4F", colour = "#1B3A2D", alpha = 0.82) +
  geom_density(colour = "#1B3A2D", linewidth = 1.0) +
  geom_vline(xintercept = mean(cv_store$r, na.rm = TRUE),
             colour = "#B85042", linetype = "dashed", linewidth = 1.0) +
  annotate("text",
           x = mean(cv_store$r, na.rm = TRUE) + 0.005, y = Inf,
           vjust = 1.5, hjust = 0,
           label = paste0("Mean r = ", round(mean(cv_store$r, na.rm = TRUE), 3)),
           colour = "#B85042", size = 3.8) +
  labs(title    = "GBLUP CV Accuracy — Deregressed BLUPs",
       subtitle = paste0(N_FOLDS, "-fold CV × ", N_REPS,
                         " reps  |  n = ", n_total, " accessions\n",
                         "Response: deregressed BLUPs (Garrick et al. 2009)"),
       x = "Pearson r (observed vs predicted)",
       y = "Density") +
  theme_classic(base_size = 12) +
  theme(plot.title    = element_text(face = "bold", colour = "#1B3A2D"),
        plot.subtitle = element_text(colour = "#5C7A65", size = 9))

# -- Predicted vs observed scatter --------------------------------------------
p_scat <- ggplot(scatter_df, aes(x = predicted, y = deBLUP)) +
  geom_point(colour = "#2E6B4F", alpha = 0.72, size = 2.5, shape = 16) +
  geom_smooth(method = "lm", colour = "#B85042", se = TRUE,
              linewidth = 1.0, fill = "#EAF4EE") +
  geom_abline(slope = 1, intercept = 0, colour = "#AAAAAA",
              linewidth = 0.8, linetype = "dashed") +
  annotate("text", x = -Inf, y = Inf, hjust = -0.1, vjust = 1.5,
           label = paste0("r = ", round(fr, 3)),
           colour = "#1B3A2D", fontface = "bold", size = 4.5) +
  labs(title    = "Predicted vs Observed (Full Model)",
       subtitle = paste0("GBLUP | Deregressed BLUPs | r = ", round(fr, 3)),
       x = "Genomic Predicted Value",
       y = "Deregressed BLUP of Oleic Acid") +
  theme_classic(base_size = 12) +
  theme(plot.title    = element_text(face = "bold", colour = "#1B3A2D"),
        plot.subtitle = element_text(colour = "#5C7A65"))

p_gp <- (p_cv | p_scat) +
  plot_annotation(
    caption = "GBLUP | Deregressed BLUPs | Garrick et al. (2009)",
    theme   = theme(plot.caption = element_text(colour = "#5C7A65", size = 8)))

pdf("results/plots/12_gblup_deregressed_cv.pdf", width = 12, height = 5)
print(p_gp)
dev.off()
cat("Plot saved: results/plots/12_gblup_deregressed_cv.pdf\n\n")


# ==============================================================================
# 9 — COMPARISON: ORIGINAL VS DEREGRESSED
# ==============================================================================
# Load the original CV results from the main script for direct comparison.

cat("----------------------------------------------------------------------\n")
cat("STEP 9 | Compare Original vs Deregressed CV Accuracy\n")
cat("----------------------------------------------------------------------\n")

orig_cv_file <- "results/gp/gblup_cv_results.csv"

if (file.exists(orig_cv_file)) {
  orig_cv <- read.csv(orig_cv_file, stringsAsFactors = FALSE)

  cat(sprintf("Original   CV mean r : %.4f  (SD: %.4f)\n",
              mean(orig_cv$r,    na.rm = TRUE),
              sd(orig_cv$r,      na.rm = TRUE)))
  cat(sprintf("Deregressed CV mean r: %.4f  (SD: %.4f)\n",
              mean(cv_store$r,   na.rm = TRUE),
              sd(cv_store$r,     na.rm = TRUE)))

  diff_r <- mean(cv_store$r, na.rm = TRUE) -
            mean(orig_cv$r,  na.rm = TRUE)
  cat(sprintf("Difference           : %+.4f\n\n", diff_r))

  if (diff_r > 0) {
    cat(">> Deregressed BLUPs yield higher PA.\n")
    cat("   Original BLUPs were double-shrunk; deregression corrects this.\n\n")
  } else {
    cat(">> Original BLUPs gave similar or higher PA.\n")
    cat("   With n=108, the PEV estimates may be noisy.\n")
    cat("   Deregression is still the theoretically correct approach.\n\n")
  }

  # Side-by-side box plot
  comp_cv <- bind_rows(
    mutate(orig_cv,  method = "Original BLUPs"),
    mutate(cv_store, method = "Deregressed BLUPs")
  )

  p_comp <- ggplot(comp_cv, aes(x = method, y = r, fill = method)) +
    geom_boxplot(alpha = 0.75, width = 0.45, colour = "#1B3A2D",
                 outlier.colour = "#B85042", outlier.size = 1.8) +
    geom_jitter(width = 0.08, size = 1.2, alpha = 0.35, colour = "#1B3A2D") +
    scale_fill_manual(
      values = c("Original BLUPs"     = "#5BAD7F",
                 "Deregressed BLUPs"  = "#2E6B4F"),
      guide = "none") +
    labs(title    = "CV Accuracy: Original vs Deregressed BLUPs",
         subtitle = paste0(N_FOLDS, "-fold CV × ", N_REPS, " reps"),
         x        = NULL,
         y        = "Pearson r (PA)",
         caption  = "") +
    theme_classic(base_size = 12) +
    theme(plot.title    = element_text(face = "bold", colour = "#1B3A2D"),
          plot.subtitle = element_text(colour = "#5C7A65"))

  pdf("results/plots/13_cv_original_vs_deregressed.pdf", width = 6, height = 5)
  print(p_comp)
  dev.off()
  cat("Plot saved: results/plots/13_cv_original_vs_deregressed.pdf\n\n")

} else {
  cat("Original CV results not found at", orig_cv_file, "\n")
  cat("Run CSES7160_Fritzner_Analysis.R Step 10 first if you want the comparison.\n\n")
}


# ==============================================================================
# 10 — SUMMARY
# ==============================================================================
cat("======================================================================\n")
cat("SUMMARY — DEREGRESSED BLUP GENOMIC PREDICTION\n")
cat("======================================================================\n")
cat(sprintf("Accessions (n)              : %d\n",  n_total))
cat(sprintf("SNPs                        : %d\n",  ncol(geno_num)))
cat(sprintf("sigma²_g                    : %.6f\n", Vu))
cat(sprintf("sigma²_e                    : %.6f\n", Ve))
cat(sprintf("h²_g                        : %.4f\n", h2g))
cat(sprintf("Mean reliability (r²)       : %.4f\n", mean(rel)))
cat(sprintf("Garrick c constant          : %.2f\n", c_const))
cat("----------------------------------------------------------------------\n")
cat(sprintf("CV mean r (deregressed)     : %.4f\n",
            mean(cv_store$r, na.rm = TRUE)))
cat(sprintf("CV 95%% CI                   : [%.4f, %.4f]\n",
            quantile(cv_store$r, 0.025, na.rm = TRUE),
            quantile(cv_store$r, 0.975, na.rm = TRUE)))
cat(sprintf("CV mean RMSE                : %.4f\n",
            mean(cv_store$RMSE, na.rm = TRUE)))
cat(sprintf("CV mean bias (slope)        : %.4f  (1.0 = unbiased)\n",
            mean(cv_store$bias, na.rm = TRUE)))
cat("----------------------------------------------------------------------\n")
cat("Output files:\n")
cat("  results/gp/deregressed_blups.csv\n")
cat("  results/gp/gblup_deregressed_cv.csv\n")
cat("  results/plots/12_gblup_deregressed_cv.pdf\n")
cat("  results/plots/13_cv_original_vs_deregressed.pdf  (if orig exists)\n")
cat("======================================================================\n")
cat("Done.\n")

tryCatch({
sink("results/gp/session_deregressed.txt")
cat("Deregressed BLUP Genomic Prediction \n")
cat("Date:", format(Sys.time(), "%Y-%m-%d %H:%M:%S"), "\n\n")
print(sessionInfo())
}, finally = {
  sink()
})
