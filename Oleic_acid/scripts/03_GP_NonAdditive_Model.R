# ==============================================================================
# Non-Additive Genomic Model — Additive + Dominance GBLUP (Independent Script)
# Genomic Dissection of Oleic Acid Concentration in Peanut
#
# Author  : Fritzner Pierre
# Course  : CSES 7160 — Genetic Data Analysis
#
# Revision applied:
#   The main script estimates only additive genomic heritability (h²_g) via
#   an additive GBLUP. When H² (broad-sense, from lmer) substantially exceeds
#   h²_g (marker-based additive), the gap suggests non-additive genetic
#   variance — dominance and/or epistasis — is contributing to the phenotype
#   but is not captured by the additive model.
#
#   For oleic acid in allotetraploid peanut, this is biologically plausible:
#   FAD2 alleles across the A- and B-subgenomes are known to act with partial
#   dominance, and complementary allele combinations can produce non-additive
#   effects on fatty acid desaturation.
#
#   This script:
#     1. Re-estimates h²_g from an additive GBLUP (replicates Step 8)
#     2. Computes the H² – h²_g gap as evidence for non-additive variance
#     3. Constructs the dominance relationship matrix G_D (Vitezica et al. 2013)
#     4. Fits additive + dominance GBLUP:
#          y = mu + g_a + g_d + e
#          g_a ~ N(0, G_A * sigma²_a)
#          g_d ~ N(0, G_D * sigma²_d)
#     5. Partitions variance into sigma²_a, sigma²_d, sigma²_e
#     6. Runs 5-fold CV × 10 reps comparing additive vs add+dom models
#     7. Produces a side-by-side PA comparison plot
#
# Requirements (produced by CSES7160_Fritzner_Analysis.R):
#   results/geno_num.rds                  — post-QC genotype matrix {0,1,2}
#   results/tables/BLUPs_OleicAcid.csv    — lmer BLUPs per accession
#   results/tables/heritability.csv       — h²_g and H² estimates
#   results/gp/deregressed_blups.csv      — deBLUPs (from GP_Deregressed_BLUPs.R)
#
# Outputs:
#   results/tables/variance_components_adddom.csv
#   results/gp/cv_additive_vs_dominance.csv
#   results/plots/14_cv_additive_vs_dominance.pdf
#   results/plots/15_variance_partition.pdf
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
# install.packages(c("rrBLUP", "sommer", "ggplot2", "dplyr", "patchwork"))
# sommer is used for the two-kernel GBLUP (additive + dominance).
# rrBLUP is kept for A.mat() and the single-kernel additive model.

pkgs <- c("rrBLUP","sommer","ggplot2","dplyr","patchwork")

missing <- pkgs[!sapply(pkgs, requireNamespace, quietly = TRUE)]

if (length(missing) > 0) {
  stop("Missing packages: ", paste(missing, collapse = ", "))
}

suppressPackageStartupMessages(
  lapply(pkgs, library, character.only = TRUE)
)

set.seed(2026)

for (d in c("results/gp", "results/plots", "results/tables"))
  dir.create(d, showWarnings = FALSE, recursive = TRUE)

cat("======================================================================\n")
cat(" Non-Additive Genomic Model \n")
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
  stop("results/geno_num.rds not found.\nRun CSES7160_Fritzner_Analysis.R first.")
geno_num <- readRDS("results/geno_num.rds")
cat(" Loaded: results/geno_num.rds —", nrow(geno_num), "taxa x",
    ncol(geno_num), "SNPs\n")

# lmer BLUPs
blup_file <- "results/tables/BLUPs_OleicAcid.csv"
if (!file.exists(blup_file))
  stop(blup_file, " not found.\nRun CSES7160_Fritzner_Analysis.R first.")
blup_df <- read.csv(blup_file, stringsAsFactors = FALSE)
cat(" Loaded:", blup_file, "—", nrow(blup_df), "accessions\n")

# Heritability table (for H²)
h2_table <- NULL
if (file.exists("results/tables/heritability.csv")) {
  h2_table <- read.csv("results/tables/heritability.csv",
                       stringsAsFactors = FALSE)
  cat(" Loaded: results/tables/heritability.csv\n")
} else {
  cat(" Warning: heritability.csv not found — H² will not be available for gap check.\n")
}

# Deregressed BLUPs (from GP_Deregressed_BLUPs.R)
deblup_df <- NULL
if (file.exists("results/gp/deregressed_blups.csv")) {
  deblup_df <- read.csv("results/gp/deregressed_blups.csv",
                        stringsAsFactors = FALSE)
  cat(" Loaded: results/gp/deregressed_blups.csv —",
      nrow(deblup_df), "accessions\n")
} else {
  cat(" Warning: deregressed_blups.csv not found.\n")
  cat("   Run GP_Deregressed_BLUPs.R first.\n")
  cat("   Falling back to raw BLUPs for CV (not ideal — see revision notes).\n")
}
cat("\n")

# Align taxa
common_taxa  <- intersect(blup_df$Taxa, rownames(geno_num))
blup_aligned <- blup_df[match(common_taxa, blup_df$Taxa), ]
geno_num     <- geno_num[common_taxa, , drop = FALSE]
y            <- blup_aligned$BLUP
cat("Common taxa (pheno ∩ geno):", length(common_taxa), "\n\n")

geno_rrblup <- geno_num - 1   # {-1, 0, 1} for rrBLUP


# ==============================================================================
# 3 — ADDITIVE RELATIONSHIP MATRIX (VanRaden 2008)
# ==============================================================================
cat("----------------------------------------------------------------------\n")
cat("STEP 3 | Build G_A (Additive)\n")
cat("----------------------------------------------------------------------\n")

G_A <- A.mat(geno_rrblup, min.MAF = 0)
cat("G_A:", nrow(G_A), "x", ncol(G_A), "\n\n")


# ==============================================================================
# 4 — DOMINANCE RELATIONSHIP MATRIX (Vitezica et al. 2013)
# ==============================================================================
# For {0,1,2} coded SNPs (0=AA, 1=AB, 2=BB), the dominance deviation
# at locus j for individual i is:
#   d_ij = { -2q²_j    if genotype = 0
#          {  2p_j q_j  if genotype = 1
#          { -2p²_j     if genotype = 2
# where p_j = freq of allele B, q_j = 1 - p_j.
#
# G_D = D D' / sum_j(2 p_j q_j)² * 4  (Vitezica eq. 3 / Su et al. 2012)
# Diagonal values ≈ 1, off-diagonal = dominance relationship coefficients.

cat("----------------------------------------------------------------------\n")
cat("STEP 4 | Build G_D (Dominance, Vitezica et al. 2013)\n")
cat("----------------------------------------------------------------------\n")

build_dominance_matrix <- function(G012) {
  # G012: taxa × SNPs matrix, coded {0, 1, 2}
  n   <- nrow(G012)
  m   <- ncol(G012)

  # Allele frequencies (freq of allele coded as 2)
  p   <- colMeans(G012, na.rm = TRUE) / 2    # freq of B allele
  q   <- 1 - p

  # Dominance coding matrix D: n × m
  D <- matrix(0, nrow = n, ncol = m)
  for (j in seq_len(m)) {
    pj <- p[j]; qj <- q[j]
    D[G012[, j] == 0, j] <- -2 * qj^2
    D[G012[, j] == 1, j] <-  2 * pj * qj
    D[G012[, j] == 2, j] <- -2 * pj^2
  }

  # Scale factor (sum of squared heterozygosity)
  scale_f <- sum(4 * (p * q)^2)

  G_D <- tcrossprod(D) / scale_f
  rownames(G_D) <- colnames(G_D) <- rownames(G012)
  G_D
}

G_D <- build_dominance_matrix(geno_num)

# Add small ridge to ensure positive definiteness
G_D_pd <- G_D + diag(1e-4, nrow(G_D))

cat("G_D:", nrow(G_D), "x", ncol(G_D), "\n")
cat("G_D diagonal mean :", round(mean(diag(G_D)), 3), "\n")
cat("G_D off-diag mean :", round(mean(G_D[lower.tri(G_D)]), 4), "\n\n")


# ==============================================================================
# 5 — EVIDENCE FOR NON-ADDITIVITY: H² VS h²_g GAP
# ==============================================================================
cat("----------------------------------------------------------------------\n")
cat("STEP 5 | Evidence for Non-Additivity (H² – h²_g Gap)\n")
cat("----------------------------------------------------------------------\n")

# Re-estimate additive h²_g via rrBLUP for reproducibility
lmm_add <- mixed.solve(y = y, K = G_A, SE = FALSE)
h2g_add  <- lmm_add$Vu / (lmm_add$Vu + lmm_add$Ve)
cat(sprintf("Additive sigma²_g : %.6f\n", lmm_add$Vu))
cat(sprintf("Additive sigma²_e : %.6f\n", lmm_add$Ve))
cat(sprintf("h²_g (additive)   : %.4f\n", h2g_add))

# Load H² (broad-sense) from lmer output — NOT from heritability.csv
# heritability.csv contains h²_g (genomic), not H² (broad-sense)
H2 <- NA_real_
h2_broad_file <- "results/tables/H2_broad_lmer.csv"
if (file.exists(h2_broad_file)) {
  h2_broad_table <- read.csv(h2_broad_file, stringsAsFactors = FALSE)
  h2_row <- h2_broad_table[h2_broad_table$Parameter == "H2_broad", ]
  if (nrow(h2_row) > 0) {
    H2 <- as.numeric(h2_row$Estimate[1])
    cat(sprintf("H²  (broad-sense) : %.4f  (from H2_broad_lmer.csv)\n", H2))
  }
} else {
  cat("Warning: H2_broad_lmer.csv not found.\n")
  cat("  Run CSES7160_Fritzner_Analysis.R Step 3 first.\n\n")
}

if (!is.na(H2)) {
  gap <- H2 - h2g_add
  cat(sprintf("Gap H² – h²_g     : %.4f\n\n", gap))
  if (gap > 0.10) {
    cat(">> Gap > 0.10 — substantial non-additive variance detected.\n")
    cat("   Fitting additive + dominance model (Step 6).\n\n")
  } else {
    cat(">> Gap <= 0.10 — additive model is sufficient.\n")
    cat("   Fitting dominance model for completeness only.\n\n")
  }
} else {
  cat("H² not available — fitting non-additive model for completeness.\n\n")
}


# ==============================================================================
# 6 — FIT ADDITIVE + DOMINANCE GBLUP (sommer::mmer)
# ==============================================================================
cat("----------------------------------------------------------------------\n")
cat("STEP 6 | Fit Additive + Dominance GBLUP (sommer::mmer)\n")
cat("----------------------------------------------------------------------\n")

# Prepare data frame for sommer
pheno_sommer <- data.frame(
  Taxa  = factor(common_taxa),
  oleic = y,
  stringsAsFactors = FALSE
)
rownames(pheno_sommer) <- common_taxa

# Ensure G_A and G_D row/col names match Taxa factor levels
rownames(G_A)    <- colnames(G_A)    <- common_taxa
rownames(G_D_pd) <- colnames(G_D_pd) <- common_taxa

cat("Fitting two-kernel model via REML (may take ~1-2 minutes)...\n")

fit_adddom <- mmer(
  fixed   = oleic ~ 1,
  random  = ~ vsr(Taxa, Gu = G_A) +
    vsr(Taxa, Gu = G_D_pd),
  rcov    = ~ units,
  data    = pheno_sommer,
  verbose = FALSE
)

# Extract variance components (with safety check for sommer version differences)
sigma_a2 <- as.numeric(fit_adddom$sigma[[1]])
sigma_d2 <- as.numeric(fit_adddom$sigma[[2]])
sigma_e2 <- if ("units" %in% names(fit_adddom$sigma)) {
  as.numeric(fit_adddom$sigma$units)
} else {
  as.numeric(fit_adddom$sigma[[3]])
}
Vp_ad        <- sigma_a2 + sigma_d2 + sigma_e2
h2_add_ad    <- sigma_a2 / Vp_ad
h2_total_ad  <- (sigma_a2 + sigma_d2) / Vp_ad
dom_ratio    <- sigma_d2 / Vp_ad

cat(sprintf("\nVariance components (Additive + Dominance model):\n"))
cat(sprintf("  sigma²_a             = %.6f\n", sigma_a2))
cat(sprintf("  sigma²_d             = %.6f\n", sigma_d2))
cat(sprintf("  sigma²_e             = %.6f\n", sigma_e2))
cat(sprintf("  Vp                   = %.6f\n", Vp_ad))
cat(sprintf("  h² (additive only)   = %.4f\n", h2_add_ad))
cat(sprintf("  h² (add + dom)       = %.4f\n", h2_total_ad))
cat(sprintf("  Dominance ratio      = %.4f  (Vd/Vp)\n", dom_ratio))

# Interpret dominance ratio
if (dom_ratio > 0.05) {
  cat(sprintf("\n>> Vd/Vp = %.4f — dominance variance is non-trivial.\n",
              dom_ratio))
  cat("   Non-additive model is warranted.\n\n")
} else {
  cat(sprintf("\n>> Vd/Vp = %.4f — dominance variance is negligible.\n",
              dom_ratio))
  cat("   Additive model is sufficient. Results retained as evidence.\n\n")
}

# Save variance components
varcomp_df <- data.frame(
  Parameter = c("sigma2_a", "sigma2_d", "sigma2_e", "Vp",
                "h2_additive", "h2_add_plus_dom", "dominance_ratio",
                "h2g_additive_only", "H2_broad_sense"),
  Estimate  = round(c(sigma_a2, sigma_d2, sigma_e2, Vp_ad,
                      h2_add_ad, h2_total_ad, dom_ratio,
                      h2g_add, H2), 6),
  Description = c(
    "Additive genomic variance (sommer)",
    "Dominance genomic variance (sommer)",
    "Residual variance (sommer)",
    "Total phenotypic variance",
    "h2 from additive component only",
    "h2 from additive + dominance",
    "Dominance share of Vp (Vd/Vp)",
    "h2_g additive GBLUP (rrBLUP)",
    "Broad-sense H2 (lmer, from H2_broad_lmer.csv)"
  ),
  stringsAsFactors = FALSE
)

write.csv(varcomp_df,
          "results/tables/variance_components_adddom.csv",
          row.names = FALSE)
cat("Table saved: results/tables/variance_components_adddom.csv\n\n")


# ==============================================================================
# 7 — 5-FOLD CV × 10 REPS: ADDITIVE vs ADDITIVE + DOMINANCE
# ==============================================================================
# Both models use deregressed BLUPs as the response (if available),
# otherwise fall back to raw BLUPs.

cat("----------------------------------------------------------------------\n")
cat("STEP 7 | 5-Fold CV × 10 Reps: Additive vs Add+Dom\n")
cat("----------------------------------------------------------------------\n")

# Select response
if (!is.null(deblup_df)) {
  deblup_aligned <- deblup_df[match(common_taxa, deblup_df$Taxa), ]
  y_cv_response  <- deblup_aligned$deBLUP
  w_i            <- deblup_aligned$weight
  W_sqrt         <- sqrt(w_i)
  cat("Using deregressed BLUPs as response (recommended).\n\n")
} else {
  y_cv_response  <- y
  w_i            <- rep(1, length(y))
  W_sqrt         <- rep(1, length(y))
  cat("Using raw BLUPs as response (deregressed_blups.csv not found).\n\n")
}

N_FOLDS <- 5
N_REPS  <- 10
n_total <- length(common_taxa)

cv_add <- data.frame(rep = integer(), fold = integer(), r = numeric(), model = character())
cv_dom <- data.frame(rep = integer(), fold = integer(), r = numeric(), model = character())

for (rep in seq_len(N_REPS)) {
  fid <- sample(rep(seq_len(N_FOLDS), length.out = n_total))

  for (fold in seq_len(N_FOLDS)) {
    test_idx <- which(fid == fold)

    # ── Model A: Additive ──────────────────────────────────────────────────
    y_a         <- y_cv_response * W_sqrt
    y_a[test_idx] <- NA
    K_sc        <- diag(W_sqrt) %*% G_A %*% diag(W_sqrt)

    fv_a <- tryCatch(
      mixed.solve(y = y_a, K = K_sc),
      error = function(e) NULL
    )

    if (!is.null(fv_a)) {
      pred_a <- (as.numeric(fv_a$u) + as.numeric(fv_a$beta))[test_idx] /
                W_sqrt[test_idx]
      obs_t  <- y_cv_response[test_idx]
      cv_add <- rbind(cv_add, data.frame(
        rep   = rep, fold = fold,
        r     = cor(obs_t, pred_a, use = "complete.obs"),
        model = "Additive"
      ))
    }

    # ── Model B: Additive + Dominance ─────────────────────────────────────
    pheno_cv              <- pheno_sommer
    pheno_cv$oleic        <- y_cv_response
    pheno_cv$oleic[test_idx] <- NA

    fv_d <- tryCatch(
      mmer(
        fixed   = oleic ~ 1,
        random  = ~ vsr(Taxa, Gu = G_A) +
                    vsr(Taxa, Gu = G_D_pd),
        rcov    = ~ units,
        data    = pheno_cv,
        verbose = FALSE
      ),
      error = function(e) NULL
    )

    if (!is.null(fv_d) && length(fv_d$U) >= 2) {
      # Sum additive and dominance GEBVs
      u_names <- names(fv_d$U)
      u_a    <- as.numeric(fv_d$U[[u_names[1]]]$oleic)
      u_d    <- as.numeric(fv_d$U[[u_names[2]]]$oleic)
      mu_hat <- as.numeric(fv_d$Beta$Estimate)
      pred_d <- (u_a + u_d + mu_hat)[test_idx]
      obs_t  <- y_cv_response[test_idx]

      cv_dom <- rbind(cv_dom, data.frame(
        rep   = rep, fold = fold,
        r     = cor(obs_t, pred_d, use = "complete.obs"),
        model = "Additive + Dominance"
      ))
    }
  }

  if (rep %% 2 == 0)
    cat(sprintf("  Rep %2d / %d  |  Add r = %.4f  |  Add+Dom r = %.4f\n",
                rep, N_REPS,
                mean(cv_add$r[cv_add$rep <= rep], na.rm = TRUE),
                mean(cv_dom$r[cv_dom$rep <= rep], na.rm = TRUE)))
}

cat("\n--- CV Summary ---\n")
cat(sprintf("Additive only  — mean r: %.4f  (SD: %.4f)\n",
            mean(cv_add$r, na.rm = TRUE), sd(cv_add$r, na.rm = TRUE)))
cat(sprintf("Add + Dom      — mean r: %.4f  (SD: %.4f)\n",
            mean(cv_dom$r, na.rm = TRUE), sd(cv_dom$r, na.rm = TRUE)))
cat(sprintf("PA gain from dominance : %+.4f\n\n",
            mean(cv_dom$r, na.rm = TRUE) - mean(cv_add$r, na.rm = TRUE)))

# Save
cv_all <- bind_rows(cv_add, cv_dom)
write.csv(cv_all, "results/gp/cv_additive_vs_dominance.csv", row.names = FALSE)
cat("Table saved: results/gp/cv_additive_vs_dominance.csv\n\n")


# ==============================================================================
# 8 — PLOTS
# ==============================================================================
cat("----------------------------------------------------------------------\n")
cat("STEP 8 | Plots\n")
cat("----------------------------------------------------------------------\n")

# -- CV comparison box plot ---------------------------------------------------
p_cv_comp <- ggplot(cv_all, aes(x = model, y = r, fill = model)) +
  geom_boxplot(alpha = 0.78, width = 0.45, colour = "#1B3A2D",
               outlier.colour = "#B85042", outlier.size = 1.8) +
  geom_jitter(width = 0.08, size = 1.2, alpha = 0.35, colour = "#1B3A2D") +
  scale_fill_manual(
    values = c("Additive"              = "#2E6B4F",
               "Additive + Dominance"  = "#5BAD7F"),
    guide = "none") +
  labs(title    = "Prediction Accuracy: Additive vs Additive + Dominance",
       subtitle = paste0(N_FOLDS, "-fold CV × ", N_REPS, " reps  |  ",
                         "Response: ",
                         ifelse(!is.null(deblup_df),
                                "deregressed BLUPs", "raw BLUPs")),
       x        = NULL,
       y        = "Pearson r (PA)",
       caption  = "") +
  theme_classic(base_size = 12) +
  theme(plot.title    = element_text(face = "bold", colour = "#1B3A2D"),
        plot.subtitle = element_text(colour = "#5C7A65"))

# -- Variance partition bar chart ---------------------------------------------
# Build Vp = Va + Vd + Ve for both models
vc_df <- data.frame(
  Component = c("Additive (Va)", "Dominance (Vd)", "Residual (Ve)"),
  Variance  = c(sigma_a2, sigma_d2, sigma_e2),
  Share     = round(c(sigma_a2, sigma_d2, sigma_e2) / Vp_ad * 100, 1)
)
vc_df$Component <- factor(vc_df$Component,
                           levels = c("Residual (Ve)", "Dominance (Vd)",
                                      "Additive (Va)"))

p_vpart <- ggplot(vc_df, aes(x = "", y = Variance, fill = Component)) +
  geom_col(width = 0.5, colour = "#1B3A2D", alpha = 0.85) +
  geom_text(aes(label = paste0(Share, "%")),
            position = position_stack(vjust = 0.5),
            colour = "white", fontface = "bold", size = 4.5) +
  scale_fill_manual(
    values = c("Additive (Va)"   = "#2E6B4F",
               "Dominance (Vd)"  = "#C7A84F",
               "Residual (Ve)"   = "#AAAAAA"),
    name = "Variance component") +
  labs(title    = "Variance Partition — Add + Dom GBLUP",
       subtitle = paste0("h²_add = ", round(h2_add_ad, 3),
                         "  |  h²_total = ", round(h2_total_ad, 3),
                         "  |  Vd/Vp = ", round(dom_ratio, 3)),
       x = NULL, y = "Variance",
       caption = "") +
  theme_classic(base_size = 12) +
  theme(plot.title    = element_text(face = "bold", colour = "#1B3A2D"),
        plot.subtitle = element_text(colour = "#5C7A65"),
        axis.text.x   = element_blank(),
        axis.ticks.x  = element_blank())

p_combined <- (p_cv_comp | p_vpart) +
  plot_layout(widths = c(2, 1)) +        # CV gets 2/3, variance gets 1/3
  plot_annotation(
    title   = "Non-Additive Genomic Model — Oleic Acid in Peanut",
    caption = "sommer::mmer | Vitezica et al. (2013)",
    theme   = theme(
      plot.title   = element_text(face = "bold", size = 13, colour = "#1B3A2D"),
      plot.caption = element_text(colour = "#5C7A65", size = 8))
  )

pdf("results/plots/14_cv_additive_vs_dominance.pdf", width = 12, height = 5.5)
print(p_combined)
dev.off()
cat("Plot saved: results/plots/14_cv_additive_vs_dominance.pdf\n\n")

pdf("results/plots/15_variance_partition.pdf", width = 6, height = 5)
print(p_vpart)
dev.off()
cat("Plot saved: results/plots/15_variance_partition.pdf\n\n")


# ==============================================================================
# 9 — SUMMARY
# ==============================================================================
cat("======================================================================\n")
cat("SUMMARY — NON-ADDITIVE GENOMIC MODEL\n")
cat("======================================================================\n")
cat(sprintf("Accessions (n)              : %d\n",  length(common_taxa)))
cat(sprintf("SNPs                        : %d\n",  ncol(geno_num)))
cat("----------------------------------------------------------------------\n")
cat(sprintf("h²_g (additive, rrBLUP)     : %.4f\n", h2g_add))
if (!is.na(H2)) {
  cat(sprintf("H²  (broad-sense, lmer)     : %.4f\n", H2))
  cat(sprintf("H² – h²_g gap               : %.4f\n", H2 - h2g_add))
} else {
  cat("H²  (broad-sense, lmer)     : not available\n")
  cat("H² – h²_g gap               : not available\n")
}
cat("----------------------------------------------------------------------\n")
cat(sprintf("sigma²_a                    : %.6f\n", sigma_a2))
cat(sprintf("sigma²_d                    : %.6f\n", sigma_d2))
cat(sprintf("sigma²_e                    : %.6f\n", sigma_e2))
cat(sprintf("h² (additive only)          : %.4f\n", h2_add_ad))
cat(sprintf("h² (additive + dominance)   : %.4f\n", h2_total_ad))
cat(sprintf("Dominance ratio (Vd/Vp)     : %.4f\n", dom_ratio))
cat("----------------------------------------------------------------------\n")
cat(sprintf("CV PA — Additive only        : %.4f  (SD: %.4f)\n",
            mean(cv_add$r, na.rm = TRUE), sd(cv_add$r, na.rm = TRUE)))
cat(sprintf("CV PA — Additive + Dominance : %.4f  (SD: %.4f)\n",
            mean(cv_dom$r, na.rm = TRUE), sd(cv_dom$r, na.rm = TRUE)))
cat(sprintf("PA gain from dominance       : %+.4f\n",
            mean(cv_dom$r, na.rm = TRUE) - mean(cv_add$r, na.rm = TRUE)))
cat("----------------------------------------------------------------------\n")
cat("CONCLUSION:\n")
if (dom_ratio <= 0.05) {
  cat("  Dominance variance is negligible (Vd/Vp <= 0.05).\n")
  cat("  PA gain from adding dominance is minimal.\n")
  cat("  >> Additive model is sufficient for this dataset.\n")
  cat("     Results above are retained as evidence.\n")
} else {
  cat("  Dominance variance is non-trivial (Vd/Vp > 0.05).\n")
  cat("  >> Additive + dominance model is warranted.\n")
}
cat("----------------------------------------------------------------------\n")
cat("Output files:\n")
cat("  results/tables/variance_components_adddom.csv\n")
cat("  results/gp/cv_additive_vs_dominance.csv\n")
cat("  results/plots/14_cv_additive_vs_dominance.pdf\n")
cat("  results/plots/15_variance_partition.pdf\n")
cat("======================================================================\n")
cat("Done.\n")

tryCatch({
  sink("results/gp/session_nonadditive.txt")
  cat("Non-Additive Genomic Model | CSES 7160 | Fritzner Pierre\n")
  cat("Date:", format(Sys.time(), "%Y-%m-%d %H:%M:%S"), "\n\n")
  print(sessionInfo())
}, finally = {
  sink()
})
