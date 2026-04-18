# ==============================================================================
# CSES 7160 — FarmCPU vs. MLM GWAS Comparison (INDEPENDENT SCRIPT)
# Oleic Acid Concentration in Peanut (Arachis hypogaea)
#
# Author  : Fritzner Pierre
# Course  : CSES 7160 — Genetic Data Analysis
# Date    : 2026
#
# Description:
#   Fully self-contained script. Reads the same three input files as the main
#   analysis script, rebuilds all required objects from scratch, then runs
#   both FarmCPU and MLM via GAPIT3 and produces a full comparison:
#
#     1. Genomic inflation (lambda_GC) — type I error control
#     2. Genome-wide and suggestive hit counts
#     3. SNP overlap between models (shared vs. model-specific)
#     4. Effect size (Beta) direction concordance for shared hits
#     5. -log10(p) correlation across all SNPs
#     6. Side-by-side Manhattan + QQ plots (PDF)
#     7. FarmCPU vs. MLM -log10(p) scatter plot (PDF)
#     8. Summary comparison table (CSV)
#
# Required input files (same as main analysis — set paths in USER SETTINGS):
#   values.csv                                     (phenotype)
#   arahy.Tifrunner.gnm1.mrk.Axiom_Arachis_58K.gff3  (SNP position map)
#   aradu1_araip1.gnm1.div.Otyama_Kulkarni_2020.main.vcf  (VCF genotype)
#
# Outputs written to:  <PROJECT_DIR>/results/gwas_comparison/
# ==============================================================================


# ==============================================================================
# SECTION 0 — USER SETTINGS
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

PHENO_FILE  <- "data/values.csv"
GFF3_FILE   <- "data/arahy.Tifrunner.gnm1.mrk.Axiom_Arachis_58K.gff3"
GENO_FILE   <- "data/aradu1_araip1.gnm1.div.Otyama_Kulkarni_2020.main.vcf"

TRAIT       <- "OleicAcid"
MAF_THOLD   <- 0.05     # minor allele frequency filter
MISS_THOLD  <- 0.10     # SNP missingness filter
N_PCS       <- 3        # PCs to use as covariates

set.seed(2026)


# ==============================================================================
# SECTION 1 — PACKAGES
# ==============================================================================

pkgs <- c(
  "vcfR", "rrBLUP", "ggplot2","dplyr","patchwork",
  "lme4"
)

missing <- pkgs[!sapply(pkgs, requireNamespace, quietly = TRUE)]

if (length(missing) > 0) {
  stop("Missing packages: ", paste(missing, collapse = ", "),
       "\nInstall them before running the script.")
}

suppressPackageStartupMessages(
  lapply(pkgs, library, character.only = TRUE)
)

# ---- Ensure devtools exists ----
if (!requireNamespace("devtools", quietly = TRUE))
  install.packages("devtools", repos = "https://cloud.r-project.org")

# ---- Install GAPIT3 from GitHub if missing ----
if (!requireNamespace("GAPIT", quietly = TRUE))
  devtools::install_github("jiabowang/GAPIT3", force = TRUE)

suppressPackageStartupMessages(library(GAPIT))


# Output directories
for (d in c("results/gwas_comparison",
            "results/gwas_comparison/farmcpu",
            "results/gwas_comparison/mlm",
            "results/gwas_comparison/plots",
            "results/gwas_comparison/tables"))
  dir.create(d, showWarnings = FALSE, recursive = TRUE)

cat("======================================================================\n")
cat(" FarmCPU vs. MLM — GWAS Comparison \n")
cat("======================================================================\n")
cat("R version :", R.version$version.string, "\n")
cat("Date/Time :", format(Sys.time(), "%Y-%m-%d %H:%M:%S"), "\n\n")


# ==============================================================================
# SECTION 2 — DATA IMPORT
# ==============================================================================

cat("----------------------------------------------------------------------\n")
cat("SECTION 2 | Data Import\n")
cat("----------------------------------------------------------------------\n")

# ── 2.1 Phenotype ─────────────────────────────────────────────────────────────
stopifnot(file.exists(PHENO_FILE))
pheno_raw <- read.csv(PHENO_FILE, stringsAsFactors = FALSE)
pheno_raw <- rename(pheno_raw, Taxa = accession_name, Value = phenotype_value)
cat("Phenotype rows      :", nrow(pheno_raw), "\n")
cat("Unique accessions   :", length(unique(pheno_raw$Taxa)), "\n")

# ── 2.2 GFF3 SNP position map ─────────────────────────────────────────────────
stopifnot(file.exists(GFF3_FILE))
gff3_raw <- readLines(GFF3_FILE)
gff3_raw <- gff3_raw[!grepl("^#", gff3_raw) & nchar(gff3_raw) > 0]

gff3_df <- do.call(rbind, lapply(gff3_raw, function(ln) {
  p   <- strsplit(ln, "\t")[[1]]
  chr <- as.integer(sub(".*Arahy\\.(\\d+)$", "\\1", p[1]))
  pos <- as.integer(p[4])
  snp <- sub(".*Name=([^;[:space:]]+).*", "\\1", p[9])
  data.frame(SNP = trimws(snp), Chr = chr, Pos = pos, stringsAsFactors = FALSE)
}))
cat("GFF3 markers loaded :", nrow(gff3_df), "on",
    length(unique(gff3_df$Chr)), "chromosomes\n")

# ── 2.3 Genotype VCF ──────────────────────────────────────────────────────────
stopifnot(file.exists(GENO_FILE))
cat("Reading VCF (may take a minute)...\n")
vcf_raw <- read.vcfR(GENO_FILE, verbose = FALSE)
cat("VCF loaded:", nrow(vcf_raw@fix), "variants x",
    ncol(vcf_raw@gt) - 1, "samples\n")

# Recode GT to allele dosage {0, 1, 2}
gt_mat <- extract.gt(vcf_raw, element = "GT",
                     as.numeric = FALSE, return.alleles = FALSE)
recode_gt <- function(x) {
  x[x == "./."] <- NA
  x[x == "0/0"] <- "0"
  x[x == "0/1" | x == "1/0"] <- "1"
  x[x == "1/1"] <- "2"
  as.numeric(x)
}
gt_num   <- apply(gt_mat, 2, recode_gt)
geno_num <- t(gt_num)
rownames(geno_num) <- colnames(gt_mat)

# Build SNP position map from VCF
chr_raw <- vcf_raw@fix[, "CHROM"]
pos_raw <- as.integer(vcf_raw@fix[, "POS"])
snp_ids <- paste0(chr_raw, "_", pos_raw)

# Convert chromosome strings to integers:
#   Aradu.A01-A10 -> 1-10  (A-subgenome, A. duranensis)
#   Araip.B01-B10 -> 11-20 (B-subgenome, A. ipaensis)
# IMPORTANT: ifelse must test the INPUT vector (ch_vec), not close over chr_raw.
chr_to_int <- function(ch_vec) {
  is_b <- grepl("Araip", ch_vec)
  num  <- suppressWarnings(
    as.integer(sub("^.*[AB]0?", "", ch_vec))  # strip prefix, keep trailing digits
  )
  ifelse(is_b, num + 10L, num)
}
chr_int <- chr_to_int(chr_raw)

geno_map <- data.frame(SNP = snp_ids, Chr = chr_int, Pos = pos_raw,
                       stringsAsFactors = FALSE)
colnames(geno_num) <- snp_ids

# Drop SNPs with NA chromosome (unmapped scaffold contigs)
na_chr <- is.na(geno_map$Chr)
if (any(na_chr)) {
  cat("Dropping", sum(na_chr), "SNPs with NA chromosome.\n")
  geno_map <- geno_map[!na_chr, ]
  geno_num <- geno_num[, !na_chr, drop = FALSE]
}
cat("Chromosomes present:", sort(unique(geno_map$Chr)), "\n")

rm(vcf_raw, gt_mat, gt_num); gc()

# Normalise VCF sample names -> "PI XXXXXX" format
vcf_taxa <- rownames(geno_num)
vcf_taxa <- sub("_[0-9]+$", "",    vcf_taxa)
vcf_taxa <- sub("_s$",       "",    vcf_taxa)
vcf_taxa <- sub("^PI([0-9]+)$", "PI \\1", vcf_taxa)
rownames(geno_num) <- vcf_taxa
cat("Taxa matched pheno/geno:",
    length(intersect(vcf_taxa, pheno_raw$Taxa)), "\n\n")


# ==============================================================================
# SECTION 3 — PHENOTYPE PRE-PROCESSING (BLUP EXTRACTION)
# ==============================================================================

cat("----------------------------------------------------------------------\n")
cat("SECTION 3 | Phenotype Pre-processing\n")
cat("----------------------------------------------------------------------\n")

# Outlier replicate removal (|z_within| > 3)
pheno_flagged <- pheno_raw |>
  group_by(Taxa) |>
  mutate(z_within = if (n() > 1) as.numeric(scale(Value)) else 0) |>
  ungroup()
pheno_clean <- filter(pheno_flagged, abs(z_within) <= 3)
cat("Outlier replicates removed:",
    nrow(pheno_raw) - nrow(pheno_clean), "\n")

# BLUP extraction: Value ~ 1 + (1 | Taxa)
cat("Extracting BLUPs via lmer...\n")
lmer_fit  <- lmer(Value ~ 1 + (1 | Taxa), data = pheno_clean, REML = TRUE)
blup_list <- ranef(lmer_fit)$Taxa
mu        <- as.numeric(fixef(lmer_fit))

blup_df <- data.frame(
  Taxa = rownames(blup_list),
  BLUP = blup_list[, "(Intercept)"],
  stringsAsFactors = FALSE
)
cat("Grand mean (mu) :", round(mu, 3), "%\n")
cat("BLUPs extracted :", nrow(blup_df), "accessions\n")

# Align phenotype and genotype
common_taxa  <- intersect(blup_df$Taxa, rownames(geno_num))
blup_aligned <- blup_df[match(common_taxa, blup_df$Taxa), ]
geno_num     <- geno_num[common_taxa, , drop = FALSE]
y            <- blup_aligned$BLUP
cat("Common taxa     :", length(common_taxa), "\n\n")


# ==============================================================================
# SECTION 4 — GENOTYPE QC
# ==============================================================================

cat("----------------------------------------------------------------------\n")
cat("SECTION 4 | Genotype QC\n")
cat("----------------------------------------------------------------------\n")
cat("Input:", nrow(geno_num), "x", ncol(geno_num), "\n")

# Missingness filter
miss_rate <- colMeans(is.na(geno_num))
keep_miss <- miss_rate <= MISS_THOLD
cat("Removed (miss >",  MISS_THOLD * 100, "%):", sum(!keep_miss), "\n")
geno_num <- geno_num[, keep_miss, drop = FALSE]
geno_map <- geno_map[keep_miss, ]

# MAF filter
maf_vals <- apply(geno_num, 2, function(x) {
  p <- mean(x, na.rm = TRUE) / 2; min(p, 1 - p)
})
keep_maf <- maf_vals >= MAF_THOLD
cat("Removed (MAF <", MAF_THOLD, ")    :", sum(!keep_maf), "\n")
geno_num <- geno_num[, keep_maf, drop = FALSE]
geno_map <- geno_map[keep_maf, ]

# Mean imputation
n_miss <- sum(is.na(geno_num))
if (n_miss > 0) {
  col_mu <- colMeans(geno_num, na.rm = TRUE)
  for (j in seq_len(ncol(geno_num))) {
    idx <- is.na(geno_num[, j])
    if (any(idx)) geno_num[idx, j] <- col_mu[j]
  }
}
cat("Imputed missing  :", n_miss, "\n")
cat("Final matrix     :", nrow(geno_num), "x", ncol(geno_num), "\n\n")

geno_rrblup <- geno_num - 1   # {-1, 0, 1} for rrBLUP


# ==============================================================================
# SECTION 5 — POPULATION STRUCTURE (PCA) AND GRM
# ==============================================================================
cat("----------------------------------------------------------------------\n")
cat("SECTION 5 | PCA + GRM\n")
cat("----------------------------------------------------------------------\n")

geno_sc <- scale(geno_rrblup, center = TRUE, scale = TRUE)
pca_res <- prcomp(geno_sc, retx = TRUE, center = FALSE, scale. = FALSE)
pve     <- (pca_res$sdev^2) / sum(pca_res$sdev^2) * 100
cat(sprintf("PC1: %.2f%%  PC2: %.2f%%  PC3: %.2f%%\n",
            pve[1], pve[2], pve[3]))

scores      <- as.data.frame(pca_res$x[, 1:10])
scores$Taxa <- rownames(geno_rrblup)

# Covariate matrix for GAPIT (PC1-PC3)
CV <- scores[, c("Taxa", paste0("PC", seq_len(N_PCS)))]

# Genomic relationship matrix (VanRaden 2008)
cat("Building GRM...\n")
# Manual VanRaden GRM — avoids A.mat() C-level error at small n
Z_grm   <- geno_rrblup
p_grm   <- (colMeans(Z_grm) + 1) / 2
Z_c_grm <- sweep(Z_grm, 2, 2 * (p_grm - 0.5))
K_mat   <- tcrossprod(Z_c_grm) / (2 * sum(p_grm * (1 - p_grm)))
rownames(K_mat) <- colnames(K_mat) <- rownames(geno_rrblup)
K_mat   <- K_mat + diag(1e-4, nrow(K_mat))
cat("GRM          :", nrow(K_mat), "x", ncol(K_mat), "\n")
cat("Diagonal mean:", round(mean(diag(K_mat)), 3), "\n\n")


# ==============================================================================
# SECTION 6 — PREPARE GAPIT INPUTS
# ==============================================================================

cat("----------------------------------------------------------------------\n")
cat("SECTION 6 | Prepare GAPIT inputs\n")
cat("----------------------------------------------------------------------\n")

# Phenotype data.frame
pheno_gapit <- data.frame(
  Taxa      = blup_aligned$Taxa,
  OleicAcid = blup_aligned$BLUP,
  stringsAsFactors = FALSE
)

# Genotype matrix: enforce 0/1/2 integers
GD <- matrix(
  as.integer(round(geno_num)),
  nrow = nrow(geno_num), ncol = ncol(geno_num),
  dimnames = dimnames(geno_num)
)
GD[GD < 0L] <- 0L
GD[GD > 2L] <- 2L

# Shorten SNP names (GAPIT can crash on long names)
# Include Chr_orig so results tables report original chromosome numbers
snp_map <- data.frame(
  ID_short = paste0("SNP", seq_len(nrow(geno_map))),
  ID_orig  = geno_map$SNP,
  Chr      = geno_map$Chr,
  Pos      = geno_map$Pos,
  stringsAsFactors = FALSE
)
GM <- geno_map
GM$SNP          <- snp_map$ID_short
colnames(GD)    <- snp_map$ID_short

# Additional MAF filter on GAPIT GD
maf_gapit <- apply(GD, 2, function(x) {
  freq <- sum(x, na.rm = TRUE) / (2 * length(x))
  pmin(freq, 1 - freq)
})
keep_gapit <- names(maf_gapit[maf_gapit >= MAF_THOLD])
GD         <- GD[, keep_gapit]
GM         <- GM[GM$SNP %in% keep_gapit, ]
snp_map    <- snp_map[snp_map$ID_short %in% keep_gapit, ]
cat("SNPs for GAPIT (MAF >=", MAF_THOLD, "):", ncol(GD), "\n")

# Remove monomorphic SNPs
valid_snps <- apply(GD, 2, function(x) length(unique(x[!is.na(x)])) > 1)
GD         <- GD[, valid_snps]
GM         <- GM[valid_snps, ]
snp_map    <- snp_map[snp_map$ID_short %in% colnames(GD), ]
cat("SNPs after monomorphic removal:", ncol(GD), "\n")

# GAPIT indexing bug workaround (remove last SNP/row)
GD <- GD[, -ncol(GD)]
GM <- GM[-nrow(GM), ]
cat("Final SNP count for GAPIT:", ncol(GD), "\n")

# ── GAPIT chromosome continuity fix ──────────────────────────────────────────
# GAPIT's internal Manhattan plotter uses seq(1, nchr, by = ncycle) which
# crashes when:
#   (a) GM$Chr contains NA values
#   (b) chromosomes are non-contiguous after QC (e.g. 1,3,5 with gaps)
#   (c) only one chromosome is present (nchr = 1 -> seq(1,1,...) issue)
#
# Fix: drop NA chrs, then remap chromosome integers to a contiguous 1:k scale.

# (a) Drop any remaining NA chromosomes in GM/GD
na_gm <- is.na(GM$Chr)
if (any(na_gm)) {
  cat("Dropping", sum(na_gm), "SNPs with NA Chr from GM.\n")
  GD <- GD[, !na_gm, drop = FALSE]
  GM <- GM[!na_gm, ]
}

# (b) Remap chromosomes to contiguous integers 1:k
chr_levels  <- sort(unique(GM$Chr))
chr_remap   <- setNames(seq_along(chr_levels), chr_levels)
# Store reverse lookup to restore original chr numbers after GAPIT finishes
chr_restore <- setNames(chr_levels, seq_along(chr_levels))
GM$Chr      <- chr_remap[as.character(GM$Chr)]
# NOTE: do NOT add Chr_orig as a column — extra columns in GM confuse GAPIT
cat(sprintf("Chromosomes remapped: %d original -> 1:%d contiguous\n",
            length(chr_levels), length(chr_levels)))
cat("Original chr values:", chr_levels, "\n\n")

# Keep GM to exactly the 3 columns GAPIT expects: SNP, Chr, Pos
GM <- GM[, c("SNP", "Chr", "Pos")]

# GD as data.frame with Taxa column (GAPIT alternate input path)
GD_df <- cbind(
  data.frame(Taxa = rownames(GD), stringsAsFactors = FALSE),
  as.data.frame(GD)
)
rownames(GD_df) <- NULL

# Kinship as data.frame with Taxa column
KI_gapit <- cbind(
  data.frame(Taxa = rownames(K_mat)),
  as.data.frame(K_mat)
)

# Verify alignment
stopifnot(
  nrow(pheno_gapit) == nrow(GD),
  all(pheno_gapit$Taxa == rownames(GD)),
  all(pheno_gapit$Taxa == CV$Taxa),
  ncol(GD) == nrow(GM),
  all(colnames(GD) == GM$SNP)
)
cat("All input alignment checks passed.\n")
write.csv(snp_map, "results/gwas_comparison/tables/snp_name_mapping.csv",
          row.names = FALSE)
cat("\n")


# ==============================================================================
# SECTION 7 — GAPIT PATCHES (required each session)
#
# Patch 1: GAPIT.Genotype — skips the broken internal MAF filter that causes
#           a dimension mismatch crash during genotype processing.
#
# Patch 2: GAPIT.Multiple.Manhattan — guards against the
#           "seq.default(1, nchr, by = ncycle): invalid '(to-from)/by'" crash.
#           This fires when nchr = 0, which happens when GAPIT's internal
#           plotter removes chr 0 / chr 99 and is left with no chromosomes.
#           ncycle = ceiling(0/5) = 0, making seq(1, 0, by=0) invalid.
#           Fix: wrap the offending seq() call with a max(1, ncycle) guard
#           and skip plotting entirely when nchr = 0.
# ==============================================================================

cat("Applying GAPIT patches...\n")

# ── GAPIT.Multiple.Manhattan patch: guard seq(1, nchr, by=ncycle) ───────
GAPIT.MM   <- GAPIT:::GAPIT.Multiple.Manhattan
src_mm     <- deparse(body(GAPIT.MM))
new_mm     <- src_mm

ncycle_ln  <- grep("ncycle=ceiling", src_mm)
if (length(ncycle_ln) > 0)
  new_mm[ncycle_ln[1]] <- sub("ncycle=ceiling\\(nchr/5\\)",
                              "ncycle=max(1L, ceiling(nchr/5))",
                              new_mm[ncycle_ln[1]])

body(GAPIT.MM) <- parse(text = paste(new_mm, collapse = "\n"))[[1]]
assignInNamespace("GAPIT.Multiple.Manhattan", GAPIT.MM, ns = "GAPIT")
cat("GAPIT.Multiple.Manhattan patch applied.\n")


# ==============================================================================
# SECTION 8 — SIGNIFICANCE THRESHOLDS
# ==============================================================================

n_snps <- nrow(GM)
gwt    <- -log10(0.05 / n_snps)   # Bonferroni genome-wide
sugt   <- -log10(1    / n_snps)   # Suggestive
cat(sprintf("SNPs in analysis     : %d\n",   n_snps))
cat(sprintf("Bonferroni threshold : %.4f (-log10p)\n", gwt))
cat(sprintf("Suggestive threshold : %.4f (-log10p)\n", sugt))
cat("\n")


# ==============================================================================
# SECTION 9 — RUN FarmCPU
# ==============================================================================

cat("----------------------------------------------------------------------\n")
cat("SECTION 9 | Running FarmCPU\n")
cat("----------------------------------------------------------------------\n")

setwd(file.path(PROJECT_DIR, "results/gwas_comparison/farmcpu"))

# file.output = FALSE disables GAPIT's internal plots entirely (avoids the
# Manhattan plotter crash). The script builds its own plots from the returned
# $GWAS data.frame instead.
farmcpu_out <- tryCatch(
  GAPIT(
    Y                = pheno_gapit,
    GD               = GD_df,
    GM               = GM,
    CV               = CV,
    KI               = KI_gapit,
    model            = "FarmCPU",
    PCA.total        = 0,
    SNP.MAF          = 0,
    file.output      = FALSE,
    Geno.View.output = FALSE,
    PCA.View.output  = FALSE
  ),
  error = function(e) {
    cat("WARNING: GAPIT threw an error after completing GWAS:\n  ", conditionMessage(e), "\n")
    cat("  This usually means the internal Manhattan plotter crashed.\n")
    cat("  Attempting to recover results from output files...\n")
    NULL
  }
)

setwd(PROJECT_DIR)

# Save results CSV — first try the returned object, then fall back to any
# CSV that GAPIT wrote before the plotter crashed.
farmcpu_res <- NULL
if (!is.null(farmcpu_out)) {
  farmcpu_res <- farmcpu_out$GWAS
}
if (is.null(farmcpu_res)) {
  # Attempt recovery from CSV written before the crash
  csv_candidates <- list.files("results/gwas_comparison/farmcpu",
                                pattern = "GWAS_Results.*\\.csv$",
                                full.names = TRUE, recursive = TRUE)
  if (length(csv_candidates) > 0) {
    farmcpu_res <- read.csv(csv_candidates[1], stringsAsFactors = FALSE)
    cat("FarmCPU: recovered", nrow(farmcpu_res), "SNPs from",
        basename(csv_candidates[1]), "\n\n")
  } else {
    stop("FarmCPU: no results found in returned object or output files.")
  }
} else {
  write.csv(farmcpu_res,
            "results/gwas_comparison/farmcpu/FarmCPU_Results.csv",
            row.names = FALSE)
  cat("FarmCPU results saved:", nrow(farmcpu_res), "SNPs\n\n")
}


# ==============================================================================
# SECTION 10 — RUN MLM
#
# MLM (Mixed Linear Model / EMMA / P3D):
#   - Uses the full kinship matrix K as a random effect
#   - Controls ALL background polygenic variation in one term
#   - PC1-PC3 included as fixed-effect covariates
#   - No pseudo-QTN selection (unlike FarmCPU)
#
# Expected behaviour at n=108:
#   At small n, the K matrix is highly correlated with individual tested
#   markers, causing over-correction (lambda_GC << 1, loss of power).
#   This is a known limitation of MLM in small panels.
# ==============================================================================

cat("----------------------------------------------------------------------\n")
cat("SECTION 10 | Running MLM\n")
cat("----------------------------------------------------------------------\n")

setwd(file.path(PROJECT_DIR, "results/gwas_comparison/mlm"))

mlm_out <- tryCatch(
  GAPIT(
    Y                = pheno_gapit,
    GD               = GD_df,
    GM               = GM,
    CV               = CV,
    KI               = KI_gapit,
    model            = "MLM",
    PCA.total        = 0,
    SNP.MAF          = 0,
    file.output      = FALSE,
    Geno.View.output = FALSE,
    PCA.View.output  = FALSE
  ),
  error = function(e) {
    cat("WARNING: GAPIT threw an error after completing GWAS:\n  ", conditionMessage(e), "\n")
    cat("  Attempting to recover results from output files...\n")
    NULL
  }
)

setwd(PROJECT_DIR)

# Save / recover MLM results
mlm_res <- NULL
if (!is.null(mlm_out)) {
  mlm_res <- mlm_out$GWAS
}
if (is.null(mlm_res)) {
  csv_candidates <- list.files("results/gwas_comparison/mlm",
                                pattern = "GWAS_Results.*\\.csv$",
                                full.names = TRUE, recursive = TRUE)
  if (length(csv_candidates) > 0) {
    mlm_res <- read.csv(csv_candidates[1], stringsAsFactors = FALSE)
    cat("MLM: recovered", nrow(mlm_res), "SNPs from",
        basename(csv_candidates[1]), "\n\n")
  } else {
    stop("MLM: no results found in returned object or output files.")
  }
} else {
  write.csv(mlm_res,
            "results/gwas_comparison/mlm/MLM_Results.csv",
            row.names = FALSE)
  cat("MLM results saved:", nrow(mlm_res), "SNPs\n\n")
}


# ==============================================================================
# SECTION 11 — LOAD AND PARSE RESULTS
# ==============================================================================

cat("----------------------------------------------------------------------\n")
cat("SECTION 11 | Loading and parsing results\n")
cat("----------------------------------------------------------------------\n")

load_gwas_results <- function(input, model_name, snp_map, chr_restore,
                              gwt, sugt) {
  # Accept either an in-memory data.frame or a file path string
  if (is.data.frame(input)) {
    df <- input
  } else {
    if (!file.exists(input))
      stop("Results CSV not found: ", input)
    df <- read.csv(input, stringsAsFactors = FALSE)
  }
  
  # ── Normalise column names (GAPIT versions differ in capitalisation) ─────────
  names(df) <- trimws(names(df))
  # P-value
  names(df)[tolower(names(df)) == "p.value"]  <- "P"
  names(df)[tolower(names(df)) == "pvalue"]   <- "P"
  # Effect / beta
  names(df)[tolower(names(df)) == "effect"]   <- "Beta"
  # Position — some versions write "pos", "Pos", "Position", "BP"
  pos_col <- which(tolower(names(df)) %in% c("pos","position","bp"))
  if (length(pos_col) > 0 && !"Pos" %in% names(df))
    names(df)[pos_col[1]] <- "Pos"
  # Chromosome
  chr_col <- which(tolower(names(df)) %in% c("chr","chrom","chromosome"))
  if (length(chr_col) > 0 && !"Chr" %in% names(df))
    names(df)[chr_col[1]] <- "Chr"
  
  # ── Basic filters ─────────────────────────────────────────────────────────────
  df <- df[!is.na(df$P) & df$P > 0, ]
  df$logP  <- -log10(df$P)
  df$model <- model_name
  
  # ── Restore original chromosome numbers ───────────────────────────────────────
  df$Chr <- chr_restore[as.character(as.integer(df$Chr))]
  df$Chr <- as.integer(df$Chr)
  
  # ── Recover original SNP IDs (join ONLY ID_orig — no Pos/Chr columns) ────────
  if ("SNP" %in% names(df)) {
    df <- left_join(df,
                    snp_map[, c("ID_short", "ID_orig")],  # only 2 cols — no conflict
                    by = c("SNP" = "ID_short"))
  }
  
  # Confirm Pos exists before continuing
  if (!"Pos" %in% names(df))
    stop(model_name, ": 'Pos' column not found after column normalisation.\n",
         "  Columns present: ", paste(names(df), collapse = ", "))
  
  df <- df[order(df$Chr, df$Pos), ]
  
  # ── Significance labels ───────────────────────────────────────────────────────
  df$sig <- ifelse(df$logP >= gwt,  "Genome-wide",
                   ifelse(df$logP >= sugt, "Suggestive", "Null"))
  
  # ── Cumulative chromosome positions for Manhattan plot ────────────────────────
  chr_meta <- df |>
    group_by(Chr) |>
    summarise(max_pos = max(Pos), .groups = "drop") |>
    mutate(cum_add = lag(cumsum(as.numeric(max_pos)), default = 0))
  
  df <- left_join(df, chr_meta, by = "Chr") |>
    mutate(pos_cum = Pos + cum_add)
  
  cat(sprintf("  %-10s : %d SNPs | Genome-wide: %d | Suggestive: %d\n",
              model_name, nrow(df),
              sum(df$sig == "Genome-wide"),
              sum(df$sig == "Suggestive")))
  df
}


# Safety check: chr_restore maps contiguous GAPIT chr integers back to the
# original chromosome numbers. It is built in Section 6. If the session is
# being resumed or Section 6 was run from another script version, rebuild it.
if (!exists("chr_restore")) {
  if (!exists("chr_levels"))
    chr_levels <- sort(unique(snp_map$Chr))
  chr_restore <- setNames(chr_levels, seq_along(chr_levels))
  cat("Note: chr_restore rebuilt from chr_levels.\n")
}

farmcpu_df <- load_gwas_results(
  if (!is.null(farmcpu_res)) farmcpu_res
  else "results/gwas_comparison/farmcpu/FarmCPU_Results.csv",
  "FarmCPU", snp_map, chr_restore, gwt, sugt)

mlm_df <- load_gwas_results(
  if (!is.null(mlm_res)) mlm_res
  else "results/gwas_comparison/mlm/MLM_Results.csv",
  "MLM", snp_map, chr_restore, gwt, sugt)

cat("\n")


# ==============================================================================
# SECTION 12 — LAMBDA GC
# ==============================================================================

compute_lambda <- function(df) {
  chi_obs <- qchisq(10^(-df$logP), df = 1, lower.tail = FALSE)
  median(chi_obs, na.rm = TRUE) / qchisq(0.5, df = 1)
}

lambda_farm <- compute_lambda(farmcpu_df)
lambda_mlm  <- compute_lambda(mlm_df)

interpret_lambda <- function(l) {
  if      (l > 1.10) "INFLATED (risk of false positives)"
  else if (l < 0.90) "CONSERVATIVE (over-corrected, potential loss of power)"
  else               "WELL-CONTROLLED"
}

cat(sprintf("Lambda GC — FarmCPU : %.4f  [%s]\n",
            lambda_farm, interpret_lambda(lambda_farm)))
cat(sprintf("Lambda GC — MLM     : %.4f  [%s]\n\n",
            lambda_mlm, interpret_lambda(lambda_mlm)))


# ==============================================================================
# SECTION 13 — HIT COMPARISON
# ==============================================================================

cat("----------------------------------------------------------------------\n")
cat("SECTION 13 | Hit comparison\n")
cat("----------------------------------------------------------------------\n")

keep_cols <- c("SNP","ID_orig","Chr","Pos","P","logP","sig","Beta","model")
top_farm  <- farmcpu_df[farmcpu_df$sig != "Null", ]
top_farm  <- top_farm[, intersect(keep_cols, names(top_farm)), drop = FALSE]
top_farm  <- top_farm[order(top_farm$logP, decreasing = TRUE), ]

top_mlm   <- mlm_df[mlm_df$sig != "Null", ]
top_mlm   <- top_mlm[, intersect(keep_cols, names(top_mlm)), drop = FALSE]
top_mlm   <- top_mlm[order(top_mlm$logP, decreasing = TRUE), ]

cat("FarmCPU significant hits:\n")
if (nrow(top_farm) > 0) {
  print(top_farm[, intersect(c("ID_orig","Chr","Pos","P","logP","sig"), names(top_farm))])
} else { cat("  None above suggestive threshold.\n") }

cat("\nMLM significant hits:\n")
if (nrow(top_mlm) > 0) {
  print(top_mlm[, intersect(c("ID_orig","Chr","Pos","P","logP","sig"), names(top_mlm))])
} else { cat("  None above suggestive threshold.\n") }

# Overlap
shared_snps  <- intersect(top_farm$SNP, top_mlm$SNP)
only_farm    <- setdiff(top_farm$SNP, top_mlm$SNP)
only_mlm     <- setdiff(top_mlm$SNP, top_farm$SNP)

cat(sprintf("\nSNP overlap (significant hits):\n"))
cat(sprintf("  Shared by both models : %d\n", length(shared_snps)))
cat(sprintf("  FarmCPU only          : %d\n", length(only_farm)))
cat(sprintf("  MLM only              : %d\n\n", length(only_mlm)))

# Effect size concordance
if ("Beta" %in% names(farmcpu_df) && "Beta" %in% names(mlm_df) &&
    length(shared_snps) > 0) {

  cat("Effect size concordance (shared significant SNPs):\n")
  beta_comp <- merge(
    farmcpu_df[farmcpu_df$SNP %in% shared_snps,
               c("SNP","ID_orig","Chr","Pos","Beta","logP")],
    mlm_df[mlm_df$SNP %in% shared_snps, c("SNP","Beta","logP")],
    by = "SNP", suffixes = c("_FarmCPU","_MLM")
  )
  beta_comp$same_direction <- sign(beta_comp$Beta_FarmCPU) ==
                              sign(beta_comp$Beta_MLM)
  print(beta_comp)
  cat("\n")
}

# -log10(p) correlation across all common SNPs
common_snps <- intersect(farmcpu_df$SNP, mlm_df$SNP)
r_lp <- NA_real_
if (length(common_snps) > 10) {
  lp_f <- farmcpu_df$logP[match(common_snps, farmcpu_df$SNP)]
  lp_m <- mlm_df$logP[match(common_snps, mlm_df$SNP)]
  r_lp <- cor(lp_f, lp_m, use = "complete.obs")
  cat(sprintf("SNPs in both outputs       : %d\n", length(common_snps)))
  cat(sprintf("Pearson r of -log10(p)     : %.4f\n", r_lp))
  cat("  r > 0.8 = models largely agree on signal ranking\n")
  cat("  r < 0.5 = models diverge substantially\n\n")
}


# ==============================================================================
# SECTION 14 — SUMMARY TABLE
# ==============================================================================

summary_tab <- data.frame(
  Metric = c(
    "SNPs analysed",
    "Genome-wide hits (Bonferroni)",
    "Suggestive hits",
    "Total hits",
    "Lambda GC",
    "Lambda interpretation",
    "Top hit chromosome",
    "Top hit -log10(p)",
    "Shared significant SNPs",
    "-log10(p) correlation (all SNPs)"
  ),
  FarmCPU = c(
    nrow(farmcpu_df),
    sum(farmcpu_df$sig == "Genome-wide"),
    sum(farmcpu_df$sig == "Suggestive"),
    sum(farmcpu_df$sig != "Null"),
    round(lambda_farm, 4),
    interpret_lambda(lambda_farm),
    if (nrow(top_farm) > 0) as.character(top_farm$Chr[1]) else "—",
    if (nrow(top_farm) > 0) round(top_farm$logP[1], 3)   else "—",
    length(shared_snps),
    if (!is.na(r_lp)) round(r_lp, 4) else "—"
  ),
  MLM = c(
    nrow(mlm_df),
    sum(mlm_df$sig == "Genome-wide"),
    sum(mlm_df$sig == "Suggestive"),
    sum(mlm_df$sig != "Null"),
    round(lambda_mlm, 4),
    interpret_lambda(lambda_mlm),
    if (nrow(top_mlm) > 0) as.character(top_mlm$Chr[1]) else "—",
    if (nrow(top_mlm) > 0) round(top_mlm$logP[1], 3)   else "—",
    length(shared_snps),
    if (!is.na(r_lp)) round(r_lp, 4) else "—"
  ),
  stringsAsFactors = FALSE
)

cat("--- Summary ---\n")
print(summary_tab, row.names = FALSE)
cat("\n")

write.csv(summary_tab,
          "results/gwas_comparison/tables/model_comparison_summary.csv",
          row.names = FALSE)
write.csv(rbind(top_farm, top_mlm),
          "results/gwas_comparison/tables/all_sig_hits_both_models.csv",
          row.names = FALSE)
cat("Tables saved to results/gwas_comparison/tables/\n\n")


# ==============================================================================
# SECTION 15 — PLOTS
# ==============================================================================

cat("----------------------------------------------------------------------\n")
cat("SECTION 15 | Plots\n")
cat("----------------------------------------------------------------------\n")

# ── Helper: Manhattan ─────────────────────────────────────────────────────────
# Alternating chromosome colours (teal / navy), suggestive = amber, gw = red
make_manhattan <- function(df, title_label, gwt, sugt) {
  
  chr_rects <- df |>
    group_by(Chr) |>
    summarise(xmin = min(pos_cum), xmax = max(pos_cum), .groups = "drop") |>
    filter(Chr %% 2 == 0)
  
  chr_ctrs <- df |>
    group_by(Chr) |>
    summarise(centre = mean(pos_cum), .groups = "drop")
  
  # Assign alternating colours to chromosomes for null SNPs
  chr_uniq  <- sort(unique(df$Chr))
  chr_cols  <- setNames(
    rep(c("#1D9E75", "#185FA5"), length.out = length(chr_uniq)),
    chr_uniq
  )
  df$chr_col <- chr_cols[as.character(df$Chr)]
  
  ggplot(df, aes(x = pos_cum, y = logP)) +
    geom_rect(data = chr_rects,
              aes(xmin = xmin, xmax = xmax, ymin = -Inf, ymax = Inf),
              inherit.aes = FALSE, fill = "#F4F4F4", alpha = 0.6) +
    geom_point(data = subset(df, sig == "Null"),
               aes(colour = chr_col), size = 0.7, alpha = 0.50, shape = 16) +
    geom_point(data = subset(df, sig == "Suggestive"),
               colour = "#D4860A", size = 2.4, alpha = 0.92, shape = 17) +
    geom_point(data = subset(df, sig == "Genome-wide"),
               colour = "#B85042", size = 3.4, alpha = 0.97, shape = 16) +
    geom_hline(yintercept = gwt,  colour = "#B85042",
               linewidth = 0.85, linetype = "dashed") +
    geom_hline(yintercept = sugt, colour = "#D4860A",
               linewidth = 0.65, linetype = "dotted") +
    annotate("text", x = Inf, y = gwt + 0.08, hjust = 1.05,
             label = paste0("Bonferroni (",
                            formatC(0.05 / nrow(df), format = "e", digits = 1),
                            ")"),
             colour = "#B85042", size = 2.8) +
    scale_colour_identity() +
    scale_x_continuous(breaks = chr_ctrs$centre, labels = chr_ctrs$Chr,
                       expand = expansion(mult = 0.01)) +
    labs(title = title_label,
         x     = "Chromosome",
         y     = expression(-log[10](italic(p)))) +
    theme_classic(base_size = 10) +
    theme(plot.title  = element_text(face = "bold", size = 11),
          axis.text.x = element_text(size = 7))
}

# ── Helper: QQ ────────────────────────────────────────────────────────────────
# CI ribbon = light teal, null = teal, suggestive = amber, gw = red
make_qq <- function(df, lambda_val) {
  
  n_s    <- nrow(df)
  obs_lp <- sort(df$logP, decreasing = TRUE)
  exp_lp <- sort(-log10(seq(1/n_s, 1, length.out = n_s)), decreasing = TRUE)
  ci_u   <- -log10(qbeta(0.025, seq_len(n_s), n_s - seq_len(n_s) + 1))
  ci_l   <- -log10(qbeta(0.975, seq_len(n_s), n_s - seq_len(n_s) + 1))
  
  qq_df <- data.frame(
    exp = sort(exp_lp, decreasing = TRUE),
    obs = obs_lp,
    cil = sort(ci_l,   decreasing = TRUE),
    ciu = sort(ci_u,   decreasing = TRUE),
    sig = ifelse(obs_lp >= gwt,  "Genome-wide",
                 ifelse(obs_lp >= sugt, "Suggestive", "Null"))
  )
  
  ggplot(qq_df, aes(x = exp, y = obs)) +
    geom_ribbon(aes(ymin = cil, ymax = ciu),
                fill = "#1D9E75", alpha = 0.15) +
    geom_abline(slope = 1, intercept = 0,
                colour = "#AAAAAA", linewidth = 0.9, linetype = "dashed") +
    geom_point(data = subset(qq_df, sig == "Null"),
               colour = "#185FA5", size = 0.9, alpha = 0.55) +
    geom_point(data = subset(qq_df, sig == "Suggestive"),
               colour = "#D4860A", size = 1.8, alpha = 0.92) +
    geom_point(data = subset(qq_df, sig == "Genome-wide"),
               colour = "#B85042", size = 2.5, alpha = 0.97) +
    annotate("text", x = 0.1, y = max(obs_lp) * 0.92, hjust = 0,
             label = paste0("Lambda_GC = ", round(lambda_val, 3)),
             colour = "#1C3A2A", size = 3.5, fontface = "bold") +
    labs(x = expression("Expected " * -log[10](italic(p))),
         y = expression("Observed "  * -log[10](italic(p)))) +
    theme_classic(base_size = 10)
}

# ── Build Manhattan + QQ panels ───────────────────────────────────────────────
p_man_farm <- make_manhattan(farmcpu_df,
                             paste0("FarmCPU | Lambda_GC = ", round(lambda_farm, 3),
                                    " | hits: ", sum(farmcpu_df$sig != "Null")),
                             gwt, sugt)

p_man_mlm  <- make_manhattan(mlm_df,
                             paste0("MLM | Lambda_GC = ", round(lambda_mlm, 3),
                                    " | hits: ", sum(mlm_df$sig != "Null")),
                             gwt, sugt)

p_qq_farm  <- make_qq(farmcpu_df, lambda_farm)
p_qq_mlm   <- make_qq(mlm_df,    lambda_mlm)

p_compare <- (p_man_farm + p_qq_farm) / (p_man_mlm + p_qq_mlm) +
  plot_layout(widths = c(3, 1)) +
  plot_annotation(
    title    = paste0("GWAS Comparison: FarmCPU vs. MLM \u2014 ",
                      TRAIT, " in Peanut"),
    subtitle = paste0("n = ", length(y), " accessions (BLUPs) | ",
                      n_snps, " SNPs | Covariates: PC1\u2013PC3 | ",
                      "Triangles = suggestive; circles = genome-wide"),
    caption  = "",
    theme    = theme(
      plot.title    = element_text(face = "bold", size = 13),
      plot.subtitle = element_text(size = 9,  colour = "#555555"),
      plot.caption  = element_text(size = 8,  colour = "#888888")))

pdf("results/gwas_comparison/plots/gwas_farmcpu_vs_mlm_manhattan_qq.pdf",
    width = 16, height = 9)
print(p_compare)
dev.off()
cat("Saved: results/gwas_comparison/plots/gwas_farmcpu_vs_mlm_manhattan_qq.pdf\n")

# ── -log10(p) scatter: FarmCPU vs MLM ────────────────────────────────────────
if (length(common_snps) > 10) {
  
  scatter_df <- data.frame(
    SNP       = common_snps,
    logP_farm = farmcpu_df$logP[match(common_snps, farmcpu_df$SNP)],
    logP_mlm  = mlm_df$logP[match(common_snps, mlm_df$SNP)],
    sig_farm  = farmcpu_df$sig[match(common_snps, farmcpu_df$SNP)],
    sig_mlm   = mlm_df$sig[match(common_snps,  mlm_df$SNP)]
  ) |>
    mutate(highlight = case_when(
      sig_farm != "Null" & sig_mlm != "Null" ~ "Both models",
      sig_farm != "Null"                     ~ "FarmCPU only",
      sig_mlm  != "Null"                     ~ "MLM only",
      TRUE                                   ~ "Null"
    ))
  
  p_scatter <- ggplot(scatter_df, aes(x = logP_farm, y = logP_mlm,
                                      colour = highlight)) +
    geom_point(data = subset(scatter_df, highlight == "Null"),
               size = 0.6, alpha = 0.30) +
    geom_point(data = subset(scatter_df, highlight != "Null"),
               size = 2.0, alpha = 0.90) +
    geom_abline(slope = 1, intercept = 0,
                colour = "#AAAAAA", linetype = "dashed", linewidth = 0.8) +
    scale_colour_manual(
      values = c("Both models"  = "#B85042",
                 "FarmCPU only" = "#D4860A",
                 "MLM only"     = "#1D9E75",
                 "Null"         = "#CCCCCC"),
      name = "Significance") +
    annotate("text", x = -Inf, y = Inf, hjust = -0.1, vjust = 1.4,
             label = paste0("r = ", round(r_lp, 3)),
             fontface = "bold", size = 4, colour = "#1C3A2A") +
    labs(title = expression("FarmCPU vs. MLM: " * -log[10](italic(p))),
         x = expression("FarmCPU  " * -log[10](italic(p))),
         y = expression("MLM  "     * -log[10](italic(p)))) +
    theme_classic(base_size = 11) +
    theme(plot.title  = element_text(face = "bold"),
          legend.position = "bottom",
          legend.text = element_text(size = 9))
  
  pdf("results/gwas_comparison/plots/gwas_logp_scatter.pdf",
      width = 6, height = 5.5)
  print(p_scatter)
  dev.off()
  cat("Saved: results/gwas_comparison/plots/gwas_logp_scatter.pdf\n")
}

# ==============================================================================
# SECTION 16 — INTERPRETATION & SESSION INFO
# ==============================================================================

cat("\n======================================================================\n")
cat(" INTERPRETATION\n")
cat("======================================================================\n")
cat(sprintf("FarmCPU lambda = %.4f  [%s]\n", lambda_farm, interpret_lambda(lambda_farm)))
cat(sprintf("MLM     lambda = %.4f  [%s]\n\n", lambda_mlm, interpret_lambda(lambda_mlm)))

cat("At n=108:\n")
cat("  FarmCPU uses pseudo-QTN fixed effects to control polygenic background.\n")
cat("  This avoids the K-matrix confounding that hurts MLM at small n.\n")
cat("  MLM with full K often over-corrects (lambda << 1) -> loss of power.\n\n")

cat("High-confidence loci: SNPs significant in BOTH models.\n")
cat("FarmCPU-only hits: may reflect power gain from pseudo-QTN approach.\n")
cat("MLM-only hits:     rare at small n; inspect carefully for artefacts.\n\n")

cat("======================================================================\n")
cat(" OUTPUT FILES\n")
cat("======================================================================\n")
cat("  results/gwas_comparison/tables/model_comparison_summary.csv\n")
cat("  results/gwas_comparison/tables/all_sig_hits_both_models.csv\n")
cat("  results/gwas_comparison/tables/snp_name_mapping.csv\n")
cat("  results/gwas_comparison/plots/gwas_farmcpu_vs_mlm_manhattan_qq.pdf\n")
cat("  results/gwas_comparison/plots/gwas_logp_scatter.pdf\n")
cat("  results/gwas_comparison/farmcpu/  (full GAPIT3 output)\n")
cat("  results/gwas_comparison/mlm/      (full GAPIT3 output)\n\n")

tryCatch({
sink("results/gwas_comparison/session_info.txt")
cat("FarmCPU vs. MLM Comparison | CSES 7160 | Fritzner Pierre\n")
cat("Date:", format(Sys.time(), "%Y-%m-%d %H:%M:%S"), "\n\n")
cat("Settings: MAF =", MAF_THOLD, "| Missingness =", MISS_THOLD,
    "| PCs =", N_PCS, "\n\n")
print(sessionInfo())
}, finally = {
  sink()
})
cat("Session info: results/gwas_comparison/session_info.txt\n")
cat("Done.\n")


# Near-threshold concordance check
# SNPs significant in FarmCPU but just below suggestive in MLM (logP > sugt - 0.3)
near_thresh <- farmcpu_df[farmcpu_df$sig != "Null", "SNP"]
near_thresh_mlm <- mlm_df[mlm_df$SNP %in% near_thresh, 
                          c("SNP", "ID_orig", "Chr", "Pos", "logP", "sig")]
near_thresh_mlm$logP_FarmCPU <- farmcpu_df$logP[match(near_thresh_mlm$SNP, farmcpu_df$SNP)]
near_thresh_mlm$above_80pct_threshold <- near_thresh_mlm$logP >= (sugt * 0.80)
cat("\nFarmCPU hits and their MLM -log10(p):\n")
print(near_thresh_mlm)

