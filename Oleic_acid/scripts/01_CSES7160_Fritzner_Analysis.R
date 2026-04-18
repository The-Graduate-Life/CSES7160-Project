# ==============================================================================
# Genomic Dissection and Prediction of Oleic Acid Concentration in Peanut
# Using High-Density SNP Markers (Axiom_Arachis_58K)
#
# Author  : Fritzner Pierre
# Course  : CSES 7160 — Genetic Data Analysis
# Date    : 2026
#
# Breeding Question:
#   What genomic regions control variation in oleic acid concentration in
#   peanut (Arachis hypogaea), and how accurately can genomic prediction
#   models predict oil quality?
#
# ── Attached data files ──────────────────────────────────────────────────────
#   values.csv
#       Source  : PeanutBase (https://arachispheno.peanutbase.org/phenotypes/)
#       Content : Oleic acid concentration (%) — U.S. peanut mini core
#                 315 observations across 108 accessions (1–4 reps, unbalanced)
#       Columns : phenotype_name, accession_id, accession_name,
#                 accession_cs_number, accession_longitude, accession_latitude,
#                 accession_country, phenotype_value, obs_unit_id
#
#   arahy_Tifrunner_gnm1_mrk_Axiom_Arachis_58K.gff3
#       Source  : Legume Information System (LIS)
#                 data.legumeinfo.org/Arachis/hypogaea/markers/
#                 Tifrunner.gnm1.mrk.Axiom_Arachis_58K/
#       Content : Chromosomal positions of 115 representative SNP markers
#                 from the Axiom_Arachis_58K array on the Tifrunner gnm1
#                 reference genome (20 chromosomes, Arahy.01–20)
#       Role    : SNP position map (GM table) for GAPIT.
#
#   aradu1_araip1.gnm1.div.Otyama_Kulkarni_2020.main.vcf.gz
#       Source  : LIS DataStore (Otyama, Kulkarni et al. 2020, G3)
#                 data.legumeinfo.org/Arachis/hypogaea/diversity/
#                 aradu1_araip1.gnm1.div.Otyama_Kulkarni_2020/
#       Content : SNP genotype calls for U.S. peanut mini core using
#                 Axiom_Arachis_58K array, mapped to diploid A. duranensis
#                 v1 + A. ipaensis v1 reference genomes. VCF (gzipped).
#       Note    : Chr names use Aradu/Araip prefixes; script converts
#                 these to integers (1-10 = A-genome, 11-20 = B-genome).
#
# ── GWAS model ───────────────────────────────────────────────────────────────
#   FarmCPU (Fixed and Random Model Circulating Probability Unification)
#   Reference : Liu X et al. (2016) PLOS Genetics 12(2): e1005767
#   Rationale : At n = 108, FarmCPU's pseudo-QTN fixed-effect approach avoids
#               the over-correction that a full kinship-matrix covariate (MLM)
#               can cause with small samples, while still controlling for
#               background polygenic effects more explicitly than a naive model.
#
# ── Pipeline ─────────────────────────────────────────────────────────────────
#   Step 1  — Packages & setup
#   Step 2  — Data import (values.csv + GFF3 + VCF genotype)
#   Step 3  — Phenotype pre-processing (outlier QC + BLUP extraction)
#   Step 4  — Genotype QC (MAF, missingness, imputation)
#   Step 5  — Exploratory phenotypic analysis (EDA)
#   Step 6  — Population structure (PCA)
#   Step 7  — Genomic relationship matrix (GRM, VanRaden 2008)
#   Step 8  — Genomic heritability (GBLUP / REML + bootstrap CI)
#   Step 9  — GWAS: FarmCPU (GAPIT3) + Manhattan & QQ plots
#   Step 10 — Genomic prediction: GBLUP, 5-fold CV (10 reps)
#   Step 11 — Results summary & session info
# ==============================================================================


# ==============================================================================
# WORKING DIRECTORY
# ==============================================================================
# Set this to your project folder before running the script.
# All input files (values.csv, GFF3, VCF.gz) must be inside this folder,
# and all results/ output will be written here too.
#
# Windows path note: use forward slashes (/) or double backslashes (\\).
# The path below is your OneDrive-synced Auburn project folder.

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
# STEP 1 — PACKAGES & SETUP
# ==============================================================================

# -- Install once, then recomment ----------------------------------------------
# install.packages(c(
#   "rrBLUP",      # A.mat(), mixed.solve()
#   "ggplot2",     # plots
#   "dplyr",       # data wrangling
#   "tidyr",       # reshaping
#   "patchwork",   # multi-panel figures
#   "lme4",        # lmer() for BLUP extraction
#   "data.table",  # fast file I/O
#   "reshape2",    # melt() for heatmap
#   "moments"      # skewness / kurtosis
# ))
# if (!requireNamespace("devtools", quietly = TRUE))
#   install.packages("devtools")
# devtools::install_github("jiabowang/GAPIT3", force = TRUE)

pkgs <- c(
  "rrBLUP","ggplot2","dplyr","tidyr","patchwork",
  "lme4","data.table","reshape2","moments"
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

# -- User settings -------------------------------------------------------------
PHENO_FILE  <- "data/values.csv"
GFF3_FILE   <- "data/arahy.Tifrunner.gnm1.mrk.Axiom_Arachis_58K.gff3"
GENO_FILE   <- "data/aradu1_araip1.gnm1.div.Otyama_Kulkarni_2020.main.vcf"
TRAIT       <- "data/OleicAcid"
MAF_THOLD   <- 0.05
MISS_THOLD  <- 0.10
N_FOLDS     <- 5
N_REPS      <- 10
N_PCS       <- 3

set.seed(2026)

for (d in c("results", "results/plots", "results/tables",
            "results/gwas", "results/gp"))
  dir.create(d, showWarnings = FALSE, recursive = TRUE)

cat("\n======================================================================\n")
cat("Peanut Oleic Acid Genomic Analysis\n")
cat("======================================================================\n")
cat("R version :", R.version$version.string, "\n")
cat("Date/Time :", format(Sys.time(), "%Y-%m-%d %H:%M:%S"), "\n\n")


# ==============================================================================
# STEP 2 — DATA IMPORT
# ==============================================================================

cat("----------------------------------------------------------------------\n")
cat("STEP 2 | Data Import\n")
cat("----------------------------------------------------------------------\n")

# -- 2.1 Phenotype (values.csv) ------------------------------------------------
stopifnot(file.exists(PHENO_FILE))
pheno_raw <- read.csv(PHENO_FILE, stringsAsFactors = FALSE)
pheno_raw <- rename(pheno_raw, Taxa = accession_name, Value = phenotype_value)
cat("Phenotype rows       :", nrow(pheno_raw), "\n")
cat("Unique accessions    :", length(unique(pheno_raw$Taxa)), "\n")
cat("Reps per accession   :",
    min(table(pheno_raw$Taxa)), "to",
    max(table(pheno_raw$Taxa)), "\n")

# -- 2.2 GFF3 marker position map ----------------------------------------------
# 115 markers from Axiom_Arachis_58K mapped to Tifrunner gnm1.
# Chromosome names: arahy.Tifrunner.gnm1.Arahy.XX -> integer XX
stopifnot(file.exists(GFF3_FILE))
gff3_raw  <- readLines(GFF3_FILE)
gff3_raw  <- gff3_raw[!grepl("^#", gff3_raw) & nchar(gff3_raw) > 0]

gff3_df <- do.call(rbind, lapply(gff3_raw, function(ln) {
  p   <- strsplit(ln, "\t")[[1]]
  chr <- as.integer(sub(".*Arahy\\.(\\d+)$", "\\1", p[1]))
  pos <- as.integer(p[4])
  snp <- sub(".*Name=([^;[:space:]]+).*", "\\1", p[9])
  data.frame(SNP = trimws(snp), Chr = chr, Pos = pos,
             stringsAsFactors = FALSE)
}))
cat("\nGFF3 markers loaded  :", nrow(gff3_df),
    "on", length(unique(gff3_df$Chr)), "chromosomes\n")

# -- 2.3 Genotype VCF (Otyama & Kulkarni 2020, LIS DataStore) -----------------
# File: aradu1_araip1.gnm1.div.Otyama_Kulkarni_2020.main.vcf.gz
# Download from:
#   https://data.legumeinfo.org/Arachis/hypogaea/diversity/
#   aradu1_araip1.gnm1.div.Otyama_Kulkarni_2020/
#
# This VCF contains SNP calls from the Axiom_Arachis_58K array for the
# U.S. peanut mini core collection, mapped to the diploid A. duranensis v1
# and A. ipaensis v1 reference genomes (Aradu/Araip chromosome prefixes).
#
# Requires the vcfR package:
#   install.packages("vcfR")
library(vcfR)

stopifnot(file.exists(GENO_FILE))
cat("Reading VCF file (this may take a minute)...\n")
vcf_raw  <- read.vcfR(GENO_FILE, verbose = FALSE)
cat("VCF loaded:", nrow(vcf_raw@fix), "variants x",
    ncol(vcf_raw@gt) - 1, "samples\n")

# -- Extract numeric genotype matrix {0, 1, 2} --------------------------------
# vcfR::extract.gt() returns character genotypes; convert to allele dosage.
gt_mat <- extract.gt(vcf_raw, element = "GT",
                     as.numeric = FALSE, return.alleles = FALSE)
# Recode: "0/0"->0, "0/1" or "1/0"->1, "1/1"->2, "./." or NA -> NA
recode_gt <- function(x) {
  x[x == "./."] <- NA
  x[x == "0/0"] <- "0"
  x[x == "0/1" | x == "1/0"] <- "1"
  x[x == "1/1"] <- "2"
  as.numeric(x)
}
gt_num <- apply(gt_mat, 2, recode_gt)   # rows = SNPs, cols = samples
geno_num <- t(gt_num)                   # transpose: rows = taxa, cols = SNPs
rownames(geno_num) <- colnames(gt_mat)
cat("Genotype matrix:", nrow(geno_num), "accessions x", ncol(geno_num), "SNPs\n")

# -- Build SNP position map (GM) from VCF CHROM/POS fields -------------------
# Chromosome names in this VCF the diploid reference is used:
#   Aradu.A01 ... Aradu.A10 (A-subgenome, A. duranensis)
#   Araip.B01 ... Araip.B10 (B-subgenome, A. ipaensis)
# GAPIT requires integer chromosome codes. We map:
#   Aradu.A01-A10 -> 1-10,  Araip.B01-B10 -> 11-20
chr_raw  <- vcf_raw@fix[, "CHROM"]
pos_raw  <- as.integer(vcf_raw@fix[, "POS"])
snp_ids  <- paste0(chr_raw, "_", pos_raw)   # unique SNP ID: CHROM_POS

chr_to_int <- function(ch) {
  ch <- sub("Aradu\\.A0?", "", ch)
  ch <- sub("Araip\\.B0?", "", ch)
  num <- suppressWarnings(as.integer(ch))
  # A01-A10 -> 1-10; B01-B10 -> 11-20
  ifelse(grepl("Araip", chr_raw),
         num + 10L,
         num)
}
chr_int <- chr_to_int(chr_raw)

geno_map <- data.frame(
  SNP  = snp_ids,
  Chr  = chr_int,
  Pos  = pos_raw,
  stringsAsFactors = FALSE
)
colnames(geno_num) <- snp_ids   # align SNP names
rm(vcf_raw, gt_mat, gt_num); gc()
cat("SNP map built:", nrow(geno_map), "SNPs on",
    length(unique(geno_map$Chr)), "chromosomes\n\n")

# -- Normalise VCF sample names to match phenotype format "PI XXXXXX" -------
vcf_taxa <- rownames(geno_num)
vcf_taxa <- sub("_[0-9]+$", "", vcf_taxa)        # strip _2 / _3
vcf_taxa <- sub("_s$",       "", vcf_taxa)        # strip _s
vcf_taxa <- sub("^PI([0-9]+)$", "PI \\1", vcf_taxa)  # "PI200441" -> "PI 200441"
rownames(geno_num) <- vcf_taxa

n_match <- length(intersect(vcf_taxa, pheno_raw$Taxa))
cat("Taxa matched after name normalisation:", n_match,
    "/ 108 expected\n")
print(head(vcf_taxa, 10))


# Document the accessions dropped due to absence from VCF
dropped <- setdiff(pheno_raw$Taxa, vcf_taxa)
if (length(dropped) > 0) {
  cat("Accessions in phenotype but absent from VCF:", length(dropped), "\n")
  cat(" ", paste(dropped, collapse = ", "), "\n")
  cat("  These will be excluded from GWAS and GP (no genotype data).\n\n")
  write.csv(
    data.frame(Taxa = dropped,
               Reason = "Absent from Otyama_Kulkarni_2020 VCF"),
    "results/tables/accessions_dropped_no_genotype.csv",
    row.names = FALSE)
}

if (n_match == 0)
  stop("Still no taxa match after normalisation. ",
       "Check rownames(geno_num) vs blup_df$Taxa manually.")

# ==============================================================================
# STEP 3 — PHENOTYPE PRE-PROCESSING
# ==============================================================================
# Rationale: 108 accessions have 1-4 unbalanced replicates. Using the raw mean
# as phenotype ignores measurement error heterogeneity across accessions.
# Solution:
#   (a) Detect and remove extreme within-accession outlier replicates
#   (b) Fit lmer: Value ~ 1 + (1|Taxa) and extract BLUPs of Taxa
#   (c) Use BLUPs as the adjusted genotype-level trait for GWAS and GP
# ==============================================================================

cat("----------------------------------------------------------------------\n")
cat("STEP 3 | Phenotype Pre-processing\n")
cat("----------------------------------------------------------------------\n")

# -- 3.1 Flag within-accession outlier replicates (|z_within| > 3) -----------
pheno_flagged <- pheno_raw |>
  group_by(Taxa) |>
  mutate(z_within = if (n() > 1) as.numeric(scale(Value)) else 0) |>
  ungroup()

outlier_reps <- filter(pheno_flagged, abs(z_within) > 3)
cat("Replicate-level outliers (|z_within| > 3):", nrow(outlier_reps), "\n")
if (nrow(outlier_reps) > 0) {
  cat("  Flagged:\n")
  print(select(outlier_reps, Taxa, Value, obs_unit_id, z_within))
}

pheno_clean <- filter(pheno_flagged, abs(z_within) <= 3)
cat("Observations retained:", nrow(pheno_clean), "/", nrow(pheno_raw), "\n")

# Report accessions with high within-rep SD (> 10%)
high_var <- pheno_clean |>
  group_by(Taxa) |>
  filter(n() > 1) |>
  summarise(SD = sd(Value), .groups = "drop") |>
  filter(SD > 10) |>
  arrange(desc(SD))
if (nrow(high_var) > 0) {
  cat("Accessions with within-rep SD > 10 (post-outlier-removal):\n")
  print(high_var)
}

# -- 3.2 BLUP extraction via lmer ----------------------------------------------
cat("\nExtracting BLUPs with lmer: Value ~ 1 + (1|Taxa)...\n")
lmer_fit  <- lmer(Value ~ 1 + (1 | Taxa), data = pheno_clean, REML = TRUE)
blup_list <- ranef(lmer_fit)$Taxa
mu        <- as.numeric(fixef(lmer_fit))

blup_df <- data.frame(
  Taxa  = rownames(blup_list),
  BLUP  = blup_list[, "(Intercept)"],
  GEBV  = blup_list[, "(Intercept)"] + mu,
  stringsAsFactors = FALSE
)
cat("Grand mean (mu)      :", round(mu, 3), "%\n")
cat("BLUPs extracted for  :", nrow(blup_df), "accessions\n")
cat("BLUP range           :", round(min(blup_df$BLUP), 3), "to",
    round(max(blup_df$BLUP), 3), "\n")
write.csv(blup_df, "results/tables/BLUPs_OleicAcid.csv", row.names = FALSE)

# -- 3.3 Broad-sense heritability from lmer ------------------------------------
vc_lmer   <- as.data.frame(VarCorr(lmer_fit))
sigma2_G  <- vc_lmer[vc_lmer$grp == "Taxa",     "vcov"]
sigma2_eL <- vc_lmer[vc_lmer$grp == "Residual", "vcov"]
H2_broad  <- sigma2_G / (sigma2_G + sigma2_eL)

cat("\n--- Broad-Sense Heritability (lmer variance components) ---\n")
cat("sigma2_G  (genetic variance) :", round(sigma2_G,  4), "\n")
cat("sigma2_e  (residual variance):", round(sigma2_eL, 4), "\n")
cat("H2 (broad-sense)             :", round(H2_broad,  4), "\n")

# Save alongside BLUPs
H2_table <- data.frame(
  Parameter = c("sigma2_G", "sigma2_e", "H2_broad"),
  Estimate  = round(c(sigma2_G, sigma2_eL, H2_broad), 6),
  Description = c("Among-accession genetic variance (lmer)",
                  "Within-accession residual variance (lmer)",
                  "Broad-sense heritability")
)
write.csv(H2_table, "results/tables/H2_broad_lmer.csv", row.names = FALSE)
cat("Table: results/tables/H2_broad_lmer.csv\n")

# -- 3.3 Align with genotype taxa ---------------------------------------------
common_taxa  <- intersect(blup_df$Taxa, rownames(geno_num))
if (length(common_taxa) == 0)
  stop("No taxa match between phenotype and genotype. Check accession ID format.")

blup_aligned <- blup_df[match(common_taxa, blup_df$Taxa), ]
geno_num     <- geno_num[common_taxa, , drop = FALSE]
y            <- blup_aligned$BLUP
cat("Taxa in common (pheno & geno):", length(common_taxa), "\n\n")


# ==============================================================================
# STEP 4 — GENOTYPE QUALITY CONTROL
# ==============================================================================

cat("----------------------------------------------------------------------\n")
cat("STEP 4 | Genotype QC\n")
cat("----------------------------------------------------------------------\n")
cat("Input:", nrow(geno_num), "x", ncol(geno_num), "\n")

# Missingness filter
miss_rate <- colMeans(is.na(geno_num))
keep_miss <- miss_rate <= MISS_THOLD
cat("Removed (missingness >", MISS_THOLD * 100, "%):", sum(!keep_miss), "\n")
geno_num  <- geno_num[, keep_miss, drop = FALSE]
geno_map  <- geno_map[keep_miss, ]

# MAF filter
maf_vals  <- apply(geno_num, 2, function(x) {
  p <- mean(x, na.rm = TRUE) / 2; min(p, 1 - p)
})
keep_maf  <- maf_vals >= MAF_THOLD
cat("Removed (MAF <", MAF_THOLD, ")         :", sum(!keep_maf), "\n")
geno_num  <- geno_num[, keep_maf, drop = FALSE]
geno_map  <- geno_map[keep_maf, ]
maf_vals  <- maf_vals[keep_maf]

# Mean imputation
n_miss <- sum(is.na(geno_num))
if (n_miss > 0) {
  col_mu <- colMeans(geno_num, na.rm = TRUE)
  for (j in seq_len(ncol(geno_num))) {
    idx <- is.na(geno_num[, j])
    if (any(idx)) geno_num[idx, j] <- col_mu[j]
  }
}
cat("Imputed missing      :", n_miss, "\n")
cat("Final matrix         :", nrow(geno_num), "x", ncol(geno_num), "\n")

geno_rrblup <- geno_num - 1   # {-1, 0, 1} for rrBLUP

# MAF plot
pdf("results/plots/00_maf_distribution.pdf", width = 7, height = 4)
print(ggplot(data.frame(MAF = maf_vals), aes(x = MAF)) +
  geom_histogram(binwidth = 0.02, fill = "#2E6B4F",
                 colour = "#1B3A2D", alpha = 0.85) +
  geom_vline(xintercept = MAF_THOLD, colour = "#B85042",
             linetype = "dashed", linewidth = 0.9) +
  labs(title    = "MAF Distribution (post-QC SNPs)",
       subtitle = paste0(ncol(geno_num), " SNPs | MAF >= ", MAF_THOLD),
       x = "Minor Allele Frequency", y = "Count") +
  theme_classic(base_size = 11) +
  theme(plot.title = element_text(face = "bold", colour = "#1B3A2D")))
dev.off()
cat("Plot: results/plots/00_maf_distribution.pdf\n\n")


# ==============================================================================
# STEP 5 — EXPLORATORY PHENOTYPIC ANALYSIS
# ==============================================================================

cat("----------------------------------------------------------------------\n")
cat("STEP 5 | Exploratory Phenotypic Analysis\n")
cat("----------------------------------------------------------------------\n")

raw_vals <- pheno_clean$Value
sw_raw   <- shapiro.test(raw_vals)
sw_blup  <- shapiro.test(y)

summ_tbl <- data.frame(
  Statistic = c("n accessions", "n replicates",
                "Mean raw (%)", "SD raw (%)", "CV raw (%)",
                "Min raw (%)", "Max raw (%)",
                "Skewness raw", "Kurtosis raw", "Shapiro-Wilk p (raw)",
                "Mean BLUP", "SD BLUP", "Shapiro-Wilk p (BLUP)"),
  Value = round(c(
    length(common_taxa), nrow(pheno_clean),
    mean(raw_vals), sd(raw_vals), 100 * sd(raw_vals) / mean(raw_vals),
    min(raw_vals), max(raw_vals),
    skewness(raw_vals), kurtosis(raw_vals), sw_raw$p.value,
    mean(y), sd(y), sw_blup$p.value), 4)
)
print(summ_tbl, row.names = FALSE)
write.csv(summ_tbl, "results/tables/phenotype_summary.csv", row.names = FALSE)

# Histogram of raw values
p_raw <- ggplot(data.frame(v = raw_vals), aes(x = v)) +
  geom_histogram(aes(y = after_stat(density)), binwidth = 3,
                 fill = "#2E6B4F", colour = "#1B3A2D", alpha = 0.80) +
  geom_density(colour = "#1B3A2D", linewidth = 1.1, adjust = 0.9) +
  geom_rug(colour = "#5BAD7F", alpha = 0.45, linewidth = 0.5) +
  geom_vline(xintercept = mean(raw_vals), colour = "#B85042",
             linetype = "dashed", linewidth = 0.9) +
  annotate("text", x = mean(raw_vals) + 1.5, y = Inf, vjust = 1.8,
           label = paste0("Mean = ", round(mean(raw_vals), 1), "%"),
           colour = "#B85042", size = 3.5) +
  labs(title = "Raw Replicates",
       subtitle = paste0("n = ", nrow(pheno_clean), " observations"),
       x = "Oleic Acid (%)", y = "Density") +
  theme_classic(base_size = 11) +
  theme(plot.title    = element_text(face = "bold", colour = "#1B3A2D"),
        plot.subtitle = element_text(colour = "#5C7A65"))

# Histogram of BLUPs
p_blup <- ggplot(data.frame(v = y), aes(x = v)) +
  geom_histogram(aes(y = after_stat(density)), binwidth = 2,
                 fill = "#5BAD7F", colour = "#1B3A2D", alpha = 0.80) +
  geom_density(colour = "#1B3A2D", linewidth = 1.1, adjust = 0.9) +
  geom_rug(colour = "#2E6B4F", alpha = 0.55, linewidth = 0.5) +
  geom_vline(xintercept = mean(y), colour = "#B85042",
             linetype = "dashed", linewidth = 0.9) +
  labs(title = "BLUPs  [used for GWAS & GP]",
       subtitle = paste0("n = ", length(y), " accessions"),
       x = "BLUP of Oleic Acid", y = "Density") +
  theme_classic(base_size = 11) +
  theme(plot.title    = element_text(face = "bold", colour = "#1B3A2D"),
        plot.subtitle = element_text(colour = "#5C7A65"))

# Boxplot comparison raw vs BLUPs
comp_df <- bind_rows(
  data.frame(Group = "Raw replicates", Value = raw_vals),
  data.frame(Group = "BLUPs",          Value = y)
)
p_cmp <- ggplot(comp_df, aes(x = Group, y = Value, fill = Group)) +
  geom_boxplot(colour = "#1B3A2D", outlier.colour = "#B85042",
               outlier.size = 2, alpha = 0.80, width = 0.5) +
  geom_jitter(width = 0.1, alpha = 0.30, size = 1.2, colour = "#1B3A2D") +
  scale_fill_manual(values = c("Raw replicates" = "#2E6B4F",
                                "BLUPs"          = "#5BAD7F"),
                    guide = "none") +
  labs(title = "Raw vs BLUPs", y = "Oleic Acid (%) / BLUP", x = NULL) +
  theme_classic(base_size = 11) +
  theme(plot.title = element_text(face = "bold", colour = "#1B3A2D"))

# Q-Q plot for BLUPs
p_qq0 <- ggplot(data.frame(v = y), aes(sample = v)) +
  stat_qq(colour = "#2E6B4F", alpha = 0.75, size = 2.2) +
  stat_qq_line(colour = "#B85042", linewidth = 1.0, linetype = "dashed") +
  labs(title = "Q-Q Plot (BLUPs)",
       x = "Theoretical Quantiles", y = "Sample Quantiles") +
  theme_classic(base_size = 11) +
  theme(plot.title = element_text(face = "bold", colour = "#1B3A2D"))

p_eda <- (p_raw | p_blup) / (p_cmp | p_qq0) +
  plot_annotation(
    title   = "Exploratory Phenotypic Analysis — Oleic Acid Concentration",
    caption = paste0("CSES 7160 | Fritzner Pierre | n=108 accessions, ",
                     "315 observations | U.S. peanut mini core"),
    theme   = theme(
      plot.title   = element_text(face = "bold", size = 13, colour = "#1B3A2D"),
      plot.caption = element_text(colour = "#5C7A65", size = 9))
  )
pdf("results/plots/01_phenotype_eda.pdf", width = 12, height = 9)
print(p_eda)
dev.off()
cat("Plot: results/plots/01_phenotype_eda.pdf\n\n")


# ==============================================================================
# STEP 6 — POPULATION STRUCTURE (PCA)
# ==============================================================================

cat("----------------------------------------------------------------------\n")
cat("STEP 6 | Population Structure (PCA)\n")
cat("----------------------------------------------------------------------\n")

geno_sc  <- scale(geno_rrblup, center = TRUE, scale = TRUE)
pca_res  <- prcomp(geno_sc, retx = TRUE, center = FALSE, scale. = FALSE)
pve      <- (pca_res$sdev^2) / sum(pca_res$sdev^2) * 100
cat("PC1:", round(pve[1], 2), "% | PC2:", round(pve[2], 2),
    "% | PC3:", round(pve[3], 2), "%\n")

scores      <- as.data.frame(pca_res$x[, 1:10])
scores$Taxa <- rownames(geno_rrblup)
write.csv(scores, "results/tables/pca_scores.csv", row.names = FALSE)

# Colour points by BLUP quartile to overlay phenotypic signal
q_df <- data.frame(
  Taxa     = blup_aligned$Taxa,
  Quartile = cut(blup_aligned$BLUP,
                 breaks = quantile(blup_aligned$BLUP, c(0,.25,.5,.75,1)),
                 labels = c("Q1 Low","Q2","Q3","Q4 High"),
                 include.lowest = TRUE)
)
sc_q <- left_join(scores, q_df, by = "Taxa")
qpal <- c("Q1 Low" = "#1B3A2D", "Q2" = "#2E6B4F",
          "Q3"     = "#C7A84F", "Q4 High" = "#B85042")

p_pc12 <- ggplot(sc_q, aes(x = PC1, y = PC2, colour = Quartile)) +
  geom_hline(yintercept = 0, colour = "#CCDDCC", linewidth = 0.4) +
  geom_vline(xintercept = 0, colour = "#CCDDCC", linewidth = 0.4) +
  geom_point(alpha = 0.80, size = 2.8) +
  stat_ellipse(aes(group = Quartile), level = 0.90,
               linewidth = 0.7, linetype = "dashed") +
  scale_colour_manual(values = qpal, name = "BLUP Quartile") +
  labs(title = "PC1 vs PC2  (coloured by BLUP quartile)",
       x = paste0("PC1 (", round(pve[1],1), "%)"),
       y = paste0("PC2 (", round(pve[2],1), "%)")) +
  theme_classic(base_size = 11) +
  theme(plot.title = element_text(face = "bold", colour = "#1B3A2D"))

scree_df <- data.frame(PC = 1:20, PVE = pve[1:20],
                        CumPVE = cumsum(pve[1:20]))
p_scree <- ggplot(scree_df, aes(x = PC)) +
  geom_col(aes(y = PVE), fill = "#2E6B4F", alpha = 0.85, width = 0.7) +
  geom_line(aes(y = CumPVE / 4), colour = "#B85042",
            linewidth = 0.9, linetype = "dashed") +
  geom_point(aes(y = CumPVE / 4), colour = "#B85042", size = 2.0) +
  scale_y_continuous(name = "Variance Explained (%)",
                     sec.axis = sec_axis(~.*4, name = "Cumulative (%)")) +
  scale_x_continuous(breaks = 1:20) +
  labs(title = "Scree Plot") +
  theme_classic(base_size = 10) +
  theme(plot.title = element_text(face = "bold", colour = "#1B3A2D"))

p_pca <- p_pc12 / p_scree +
  plot_annotation(
    title   = "Population Structure — PCA of Axiom_Arachis_58K SNPs",
    caption = "CSES 7160 | Fritzner Pierre | n = 108 accessions",
    theme   = theme(
      plot.title   = element_text(face = "bold", size = 13, colour = "#1B3A2D"),
      plot.caption = element_text(colour = "#5C7A65", size = 9))
  )
pdf("results/plots/02_pca.pdf", width = 10, height = 10)
print(p_pca)
dev.off()
cat("Plot : results/plots/02_pca.pdf\n")
cat("Table: results/tables/pca_scores.csv\n\n")

cv_farmcpu <- scores[, c("Taxa", "PC1", "PC2", "PC3")]


# ==============================================================================
# STEP 7 — GENOMIC RELATIONSHIP MATRIX (VanRaden 2008)
# ==============================================================================

cat("----------------------------------------------------------------------\n")
cat("STEP 7 | Genomic Relationship Matrix\n")
cat("----------------------------------------------------------------------\n")

K_mat     <- A.mat(geno_rrblup, min.MAF = NULL,
                   impute.method = "mean", return.imputed = FALSE)
diag_v    <- diag(K_mat)
off_v     <- K_mat[lower.tri(K_mat)]
cat("GRM dim    :", nrow(K_mat), "x", ncol(K_mat), "\n")
cat("Diag mean  :", round(mean(diag_v), 3), " range [",
    round(min(diag_v), 3), ",", round(max(diag_v), 3), "]\n")
cat("Off-diag   :", round(mean(off_v), 4), " range [",
    round(min(off_v),  3), ",", round(max(off_v),  3), "]\n")

pdf("results/plots/03_grm_heatmap.pdf", width = 7, height = 6)
print(ggplot(melt(K_mat), aes(x = Var1, y = Var2, fill = value)) +
  geom_raster() +
  scale_fill_gradientn(
    colours = c("#1B3A2D","#2E6B4F","#5BAD7F","#EAF4EE","#FFFFFF"),
    name = "Relatedness") +
  labs(title    = "Genomic Relationship Matrix (GRM)",
       subtitle = "VanRaden (2008) | Axiom_Arachis_58K",
       x = "Accession", y = "Accession") +
  theme_classic(base_size = 10) +
  theme(axis.text  = element_blank(), axis.ticks = element_blank(),
        plot.title = element_text(face = "bold", colour = "#1B3A2D"),
        plot.subtitle = element_text(colour = "#5C7A65"),
        legend.key.width = unit(0.8,"cm")))
dev.off()
cat("Plot: results/plots/03_grm_heatmap.pdf\n\n")


# ==============================================================================
# STEP 8 — GENOMIC HERITABILITY
# ==============================================================================
# GBLUP model:  y_BLUP = mu + g + e
#   g ~ N(0, K*sigma2_g),  e ~ N(0, I*sigma2_e)
# h2g = sigma2_g / (sigma2_g + sigma2_e)
# 95% CI via 500 bootstrap resamples
# ==============================================================================

cat("----------------------------------------------------------------------\n")
cat("STEP 8 | Genomic Heritability\n")
cat("----------------------------------------------------------------------\n")

lmm   <- mixed.solve(y = y, K = K_mat, SE = TRUE, return.Hinv = FALSE)
Vu    <- lmm$Vu;  Ve <- lmm$Ve;  Vp <- Vu + Ve
h2g   <- Vu / Vp
cat("sigma2_g :", round(Vu, 6), "\n")
cat("sigma2_e :", round(Ve, 6), "\n")
cat("h2g      :", round(h2g, 4), "\n")

cat("Bootstrap CI (500 resamples)...\n")
boot_h <- replicate(500, {
  idx <- sample(length(y), replace = TRUE)
  fb  <- tryCatch(mixed.solve(y = y[idx], K = K_mat[idx, idx]),
                  error = function(e) NULL)
  if (!is.null(fb)) fb$Vu / (fb$Vu + fb$Ve) else NA_real_
})
ci <- quantile(boot_h, c(0.025, 0.975), na.rm = TRUE)
cat("Bootstrap 95% CI: [", round(ci[1], 4), ",", round(ci[2], 4), "]\n\n")

write.csv(
  data.frame(
    Parameter   = c("sigma2_g","sigma2_e","sigma2_P",
                    "h2g","h2g_CI_lower","h2g_CI_upper"),
    Estimate    = round(c(Vu, Ve, Vp, h2g, ci[1], ci[2]), 6),
    Description = c("Genomic variance","Residual variance",
                    "Phenotypic variance of BLUPs","Genomic heritability",
                    "Bootstrap 95% CI lower","Bootstrap 95% CI upper")),
  "results/tables/heritability.csv", row.names = FALSE)
cat("Table: results/tables/heritability.csv\n\n")


# ==============================================================================
# STEP 9 — GWAS: FarmCPU (GAPIT3) — SAFE & EFFICIENT
# ==============================================================================
cat("----------------------------------------------------------------------\n")
cat("STEP 9 | GWAS — FarmCPU (GAPIT3) SAFE & EFFICIENT\n")
cat("----------------------------------------------------------------------\n")

# -- Prepare phenotype --------------------------------------------------------
pheno_gapit <- data.frame(
  Taxa      = blup_aligned$Taxa,
  OleicAcid = blup_aligned$BLUP,
  stringsAsFactors = FALSE
)

# -- Prepare genotype ---------------------------------------------------------
GD <- matrix(
  as.integer(round(geno_num)),  # enforce 0/1/2 integers
  nrow = nrow(geno_num),
  ncol = ncol(geno_num),
  dimnames = dimnames(geno_num)
)
GD[GD < 0L] <- 0L
GD[GD > 2L] <- 2L

cat("GD unique values  :", sort(unique(as.vector(GD))), "\n")
cat("GD dimensions     :", nrow(GD), "x", ncol(GD), "\n")

# -- Clean SNP map ------------------------------------------------------------
chr_clean     <- sub("_[0-9]+$", "", geno_map$SNP)
chr_last      <- sub(".*\\.", "", chr_clean)
is_b          <- grepl("^B", chr_last)
chr_num       <- suppressWarnings(as.integer(sub("^[AB]0?", "", chr_last)))
geno_map$Chr  <- ifelse(is_b, chr_num + 10L, chr_num)

# Shorten SNP names to avoid GAPIT errors
snp_map <- data.frame(
  ID_short = paste0("SNP", seq_len(nrow(geno_map))),
  ID_orig  = geno_map$SNP,
  Chr      = geno_map$Chr,
  Pos      = geno_map$Pos,
  stringsAsFactors = FALSE
)
GM <- geno_map
GM$SNP <- snp_map$ID_short
colnames(GD) <- snp_map$ID_short

# -- Optional: Manual MAF filter ---------------------------------------------
# Compute minor allele frequency for each SNP
maf <- apply(GD, 2, function(x) {
  freq <- sum(x, na.rm = TRUE) / (2 * length(x))  # allele frequency
  pmin(freq, 1 - freq)                             # minor allele frequency
})

# Keep only SNPs with MAF >= 0.05 (adjust threshold as needed)
keep_snps <- names(maf[maf >= 0.05])
GD <- GD[, keep_snps]
GM <- GM[GM$SNP %in% keep_snps, ]
cat("SNPs kept after MAF >= 0.05:", ncol(GD), "\n")

# -- Save SNP mapping ---------------------------------------------------------
#dir.create("results/tables", showWarnings = FALSE, recursive = TRUE)
write.csv(snp_map, "results/tables/snp_name_mapping.csv", row.names = FALSE)
cat("SNP names mapping saved -> results/tables/snp_name_mapping.csv\n")

# -- Covariates ---------------------------------------------------------------
CV <- cv_farmcpu

# Remove monomorphic SNPs (all 0s, 1s, or 2s) or SNPs with any NA
valid_snps <- apply(GD, 2, function(x) length(unique(x[!is.na(x)])) > 1)
cat("Removing", ncol(GD) - sum(valid_snps), "monomorphic or NA SNPs\n")
GD <- GD[, valid_snps]
GM <- GM[valid_snps, ]


# GAPIT FarmCPU indexing bug workaround
GD <- GD[, -ncol(GD)]
GM <- GM[-nrow(GM), ]

cat("Final SNP count:", ncol(GD), "\n")


# -- Final alignment checks ---------------------------------------------------
stopifnot(
  nrow(pheno_gapit) == nrow(GD),
  all(pheno_gapit$Taxa == rownames(GD)),
  all(pheno_gapit$Taxa == CV$Taxa),
  ncol(GD) == nrow(GM),
  all(colnames(GD) == GM$SNP)
)
cat("\n--- All input checks passed ---\n")

# -- Run FarmCPU ---------------------------------------------------------------
dir.create("results/gwas", showWarnings = FALSE, recursive = TRUE)
old_wd <- getwd()
setwd("results/gwas")

# Convert the existing K_mat to GAPIT's expected format
# GAPIT wants kinship as a data.frame with Taxa as first column
KI_gapit <- cbind(
  data.frame(Taxa = rownames(K_mat)),
  as.data.frame(K_mat)
)

# Convert GD matrix to data.frame with Taxa as first column
# This routes through GAPIT's alternate input path and avoids the
# GD/GI dimension mismatch in myGenotype$GD[, !is.na(myGenotype$GI[,1])]
GD_df <- cbind(
  data.frame(Taxa = rownames(GD), stringsAsFactors = FALSE),
  as.data.frame(GD)
)
rownames(GD_df) <- NULL

cat("GD_df dimensions:", nrow(GD_df), "rows x", ncol(GD_df), "cols\n")
cat("First column name:", names(GD_df)[1], "\n")
cat("First few column names:", names(GD_df)[1:5], "\n")

setwd(PROJECT_DIR); setwd("results/gwas")
gwas_out <- GAPIT(
  Y                = pheno_gapit,
  GD               = GD_df,
  GM               = GM,
  CV               = CV,
  KI               = KI_gapit,
  model            = "FarmCPU",
  PCA.total        = 0,
  SNP.MAF          = 0,
  file.output      = TRUE,
  Geno.View.output = FALSE,
  PCA.View.output  = FALSE
)
setwd(PROJECT_DIR)


cat("FarmCPU complete. Output written to results/gwas/\n")

# -- Load results and plot with ggplot2 ----------------------------------------
res_files <- list.files("results/gwas",
                        pattern = "FarmCPU.*Results.csv",
                        full.names = TRUE)

if (length(res_files) == 0) {
  warning("FarmCPU results CSV not found in results/gwas/ — check GAPIT3 output.")
} else {
  
  gwas_df <- read.csv(res_files[1], stringsAsFactors = FALSE)
  
  # Standardise column name
  names(gwas_df)[names(gwas_df) == "P.value"] <- "P"
  gwas_df <- gwas_df[!is.na(gwas_df$P) & gwas_df$P > 0, ]
  gwas_df$logP <- -log10(gwas_df$P)
  gwas_df$Chr  <- as.integer(gwas_df$Chr)
  gwas_df      <- gwas_df[order(gwas_df$Chr, gwas_df$Pos), ]
  
  # Recover original SNP names for the significant hits table
  gwas_df <- left_join(gwas_df,
                       snp_map[, c("ID_short", "ID_orig")],
                       by = c("SNP" = "ID_short"))
  
  # Significance thresholds (Bonferroni)
  gwt  <- -log10(0.05 / nrow(gwas_df))
  sugt <- -log10(1    / nrow(gwas_df))
  cat("Bonferroni threshold : -log10(p) =", round(gwt,  2), "\n")
  cat("Suggestive threshold : -log10(p) =", round(sugt, 2), "\n")
  
  # Cumulative x-axis positions
  chr_meta <- gwas_df |>
    group_by(Chr) |>
    summarise(max_pos = max(Pos), .groups = "drop") |>
    mutate(cum_add = lag(cumsum(as.numeric(max_pos)), default = 0))
  
  gwas_df <- left_join(gwas_df, chr_meta, by = "Chr") |>
    mutate(pos_cum = Pos + cum_add)
  
  chr_ctrs <- gwas_df |>
    group_by(Chr) |>
    summarise(centre = mean(pos_cum), .groups = "drop")
  
  gwas_df$sig <- with(gwas_df,
                      ifelse(logP >= gwt,  "Genome-wide",
                             ifelse(logP >= sugt, "Suggestive", "Null")))
  
  cat("Genome-wide hits:", sum(gwas_df$sig == "Genome-wide"), "\n")
  cat("Suggestive hits :", sum(gwas_df$sig == "Suggestive"),  "\n")
  
  chr_rects <- gwas_df |>
    group_by(Chr) |>
    summarise(xmin = min(pos_cum), xmax = max(pos_cum), .groups = "drop") |>
    filter(Chr %% 2 == 0)
  
  # Manhattan plot
  
  # Save the Manhattan + QQ plots to PDF
  dir.create("results/plots", showWarnings = FALSE, recursive = TRUE)
  
  chr_rects <- gwas_df |>
    group_by(Chr) |>
    summarise(xmin = min(pos_cum), xmax = max(pos_cum), .groups = "drop") |>
    filter(Chr %% 2 == 0)
  
  p_man <- ggplot(gwas_df, aes(x = pos_cum, y = logP, colour = sig)) +
    geom_rect(data = chr_rects,
              aes(xmin = xmin, xmax = xmax, ymin = -Inf, ymax = Inf),
              inherit.aes = FALSE, fill = "#EAF4EE", alpha = 0.55) +
    geom_point(data = subset(gwas_df, sig == "Null"),
               size = 0.8, alpha = 0.50, shape = 16) +
    geom_point(data = subset(gwas_df, sig == "Suggestive"),
               size = 2.2, alpha = 0.90, shape = 16) +
    geom_point(data = subset(gwas_df, sig == "Genome-wide"),
               size = 3.2, alpha = 0.95, shape = 16) +
    geom_hline(yintercept = gwt,  colour = "#B85042",
               linewidth = 0.9, linetype = "dashed") +
    geom_hline(yintercept = sugt, colour = "#C7A84F",
               linewidth = 0.7, linetype = "dotted") +
    annotate("text", x = Inf, y = gwt + 0.1, hjust = 1.05,
             label = paste0("Bonferroni (p < ",
                            formatC(0.05/nrow(gwas_df), format="e", digits=1), ")"),
             colour = "#B85042", size = 3.0) +
    scale_x_continuous(breaks = chr_ctrs$centre, labels = chr_ctrs$Chr,
                       expand = expansion(mult = 0.01)) +
    scale_colour_manual(
      values = c("Null"="#2E6B4F", "Suggestive"="#C7A84F", "Genome-wide"="#B85042"),
      guide = "none") +
    labs(title    = paste0("Manhattan Plot — FarmCPU | Trait: ", TRAIT),
         subtitle = paste0("n = ", length(y), " accessions (BLUPs) | ",
                           nrow(gwas_df), " SNPs | Covariates: PC1-PC3"),
         x = "Chromosome",
         y = expression(-log[10](italic(p)))) +
    theme_classic(base_size = 11) +
    theme(plot.title    = element_text(face = "bold", colour = "#1B3A2D"),
          plot.subtitle = element_text(colour = "#5C7A65", size = 9),
          axis.text.x   = element_text(size = 8))
  
  n_s     <- nrow(gwas_df)
  obs_lp  <- sort(gwas_df$logP, decreasing = TRUE)
  exp_lp  <- sort(-log10(seq(1/n_s, 1, length.out = n_s)), decreasing = TRUE)
  chi_obs <- qchisq(10^(-gwas_df$logP), df = 1, lower.tail = FALSE)
  lambda  <- median(chi_obs, na.rm = TRUE) / qchisq(0.5, df = 1)
  cat("Lambda GC:", round(lambda, 4), "\n")
  
  ci_u <- -log10(qbeta(0.025, seq_len(n_s), n_s - seq_len(n_s) + 1))
  ci_l <- -log10(qbeta(0.975, seq_len(n_s), n_s - seq_len(n_s) + 1))
  
  qq_df <- data.frame(
    exp = sort(exp_lp, decreasing = TRUE),
    obs = obs_lp,
    cil = sort(ci_l, decreasing = TRUE),
    ciu = sort(ci_u, decreasing = TRUE),
    sig = ifelse(obs_lp >= gwt,  "Genome-wide",
                 ifelse(obs_lp >= sugt, "Suggestive", "Null"))
  )
  
  p_qq <- ggplot(qq_df, aes(x = exp, y = obs)) +
    geom_ribbon(aes(ymin = cil, ymax = ciu), fill = "#5BAD7F", alpha = 0.18) +
    geom_abline(slope = 1, intercept = 0, colour = "#AAAAAA",
                linewidth = 1.0, linetype = "dashed") +
    geom_point(data = subset(qq_df, sig == "Null"),
               colour = "#2E6B4F", size = 1.0, alpha = 0.55) +
    geom_point(data = subset(qq_df, sig == "Suggestive"),
               colour = "#C7A84F", size = 2.0, alpha = 0.90) +
    geom_point(data = subset(qq_df, sig == "Genome-wide"),
               colour = "#B85042", size = 2.8, alpha = 0.95) +
    annotate("text", x = 0.1, y = max(obs_lp) * 0.93, hjust = 0, vjust = 1,
             label = paste0("lambda_GC = ", round(lambda, 3)),
             colour = "#1B3A2D", size = 4.0, fontface = "bold") +
    labs(title = "QQ Plot",
         x = expression("Expected " * -log[10](italic(p))),
         y = expression("Observed "  * -log[10](italic(p)))) +
    theme_classic(base_size = 11) +
    theme(plot.title = element_text(face = "bold", colour = "#1B3A2D"))
  
  p_gwas <- (p_man + p_qq) +
    plot_layout(widths = c(3, 1)) +
    plot_annotation(
      caption = "FarmCPU (GAPIT) | CSES 7160 | Fritzner Pierre",
      theme   = theme(plot.caption = element_text(colour = "#5C7A65", size = 8)))
  
  pdf("results/plots/04_gwas_manhattan_qq.pdf", width = 16, height = 5.5)
  print(p_gwas)
  dev.off()
  cat("Plot saved: results/plots/04_gwas_manhattan_qq.pdf\n")
  
  # Save tables
  dir.create("results/tables", showWarnings = FALSE, recursive = TRUE)
  write.csv(gwas_df[, c("SNP","ID_orig","Chr","Pos","P","logP","sig")],
            "results/tables/gwas_all_results_FarmCPU.csv", row.names = FALSE)
  
  sig_hits <- gwas_df[gwas_df$sig != "Null",
                      c("SNP","ID_orig","Chr","Pos","P","logP","sig")] |>
    rename(SNP_short = SNP, SNP_orig = ID_orig) |>
    arrange(desc(logP))
  write.csv(sig_hits,
            "results/tables/gwas_significant_hits_FarmCPU.csv", row.names = FALSE)
  
  cat("Genome-wide hits:", sum(gwas_df$sig == "Genome-wide"), "\n")
  cat("Suggestive hits :", sum(gwas_df$sig == "Suggestive"), "\n")
  cat("Top suggestive SNPs:\n")
  print(sig_hits[, c("SNP_orig","Chr","Pos","P","logP")])
}

# ==============================================================================
# STEP 10 — GENOMIC PREDICTION: GBLUP + 5-FOLD CV (10 REPS)
# ==============================================================================

cat("----------------------------------------------------------------------\n")
cat("STEP 10 | Genomic Prediction  (GBLUP, 5-fold CV x 10 reps)\n")
cat("----------------------------------------------------------------------\n")

n_total  <- length(y)
cv_store <- data.frame(rep=integer(), fold=integer(),
                       r=numeric(), RMSE=numeric(), bias=numeric())

for (rep in seq_len(N_REPS)) {
  fid <- sample(rep(seq_len(N_FOLDS), length.out=n_total))

  for (fold in seq_len(N_FOLDS)) {
    idx_t  <- which(fid == fold)
    y_cv   <- y;  y_cv[idx_t] <- NA

    fv <- tryCatch(mixed.solve(y=y_cv, K=K_mat), error=function(e) NULL)
    if (!is.null(fv)) {
      yp <- as.numeric(fv$u)[idx_t] + as.numeric(fv$beta)
      yo <- y[idx_t]
      cv_store <- rbind(cv_store, data.frame(
        rep  = rep, fold = fold,
        r    = cor(yo, yp, use="complete.obs"),
        RMSE = sqrt(mean((yo-yp)^2, na.rm=TRUE)),
        bias = coef(lm(yo ~ yp))[2]))
    }
  }
  if (rep %% 2 == 0)
    cat("  Rep", rep, "/", N_REPS, "| mean r =",
        round(mean(cv_store$r, na.rm=TRUE), 4), "\n")
}

cat("\n--- Cross-Validation Summary ---\n")
cat("Mean r  :", round(mean(cv_store$r,    na.rm=TRUE), 4), "\n")
cat("SD r    :", round(sd(cv_store$r,      na.rm=TRUE), 4), "\n")
cat("95% CI  : [",
    round(quantile(cv_store$r, 0.025, na.rm=TRUE), 4), ",",
    round(quantile(cv_store$r, 0.975, na.rm=TRUE), 4), "]\n")
cat("Mean RMSE:", round(mean(cv_store$RMSE, na.rm=TRUE), 4), "\n")
cat("Bias (slope):", round(mean(cv_store$bias, na.rm=TRUE), 4),
    "(1.0 = unbiased)\n\n")

p_cv <- ggplot(cv_store, aes(x=r)) +
  geom_histogram(aes(y=after_stat(density)), binwidth=0.03,
                 fill="#2E6B4F", colour="#1B3A2D", alpha=0.82) +
  geom_density(colour="#1B3A2D", linewidth=1.0) +
  geom_vline(xintercept=mean(cv_store$r, na.rm=TRUE),
             colour="#B85042", linetype="dashed", linewidth=1.0) +
  annotate("text",
           x=mean(cv_store$r, na.rm=TRUE)+0.005, y=Inf,
           vjust=1.5, hjust=0,
           label=paste0("Mean r = ", round(mean(cv_store$r, na.rm=TRUE),3)),
           colour="#B85042", size=3.8) +
  labs(title    = "GBLUP Cross-Validation Accuracy",
       subtitle = paste0(N_FOLDS, "-fold CV x ", N_REPS,
                         " reps | n=", n_total, " accessions"),
       x = "Pearson r (observed vs predicted)", y = "Density") +
  theme_classic(base_size=12) +
  theme(plot.title    = element_text(face="bold", colour="#1B3A2D"),
        plot.subtitle = element_text(colour="#5C7A65"))

fit_full   <- mixed.solve(y=y, K=K_mat)
scatter_df <- data.frame(
  Taxa      = common_taxa,
  observed  = y,
  predicted = as.numeric(fit_full$u) + as.numeric(fit_full$beta)
)
fr <- cor(scatter_df$predicted, scatter_df$observed)

p_scat <- ggplot(scatter_df, aes(x=predicted, y=observed)) +
  geom_point(colour="#2E6B4F", alpha=0.72, size=2.5, shape=16) +
  geom_smooth(method="lm", colour="#B85042", se=TRUE,
              linewidth=1.0, fill="#EAF4EE") +
  geom_abline(slope=1, intercept=0, colour="#AAAAAA",
              linewidth=0.8, linetype="dashed") +
  annotate("text", x=-Inf, y=Inf, hjust=-0.1, vjust=1.5,
           label=paste0("r = ", round(fr,3)),
           colour="#1B3A2D", fontface="bold", size=4.5) +
  labs(title    = "Predicted vs Observed (Full Model)",
       subtitle = paste0("GBLUP | r = ", round(fr,3)),
       x = "Genomic Predicted Value (GBLUP)",
       y = "Observed BLUP of Oleic Acid") +
  theme_classic(base_size=12) +
  theme(plot.title    = element_text(face="bold", colour="#1B3A2D"),
        plot.subtitle = element_text(colour="#5C7A65"))

p_gp <- (p_cv | p_scat) +
  plot_annotation(
    caption = "GBLUP | rrBLUP | CSES 7160 | Fritzner Pierre",
    theme   = theme(plot.caption=element_text(colour="#5C7A65", size=8)))

setwd(PROJECT_DIR)

pdf("results/plots/05_gblup_cv_and_fit.pdf", width=12, height=5)
print(p_gp); dev.off()

write.csv(cv_store,   "results/gp/gblup_cv_results.csv",      row.names=FALSE)
write.csv(scatter_df, "results/gp/gblup_full_predictions.csv", row.names=FALSE)
cat("Plot : results/plots/05_gblup_cv_and_fit.pdf\n")
cat("Table: results/gp/gblup_cv_results.csv\n")
cat("Table: results/gp/gblup_full_predictions.csv\n\n")


# ==============================================================================
# STEP 11 — RESULTS SUMMARY & SESSION INFO
# ==============================================================================

cat("======================================================================\n")
cat("FINAL RESULTS SUMMARY\n")
cat("======================================================================\n")
cat(sprintf("Accessions (n)            : %d\n", length(y)))
cat(sprintf("Replicate observations    : %d\n", nrow(pheno_clean)))
cat(sprintf("SNPs after QC             : %d\n", ncol(geno_num)))
cat(sprintf("Mean oleic acid (raw)     : %.2f %%\n", mean(pheno_clean$Value)))
cat(sprintf("SD oleic acid  (raw)      : %.2f %%\n", sd(pheno_clean$Value)))
cat(sprintf("Genomic variance (sigma2g): %.6f\n", Vu))
cat(sprintf("Residual var   (sigma2e)  : %.6f\n", Ve))
cat(sprintf("Genomic heritability h2g  : %.4f\n", h2g))
cat(sprintf("Bootstrap 95%% CI          : [%.4f, %.4f]\n", ci[1], ci[2]))
if (exists("gwas_df"))
  cat(sprintf("Genome-wide GWAS hits     : %d\n",
              sum(gwas_df$sig=="Genome-wide")))
cat(sprintf("GBLUP CV mean r           : %.4f\n",
            mean(cv_store$r, na.rm=TRUE)))
cat(sprintf("GBLUP CV 95%% CI           : [%.4f, %.4f]\n",
            quantile(cv_store$r, 0.025, na.rm=TRUE),
            quantile(cv_store$r, 0.975, na.rm=TRUE)))
cat("======================================================================\n\n")

tryCatch({
  sink("results/session_info.txt")
  cat("Peanut Oleic Acid Analysis\n")
  cat("Date:", format(Sys.time(), "%Y-%m-%d %H:%M:%S"), "\n\n")
  cat("Files:\n  Pheno:", PHENO_FILE, "\n  GFF3:", GFF3_FILE,
      "\n  Geno (VCF):", GENO_FILE, "\n\n")
  cat("Settings: MAF=", MAF_THOLD, "| Miss=", MISS_THOLD,
      "| GWAS=FarmCPU | GP=GBLUP | CV=", N_FOLDS, "fold x", N_REPS, "reps\n\n")
  print(sessionInfo())
}, finally = {
  sink()
})
cat("Session info: results/session_info.txt\nDone.\n")


# ==============================================================================
# SAVE OBJECTS FOR INDEPENDENT LD ANALYSIS
# ==============================================================================
saveRDS(geno_num,  "results/geno_num.rds")
saveRDS(geno_map,  "results/geno_map.rds")
saveRDS(gwas_df,   "results/gwas_df.rds")
cat("Objects saved for independent LD analysis.\n")



# Extended accuracy diagnostics
r_full   <- cor(scatter_df$observed, scatter_df$predicted)
mae      <- mean(abs(scatter_df$observed - scatter_df$predicted))
slope_op <- coef(lm(observed ~ predicted, data = scatter_df))[2]
sd_ratio <- sd(scatter_df$predicted) / sd(scatter_df$observed)
spearman <- cor(scatter_df$observed, scatter_df$predicted, method = "spearman")
overfit_gap <- r_full - mean(cv_store$r)

cat(sprintf("Full r: %.4f | CV r: %.4f | Gap: %.4f (%.1f%%)\n",
            r_full, mean(cv_store$r), overfit_gap, overfit_gap/r_full*100))
cat(sprintf("MAE: %.4f | Spearman: %.4f\n", mae, spearman))
cat(sprintf("Regression slope (obs~pred): %.4f  (ideal=1.0)\n", slope_op))
cat(sprintf("SD ratio pred/obs: %.4f  (GBLUP shrinkage)\n", sd_ratio))

