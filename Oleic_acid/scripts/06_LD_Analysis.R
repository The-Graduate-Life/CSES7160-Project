# ==============================================================================
# LD Analysis — Top GWAS Hits (Independent Script)
# Genomic Dissection of Oleic Acid Concentration in Peanut
#
# Author  : Fritzner Pierre
# Course  : CSES 7160 — Genetic Data Analysis
#
# Requirements:
#   results/geno_num.rds   — post-QC genotype matrix {0,1,2}
#   results/geno_map.rds   — SNP position map (SNP, Chr, Pos)
#   results/gwas_df.rds    — FarmCPU GWAS results
#
#   These are saved at the end of CSES7160_Fritzner_Analysis.R.
#   Run that script first if the .rds files do not exist yet.
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
# Install once if needed:
#install.packages("LDheatmap")
#BiocManager::install("snpStats")
#devtools::install_github("SFUStatgen/LDheatmap", dependencies = FALSE)

pkgs <- c("ggplot2","dplyr","LDheatmap")

for (p in pkgs) {
  if (!requireNamespace(p, quietly = TRUE)) {
    stop("Package missing: ", p)
  }
}

suppressPackageStartupMessages(
  lapply(pkgs, library, character.only = TRUE)
)


# ==============================================================================
# 2 — SETTINGS
# ==============================================================================
WINDOW_BP  <- 500000   # ±500 kb around the top hit; widen if block looks truncated
R2_STRONG  <- 0.80     # strong LD threshold
R2_MOD     <- 0.50     # moderate LD threshold
TOP_N_FALLBACK <- 5    # how many top SNPs to use if no significant hits exist

set.seed(2026)

for (d in c("results/tables", "results/plots"))
  dir.create(d, showWarnings = FALSE, recursive = TRUE)

cat("======================================================================\n")
cat(" LD Analysis | Peanut Oleic Acid\n")
cat("======================================================================\n\n")


# ==============================================================================
# 3 — LOAD SAVED OBJECTS
# ==============================================================================
cat("Loading saved objects...\n")

rds_files <- c(
  geno_num = "results/geno_num.rds",
  geno_map = "results/geno_map.rds",
  gwas_df  = "results/gwas_df.rds"
)

for (nm in names(rds_files)) {
  if (!file.exists(rds_files[nm]))
    stop("File not found: ", rds_files[nm],
         "\nRun CSES7160_Fritzner_Analysis.R first to generate it.")
  assign(nm, readRDS(rds_files[nm]))
  cat(" Loaded:", rds_files[nm], "\n")
}

cat("\ngeno_num  :", nrow(geno_num), "taxa x", ncol(geno_num), "SNPs\n")
cat("geno_map  :", nrow(geno_map), "SNPs\n")
cat("gwas_df   :", nrow(gwas_df),  "SNPs\n\n")


# ==============================================================================
# 4 — IDENTIFY TOP GWAS HITS
# ==============================================================================
cat("----------------------------------------------------------------------\n")
cat("STEP 4 | Identify Top GWAS Hits\n")
cat("----------------------------------------------------------------------\n")

sig_snps <- gwas_df[gwas_df$sig != "Null", ]

if (nrow(sig_snps) == 0) {
  cat("No significant/suggestive hits found.\n")
  cat("Falling back to top", TOP_N_FALLBACK, "SNPs by -log10(p).\n")
  sig_snps <- gwas_df[order(gwas_df$logP, decreasing = TRUE), ][1:TOP_N_FALLBACK, ]
}

cat("Seed SNPs for LD analysis:\n")
print(sig_snps[, c("SNP", "Chr", "Pos", "logP", "sig")])

# Top hit = SNP with the highest -log10(p)
top_hit   <- sig_snps[which.max(sig_snps$logP), ]
focal_chr <- top_hit$Chr
focal_pos <- top_hit$Pos

cat(sprintf("\nFocal SNP : %s  (Chr %d, Pos %d)\n",
            top_hit$SNP, focal_chr, focal_pos))
cat(sprintf("Window    : ±%d kb\n\n", WINDOW_BP / 1000))


# ==============================================================================
# 5 — EXTRACT REGIONAL GENOTYPE SUBMATRIX
# ==============================================================================
cat("----------------------------------------------------------------------\n")
cat("STEP 5 | Extract Regional Genotype Submatrix\n")
cat("----------------------------------------------------------------------\n")

region_idx <- which(
  geno_map$Chr == focal_chr &
    geno_map$Pos >= (focal_pos - WINDOW_BP) &
    geno_map$Pos <= (focal_pos + WINDOW_BP)
)

cat("SNPs in window:", length(region_idx), "\n")

if (length(region_idx) < 2)
  stop(
    "Fewer than 2 SNPs found in the ±", WINDOW_BP / 1000, " kb window.\n",
    "Try increasing WINDOW_BP, or check that Chr/Pos coding matches between\n",
    "geno_map and gwas_df (both use integer chromosomes 1-20)."
  )

geno_region <- geno_num[, region_idx, drop = FALSE]   # taxa × region SNPs
map_region  <- geno_map[region_idx, ]

cat("Regional matrix:", nrow(geno_region), "taxa x",
    ncol(geno_region), "SNPs\n\n")


# ==============================================================================
# 6 — COMPUTE PAIRWISE r²
# ==============================================================================
cat("----------------------------------------------------------------------\n")
cat("STEP 6 | Compute Pairwise r²\n")
cat("----------------------------------------------------------------------\n")

compute_r2_matrix <- function(G) {
  # G : n_taxa × n_snp matrix {0, 1, 2}
  # Returns an n_snp × n_snp r² matrix
  G_sc <- scale(G, center = TRUE, scale = TRUE)
  G_sc[is.na(G_sc)] <- 0        # residual NAs → 0 after mean imputation
  n  <- nrow(G_sc)
  R  <- crossprod(G_sc) / (n - 1)   # correlation matrix
  R^2
}

cat("Computing pairwise r² for", ncol(geno_region), "SNPs...\n")
r2_mat <- compute_r2_matrix(geno_region)
rownames(r2_mat) <- colnames(r2_mat) <- map_region$SNP
cat("r² matrix:", nrow(r2_mat), "x", ncol(r2_mat), "\n\n")


# ==============================================================================
# STEP 7 — FLAG SNPs IN LD WITH THE TOP HIT
# ==============================================================================

cat("----------------------------------------------------------------------\n")
cat("STEP 7 | Flag SNPs in LD with Top Hit\n")
cat("----------------------------------------------------------------------\n")

# -- Load SNP name mapping (short "SNP####" -> original "CHROM_POS") ----------
top_snp_name <- gwas_df$SNP[which.min(gwas_df$P)]

snp_map_file <- "results/tables/snp_name_mapping.csv"

if (!file.exists(snp_map_file))
  stop("snp_name_mapping.csv not found.\n",
       "Make sure CSES7160_Fritzner_Analysis.R has completed Step 9.")

snp_map <- read.csv(snp_map_file, stringsAsFactors = FALSE)

# -- Resolve top hit name to original CHROM_POS name -------------------------
# gwas_df uses short IDs (e.g. "SNP1042"); r2_mat uses original CHROM_POS IDs
top_snp_orig <- snp_map$ID_orig[snp_map$ID_short == top_snp_name]

if (length(top_snp_orig) == 0)
  stop("Could not resolve original name for: ", top_snp_name,
       "\nCheck snp_name_mapping.csv.")

cat("Short name :", top_snp_name, "\n")
cat("Orig  name :", top_snp_orig, "\n")

# -- Match in r2_mat ----------------------------------------------------------
match_row <- which(rownames(r2_mat) == top_snp_orig)

if (length(match_row) == 0)
  stop("'", top_snp_orig, "' not found in r2_mat.\n",
       "First 5 rownames: ",
       paste(head(rownames(r2_mat), 5), collapse = ", "))

r2_with_top <- r2_mat[match_row, ]
cat("r2_with_top length:", length(r2_with_top), "\n\n")

# -- Build LD table -----------------------------------------------------------
ld_table <- data.frame(
  SNP         = names(r2_with_top),
  Chr         = map_region$Chr,
  Pos         = map_region$Pos,
  r2_with_top = as.numeric(r2_with_top),
  stringsAsFactors = FALSE
) |> arrange(desc(r2_with_top))

# -- Classify by LD strength --------------------------------------------------
strong_ld <- ld_table[ld_table$r2_with_top >= R2_STRONG, ]
mod_ld    <- ld_table[ld_table$r2_with_top >= R2_MOD &
                        ld_table$r2_with_top <  R2_STRONG, ]

cat(sprintf("Strong LD (r² ≥ %.2f)              : %d SNPs\n",
            R2_STRONG, nrow(strong_ld)))
print(strong_ld)

cat(sprintf("\nModerate LD (%.2f ≤ r² < %.2f) : %d SNPs\n",
            R2_MOD, R2_STRONG, nrow(mod_ld)))
print(mod_ld)

# -- Save ---------------------------------------------------------------------
write.csv(ld_table,
          "results/tables/ld_with_top_hit.csv",
          row.names = FALSE)
cat("\nTable saved: results/tables/ld_with_top_hit.csv\n\n")

# ==============================================================================
# STEP 8 | LD Heatmap
# ==============================================================================
cat("----------------------------------------------------------------------\n")
cat("STEP 8 | LD Heatmap\n")
cat("----------------------------------------------------------------------\n")

positions_mb <- map_region$Pos / 1e6   # bp → Mb for axis labels

# LDheatmap requires:
#   (1) a square symmetric matrix of LD values (r²)
#   (2) genetic.distances of the same length as nrow/ncol of the matrix
# r2_mat is already n_snp × n_snp — use it directly

cat("r2_mat dimensions  :", nrow(r2_mat), "x", ncol(r2_mat), "\n")
cat("positions_mb length:", length(positions_mb), "\n")

# Confirm dimensions match before calling LDheatmap
stopifnot(
  nrow(r2_mat) == ncol(r2_mat),
  nrow(r2_mat) == length(positions_mb)
)

cat("Rendering LD heatmap...\n")
pdf("results/plots/06_ld_heatmap_top_region.pdf", width = 8, height = 7)

LDheatmap(
  r2_mat,
  genetic.distances = positions_mb,
  distances         = "physical",
  color             = heat.colors(20),
  title             = paste0("LD Heatmap — Chr ", focal_chr,
                             "  ±", WINDOW_BP / 1000, " kb | Top hit: ",
                             top_snp_name),
  add.map           = TRUE,
  name              = "LD_top_region"
)

dev.off()
cat("Plot saved: results/plots/06_ld_heatmap_top_region.pdf\n\n")


# ==============================================================================
# 9 — LD DECAY BAR PLOT
# ==============================================================================
cat("----------------------------------------------------------------------\n")
cat("STEP 9 | LD Decay Bar Plot\n")
cat("----------------------------------------------------------------------\n")

ld_ordered <- ld_table[order(ld_table$Pos), ]

p_ld_bar <- ggplot(ld_ordered,
                   aes(x     = Pos / 1e6,
                       y     = r2_with_top,
                       fill  = r2_with_top >= R2_STRONG)) +
  geom_col(width = 0.01, alpha = 0.85) +
  geom_hline(yintercept = R2_STRONG, colour = "#B85042",
             linetype = "dashed", linewidth = 0.9) +
  geom_hline(yintercept = R2_MOD, colour = "#C7A84F",
             linetype = "dotted", linewidth = 0.7) +
  annotate("text",
           x     = min(ld_ordered$Pos) / 1e6,
           y     = R2_STRONG + 0.02,
           label = paste0("r² = ", R2_STRONG, " (strong LD)"),
           hjust = 0, colour = "#B85042", size = 3.2) +
  annotate("text",
           x     = min(ld_ordered$Pos) / 1e6,
           y     = R2_MOD + 0.02,
           label = paste0("r² = ", R2_MOD, " (moderate LD)"),
           hjust = 0, colour = "#C7A84F", size = 3.2) +
  scale_fill_manual(
    values = c("TRUE" = "#B85042", "FALSE" = "#2E6B4F"),
    guide  = "none") +
  labs(
    title    = paste0("LD Decay from Top Hit — Chr ", focal_chr),
    subtitle = paste0("Focal SNP: ", top_snp_name,
                      "  |  Window ±", WINDOW_BP / 1000, " kb",
                      "  |  n = ", nrow(ld_ordered), " SNPs"),
    x = "Position (Mb)",
    y = expression(r^2 ~ "with top hit"),
    caption = ""
  ) +
  theme_classic(base_size = 11) +
  theme(
    plot.title    = element_text(face = "bold", colour = "#1B3A2D"),
    plot.subtitle = element_text(colour = "#5C7A65"),
    plot.caption  = element_text(colour = "#5C7A65", size = 8)
  )

pdf("results/plots/07_ld_decay_top_hit.pdf", width = 9, height = 4)
print(p_ld_bar)
dev.off()
cat("Plot saved: results/plots/07_ld_decay_top_hit.pdf\n\n")


# ==============================================================================
# 10 — SUMMARY
# ==============================================================================
cat("======================================================================\n")
cat("LD ANALYSIS SUMMARY\n")
cat("======================================================================\n")
cat(sprintf("Focal SNP              : %s\n", top_snp_name))
cat(sprintf("Chromosome             : %d\n", focal_chr))
cat(sprintf("Position               : %d bp\n", focal_pos))
cat(sprintf("Window                 : ±%d kb\n", WINDOW_BP / 1000))
cat(sprintf("SNPs in window         : %d\n", nrow(ld_ordered)))
cat(sprintf("Strong LD (r²≥%.2f)   : %d SNPs\n", R2_STRONG, nrow(strong_ld)))
cat(sprintf("Moderate LD (r²≥%.2f) : %d SNPs\n", R2_MOD,    nrow(mod_ld)))
cat("----------------------------------------------------------------------\n")
cat("Output files:\n")
cat("  results/tables/ld_with_top_hit.csv\n")
cat("  results/plots/06_ld_heatmap_top_region.pdf\n")
cat("  results/plots/07_ld_decay_top_hit.pdf\n")
cat("======================================================================\n")
cat("Done.\n")


