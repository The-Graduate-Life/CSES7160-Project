# ==============================================================================
# Manhattan Plot Confirmation & LD Verification (Independent Script)
# Genomic Dissection of Oleic Acid Concentration in Peanut
#
# Author  : Fritzner Pierre
# Course  : CSES 7160 — Genetic Data Analysis
#
# Requirements:
#   results/gwas_df.rds                    — FarmCPU GWAS results
#   results/tables/snp_name_mapping.csv    — short ID <-> original ID mapping
#
# Output:
#   results/plots/08_manhattan_confirmed_top_pair.pdf
#   results/plots/09_chr16_zoom.pdf
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
# install.packages(c("ggplot2", "dplyr", "ggrepel"))

pkgs <- c("ggplot2","dplyr","ggrepel")

missing <- pkgs[!sapply(pkgs, requireNamespace, quietly = TRUE)]
if (length(missing) > 0)
  stop("Missing packages: ", paste(missing, collapse = ", "))

suppressPackageStartupMessages(
  lapply(pkgs, library, character.only = TRUE)
)


# ==============================================================================
# 2 — SETTINGS
# ==============================================================================
# The two perfectly linked SNPs (r² = 1.0) identified in the LD analysis
TOP_SNP_ORIGINALS <- c(
  "araip.K30076.gnm1.Araip.B06_2010906",
  "araip.K30076.gnm1.Araip.B06_2335571"
)

# Chr 16 block boundaries (from LD analysis) ± buffer
BLOCK_CHR    <- 16
BLOCK_START  <- 2010906 - 500000
BLOCK_END    <- 2335571 + 500000

for (d in c("results/plots", "results/tables"))
  dir.create(d, showWarnings = FALSE, recursive = TRUE)

cat("======================================================================\n")
cat(" Manhattan Confirmation | Fritzner Pierre | Peanut Oleic Acid\n")
cat("======================================================================\n\n")


# ==============================================================================
# 3 — LOAD DATA
# ==============================================================================
cat("Loading saved objects...\n")

# GWAS results
if (!file.exists("results/gwas_df.rds"))
  stop("results/gwas_df.rds not found.\n",
       "Run CSES7160_Fritzner_Analysis.R first.")
gwas_df <- readRDS("results/gwas_df.rds")
cat(" Loaded: results/gwas_df.rds —", nrow(gwas_df), "SNPs\n")

# SNP name mapping
if (!file.exists("results/tables/snp_name_mapping.csv"))
  stop("snp_name_mapping.csv not found.\n",
       "Run CSES7160_Fritzner_Analysis.R first.")
snp_map <- read.csv("results/tables/snp_name_mapping.csv",
                    stringsAsFactors = FALSE)
cat(" Loaded: results/tables/snp_name_mapping.csv —",
    nrow(snp_map), "entries\n\n")


# ==============================================================================
# 4 — RESOLVE SNP NAMES & TAG TOP PAIR
# ==============================================================================
cat("----------------------------------------------------------------------\n")
cat("STEP 4 | Resolve SNP Names and Tag Top Pair\n")
cat("----------------------------------------------------------------------\n")

# Resolve original names -> short names used in gwas_df
top_snp_shorts <- snp_map$ID_short[snp_map$ID_orig %in% TOP_SNP_ORIGINALS]

if (length(top_snp_shorts) == 0)
  stop("Could not resolve short names for top SNP pair.\n",
       "Check TOP_SNP_ORIGINALS match entries in snp_name_mapping.csv.")

cat("Top SNP pair resolved:\n")
print(data.frame(
  Original = TOP_SNP_ORIGINALS,
  Short    = top_snp_shorts
))

# Tag SNPs in gwas_df by class
gwas_df <- gwas_df |>
  mutate(highlight = case_when(
    SNP %in% top_snp_shorts      ~ "Top pair (r²=1.0)",
    sig == "Genome-wide"         ~ "Genome-wide",
    sig == "Suggestive"          ~ "Suggestive",
    TRUE                         ~ "Null"
  ))

cat("\nSNP class counts:\n")
print(table(gwas_df$highlight))
cat("\n")


# ==============================================================================
# 5 — RECOMPUTE CUMULATIVE POSITIONS (CLEAN)
# ==============================================================================
cat("----------------------------------------------------------------------\n")
cat("STEP 5 | Recompute Cumulative Positions\n")
cat("----------------------------------------------------------------------\n")

# Drop any pre-existing cumulative columns to avoid join conflicts
# Drop pre-existing cumulative columns (base R — compatible with all dplyr versions)
drop_cols <- c("cum_add", "pos_cum", "max_pos")
gwas_df   <- gwas_df[, !names(gwas_df) %in% drop_cols, drop = FALSE]

# Recompute from scratch
chr_meta <- gwas_df |>
  group_by(Chr) |>
  summarise(max_pos = max(Pos), .groups = "drop") |>
  mutate(cum_add = lag(cumsum(as.numeric(max_pos)), default = 0))

gwas_df <- left_join(gwas_df, chr_meta, by = "Chr") |>
  mutate(pos_cum = Pos + cum_add)

chr_ctrs <- gwas_df |>
  group_by(Chr) |>
  summarise(centre = mean(pos_cum), .groups = "drop")

chr_rects <- gwas_df |>
  group_by(Chr) |>
  summarise(xmin = min(pos_cum), xmax = max(pos_cum), .groups = "drop") |>
  filter(Chr %% 2 == 0)

cat("Chromosomes found     :", n_distinct(gwas_df$Chr), "\n")
cat("Cumulative positions  : OK\n\n")


# ==============================================================================
# 6 — SIGNIFICANCE THRESHOLDS
# ==============================================================================
gwt  <- -log10(0.05 / nrow(gwas_df))   # Bonferroni genome-wide
sugt <- -log10(1    / nrow(gwas_df))   # Suggestive

cat("Bonferroni threshold  : -log10(p) =", round(gwt,  3), "\n")
cat("Suggestive threshold  : -log10(p) =", round(sugt, 3), "\n\n")


# ==============================================================================
# 7 — COLOUR / SIZE / SHAPE SCHEME
# ==============================================================================
col_vals <- c(
  "Null"              = "#2E6B4F",
  "Suggestive"        = "#C7A84F",
  "Genome-wide"       = "#B85042",
  "Top pair (r²=1.0)" = "#7B2FBE"
)
size_vals <- c(
  "Null"              = 0.8,
  "Suggestive"        = 2.2,
  "Genome-wide"       = 3.2,
  "Top pair (r²=1.0)" = 4.5
)
shape_vals <- c(
  "Null"              = 16,
  "Suggestive"        = 16,
  "Genome-wide"       = 16,
  "Top pair (r²=1.0)" = 18
)


# ==============================================================================
# 8 — GENOME-WIDE MANHATTAN PLOT WITH HIGHLIGHTED TOP PAIR
# ==============================================================================
cat("----------------------------------------------------------------------\n")
cat("STEP 8 | Genome-wide Manhattan Plot\n")
cat("----------------------------------------------------------------------\n")

p_confirm <- ggplot(gwas_df,
                    aes(x      = pos_cum,
                        y      = logP,
                        colour = highlight,
                        size   = highlight,
                        shape  = highlight)) +
  
  # Alternating chromosome bands
  geom_rect(data        = chr_rects,
            aes(xmin = xmin, xmax = xmax,
                ymin = -Inf, ymax = Inf),
            inherit.aes = FALSE,
            fill        = "#EAF4EE",
            alpha       = 0.55) +
  
  # Null SNPs (bottom layer)
  geom_point(data  = subset(gwas_df, highlight == "Null"),
             alpha = 0.45) +
  
  # Suggestive and genome-wide
  geom_point(data  = subset(gwas_df,
                            highlight %in% c("Suggestive", "Genome-wide")),
             alpha = 0.90) +
  
  # Top pair — uppermost layer
  geom_point(data  = subset(gwas_df, highlight == "Top pair (r²=1.0)"),
             alpha = 1.00) +
  
  # Threshold lines
  geom_hline(yintercept = gwt,  colour = "#B85042",
             linewidth  = 0.9,  linetype = "dashed") +
  geom_hline(yintercept = sugt, colour = "#C7A84F",
             linewidth  = 0.7,  linetype = "dotted") +
  
  # Labels for the two highlighted SNPs
  geom_label_repel(
    data        = subset(gwas_df, highlight == "Top pair (r²=1.0)"),
    aes(label   = paste0("Chr16: ", round(Pos / 1e6, 2), " Mb")),
    colour      = "#7B2FBE",
    size        = 3.2,
    fontface    = "bold",
    box.padding = 0.5,
    nudge_y     = 0.4,
    show.legend = FALSE
  ) +
  
  scale_x_continuous(
    breaks = chr_ctrs$centre,
    labels = chr_ctrs$Chr,
    expand = expansion(mult = 0.01)) +
  scale_colour_manual(values = col_vals,  name = "SNP class") +
  scale_size_manual(  values = size_vals, name = "SNP class") +
  scale_shape_manual( values = shape_vals,name = "SNP class") +
  
  labs(
    title    = "Manhattan Plot — Confirming Top Signal on Chr 16",
    subtitle = paste0("Purple diamonds = r²=1.0 pair  |  n = ",
                      nrow(gwas_df), " SNPs  |  FarmCPU"),
    x        = "Chromosome",
    y        = expression(-log[10](italic(p))),
    caption  = "CSES 7160 | Fritzner Pierre"
  ) +
  theme_classic(base_size = 11) +
  theme(
    plot.title      = element_text(face = "bold", colour = "#1B3A2D"),
    plot.subtitle   = element_text(colour = "#5C7A65", size = 9),
    plot.caption    = element_text(colour = "#5C7A65", size = 8),
    axis.text.x     = element_text(size = 8),
    legend.position = "right"
  )

pdf("results/plots/08_manhattan_confirmed_top_pair.pdf",
    width = 16, height = 5.5)
print(p_confirm)
dev.off()
cat("Plot saved: results/plots/08_manhattan_confirmed_top_pair.pdf\n\n")


# ==============================================================================
# 9 — CHECK FOR SECOND INDEPENDENT PEAKS
# ==============================================================================
cat("----------------------------------------------------------------------\n")
cat("STEP 9 | Check for Second Independent Peaks\n")
cat("----------------------------------------------------------------------\n")

other_peaks <- gwas_df |>
  filter(sig != "Null") |>
  mutate(in_chr16_block = (Chr  == BLOCK_CHR  &
                             Pos  >= BLOCK_START &
                             Pos  <= BLOCK_END)) |>
  filter(!in_chr16_block) |>
  arrange(desc(logP))

if (nrow(other_peaks) == 0) {
  cat("No significant or suggestive peaks outside the Chr 16 block.\n")
  cat("Conclusion: Chr 16 carries the ONLY major signal genome-wide.\n\n")
} else {
  cat("Potential second independent peaks detected:\n")
  print(other_peaks[, c("SNP", "Chr", "Pos", "logP", "sig")])
  cat("\nThese chromosomes may warrant a separate LD analysis.\n\n")
  
  # Save for reference
  write.csv(other_peaks[, c("SNP", "Chr", "Pos", "logP", "sig")],
            "results/tables/second_peaks_candidates.csv",
            row.names = FALSE)
  cat("Table saved: results/tables/second_peaks_candidates.csv\n\n")
}


# ==============================================================================
# 10 — CHR 16 ZOOM PLOT
# ==============================================================================
cat("----------------------------------------------------------------------\n")
cat("STEP 10 | Chr 16 Zoom Plot\n")
cat("----------------------------------------------------------------------\n")

chr16_df <- gwas_df[gwas_df$Chr == 16, ]
cat("SNPs on Chr 16:", nrow(chr16_df), "\n")

p_zoom <- ggplot(chr16_df,
                 aes(x      = Pos / 1e6,
                     y      = logP,
                     colour = highlight,
                     size   = highlight,
                     shape  = highlight)) +
  
  geom_point(alpha = 0.80) +
  
  geom_hline(yintercept = gwt,  colour = "#B85042",
             linewidth  = 0.9,  linetype = "dashed") +
  geom_hline(yintercept = sugt, colour = "#C7A84F",
             linewidth  = 0.7,  linetype = "dotted") +
  
  geom_label_repel(
    data        = subset(chr16_df, highlight == "Top pair (r²=1.0)"),
    aes(label   = paste0(round(Pos / 1e6, 3), " Mb")),
    colour      = "#7B2FBE",
    size        = 3.5,
    fontface    = "bold",
    box.padding = 0.6,
    nudge_y     = 0.3,
    show.legend = FALSE
  ) +
  
  scale_colour_manual(values = col_vals,  name = "SNP class") +
  scale_size_manual(  values = size_vals, name = "SNP class") +
  scale_shape_manual( values = shape_vals,name = "SNP class") +
  
  labs(
    title    = "Chr 16 (Araip.B06) — Zoomed Signal",
    subtitle = paste0("Purple diamonds = r²=1.0 pair  |  ",
                      nrow(chr16_df), " SNPs on Chr 16"),
    x        = "Position (Mb)",
    y        = expression(-log[10](italic(p))),
    caption  = "CSES 7160 | Fritzner Pierre"
  ) +
  theme_classic(base_size = 11) +
  theme(
    plot.title      = element_text(face = "bold", colour = "#1B3A2D"),
    plot.subtitle   = element_text(colour = "#5C7A65", size = 9),
    plot.caption    = element_text(colour = "#5C7A65", size = 8),
    legend.position = "right"
  )

pdf("results/plots/09_chr16_zoom.pdf", width = 10, height = 5)
print(p_zoom)
dev.off()
cat("Plot saved: results/plots/09_chr16_zoom.pdf\n\n")


# ==============================================================================
# 11 — SUMMARY
# ==============================================================================
cat("======================================================================\n")
cat("SUMMARY\n")
cat("======================================================================\n")
cat(sprintf("Total SNPs in GWAS        : %d\n",  nrow(gwas_df)))
cat(sprintf("Genome-wide hits          : %d\n",
            sum(gwas_df$sig == "Genome-wide")))
cat(sprintf("Suggestive hits           : %d\n",
            sum(gwas_df$sig == "Suggestive")))
cat(sprintf("Top pair (r²=1.0) on Chr  : %d\n",  BLOCK_CHR))
cat(sprintf("Block boundaries          : %.2f – %.2f Mb\n",
            BLOCK_START / 1e6, BLOCK_END / 1e6))
cat(sprintf("Second independent peaks  : %d\n",  nrow(other_peaks)))
cat("----------------------------------------------------------------------\n")
cat("Output files:\n")
cat("  results/plots/08_manhattan_confirmed_top_pair.pdf\n")
cat("  results/plots/09_chr16_zoom.pdf\n")
if (nrow(other_peaks) > 0)
  cat("  results/tables/second_peaks_candidates.csv\n")
cat("======================================================================\n")
cat("Done.\n")

