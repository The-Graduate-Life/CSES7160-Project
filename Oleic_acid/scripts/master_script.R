# ==============================================================================
# MASTER SCRIPT — Peanut Oleic Acid Genomic Analysis
#
# Genomic Dissection and Prediction of Oleic Acid Concentration in Peanut
#
# Author  : Fritzner Pierre
# Course  : CSES 7160 — Genetic Data Analysis
#
# ── Directory structure ───────────────────────────────────────────────────────
#   <root>/                          ← project root (.here lives here)
#     .here                          ← marks project root for here::here()
#     data/                          ← raw input files (not tracked by git)
#       values.csv
#       arahy.Tifrunner.gnm1.mrk.Axiom_Arachis_58K.gff3
#       aradu1_araip1.gnm1.div.Otyama_Kulkarni_2020.main.vcf
#     scripts/                       ← all R scripts and .qmd
#       master_script.R                ← this file
#       01_CSES7160_Fritzner_Analysis.R
#       02_GP_Deregressed_BLUPs.R
#       03_GP_NonAdditive_Model.R
#       04_CSES7160_GWAS_Comparison.R
#       05_GWAS_Confirmation.R
#       06_LD_Analysis.R
#       07_prediction_diagnostics.R
#       CSES7160_Fritzner.qmd
#     results/                       ← all outputs (auto-created by scripts)
#
# ── Execution order ───────────────────────────────────────────────────────────
#   Script 01 must run first — it generates all .rds and .csv files that the
#   downstream scripts depend on. Scripts 02–07 can be run independently after
#   Script 01 has completed successfully.
#
# ── How to run ────────────────────────────────────────────────────────────────
#   Option A — Source this file from RStudio (File > Open > master_script.R,
#              then click Source).
#   Option B — From the R console:
#              source("master_script.R")
#   Option C — From terminal (from the project root):
#              Rscript master_script.R
#
# ── Notes ─────────────────────────────────────────────────────────────────────
#   • Each script is sourced with chdir = FALSE so that here::here() continues
#     to resolve from the project root, not from the scripts/ directory.
#   • setwd(here::here()) inside each script ensures all relative paths
#     (data/, results/) resolve correctly regardless of where R was launched.
#   • Each script is wrapped in tryCatch so a failure in one does not abort
#     the rest of the pipeline. Check the log for error messages.
#   • To render the .qmd report, run the last block (quarto::quarto_render).
# ==============================================================================


# ==============================================================================
# 0 — SETUP
# ==============================================================================

if (!requireNamespace("here", quietly = TRUE))
  install.packages("here")
suppressPackageStartupMessages(library(here))

PROJECT_DIR <- here::here()
setwd(PROJECT_DIR)

cat("==============================================================\n")
cat(" MASTER SCRIPT — Peanut Oleic Acid Genomic Analysis\n")
cat("==============================================================\n")
cat("Project root :", PROJECT_DIR, "\n")
cat("Start time   :", format(Sys.time(), "%Y-%m-%d %H:%M:%S"), "\n\n")

# Verify .here file exists at root
if (!file.exists(file.path(PROJECT_DIR, ".here")))
  warning("'.here' file not found in project root. ",
          "Run `here::set_here()` from the project root to create it.")

# Verify scripts/ directory exists
SCRIPTS_DIR <- file.path(PROJECT_DIR, "scripts")
if (!dir.exists(SCRIPTS_DIR))
  stop("scripts/ directory not found at: ", SCRIPTS_DIR)

# Verify data/ directory and required files
DATA_DIR <- file.path(PROJECT_DIR, "data")
required_data <- c(
  "values.csv",
  "arahy.Tifrunner.gnm1.mrk.Axiom_Arachis_58K.gff3",
  "aradu1_araip1.gnm1.div.Otyama_Kulkarni_2020.main.vcf"
)
missing_data <- required_data[
  !file.exists(file.path(DATA_DIR, required_data))
]
if (length(missing_data) > 0) {
  cat("WARNING: The following data files are missing from data/:\n")
  for (f in missing_data) cat("  -", f, "\n")
  cat("Download them from PeanutBase / LIS DataStore before running.\n\n")
}


# ==============================================================================
# 1 — HELPER: run_script()
# ==============================================================================
# Wraps source() in tryCatch, logs timing, and resets wd after each script.

run_script <- function(script_name, description) {
  script_path <- file.path(SCRIPTS_DIR, script_name)
  
  if (!file.exists(script_path)) {
    cat(sprintf("  [SKIP] %s not found: %s\n\n", script_name, script_path))
    return(invisible(FALSE))
  }
  
  cat(sprintf("--------------------------------------------------------------\n"))
  cat(sprintf(" Running : %s\n", script_name))
  cat(sprintf(" Purpose : %s\n", description))
  cat(sprintf(" Start   : %s\n", format(Sys.time(), "%H:%M:%S")))
  cat(sprintf("--------------------------------------------------------------\n"))
  
  t0 <- proc.time()
  
  result <- tryCatch(
    withCallingHandlers(
      {
        # chdir = FALSE: keep wd at project root; script calls setwd(here::here())
        source(script_path, chdir = FALSE, echo = FALSE)
        TRUE
      },
      warning = function(w) {
        # Log warnings but let the script continue running
        cat(sprintf("\n  [WARNING] %s: %s\n", script_name, conditionMessage(w)))
        invokeRestart("muffleWarning")
      }
    ),
    error = function(e) {
      cat(sprintf("\n  [ERROR] %s failed:\n  %s\n\n",
                  script_name, conditionMessage(e)))
      FALSE
    }
  )
  
  elapsed <- (proc.time() - t0)[["elapsed"]]
  setwd(PROJECT_DIR)   # always reset wd after each script
  
  status <- if (isTRUE(result)) "OK" else "FAILED"
  cat(sprintf("\n  [%s] %s  (%.1f sec)\n\n", status, script_name, elapsed))
  invisible(isTRUE(result))
}


# ==============================================================================
# 2 — RUN PIPELINE
# ==============================================================================

pipeline_start <- proc.time()

# ── Script 01: Main analysis ──────────────────────────────────────────────────
# Phenotype QC → BLUPs → H² → Genotype QC → GRM → Genomic heritability →
# GWAS (FarmCPU) → Genomic prediction (GBLUP, 5-fold CV)
# Outputs: results/geno_num.rds, results/geno_map.rds, results/gwas_df.rds,
#          results/tables/*.csv, results/plots/01–05_*.pdf
ok1 <- run_script(
  "01_CSES7160_Fritzner_Analysis.R",
  "Main pipeline: BLUP extraction, QC, heritability, GWAS, GBLUP prediction"
)

if (!ok1) {
  cat("==============================================================\n")
  cat(" PIPELINE HALTED: Script 01 failed.\n")
  cat(" All downstream scripts depend on Script 01 outputs.\n")
  cat(" Fix the error above and re-run master_script.R.\n")
  cat("==============================================================\n")
  stop("Pipeline aborted after Script 01 failure.", call. = FALSE)
}

# ── Script 02: Deregressed BLUPs ─────────────────────────────────────────────
# Garrick (2009) deregression: removes double-shrinkage bias from raw BLUPs
# Outputs: results/gp/deregressed_blups.csv,
#          results/plots/12_gblup_deregressed_cv.pdf
run_script(
  "02_GP_Deregressed_BLUPs.R",
  "Garrick deregression of BLUPs + deregressed GBLUP cross-validation"
)

# ── Script 03: Non-additive model ─────────────────────────────────────────────
# Additive + dominance GBLUP (Vitezica 2013 G_D matrix via sommer::mmer)
# Tests whether dominance variance is non-trivial vs. additive model
# Outputs: results/tables/variance_components_adddom.csv,
#          results/gp/cv_additive_vs_dominance.csv,
#          results/plots/14_cv_additive_vs_dominance.pdf
run_script(
  "03_GP_NonAdditive_Model.R",
  "Additive + dominance GBLUP: dominance variance test and CV comparison"
)

# ── Script 04: GWAS model comparison ─────────────────────────────────────────
# FarmCPU vs MLM side-by-side: lambda_GC, hit overlap, -log10(p) correlation
# Outputs: results/gwas_comparison/tables/*.csv,
#          results/gwas_comparison/plots/*.pdf
run_script(
  "04_CSES7160_GWAS_Comparison.R",
  "FarmCPU vs MLM GWAS comparison: inflation, power, hit concordance"
)

# ── Script 05: GWAS confirmation ──────────────────────────────────────────────
# Confirms top GWAS pair (Chr 16 r²=1.0), zoomed Chr 16 Manhattan
# Outputs: results/plots/08_manhattan_confirmed_top_pair.pdf,
#          results/plots/09_chr16_zoom.pdf
run_script(
  "05_GWAS_Confirmation.R",
  "Manhattan confirmation plot and Chr 16 zoom for top GWAS signal"
)

# ── Script 06: LD analysis ────────────────────────────────────────────────────
# Pairwise r² in ±500 kb window around top hit, LD heatmap, decay bar plot
# Outputs: results/tables/ld_with_top_hit.csv,
#          results/plots/06_ld_heatmap_top_region.pdf,
#          results/plots/07_ld_decay_top_hit.pdf
run_script(
  "06_LD_Analysis.R",
  "LD analysis of top GWAS hit: pairwise r², heatmap, decay plot"
)


# ── Script 07: Prediction Diagnostic ────────────────────────────────────────────────────
# Prediction Diagnostics — Extended Accuracy Metrics (Independent Script)
# Outputs:
#   results/gp/prediction_diagnostics_raw.csv
#   results/gp/prediction_diagnostics_deregressed.csv
#   results/gp/prediction_diagnostics_combined.csv
run_script(
  "07_prediction_diagnostics.R",
  "Extended prediction accuracy diagnostics (Spearman, MAE, classification, group r)"
)

# ==============================================================================
# 3 — OPTIONAL: RENDER QMD REPORT
# ==============================================================================
# Uncomment to render the Quarto report after all scripts have run.
# Requires: install.packages("quarto") and Quarto CLI installed.
#
# qmd_path <- file.path(SCRIPTS_DIR, "CSES7160_Fritzner.qmd")
# if (file.exists(qmd_path)) {
#   if (!requireNamespace("quarto", quietly = TRUE))
#     install.packages("quarto")
#   cat("--------------------------------------------------------------\n")
#   cat(" Rendering Quarto report...\n")
#   cat("--------------------------------------------------------------\n")
#   tryCatch(
#     quarto::quarto_render(
#       input        = qmd_path,
#       output_file  = "CSES7160_Fritzner_Report.html",
#       output_dir   = PROJECT_DIR,
#       execute_dir  = PROJECT_DIR
#     ),
#     error = function(e) cat("[ERROR] Quarto render failed:", conditionMessage(e), "\n")
#   )
# } else {
#   cat("QMD file not found — skipping report render.\n")
# }


# ==============================================================================
# 4 — PIPELINE SUMMARY
# ==============================================================================

total_elapsed <- (proc.time() - pipeline_start)[["elapsed"]]

cat("==============================================================\n")
cat(" PIPELINE COMPLETE\n")
cat("==============================================================\n")
cat("Total time  :", round(total_elapsed / 60, 1), "minutes\n")
cat("End time    :", format(Sys.time(), "%Y-%m-%d %H:%M:%S"), "\n\n")
cat("Key outputs:\n")
cat("  results/tables/BLUPs_OleicAcid.csv\n")
cat("  results/tables/H2_broad_lmer.csv\n")
cat("  results/tables/heritability.csv\n")
cat("  results/tables/variance_components_adddom.csv\n")
cat("  results/tables/gwas_significant_hits_FarmCPU.csv\n")
cat("  results/tables/ld_with_top_hit.csv\n")
cat("  results/gp/deregressed_blups.csv\n")
cat("  results/gp/cv_additive_vs_dominance.csv\n")
cat("  results/gwas_comparison/tables/model_comparison_summary.csv\n")
cat("  results/plots/  (all figures)\n")
cat("  results/gp/prediction_diagnostics_raw.csv\n")
cat("  results/gp/prediction_diagnostics_deregressed.csv\n")
cat("  results/gp/prediction_diagnostics_combined.csv\n")
cat("==============================================================\n")
