# Genomic Dissection and Prediction of Oleic Acid Concentration in Peanut

**Course:** CSES 7160 — Genetic Data Analysis  
**Author:** Fritzner Pierre  
**Institution:** Auburn University  
**Semester:** Spring 2026

---

## Overview

This project investigates the genetic architecture of oleic acid concentration in cultivated peanut (*Arachis hypogaea*) using high-density SNP markers from the Axiom_Arachis_58K genotyping array. Oleic acid is a key oil-quality trait in peanut: high-oleic varieties have superior shelf life and nutritional profile.

The analysis combines genome-wide association studies (GWAS), genomic prediction (GBLUP), heritability estimation, and linkage disequilibrium (LD) analysis to answer:

> **What genomic regions control variation in oleic acid, and how accurately can marker-based models predict oil quality across the U.S. peanut mini core collection?**


## Key Findings

| Finding | Result |
|---------|--------|
| Broad-sense heritability (H²) | 0.781 |
| SNP-based genomic heritability (h²_g) | 0.309 |
| Top GWAS signal | Chr 16 (B-genome), p = 7.05 × 10⁻⁵ |
| Perfectly linked SNP pair (r² = 1.0) | 2 SNPs at Chr 16, positions 2.01 Mb & 2.34 Mb |
| Secondary signal | Chr 6 (A-genome), p < 1 × 10⁻⁵ |
| GBLUP cross-validation accuracy (r) | 0.68–0.71 |
| Dominance variance | 0.46% of V_p (minimal) |

The large gap between H² and h²_g (0.47) suggests epistatic effects or polygenic architecture not fully captured by common SNPs, rather than classical dominance.


## Data Sources

| File | Description | Source |
|------|-------------|--------|
| `Oleic_acid/data/values.csv` | Oleic acid phenotypes (315 obs., 104–108 accessions) | [PeanutBase](https://arachispheno.peanutbase.org) |
| `Oleic_acid/data/arahy.Tifrunner.gnm1.mrk.Axiom_Arachis_58K.gff3` | SNP position map | Legume Information System (LIS) |
| `Oleic_acid/data/aradu1_araip1.gnm1.div.Otyama_Kulkarni_2020.main.vcf` | Genotype calls (58 MB VCF) | LIS / Otyama & Kulkarni et al. 2020 |


## Repository Structure

```
CSES7160-Project/
├── README.md
├── FOLLOWME.md                  # Beginner-friendly reproduction guide
└── Oleic_acid/
    ├── .here                    # Project root marker (for reproducibility)
    ├── data/                    # Raw input files (DO NOT modify)
    │   ├── values.csv
    │   ├── arahy.Tifrunner.gnm1.mrk.Axiom_Arachis_58K.gff3
    │   └── aradu1_araip1.gnm1.div.Otyama_Kulkarni_2020.main.vcf
    ├── scripts/
    │   ├── master_script.R              # Run this to execute the full pipeline
    │   ├── 01_CSES7160_Fritzner_Analysis.R   # Core analysis (REQUIRED first)
    │   ├── 02_GP_Deregressed_BLUPs.R         # Garrick deregression for GP
    │   ├── 03_GP_NonAdditive_Model.R          # Additive + dominance GBLUP
    │   ├── 04_CSES7160_GWAS_Comparison.R      # FarmCPU vs MLM comparison
    │   ├── 05_GWAS_Confirmation.R             # Manhattan plots & validation
    │   ├── 06_LD_Analysis.R                   # LD heatmaps & decay
    │   └── CSES7160_Fritzner.qmd              # Quarto HTML report
    └── results/
        ├── plots/               # 16 PDF visualizations (numbered 00–15)
        ├── tables/              # 16 CSV output tables
        ├── gp/                  # Genomic prediction results
        ├── gwas/                # Raw GAPIT GWAS outputs
        ├── gwas_comparison/     # FarmCPU vs MLM comparison outputs
        ├── geno_num.rds         # Processed genotype matrix (post-QC)
        ├── geno_map.rds         # SNP position map
        ├── gwas_df.rds          # GWAS results object
        └── session_info.txt     # R session/package versions
```


## Analysis Pipeline

The six numbered scripts form a sequential pipeline. Script 01 must complete first; the rest are independent of each other.

```
Script 01 (Required) → generates .rds files and core outputs
    ├── Script 02: Deregressed BLUPs genomic prediction
    ├── Script 03: Non-additive (additive + dominance) GBLUP
    ├── Script 04: GWAS model comparison (FarmCPU vs MLM)
    ├── Script 05: GWAS confirmation & zoomed chromosome plots
    └── Script 06: LD analysis around top GWAS hit
```

| Script | Purpose | Runtime (approx.) |
|--------|---------|------------------|
| 01 | Phenotype QC, BLUP extraction, genotype QC, PCA, GRM, heritability, FarmCPU GWAS, GBLUP CV | ~15–30 min |
| 02 | Garrick (2009) deregression; removes double-shrinkage bias | ~5–10 min |
| 03 | Dominance relationship matrix; additive vs. add+dom CV | ~10–20 min |
| 04 | MLM GWAS; side-by-side comparison with FarmCPU | ~10–20 min |
| 05 | Top-hit validation; Chr 16 zoomed Manhattan | ~2–5 min |
| 06 | LD r² heatmap and decay around Chr 16 top SNP | ~2–5 min |


## Methods Summary

### Phenotype Processing
- Raw multi-environment data fitted with `lme4::lmer` (random accession effect)
- Best Linear Unbiased Predictors (BLUPs) extracted as adjusted phenotypes

### Genotype QC
- VCF parsed with `vcfR`; converted to 0/1/2 numeric matrix
- Filters: MAF >= 5%, missingness < 10%

### Population Structure
- Genomic Relationship Matrix (GRM) constructed via VanRaden (2008) method
- Principal component analysis (PCA) on GRM

### Heritability
- Narrow-sense genomic h²_g from GBLUP (rrBLUP)
- Bootstrap 95% confidence intervals
- Broad-sense H² from lmer variance components

### GWAS
- **FarmCPU** (Liu et al., 2016): chosen for its power under small n (108 accessions) via pseudo-QTN fixed effects; controls inflation better than MLM
- Significance thresholds: genome-wide (p < 5 × 10⁻⁸), suggestive (p < 1 × 10⁻⁵)

### Genomic Prediction
- GBLUP via `rrBLUP::mixed.solve`
- 5-fold cross-validation, 10 replicates
- Deregressed BLUPs (Garrick 2009) to eliminate double-shrinkage bias
- Additive + dominance model via `sommer::mmer` (Vitezica et al., 2013)

### LD Analysis
- Pairwise r² in ±500 kb window around top GWAS hit
- LD heatmap via `LDheatmap`; decay plot via ggplot2


## Software Requirements

- **R >= 4.4** (developed with R 4.5.2)
- **Quarto** (for rendering the HTML report only)

### R Packages

```r
install.packages(c(
  "lme4", "rrBLUP", "sommer", "vcfR",
  "ggplot2", "dplyr", "tidyr", "patchwork",
  "ggrepel", "LDheatmap", "boot", "here",
  "knitr", "kableExtra"
))

# GAPIT must be installed from GitHub
install.packages("devtools")
devtools::install_github("jiabowang/GAPIT3")
```

---

## **AI Assistance Disclosure**

Portions of this analysis were developed with assistance from [Claude](https://claude.ai/) (Anthropic, claude.sonnet-4-6, 2026). AI assistance was used for the following tasks:

- Debugging R code and resolving package compatibility issues (rrBLUP, GAPIT3, sommer)
- Reviewing and refining statistical model implementations (GBLUP, deregressed BLUPs, dominance relationship matrix)
- Organizing the master pipeline script

All statistical decisions, biological interpretations, and analytical choices were made by the author. AI-generated code and text were reviewed, tested, and modified as needed before inclusion. The underlying data, methods, and conclusions are the author's own work.

---


## References

- Garrick, D. J., Taylor, J. F., & Fernando, R. L. (2009). Deregressing estimated breeding values and weighting information for genomic regression analyses. *Genetics Selection Evolution, 41(1)*, 55. [https://doi.org/10.1186/1297-9686-41-55](https://doi.org/10.1186/1297-9686-41-55).

- Liu, X., Huang, M., Fan, B., Buckler, E.S., & Z. Zhang. 2016. “Iterative Usage of Fixed and Random Effect Models for Powerful and Efficient Genome-Wide Association Studies.” *PLOS Genetics 12 (2)*: e1005767. [https://doi.org/10.1371/journal.pgen.1005767](https://doi.org/10.1371/journal.pgen.1005767).

- Otyama, P. I., Wilkey, A., Kulkarni, R., Assefa, T., Chu, Y., Clevenger, J., ... & Cannon, S. B. (2019). Evaluation of linkage disequilibrium, population structure, and genetic diversity in the US peanut mini core collection. *BMC genomics, 20(1)*, 481. [https://doi.org/10.1186/s12864-019-5824-9](https://doi.org/10.1186/s12864-019-5824-9).

- VanRaden PM (2008) Efficient methods to compute genomic predictions. *Journal of Dairy Science* 91:4414–4423. [https://doi.org/10.3168/jds.2007-0980](https://doi.org/10.3168/jds.2007-0980).

- Vitezica, Z. G., Varona, L., & Legarra, A. (2013). On the additive and dominant variance and covariance of individuals within the genomic selection scope. Genetics, 195(4), 1223-1230. [https://doi.org/10.1534/genetics.113.155176](https://doi.org/10.1534/genetics.113.155176).
