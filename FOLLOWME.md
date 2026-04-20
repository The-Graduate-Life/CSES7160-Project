# A Step-by-Step Guide to Successfully Reproduce this Pipeline Analysis

This guide walks you through reproducing the entire analysis from scratch, even if you have never used R or genomics tools before. Follow each step in order.

---

## What You Will Need

- A computer with internet access
- About 1 GB of free disk space
- 1–2 hours (most of which is waiting for R to finish)


> ***You might alraedy have `R` and `RStudio` installed. If so, start directly with step 3.***
## Step 1 — Install R

R is the programming language used for all analyses.

1. Go to https://cran.r-project.org
2. Click **Download R for Windows** (or Mac / Linux)
3. Click **base**, then download the installer
4. Run the installer and accept all defaults

> **Check it worked:** Open the app called **R** or **RGui** — you should see a `>` prompt.



## Step 2 — Install RStudio (recommended)

RStudio is a friendlier interface for R.

1. Go to https://posit.co/download/rstudio-desktop/
2. Download the free **RStudio Desktop** version for your OS
3. Install it with default settings

> From now on, open **RStudio** instead of R directly.


## Step 3 — Install Quarto (needed for the HTML report)

Quarto renders the final interactive report. Check if Quarto is Already Installed. `RStudio v2022.07` and above comes with Quarto built-in. Check your version via: *`Help → About RStudio`*
If you're on a recent version, Quarto is already there. Confirm in the RStudio Terminal:
```bash
quarto --version
```
**Install Quarto (if not already present)**
+ Install the Required R Packages
In the RStudio Console, run:
```r
# Core packages
install.packages("knitr")       # Rendering engine for R
install.packages("rmarkdown")   # Needed for compatibility
install.packages("quarto")      # Quarto R package (optional helper)
```


## Step 4 — Get the Project Files

### Option A — Clone with Git (if you have Git installed)

Open a terminal (or RStudio's **Terminal** tab) and run:

```bash
git clone https://github.com/The-Graduate-Life/CSES7160-Project.git
# if you have SSH set up on your GitHub, it will smoothen the process.
# On the GitHub repository page, click the green **Code** button and copy the link under the SSH. Use:
# git clone git@github.com:The-Graduate-Life/CSES7160-Project.git

# And then, navigate to the cloned repository:

cd CSES7160-Project
```

### Option B — Download as ZIP

1. On the GitHub repository page, click the green **Code** button
2. Click **Download ZIP**
3. Unzip the folder somewhere easy to find (e.g., your Desktop)



## Step 5 — Open the Project in RStudio

1. In RStudio, go to **File > Open File...**
2. Navigate into `CSES7160-Project/Oleic_acid/scripts/`
3. Open **`master_script.R`**

You should now see the master script in the editor pane.



## Step 6 — Install Required R Packages

You only need to do this once. In the RStudio **Console** pane (bottom-left), paste and run each block. You can ignore this step since the scripts will automatically install the required packages.

```r
# Install packages from CRAN
install.packages(c(
  "lme4",
  "rrBLUP",
  "sommer",
  "vcfR",
  "ggplot2",
  "dplyr",
  "tidyr",
  "patchwork",
  "ggrepel",
  "LDheatmap",
  "boot",
  "here",
  "knitr",
  "kableExtra"
))
```

```r
# Install GAPIT from GitHub (requires devtools)
install.packages("devtools")
devtools::install_github("jiabowang/GAPIT3")
```

When asked "Do you want to install from sources?" type `n` and press Enter.  
Installation may take 5–15 minutes. You will see many lines of output — this is normal.

> **Troubleshooting:** If a package fails, try installing it individually:
> ```r
> install.packages("package_name")
> ```



## Step 7 — Verify Your Data Files

Make sure these three files exist inside `Oleic_acid/data/`:

| File | Size |
|------|------|
| `values.csv` | ~50 KB |
| `arahy.Tifrunner.gnm1.mrk.Axiom_Arachis_58K.gff3` | ~10 KB |
| `aradu1_araip1.gnm1.div.Otyama_Kulkarni_2020.main.vcf` | ~58 MB |

In RStudio's **Files** panel (bottom-right), navigate to `Oleic_acid/data/` and confirm all three files are present.

If the VCF file is missing, download it from the Legume Information System (LIS) DataStore — see [README](README.md) for the citation.



## Step 8 — Run the Full Pipeline

How to run the pipeline:
```
Option A:
Source this file from RStudio `(File > Open > master_script.R`, then click Source).

Option B:
From the R console: `source("master_script.R")`

Option C:
From terminal (from the project root): `Rscript master_script.R`
```

The master script will run all seven analysis scripts in order. You will see progress messages printed in the Console.

**Expected total runtime: 45–90 minutes** depending on your computer.

> The scripts use `tryCatch`, so if one script fails the others will still attempt to run. Any errors are printed in the Console with the script name.



## Step 9 — Monitor Progress

While the pipeline runs, the Console will print messages like:

```
========================================
Running Script 01: Main Analysis
========================================
[1] "Loading phenotype data..."
[1] "Fitting mixed model for BLUPs..."
...
Script 01 completed in 18.4 minutes.
```

Do not close RStudio while it is running. The pipeline is finished when you see:

```
========================================
All scripts completed.
========================================
```



## Step 10 — Inspect the Outputs

All results are saved to `Oleic_acid/results/`.

### Plots (PDF files)

Navigate to `Oleic_acid/results/plots/` in the Files panel. Key plots:

| File | What it shows |
|------|---------------|
| `01_phenotype_eda.pdf` | Distribution of oleic acid values |
| `02_pca.pdf` | Population structure (PC1 vs PC2) |
| `03_grm_heatmap.pdf` | Genomic relationship between accessions |
| `04_gwas_manhattan_qq.pdf` | GWAS results — Manhattan & QQ plots |
| `05_gblup_cv_and_fit.pdf` | Genomic prediction accuracy |
| `06_ld_heatmap_top_region.pdf` | Linkage disequilibrium around top SNP |
| `08_manhattan_confirmed_top_pair.pdf` | Confirmed top GWAS signal on Chr 16 |

Double-click any PDF to open it.

### Tables (CSV files)

Navigate to `Oleic_acid/results/tables/`. Key tables:

| File | Content |
|------|---------|
| `BLUPs_OleicAcid.csv` | Adjusted means (BLUPs) per accession |
| `heritability.csv` | h² estimate with 95% confidence interval |
| `gwas_all_results.csv` | All SNP p-values and effect sizes |
| `gwas_significant_hits.csv` | SNPs passing the suggestive threshold |
| `gblup_cv_results.csv` | Cross-validation accuracy per fold |

Open CSV files in Excel or any spreadsheet program.



## Step 11 — Render the HTML Report (optional)

The Quarto report combines all analyses into a single interactive HTML document.

1. Open `Oleic_acid/scripts/CSES7160_Fritzner.qmd` in RStudio
2. Click the **Render** button at the top of the editor

The report will open in your browser automatically. It includes:
- Collapsible code blocks
- Embedded plots
- Summary tables with heritability, GWAS hits, and prediction accuracy



## Run Individual Scripts (optional)

If you want to re-run only one analysis after making changes, open and source any script individually:

```r
# Example: re-run only the LD analysis
source("Oleic_acid/scripts/06_LD_Analysis.R")
```

**Note:** Scripts 02–06 depend on intermediate files created by Script 01 (`.rds` files in `results/`). Script 01 must have run successfully at least once before running any of the others.



## Common Problems & Fixes

### "there is no package called 'X'"

The package was not installed. Run `install.packages("X")` in the Console.

### "could not find function 'GAPIT'"

GAPIT was not installed from GitHub. Re-run:

```r
devtools::install_github("jiabowang/GAPIT3")
```

### Script 01 fails with "file not found"

The data files are missing or in the wrong folder. Check Step 7.

### Script runs but produces no output

Check that your working directory is set to the project root. In the Console:

```r
library(here)
here::here()  # Should print the path ending in "Oleic_acid"
```

If it does not, run:

```r
setwd("path/to/CSES7160-Project/Oleic_acid")
```

### "package 'sommer' is not available"

Try installing from a specific CRAN mirror:

```r
install.packages("sommer", repos = "https://cloud.r-project.org")
```

### VCF file is 0 bytes or missing

Re-download the VCF from LIS DataStore (see README.md for the Otyama & Kulkarni 2020 citation). The file should be ~58 MB.



## What Each Script Does (Plain Language)

| Script | What it does |
|--------|-------------|
| **01** | Loads the raw data, calculates adjusted trait values (BLUPs), filters low-quality SNPs, builds a genomic relationship matrix, estimates heritability, runs the GWAS, and tests how well the model predicts unseen data |
| **02** | Re-runs the genomic prediction using a corrected version of the trait values that removes statistical bias (Garrick deregression) |
| **03** | Tests whether dominance (non-additive) gene action improves prediction beyond the standard additive model |
| **04** | Runs a second GWAS method (MLM) and compares it with FarmCPU to justify the method choice |
| **05** | Zooms into the top-associated chromosome region and confirms the two most significant SNPs are in perfect linkage disequilibrium |
| **06** | Measures how correlated neighboring SNPs are around the top GWAS hit (LD analysis) |



## Glossary for Beginners

| Term | Plain-language definition |
|------|--------------------------|
| **SNP** | A single DNA position where individuals differ |
| **BLUP** | A statistically adjusted trait value that accounts for unequal replication |
| **GWAS** | A scan of the entire genome to find SNPs associated with a trait |
| **GBLUP** | A model that uses all SNPs simultaneously to predict trait values |
| **Heritability (h²)** | The proportion of trait variation explained by genetics |
| **Cross-validation** | A technique to estimate how well a model predicts new, unseen individuals |
| **LD (Linkage Disequilibrium)** | Non-random association between nearby SNPs on the same chromosome |
| **VCF** | A standard file format storing genotype calls for many individuals |
| **GRM** | A matrix summarizing how genetically similar each pair of individuals is |
| **FarmCPU** | A GWAS method designed for small samples; avoids false positives by accounting for relatedness |
| **Manhattan plot** | A plot showing SNP association strength across all chromosomes |

---

## How to cite


> Pierre, F. (2026). *Genomic dissection and prediction of oleic acid concentration in peanut using high-density SNP markers*. GitHub. https://github.com/The-Graduate-Life/CSES7160-Project

**Data sources you used (for reference):**
> Pandey, M. K., Agarwal, G., Kale, S. M., Clevenger, J., Nayak, S. N., Sriswathi, M., & Varshney, R. K. (2017). Development and evaluation of a high density genotyping 'Axiom_Arachis' array with 58K SNPs for accelerating genetics and breeding in groundnut. *Scientific Reports*, *7*, 40577.[https://doi.org/10.1038/sep40577](https://doi.org/10.1038/srep40577).

> Otyama, P. I., & Kulkarni, R. (2020). *Arachis hypogaea diversity: Otyama & Kulkarni 2020* [Data set]. Legume Information System DataStore.[https://data.legumeinfo.org/Arachis/hypogaea/diversity/aradu1_araip1.gnm1.div.Otyama_Kulkarni_2020/](https://data.legumeinfo.org/Arachis/hypogaea/diversity/aradu1_araip1.gnm1.div.Otyama_Kulkarni_2020/).
