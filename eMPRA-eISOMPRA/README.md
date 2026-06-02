# eMPRA-eISOMPRA QC and Analysis

## R Markdown notebooks for the quality control, organization, and downstream analysis of an in vivo / in vitro massively parallel reporter assay (MPRA) experiment. 
### The pipeline covers two assay formats:

ISOMPRA — a splicing-aware assay, with pre- and post-splice RNA fractions.

eMPRA — a barcode-level expression MPRA.

## Repository contents
### FilePurposeLibrary_QC_and_Organizing.Rmd:

Reads the raw per-sample alignment count tables, runs alignment and replicate-correlation QC, and reshapes the data into the long "UTR" format used downstream. 
Writes the formatted .rds/.csv files.

### MPRA_QC_and_Downstream_Analysis.Rmd:

End-to-end notebook. Part 1 reproduces the QC/reformatting step (the 5/23/25 re-analysis) and writes the long-format .rds files; Part 2 reads them back in, normalizes and computes expression, and runs the full downstream analysis (enhancer/repressor and allelic-effect models, comparison plots, Venn diagrams, GC-content checks, burden testing).

Library_QC_and_Organizing.Rmd is the standalone QC step. 
MPRA_QC_and_Downstream_Analysis.Rmd is self-contained: because Part 1 writes the exact files Part 2 reads, knitting it top to bottom reproduces the full pipeline.

How the pieces fit together:
raw per-sample counts ──▶ QC + reshape ──▶ *.utr.format.*.rds ──▶ downstream analysis
   (ms*.txt,                 (Part 1 /          (long format)         (Part 2)
    *_foundBCs.txt,           Library QC)
    trevino.names.csv, ...)

## Requirements:

- R (≥ 4.0 recommended) and a LaTeX distribution if you want to knit the Library QC notebook to PDF (e.g. tinytex::install_tinytex()).
- mpra_tools.R — an external helper script that provides aggregate.counts(), normalize.counts.edgeR(), compute.expression() and related functions used in Part 2. It is sourced at the top of the merged notebook and is not yet included in this repository (see Notes below).
- R packages used across the notebooks:

install.packages(c(
    "data.table", "tidyr", "dplyr", "readr", "stringr", "reshape2", "mgsub",
    "rlang", "broom", "psych", "Hmisc", "corrplot", "ggplot2", "ggpubr",
    "ggpmisc", "ggforce", "gghalves", "fplot", "umap", "RColorBrewer",
    "VennDiagram", "gridExtra"
  ))

  Bioconductor packages
  if (!requireNamespace("BiocManager", quietly = TRUE)) install.packages("BiocManager")
  BiocManager::install(c("Biostrings", "limma"))
  
Some downstream models also rely on a mixed-model package (e.g. lme4 / lmerTest for lmer()); install these if they are not already pulled in by mpra_tools.R.

## Expected input files

These are read from the working directory and are not included here:

- ms*.txt :: ISOMPRA per-sample count tables (ID + count)
- *_foundBCs.txt :: eMPRA found-barcode lists, one per sample
- trevino.names.csv :: Barcode ↔ rsID map for the eMPRA library
- nreads.txt  :: Total read counts per ISO-MPRA sample
- nreads.empra.txt :: Total read counts per eMPRA - sample
- corr.values.trevino*.csv :: Pre-computed correlation summaries for plotting
- meta.file3.trev.rds :: Library metadata (barcode_allele → ref/alt)

## Outputs

- reads.table.trevino.isompra.csv, reads.table.trevino.empra.csv — alignment QC tables
- invivo.ps1.utr.format.*.rds / .csv — ISOMPRA long format
- empratrev.utr.format.*.rds / .csv — eMPRA long format
- corr.*.csv — replicate correlation matrices
- Various .png figures (QC plots, correlograms, comparison plots, Venn diagrams)

## Running
From an R session (or RStudio) with the input files and mpra_tools.R in the working directory:

rrmarkdown::render("Library_QC_and_Organizing.Rmd")
rmarkdown::render("MPRA_QC_and_Downstream_Analysis.Rmd")

## Notes

- Absolute paths: a few chunks read meta.file3.trev.rds from a hard-coded Windows path (C:/Users/dinse/...). Update these to a relative path before running on another machine.
- External dependency: the notebooks assume mpra_tools.R and the input data files above are present. Add them to the repo (or provide a download link / small example dataset) if you want the analysis to be reproducible by others.
- The merged notebook was assembled from two previously separate notebooks ("Re_Analysis" and "Initial_Analysis"), with duplicated setup, library, and helper-function content removed. The analysis logic is unchanged.

## License
No license has been specified yet. Until one is added, the default is that the
work is copyrighted and others have no explicit permission to reuse it. Consider
adding a LICENSE file (e.g. MIT for permissive reuse) if you want to allow
others to use or build on this code.
