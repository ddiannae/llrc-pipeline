# Loss of long-range co-expression is a common trait in cancer

## RNA-Seq pipeline for TCGA and USCS Xena datasets

This pipeline is designed to generate co-expression matrices, gene expression - copy number correlation matrices and  differential expression, for normal and cancer tissues using data from the TCGA dataset (downloaded from GDC) and the USCS Xena dataset [[1]](#1), which integrates TCGA and GTEx data.

## Pipeline Overview
The pipeline consists of the following main steps:

**1. Data Acquisition:**
  - **From Xena**: Downloads counts, sample information, and annotations directly from the Xena-Toil S3 bucket. Tissues: brain, colon, esophagus, liver,  ovary, pancreas, prostate,  skin, and testis.
    - Rules file: `xena.smk`
  - **From GDC**: Queries GDC to create a manifest with STAR Counts for gene expression and ASCAT2 files for gene level copy number. It and uses the `gdc-client` tool to download files. Tissues: bladder, breast, kidney, lung, thyroid, and uterus.
    - Rules file: `gdc.smk`
    - Script: `queryGDC.py` and `queryGDCascat.py`

**2. Raw Matrix Generation:** Integrates raw counts with their respective annotations for each tissue to produce a raw matrix.

  - Rules file: `raw.smk`
  - Scripts: `addAnnotations.R`, `getRawMatrix.R`, `getAscatMatrix.R`

**3. Quality Control**: Performs the following QC steps:
  - Generates NOISeq plots.
  - Filters genes with low expression (mean raw counts < 10).
  - Produces PCA and density plots (expression per sample).
  - Removes samples with mean expression greater or lower than 2 standard deviations from the total mean.
  - Rules: `qc.smk`
  - Scripts: `NOISeqPlots.R`, `filterLowExpression.R`, `PCA.R`, `densityPlot.R` and `filterOutliers.R`

**4. Normalization**: Runs ARSyN [[2]](#2) for batch correction and gets matrices from DESeq2 and TPM normalization.
    - Rules: `normalization.smk`
    - Scripts: `runArsyn.R`, `tpmNormalization.R`, `deseq2Normalization.R` and `getNormMatrix.R`

**5. Differential gene expression**: Runs  differential gene expression analysis with DESeq2.
  - Rules: `deg.smk`
  - Scripts: `deg.R`

**6. Gene Copy Number analysis:** Integrates ASCAT files and outputs a gene copy number matrix per tissue and sample type. 
  - Rules: `ascat.smk`
  - Scripts: `getAscatHeatmap.R`, `getAscatMatrix.R`

**7. Correlation matrices**: Runs ARACNE [[3]](#3) to get co-expression matrices for each tissue and sample type.
  - Rules: `correlation.smk`
  - Scripts: `aracne_matrix.py`, `cnv_aracne_matrix.py`
   
## Repository Structure
- `bin`:  GDC Data Transfer Tool command-line utility for Linux.
- `config/config.yaml` — Example configuration file for Snakemake pipeline.
- `input`: Contains:
  - GENECODE annotation files originally used for TCGA data (v36 GTF, as described [here](https://gdc.cancer.gov/about-data/gdc-data-processing/gdc-reference-files)). 
  - GENECODE v44 (July, 2023) for gene filtering.
  - BioMart Ensembl Genes 108 to add GC content for QC purposes. 
- `workflow`: Snakemake rules and scripts used in the workflow.

<a id="1">[1]</a> :Goldman, M.J., Craft, B., Hastie, M. et al. Visualizing and interpreting cancer genomics data via the Xena platform. Nat Biotechnol (2020). https://doi.org/10.1038/s41587-020-0546-8

<a id="2">[2]</a> : Maria j. Nueda, Alberto Ferrer, Ana Conesa, ARSyN: a method for the identification and removal of systematic noise in multifactorial time course microarray experiments, Biostatistics, Volume 13, Issue 3, July 2012, Pages 553–566, https://doi.org/10.1093/biostatistics/kxr042

<a id="3">[3]</a> : Margolin, A.A., Nemenman, I., Basso, K. et al. ARACNE: An Algorithm for the Reconstruction of Gene Regulatory Networks in a Mammalian Cellular Context. BMC Bioinformatics 7 (Suppl 1), S7 (2006). https://doi.org/10.1186/1471-2105-7-S1-S7
