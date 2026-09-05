# -GEO-RNAseq-Pipeline

# Comprehensive RNA-Seq Differential Expression & Biomarker Pipeline 🧬📊

## 📌 Project Overview
This repository contains a highly robust and automated R pipeline for performing Differential Gene Expression (DGE) analysis on varied public GEO datasets. Unlike standard pipelines tailored only for raw counts (e.g., standard DESeq2 workflows), this script is customized to handle pre-normalized data (logCPM, FPKM), varied file structures, and complex metadata typical in toxicology and exposure studies.

## 🚀 Features
- **Adaptive Data Parsing:** Automatically fetches metadata via `GEOquery` and extracts expression matrices from complex supplementary files (CSV, TXT, Excel, Cuffdiff formats).
- **Limma-based DEA:** Utilizes `Limma` to accurately model pre-normalized continuous data where raw counts are unavailable.
- **Dynamic Group Assignment:** Intelligently identifies Control/Vehicle and Treatment groups from sample characteristics without hardcoding.
- **Automated Visualization:** Generates publication-ready **Volcano Plots** and **Heatmaps** (via `ggplot2` and `pheatmap`).
- **Biomarker Evaluation:** Automatically calculates Area Under the Curve (AUC) and plots **ROC curves** for the top DEGs using the `pROC` package.

## 📂 Datasets Successfully Analyzed
This pipeline was successfully tested and validated on the following GEO datasets:
- GSE158266, GSE107839, GSE119785, GSE141134, GSE21642

## 🛠️ Prerequisites & Libraries
Make sure you have R (>= 4.0) installed along with the following packages:
```R
install.packages(c("tidyverse", "pheatmap", "pROC"))
BiocManager::install(c("GEOquery", "limma", "Biobase"))
