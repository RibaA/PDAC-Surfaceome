# PDAC Surfaceome Mapping

**Pancreatic ductal adenocarcinoma (PDAC)** is one of the most lethal and heterogeneous cancers, with limited therapeutic options and few reliable biomarkers. This project explores the **surfaceome** proteins on the cell surface as a potential source of clinically relevant diagnostic, prognostic, and therapeutic targets.

## Objective

To systematically map surface protein expression across PDAC tumors using TCGA data, stratify tumor subtypes, and identify candidate surface biomarkers linked to clinical outcomes.

## Data Sources

- **Transcriptomic & clinical data**: TCGA-PAAD (via GDC or UCSC Xena)
- **Surfaceome gene sets**: Curated from literature and surface protein databases

## Methods

- Preprocessing of RNA-seq and clinical metadata
- Filtering for known surface proteins
- Signatures (IO and TME) analysis
- Association with time-to-event outcomes (OS and PFI)
- Tumor clustering based on surfaceome profiles
- Visualization and biomarker prioritization

## Folder Structure

```md
.
├── data/                   # All data-related files
│   ├── raw/                # Unprocessed input data (e.g., TCGA downloads)
│   ├── proc/               # Cleaned and filtered data (e.g., expression matrices, clinical tables)
│   └── results/            # Data outputs from analysis (e.g., signature scores, survival stats, figures)
│
├── scripts/                # R scripts for processing and analysis
│   ├── runProcData.r               # Preprocessing TCGA data (expression & clinical)
│   ├── runSurfaceomeAnalysis.r     # Surfaceome gene filtering, expression, and clustering
│   ├── runDeconAnalysis.r          # General immune deconvolution analysis
│   └── runDeconAnalysis_tme.r      # Tumor microenvironment-specific deconvolution
│
└── README.md               # Project overview and usage