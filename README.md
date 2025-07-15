# 🧬 PDAC Surfaceome Mapping

**Pancreatic ductal adenocarcinoma (PDAC)** is one of the deadliest cancers, with high heterogeneity and limited treatment options. This project focuses on the **surfaceome**—proteins on the cell surface—as a source of biomarkers and therapeutic targets. By analyzing TCGA transcriptomic and clinical data, we aim to identify surface proteins that distinguish tumor subtypes and correlate with clinical outcomes.

---

## 🎯 Objectives

- Map surface protein expression profiles across PDAC tumors
- Stratify tumors based on surfaceome signatures
- Identify surface markers associated with survival (OS, DSS, DFI, PFI)
- Explore therapeutic and diagnostic implications of surfaceome heterogeneity

---

## 📦 Environment Setup

This project uses [**Pixi**](https://prefix.dev/docs/pixi/) for reproducible environment and dependency management.

### Requirements

- [Pixi](https://prefix.dev/docs/pixi/getting-started/installation/)
- Git

### Setup Instructions

```bash
# Clone the repo
git clone https://github.com/your-username/pdac-surfaceome.git
cd pdac-surfaceome

# Create and activate the environment
pixi install
pixi run R         # Start an R session inside the environment
```

---

## 📁 Folder Structure

```md
.
├── data/                      # All data files
│   ├── raw/                   # Unprocessed TCGA data (clinical, RNA-seq)
│   ├── proc/                  # Cleaned and formatted data for analysis
│   └── results/               # Analysis output (signature scores, survival stats, figures)
│
├── scripts/                   # R scripts for pipeline execution
│   ├── runProcData.r                  # Preprocessing TCGA data
│   ├── runSurfaceomeAnalysis.r        # Surfaceome filtering & clustering
│   ├── runDeconAnalysis.r             # General deconvolution
│   └── runDeconAnalysis_tme.r         # TME-specific deconvolution
│
├── pixi.toml                  # Pixi environment configuration
├── pixi.lock                  # Locked dependencies
└── README.md                  # Project overview and setup guide
```
