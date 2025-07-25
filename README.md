# PDAC Surfaceome Mapping

**Pancreatic ductal adenocarcinoma (PDAC)** is one of the deadliest cancers, with high heterogeneity and limited treatment options. This project focuses on the **surfaceome** proteins on the cell surface as a source of biomarkers and therapeutic targets. By analyzing TCGA transcriptomic, proteomics, and clinical data, we aim to identify surface proteins that distinguish tumor subtypes and correlate with clinical outcomes.

---

## Objectives

- Map surface protein expression profiles across PDAC tumors
- Stratify tumors based on surfaceome signatures
- Identify surface markers associated with survival (OS, DSS, DFI, PFI)
- Explore therapeutic and diagnostic implications of surfaceome heterogeneity

---

## Environment Setup

This project uses [**Pixi**](https://prefix.dev/docs/pixi/) for reproducible environment and dependency management.

### Requirements

- [Pixi](https://prefix.dev/docs/pixi/getting-started/installation/)
- Git

### Setup Instructions

```bash
# Clone the repo
git clone https://github.com/your-username/PDAC-Surfaceome.git
cd PDAC-Surfaceome

# Create and activate the environment
pixi install
pixi run R         # Start an R session inside the environment
```

---

## Folder Structure

```md
.
├── data/
│   ├── raw/                     # Raw input data (TCGA RNA-seq, RPPA, clinical)
│   ├── proc/                    # Processed and curated data
│   └── results/                 # Output from all analyses
│
├── scripts/
│   ├── runProcData.r                       # TCGA data preprocessing & clinical curation
│   ├── runSurfaceomeAnalysis.r            # Gene filtering and surfaceome clustering
│   ├── runSurfaceomeCorrelationAnalysis.r # Correlation between surfaceome and TME/cell types
│   ├── runDeconCibersortAnalysis.r        # CIBERSORT-based deconvolution (removes low-expressed genes)
│   ├── runDeconxCellAnalysis.r            # xCell-based deconvolution (removes low-expressed genes)
│
├── pixi.toml                   # Pixi environment spec
├── pixi.lock                   # Locked package versions
└── README.md                   # This file

```
