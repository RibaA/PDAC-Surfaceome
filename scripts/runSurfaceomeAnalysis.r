#######################################################
## load libraries
#######################################################
library(dplyr)
library(data.table)
library(readxl)
library(qs)
library(MultiAssayExperiment)

dir_in <- 'data/raw'
dir_out <- 'data/proc'

###############################################################
## load curated TCGA PAAD data (RNA, RPPA, and clinical)
##############################################################
# Load MAE object
mae <- qread(file.path(dir_out, 'TCGA_PAAD_mae.qs'))

# Extract RNA-seq and RPPA assays
rna <- assay(mae[["RNAseq"]])
rppa <- assay(mae[["RPPA"]])
clin <- as.data.frame(colData(mae[["RNAseq"]]))
annot_rppa <- as.data.frame(rowData(mae[["RPPA"]]))
annot_rna <- as.data.frame(rowData(mae[["RNAseq"]]))
rownames(rna) <- annot_rna$hgnc_symbol

# Transfer protein's name to gene name
#sub_annot <- annot[grepl("///", annot$gene_name), ]

# Get common sample barcodes
rna_samples <- colnames(rna)
rppa_samples <- colnames(rppa)
clin_samples <- clin$sample

common_samples <- Reduce(intersect, list(rna_samples, rppa_samples, clin_samples)) # 116 samples

# Subset datasets
rna <- rna[, common_samples]
rppa <- rppa[, common_samples]
clin <- clin[clin$sample %in% common_samples, ]

# --------------------------------------------------
# Aim 1: Filter RNA and RPPA to surfaceome genes
# --------------------------------------------------
surfaceome_dat <- read_excel(file.path(dir_in, 'table_S3_surfaceome.xlsx'), sheet = "SurfaceomeMasterTable")
surfaceome_dat <- surfaceome_dat[surfaceome_dat$'Surfaceome Label' == 'surface', ]
rna_surface <- rna[rownames(rna) %in% surfaceome_dat$'UniProt gene', ] # 18,246 genes in common --> 2674 survace genes 

# Save filtered matrices
dat <- list('rna' = rna_surface,
            'rppa' = rppa,
            'annot_rna' = annot_rna,
            'annot_rppa' = annot_rppa,
            'clin' = clin)

qsave(dat, file = file.path(dir_out, 'rna_rppa_surface_filtered.qs'))
