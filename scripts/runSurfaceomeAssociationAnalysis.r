#######################################################
## load libraries
#######################################################
library(dplyr)
library(data.table)
library(readxl)
library(PredictioR)
library(qs)
library(MultiAssayExperiment)

dir_in <- 'data/proc'
dir_out <- 'data/results'

###############################################################
## load curated TCGA PAAD data (RNA, RPPA, and clinical)
##############################################################
# Load MAE object
mae <- qread(file.path(dir_in, 'TCGA_PAAD_mae.qs'))

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
# Aim 1: Filter RPPA to surfaceome genes
# --------------------------------------------------
surfaceome_dat <- read_excel(file.path('data/raw', 'table_S3_surfaceome.xlsx'), sheet = "SurfaceomeMasterTable")
surfaceome_dat <- surfaceome_dat[surfaceome_dat$'Surfaceome Label' == 'surface', ]
annot_surfaceome <- annot_rppa[annot_rppa$gene_name %in% surfaceome_dat$'UniProt gene', ]
#rna_surface <- rna[rownames(rna) %in% surfaceome_dat$'UniProt gene', ] 
rppa_surface <- rppa[rownames(rppa) %in% annot_surfaceome$updateProtein, ]  # 22 proteins

# --------------------------------------------------
# Association of RPPA with cell types
# --------------------------------------------------
# load xCell outputs
load(file.path(dir_out, 'xcell', 'xcell_tme.RData'))
cell <- as.data.frame(xcell_tme$cell[, 2:65])

cor_cell_protein <- lapply(1:ncol(cell), function(k){

  res_per_cell <- lapply(1:nrow(rppa_surface), function(i){

    fit  <- cor.test(as.numeric(cell[, k]), as.numeric(rppa_surface[i, ]))
    data.frame(cellType = substr(colnames(cell)[k], 1, nchar(colnames(cell)[k])-6),
               protein = rownames(rppa_surface)[i],
               r = fit$estimate,
               pval= fit$p.value)   

  })

     res <- do.call(rbind, res_per_cell)
     res <- res[!is.na(res$r), ]
     res$FDR <- p.adjust(res$pval, method = "BH")
     res
})

cor_xcell <- do.call(rbind, cor_cell_protein)
# cor_xcell [cor_xcell$pval < 0.05, ]

write.csv(cor_xcell , file = file.path(dir_out, 'xcell', 'surface_rppa_xcell_correlation.csv'), row.names = FALSE)
