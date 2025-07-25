#######################################################
## load libraries
#######################################################
library(dplyr)
library(data.table)
library(readxl)
library(PredictioR)
library(qs)
library(ComplexHeatmap)
library(colorRamp2)
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
rppa_surface <- rppa[rownames(rppa) %in% annot_surfaceome$updateProtein, ]  # 20 proteins --> removing NAs
rppa_surface <- na.omit(rppa_surface)
# --------------------------------------------------
# Association of RPPA with cell types (using xCell)
# --------------------------------------------------
# load xCell outputs
load(file.path(dir_out, 'xcell', 'xcell_tme.RData'))
cell <- as.data.frame(xcell_tme$cell)
rownames(cell) <- cell$ID
cell <- as.data.frame(cell[, 2:65])
cell <- cell[order(rownames(cell)), ]
rppa_surface <- rppa_surface[, order(colnames(rppa_surface))]

cor_cell_protein <- lapply(1:ncol(cell), function(k){

  res_per_cell <- lapply(1:nrow(rppa_surface), function(i){

    fit  <- cor.test(as.numeric(cell[, k]), as.numeric(rppa_surface[i, ]))
    data.frame(cellType = substr(colnames(cell)[k], 1, nchar(colnames(cell)[k])-6),
               protein = rownames(rppa_surface)[i],
               gene_name = annot_surfaceome[rownames(annot_surfaceome) == rownames(rppa_surface)[i], 'gene_name'], 
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


# -------------------------------------------------------
# Association of RPPA with cell types (using CIBERSORT)
# -------------------------------------------------------
# load cibersort outputs
load(file.path(dir_out, 'cibersort', 'cibersort_tme.RData'))
cell <- as.data.frame(cibersort_tme$cell)
rownames(cell) <- cell$ID
cell <- as.data.frame(cell[, 2:23])
cell <- cell[order(rownames(cell)), ]
rppa_surface <- rppa_surface[, order(colnames(rppa_surface))]

cor_cell_protein <- lapply(1:ncol(cell), function(k){

  res_per_cell <- lapply(1:nrow(rppa_surface), function(i){

    fit  <- cor.test(as.numeric(cell[, k]), as.numeric(rppa_surface[i, ]))
    data.frame(cellType = substr(colnames(cell)[k], 1, nchar(colnames(cell)[k])-10),
               protein = rownames(rppa_surface)[i],
               gene_name = annot_surfaceome[rownames(annot_surfaceome) == rownames(rppa_surface)[i], 'gene_name'], 
               r = fit$estimate,
               pval= fit$p.value)   

  })

     res <- do.call(rbind, res_per_cell)
     res <- res[!is.na(res$r), ]
     res$FDR <- p.adjust(res$pval, method = "BH")
     res
})

cor_cibersort <- do.call(rbind, cor_cell_protein)
# cor_cibersort [cor_cibersort$pval < 0.05, ]

write.csv(cor_cibersort , file = file.path(dir_out, 'cibersort', 'surface_rppa_cibersort_correlation.csv'), row.names = FALSE)

########################################################
# Heatmap 22 surface proteins
########################################################
df <- na.omit(rppa_surface)
df <- t(scale(t(df)))
df <- df[, order(colnames(df))]
clin <- clin[order(rownames(clin)), ]

col = list( Histo = c( "Other" = "#67A9CF" , "PDAC" = "#EF8A62" ) ,
            Sex = c( "F" = "#fee090" , "M" = "#4477AA" ),
            Stage = c("I" = '#525252', "II" = "#35978f" , "III" = "#bf812d" , "IV" = "#8073ac"))

# Create the heatmap annotation
ha <- HeatmapAnnotation(
  Histo = clin$histo,
  Sex = clin$sex,
  Stage = clin$stage,
  col = col,
  annotation_name_gp = gpar(fontsize = 8),
  simple_anno_size = unit(0.5, "cm")
)

# Combine the heatmap and the annotation
pdf(file=file.path(dir_out, "heatmap_surface_rppa.pdf"),
     width = 8, height = 6)

ht_data <- Heatmap(df, name = "expr",
                   top_annotation = ha,
                   show_row_names = TRUE, 
                   show_column_names = FALSE,
                   row_names_gp = gpar(fontsize = 8),
                   column_names_gp = gpar(fontsize = 8),
                   colorRamp2(c(-5, median(df), 9), 
                              c("#045a8d", "white", "#800026")))

draw(ht_data,
     column_title_gp = gpar(fontsize = 8, fontface = "bold"),
     #merge_legends = TRUE,
     heatmap_legend_side = "right",
     annotation_legend_side="bottom")

dev.off()
