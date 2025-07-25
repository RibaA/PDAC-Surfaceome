#######################################################
## load libraries
#######################################################
library(dplyr)
library(qs)
library(IOBR)
library(PredictioR)
library(NbClust)
library(ComplexHeatmap)
library(RColorBrewer)
library(colorRamp2)
library(MultiAssayExperiment)

dir_in <- 'data/proc'
dir_out <- 'data/results'

# --------------------------------------------------
# Aim 2: Stratify tumors by immune and molecular features
# --------------------------------------------------
## load MAE object data 
mae <- qread(file.path(dir_in, 'TCGA_PAAD_mae.qs'))

# Extract RNA-seq and RPPA assays
rna <- assay(mae[["RNAseq"]])
rppa <- assay(mae[["RPPA"]])
clin <- as.data.frame(colData(mae[["RNAseq"]]))
annot_rppa <- as.data.frame(rowData(mae[["RPPA"]]))
annot_rna <- as.data.frame(rowData(mae[["RNAseq"]]))
rownames(rna) <- annot_rna$hgnc_symbol

# Get common sample barcodes
rna_samples <- colnames(rna)
rppa_samples <- colnames(rppa)
clin_samples <- clin$sample

common_samples <- Reduce(intersect, list(rna_samples, rppa_samples, clin_samples)) # 116 samples

# Subset datasets
expr <- rna[, common_samples]
rppa <- rppa[, common_samples]
clin <- clin[clin$sample %in% common_samples, ]

clin <- clin[order(clin$sample), ]
expr <- expr[, order(colnames(expr))]

## remove low/zero-expressed genes
remove <- rem(expr, missing.perc = 0.5, const.int = 1)
if(length(remove) > 0){
  expr <- expr[-remove, ]
}
dim(expr) # ~ 25K genes

################################################################
################################################################
# Part 2: Run CIBERSORT to infer immune cell types
################################################################
################################################################
cell <- deconvo_tme(eset = expr, method = "cibersort", perm = 500, arrays = FALSE)

# Identifying TME patterns
tme <- tme_cluster(input = cell, features = colnames(cell)[2:23], 
                   id = "ID", scale = TRUE, method = "kmeans", 
                   min_nc =2, max.nc = 2)

colnames(tme) <- gsub(colnames(tme), pattern = "_CIBERSORT", replacement = "")

cibersort_tme <- list('cell' = cell,
                  'tme' = tme)

save(cibersort_tme, file = file.path(dir_out, 'cibersort', 'cibersort_tme.RData'))
################################################################
# Heatmap
################################################################
df <- tme[, -c(1,2)]
df <- t(df)
colnames(df) <- tme$ID
df <- t(scale(t(df)))

col = list( Histo = c( "Other" = "#67A9CF" , "PDAC" = "#EF8A62" ) ,
            Cluster = c("TME2" = "#4d004b",  "TME1" = "#4A4D00", "TME3" = "#004A4D"), 
            Sex = c( "F" = "#fee090" , "M" = "#4477AA" ),
            Stage = c("I" = '#525252', "II" = "#35978f" , "III" = "#bf812d" , "IV" = "#8073ac"))

# Create the heatmap annotation
ha <- HeatmapAnnotation(
  Histo = clin$histo,
  Cluster = tme$cluster,
  Sex = clin$sex,
  Stage = clin$stage,
  col = col,
  annotation_name_gp = gpar(fontsize = 8),
  simple_anno_size = unit(0.5, "cm")
)

# Combine the heatmap and the annotation
pdf(file=file.path(dir_out, "cibersort", "heatmap_cibersort.pdf"),
     width = 8, height = 8)

ht_data <- Heatmap(df, name = "expr",
                   top_annotation = ha,
                   show_row_names = TRUE, 
                   show_column_names = FALSE,
                   row_names_gp = gpar(fontsize = 8),
                   column_names_gp = gpar(fontsize = 8),
                   colorRamp2(c(-3, median(df), 11), 
                              c("#045a8d", "white", "#800026")))

draw(ht_data,
     column_title_gp = gpar(fontsize = 8, fontface = "bold"),
     #merge_legends = TRUE,
     heatmap_legend_side = "right",
     annotation_legend_side="bottom")

dev.off()

######################################################################
## boxplot
######################################################################
cols <- c('#2692a4', '#ffbe0b')
signature_name <- colnames(tme)[-c(1,2)]

for(i in 1:length(signature_name)){

df <- data.frame( y = as.numeric(tme[, signature_name[i]]),
                  group = tme$cluster)

pdf(file=file.path(dir_out, "cibersort/boxplot", paste(signature_name[i], "heatmap_cibersort.pdf", sep="")),
     width = 4, height = 4)

p <- ggplot(df, aes(x=group, y=y, fill = group)) + 
    geom_violin(trim=FALSE) +
    geom_boxplot(width=0.07, fill = "white")+
    xlab("")+
    ylab("")+
    ggtitle(signature_name[i]) +
    scale_fill_manual(values=c("#67A9CF", "#EF8A62"))+
    theme(axis.text.x=element_text(size=10,  face="bold"),
          axis.title=element_text(size=10,face="bold"),
          axis.text.y=element_text(size=10),
          panel.grid.major = element_blank(), 
          panel.grid.minor = element_blank(),
          panel.background = element_blank(),
          plot.background = element_blank(), 
          axis.line = element_line(colour = "black")) +
    guides(fill=FALSE)
  
  print(p)

dev.off()

}

