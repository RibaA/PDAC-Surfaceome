#######################################################
## load libraries
#######################################################
# devtools::install_github('dviraran/xCell')
library(dplyr)
library(qs)
library(IOBR)
library(NbClust)
library(survival)
library(survminer)
library(ComplexHeatmap)
library(PredictioR)
library(RColorBrewer)
library(colorRamp2)
library(MultiAssayExperiment)

dir_in <- 'data/proc'
dir_out <- 'data/results'

# --------------------------------------------------
# Aim 2: Stratify tumors by immune and molecular features
# --------------------------------------------------
## load data (from Aim 1)
dat <- qread(file.path(dir_in, 'rna_rppa_surface_filtered.qs'))
expr <- dat$rna
clin <- dat$clin

# Overall survival
clin$os_time <- round(clin$OS.time / 30.44)
clin$os_status <- clin$OS

# Progression-Free Interval
clin$pfi_time <- round(clin$PFI.time / 30.44)
clin$pfi_status <- clin$PFI

# stage
table(clin$pathologic_stage)
clin$stage <- clin$pathologic_stage
clin$stage <- ifelse(clin$stage %in% c('Stage IA', 'Stage IB'), 'I', clin$stage)
clin$stage <- ifelse(clin$stage %in% c('Stage IIA', 'Stage IIB'), 'II', clin$stage)
clin$stage <- ifelse(clin$stage %in% c('Stage III'), 'III', clin$stage)
clin$stage <- ifelse(clin$stage %in% c('Stage IV'), 'IV', clin$stage)

# sex
table(clin$gender)
clin$sex <- ifelse(clin$gender == 'MALE', 'M', 'F')

# histology
table(clin$histological_type)
clin$histo <- clin$histological_type
clin$histo <- ifelse(clin$histo == 'Pancreas-Adenocarcinoma Ductal Type', 'PDAC', 'Other')

clin <- clin[order(clin$sample), ]
expr <- expr[, order(colnames(expr))]

## remove low/zero-expressed genes
remove <- rem(expr, missing.perc = 0.5, const.int = 1)
if(length(remove) > 0){
  expr <- expr[-remove, ]
}
dim(expr) # 2163 genes

################################################################
# Run CIBERSORT to infer immune cell types
################################################################
cell <- deconvo_tme(eset = expr, method = "cibersort", perm = 1000, arrays = FALSE)

# Identifying TME patterns
tme <- tme_cluster(input = cell, features = colnames(cell)[2:23], 
                   id = "ID", scale = TRUE, method = "kmeans", 
                   min_nc =2, max.nc = 2)

colnames(tme) <- gsub(colnames(tme), pattern = "_CIBERSORT", replacement = "")

################################################################
# Heatmap
################################################################
df <- tme[, -c(1,2)]
df <- t(df)
colnames(df) <- tme$ID
df <- t(scale(t(df)))

col = list( Histo = c( "Other" = "#67A9CF" , "PDAC" = "#EF8A62" ) ,
            Cluster = c("TME2" = "#4d004b",  "TME1" = "#fa9fb5"),
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
     width = 8, height = 7)

ht_data <- Heatmap(df, name = "expr",
                   top_annotation = ha,
                   show_row_names = TRUE, 
                   show_column_names = FALSE,
                   row_names_gp = gpar(fontsize = 8),
                   column_names_gp = gpar(fontsize = 8),
                   colorRamp2(c(-3, 0, 9), 
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


######################################################################
## cell type association with cluster groups 
######################################################################
signature_name <- colnames(tme)[-c(1,2)]

wilcox_res <- lapply(1:length(signature_name), function(i){
  
  df <- data.frame( y = as.numeric(tme[, signature_name[i]]),
                  group = tme$cluster)
  fit <- wilcox.test(df$y ~ df$group, exact = FALSE)
  data.frame(signature = signature_name[i],
             pval = round(fit$p.value, 3))
  
})

wilcox_res <- do.call(rbind, wilcox_res)
wilcox_res$fdr <- p.adjust(wilcox_res$pval, method = "BH")
write.csv(wilcox_res, file = file.path(dir_out, 'cibersort', 'wilcox_res.csv'), row.names = FALSE)


######################################################################
## cell type association with OS
######################################################################
signature_name <- colnames(tme)[-c(1,2)]

os_res <- lapply(1:length(signature_name), function(i){
  
  fit <- coxph(Surv(os_time, os_status) ~ as.numeric(tme[, signature_name[i]]), clin) 
  data.frame(signature = signature_name[i],
             coef = summary(fit)$coefficients[ , 'coef'],
             se = summary(fit)$coefficients[ , 'se(coef)'],
             pval = summary(fit)$coefficients[ , 'Pr(>|z|)'])
  
})

os_res  <- do.call(rbind, os_res )
os_res $fdr <- p.adjust(os_res $pval, method = "BH")
write.csv(os_res , file = file.path(dir_out, 'cibersort', 'os_res .csv'), row.names = FALSE)

######################################################################
## cell type association with PFI
######################################################################
signature_name <- colnames(tme)[-c(1,2)]

pfi_res <- lapply(1:length(signature_name), function(i){
  
  fit <- coxph(Surv(pfi_time, pfi_status) ~ as.numeric(tme[, signature_name[i]]), clin) 
  data.frame(signature = signature_name[i],
             coef = summary(fit)$coefficients[ , 'coef'],
             se = summary(fit)$coefficients[ , 'se(coef)'],
             pval = summary(fit)$coefficients[ , 'Pr(>|z|)'])
  
})

pfi_res  <- do.call(rbind, pfi_res )
pfi_res $fdr <- p.adjust(pfi_res $pval, method = "BH")
write.csv(pfi_res , file = file.path(dir_out, 'cibersort', 'pfi_res .csv'), row.names = FALSE)













