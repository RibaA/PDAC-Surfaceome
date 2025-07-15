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
# Run TME signatures (186 ---> IOBR)
################################################################
sig_tme <- calculate_sig_score(pdata           = NULL,
                               eset            = expr,
                               signature       = signature_tme,
                               method          = "ssgsea")

tme <- t(sig_tme[, -1])
colnames(tme) <- sig_tme$ID

################################################################
# Heatmap
################################################################
df <- tme
df <- t(scale(t(df)))

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
pdf(file=file.path(dir_out, "tme", "heatmap_tme.pdf"),
     width = 8, height = 10)

ht_data <- Heatmap(df, name = "expr",
                   top_annotation = ha,
                   show_row_names = TRUE, 
                   show_column_names = FALSE,
                   row_names_gp = gpar(fontsize = 7),
                   column_names_gp = gpar(fontsize = 8),
                   colorRamp2(c(min(df), 0, max(df)), 
                              c("#045a8d", "white", "#800026")))

draw(ht_data,
     column_title_gp = gpar(fontsize = 8, fontface = "bold"),
     #merge_legends = TRUE,
     heatmap_legend_side = "right",
     annotation_legend_side="right")

dev.off()

######################################################################
## cell type association with OS
######################################################################
signature_name <- rownames(tme)

os_res <- lapply(1:length(signature_name), function(i){
  
  fit <- coxph(Surv(os_time, os_status) ~ as.numeric(tme[signature_name[i], ]), clin) 
  data.frame(signature = signature_name[i],
             coef = summary(fit)$coefficients[ , 'coef'],
             se = summary(fit)$coefficients[ , 'se(coef)'],
             pval = summary(fit)$coefficients[ , 'Pr(>|z|)'])
  
})

os_res  <- do.call(rbind, os_res )
os_res $fdr <- p.adjust(os_res $pval, method = "BH")
write.csv(os_res , file = file.path(dir_out, 'tme', 'os_res .csv'), row.names = FALSE)

######################################################################
## cell type association with PFI
######################################################################
signature_name <- rownames(tme)

pfi_res <- lapply(1:length(signature_name), function(i){
  
  fit <- coxph(Surv(pfi_time, pfi_status) ~ as.numeric(tme[signature_name[i], ]), clin) 
  data.frame(signature = signature_name[i],
             coef = summary(fit)$coefficients[ , 'coef'],
             se = summary(fit)$coefficients[ , 'se(coef)'],
             pval = summary(fit)$coefficients[ , 'Pr(>|z|)'])
  
})

pfi_res  <- do.call(rbind, pfi_res )
pfi_res $fdr <- p.adjust(pfi_res $pval, method = "BH")
write.csv(pfi_res , file = file.path(dir_out, 'tme', 'pfi_res .csv'), row.names = FALSE)













