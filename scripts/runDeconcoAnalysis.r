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
clin$stage <- ifelse(clin$stage %in% c('Stage IA', 'Stage IB'), 'Stage I', clin$stage)
clin$stage <- ifelse(clin$stage %in% c('Stage IIA', 'Stage IIB'), 'Stage II', clin$stage)

## remove low/zero-expressed genes
remove <- rem(expr, missing.perc = 0.5, const.int = 1)
if(length(remove) > 0){
  expr <- expr[-remove, ]
}
dim(expr) # 2163 genes

# Run CIBERSORT to infer immune cell types
cell <- deconvo_tme(eset = expr, method = "cibersort", perm = 1000, arrays = FALSE)

# Identifying TME patterns
tme <- tme_cluster(input = cell, features = colnames(cell)[2:23], 
                   id = "ID", scale = TRUE, method = "kmeans", 
                   min_nc =2, max.nc = 2)

colnames(tme) <- gsub(colnames(tme), pattern = "_CIBERSORT", replacement = "")
res <- sig_heatmap(input = tme, features = colnames(tme)[3:ncol(tme)], 
                   group = "cluster", path = "result", palette = 6)


# Genes associated with TME clusters
group <- ifelse( tme$cluster == "TME1", 1, 0 )

res <- lapply(1:nrow(expr), function(k){
  
  print(k)
  fit <- glm( group ~ as.numeric(expr[k, ]), family = binomial)
  data.frame(gene_name = rownames(expr)[k],
             coef = summary(fit)$coefficients[2, "Estimate"],
             se = summary(fit)$coefficients[2, "Std. Error"],
             pval = summary(fit)$coefficients[2, "Pr(>|z|)"])
  
})

res <- do.call(rbind, res)
res$fdr <- p.adjust(res$pval, method = "BH")


# cell abundance 

cols <- c('#2692a4', '#ffbe0b')

p1 <- sig_box(tme, variable = "cluster", 
              signature = "Macrophages_M1", jitter = TRUE,
              cols =  cols, show_pvalue = TRUE, size_of_pvalue = 4)

p2 <- sig_box(tme, variable = "cluster", 
              signature = "T_cells_CD8", jitter = TRUE,
              cols =  cols, show_pvalue = TRUE, size_of_pvalue = 4)

















































































sig_tme <- calculate_sig_score(pdata           = NULL,
                               eset            = expr,
                               signature       = signature_tme,
                               method          = "ssgsea")

# Survival analysis on high vs low expression for top surface genes
gene_list <- rownames(rna_surface)
dir.create("results/survival_curves", showWarnings = FALSE)


for (g in gene_list) {
    expr <- rna_surface[g, ]
    #group <- ifelse(expr > median(expr, na.rm=TRUE), "High", "Low")
    surv_obj <- Surv(merged_data$OS, merged_data$OS_event)
    fit <- survfit(surv_obj ~ group)
    ggsurvplot(fit, data = merged_data, title = g)
    ggsave(paste0("results/survival_curves/survival_", g, ".pdf"))
}

# Heatmap of top variable surface genes
var_genes <- rownames(rna_surface)[order(apply(rna_surface, 1, var), decreasing = TRUE)][1:50]
Heatmap(rna_surface[var_genes, ], cluster_rows = TRUE, cluster_columns = TRUE,
        name = "Expression", col = colorRampPalette(rev(brewer.pal(9, "RdBu")))(255))

# Save merged data
write.csv(merged_data, file = "results/merged_clinical_xcell.csv", row.names = TRUE)
