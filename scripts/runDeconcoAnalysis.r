#######################################################
## load libraries
#######################################################
# devtools::install_github('dviraran/xCell')
library(dplyr)
library(qs)
library(xCell)
library(survival)
library(survminer)
library(ComplexHeatmap)
library(RColorBrewer)
library(MultiAssayExperiment)

dir_in <- 'data/proc'
dir_out <- 'data/results'

###############################################################
## load data (from Aim 1)
##############################################################
dat <- qread(file.path(dir_in, 'rna_rppa_surface_filtered.qs'))
rna <- dat$rna

# --------------------------------------------------
# Aim 2: Stratify tumors by immune and molecular features
# --------------------------------------------------

# Run xCell to infer TME cell types
xcell_scores <- xCellAnalysis(rna)
xcell_scores <- t(xcell_scores)

# Merge xCell scores with clinical metadata
xcell_scores <- xcell_scores[match(rownames(clin), rownames(xcell_scores)), ]
merged_data <- cbind(clin, xcell_scores)

# Survival analysis on high vs low expression for top surface genes
gene_list <- rownames(rna_surface)
dir.create("results/survival_curves", showWarnings = FALSE)
for (g in gene_list) {
    expr <- rna_surface[g, ]
    group <- ifelse(expr > median(expr, na.rm=TRUE), "High", "Low")
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
