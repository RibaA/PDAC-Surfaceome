#######################################################
## load libraries
#######################################################
library(UCSCXenaTools)
library(dplyr)
library(data.table)
library(MultiAssayExperiment)
library(S4Vectors)
library(qs)

dir_in <- 'data/raw'
dir_out <- 'data/proc'

###############################################################
## load TCGA PAAD data 
##############################################################
# -- Load UCSCXenaTools metadata
data(XenaData)
xd <- as.data.table(XenaData)
cohort <- "TCGA Pancreatic Cancer (PAAD)"

# -- Get RNA-seq (TPM) and RPPA dataset identifiers
rna_datasets <- xd[XenaCohorts == cohort & grepl("HiSeqV2", XenaDatasets), XenaDatasets]
rppa_datasets <- xd[XenaCohorts == cohort & grepl("RPPA", DataSubtype), XenaDatasets]

# -- Get clinical matrix
clinical_datasets <- xd[XenaCohorts == cohort & Type == "clinicalMatrix", XenaDatasets]

# -- Download clinical data
clinical_data <- XenaGenerate(subset = XenaCohorts == cohort & Type == "clinicalMatrix") |>
  XenaQuery() |>
  XenaDownload(destdir = "local_data", download_probeMap = TRUE) |>
  XenaPrepare()

class(clinical_data$PAAD_clinicalMatrix)

# -- Download RNA-seq data
rna_data <- XenaGenerate(subset = XenaDatasets %in% rna_datasets) |>
  XenaQuery() |>
  XenaDownload(destdir = "local_data", download_probeMap = TRUE) |>
  XenaPrepare()

# -- Download RPPA data
rppa_data <- XenaGenerate(subset = XenaDatasets %in% rppa_datasets) |>
  XenaQuery() |>
  XenaDownload(destdir = "local_data", download_probeMap = TRUE) |>
  XenaPrepare()

# Merge clinical datasets
cData <- copy(clinical_data[["PAAD_clinicalMatrix"]])
cData <- as.data.table(cData)
setnames(cData, "sampleID", "sample")
setkeyv(cData, "sample")

sData <- copy(clinical_data[["PAAD_survival.txt"]])
sData <- as.data.table(sData)
setkeyv(sData, "sample")

merged_cData <- sData[cData, on = "sample"] # 196 samples with clinical data

# Load only RNA-seq expression and gene annotation
rna_expr <- as.data.frame(rna_data[["HiSeqV2.gz"]])
rownames(rna_expr) <- rna_expr$sample
rna_mat <- rna_expr[, -1] # 183 samples

# -- Build SummarizedExperiment objects (RNA)
rna_se <- SummarizedExperiment(
  assays = list(expr = rna_mat),
  colData = merged_cData[colnames(rna_mat), , drop = FALSE]
)

# -- Prepare RPPA matrix
rppa_data <- as.data.frame(rppa_data)
rownames(rppa_data) <- rppa_data$sample
rppa_mat <- rppa_data[, -1]

# -- Create annotation file for RPPA data
annot <- read.csv(file = file.path(dir_in, 'MCLP_Proteins.csv'))
df <- data.frame(protein = rownames(rppa_mat),
                 updateProtein = substr(rownames(rppa_mat), 1, nchar(rownames(rppa_mat))-4))

df <- df[!duplicated(df$updateProtein), ]

int <- intersect(annot$updateProtein, df$updateProtein) # 219 proteins in common
df <- df[df$updateProtein %in%int, ]
df$gene_name <- sapply(1:nrow(df), function(k){
    annot[annot$updateProtein == df$updateProtein[k], 'Gene']
})

write.csv(df, file = file.path(dir_out, 'protein_name.csv'), row.names= FALSE)

# -- Update RPPA matrix
rppa_mat <- rppa_mat[rownames(rppa_mat) %in% df$protein, ]
rppa_mat <- rppa_mat[order(rownames(rppa_mat)), ]
df <- df[order(df$protein), ]
rownames(rppa_mat) <- df$updateProtein

# -- Build SummarizedExperiment objects (RPPA)
rppa_se <- SummarizedExperiment(
  assays = list(expr = rppa_mat),
  colData = merged_cData[colnames(rppa_mat), , drop = FALSE],
  rowData = df
)

# -- Create MultiAssayExperiment
mae <- MultiAssayExperiment(
  experiments = list(RNAseq = rna_se, RPPA = rppa_se),
  metadata = list(
    source = "UCSC Xena (TCGA PAAD)",
    download_date = Sys.Date(),
    sessionInfo = sessionInfo()
  )
)

# -- Save as serialized object
qsave(mae, file = file.path(dir_out, "TCGA_PAAD_mae.qs"), nthread = getDTthreads())
