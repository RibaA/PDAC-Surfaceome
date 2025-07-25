#######################################################
## load libraries
#######################################################
library(UCSCXenaTools)
library(dplyr)
library(data.table)
library(MultiAssayExperiment)
library(S4Vectors)
library(biomaRt)
library(qs)

dir_in <- 'data/raw'
dir_out <- 'data/proc'

###############################################################
## load TCGA PAAD and clinical data  
##############################################################
# -- Load UCSCXenaTools metadata
data(XenaData)
xd <- as.data.table(XenaData)
cohort <- "TCGA Pancreatic Cancer (PAAD)"

# -- Get RPPA dataset identifiers
rppa_datasets <- xd[XenaCohorts == cohort & grepl("RPPA", DataSubtype), XenaDatasets]

# -- Get clinical matrix
clinical_datasets <- xd[XenaCohorts == cohort & Type == "clinicalMatrix", XenaDatasets]

# -- Download clinical data
clinical_data <- XenaGenerate(subset = XenaCohorts == cohort & Type == "clinicalMatrix") |>
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

clin <- sData[cData, on = "sample"] # 196 samples with clinical data

# -- curate and standarized clinical variables
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

# -- Prepare RPPA matrix
rppa_data <- as.data.frame(rppa_data)
rownames(rppa_data) <- rppa_data$sample
rppa_mat <- rppa_data[, -1] # 229 proteins and 123 samples

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
rownames(rppa_mat) <- df$updateProtein # 219 proteins and 123 samples

# -- Build SummarizedExperiment objects (RPPA)
rppa_se <- SummarizedExperiment(
  assays = list(expr = rppa_mat),
  colData = clin[colnames(rppa_mat), , drop = FALSE],
  rowData = df
)

################################################################
## load RNA-seq count data using TCGA HUB
## download manually using 'TCGA Hub' using 'https://xenabrowser.net/datapages/?cohort=GDC%20TCGA%20Pancreatic%20Cancer%20(PAAD)&removeHub=https%3A%2F%2Fxena.treehouse.gi.ucsc.edu%3A443'
################################################################ 
# -- Download RNA-seq TPM data
rna_data <- read.delim(file = file.path(dir_in, 'TCGA-PAAD.star_tpm.tsv'), stringsAsFactors = FALSE)

# Remove the version numbers in Ensembl ID.
rna_data$Ensembl_ID <- substring(rna_data$Ensembl_ID, 1, 15)
rna_data <- rna_data[!duplicated(rna_data$Ensembl_ID), ]
rownames(rna_data) <- rna_data$Ensembl_ID
rna_data <- rna_data[, -1]
colnames(rna_data) <- gsub("\\.", "-", colnames(rna_data))
colnames(rna_data) <- substr(colnames(rna_data), 1, nchar(colnames(rna_data))-1)

# gene annotation
ensembl <- useEnsembl(biomart = "genes", dataset = "hsapiens_gene_ensembl")
gene_info <- getBM(
  attributes = c("ensembl_gene_id", "hgnc_symbol"),
  filters = "ensembl_gene_id",
  values = rownames(rna_data),
  mart = ensembl
)

gene_info <- gene_info[!duplicated(gene_info$hgnc_symbol), ]
annot <- gene_info[!duplicated(gene_info$ensembl_gene_id), ]
rna_mat <- rna_data[rownames(rna_data) %in% annot$ensembl_gene_id, ]
# rownames(rna_mat) <- annot$hgnc_symbol

# -- Build SummarizedExperiment objects (RNA)
rna_se <- SummarizedExperiment(
  assays = list(expr = rna_mat),
  colData = clin[colnames(rna_mat), , drop = FALSE],
  rowData = annot
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
