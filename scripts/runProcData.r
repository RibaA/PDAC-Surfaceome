#######################################################
## load libraries
#######################################################
library(TCGAbiolinks) 
library(MultiAssayExperiment)

dir_in <- 'data/raw'
dir_out <- 'data/proc'

###############################################################
## load TCGA PDAC data 
## https://www.bioconductor.org/packages/devel/bioc/vignettes/TCGAbiolinks/inst/doc/download_prepare.html#Transcriptome_Profiling 
##############################################################
# Load TCGA RNA-seq data
TCGAbiolinks:::getProjectSummary('TCGA-PAAD')

query_rna <- GDCquery(
  project = "TCGA-PAAD",
  data.category = "Transcriptome Profiling", 
  data.type = "Gene Expression Quantification", 
  workflow.type = "STAR - Counts"
)
GDCdownload(query_rna)
expdat <- GDCprepare(
    query = query_rna,
    save = TRUE, 
    save.filename = file.path(dir_out, 'exp.rda')
)


# Load TCGA RPPA data
query_rppa <- GDCquery(
    project = "TCGA-PAAD", 
    data.category = "Proteome Profiling",
    data.type = "Protein Expression Quantification"
)
GDCdownload(query_rppa) 
rppadat <- GDCprepare(
    query = query_rppa,
    save = TRUE, 
    save.filename = file.path(dir_out, 'rppa.rda')
)

###############################################################
## Preprocessing data
###############################################################
# 1- common samples between RPPA and RNA-seq
## RNA-seq data
expr <- assay(expdat)
clin <- data.frame(colData(expdat))
expr_annot <- data.frame(rowData(expdat))

patientid <- sapply(1:ncol(expr), function(k){
  id <- strsplit(colnames(expr)[k], "-")[[1]]
  paste(id[1], id[2], id[3], id[4], sep='-')
})
colnames(expr) <- patientid
expr <- log2(expr+1)

clin$patientid <- sapply(1:nrow(clin), function(k){
  id <- strsplit(clin$barcode[k], "-")[[1]]
  paste(id[1], id[2], id[3], id[4], sep='-')
})

## RPPA data
rppa <- rppadat[, grep('TCGA-', colnames(rppadat), ignore.case = FALSE, perl = FALSE)]
rppa_annot <- rppadat[, !colnames(rppadat) %in% colnames(rppa)]

int <- intersect(colnames(expr), colnames(rppa))
expr <- expr[, int]
rppa <- rppa[, int]
clin <- clin[clin$patientid %in% int, ]

# 2- Protein-coding genes 
expr_annot <- expr_annot[expr_annot$gene_type == 'protein_coding', ]
expr_annot <- expr_annot[!duplicated(expr_annot$gene_name), ]
expr <- expr[rownames(expr) %in% rownames(expr_annot), ]
rownames(expr) <- expr_annot$gene_name

# save RPPA, RNA-seq and clinical data
paadData <- list('expr' = expr,
                 'rppa' = rppa,
                 'clin' = clin,
                 'annot_expr' = expr_annot,
                 'annot_rppa' = rppa_annot)
save(paadData, file = file.path(dir_out, 'paad_tcga.rda'))

