## ---- include = FALSE---------------------------------------------------------
knitr::opts_chunk$set(
    collapse = TRUE,tidy = TRUE,
    warning=FALSE, message=FALSE,
    comment = "##>"
)

## ---- eval = FALSE------------------------------------------------------------
#  if (!requireNamespace("BiocManager", quietly = TRUE))
#      install.packages("BiocManager")
#  BiocManager::install("bambu")
#  BiocManager::install("NanoporeRNASeq")

## -----------------------------------------------------------------------------
library(bambu)
test.bam <- system.file("extdata",
    "SGNex_A549_directRNA_replicate5_run1_chr9_1_1000000.bam",
    package = "bambu")
fa.file <- system.file("extdata",
    "Homo_sapiens.GRCh38.dna_sm.primary_assembly_chr9_1_1000000.fa",
    package = "bambu")
gtf.file <- system.file("extdata", "Homo_sapiens.GRCh38.91_chr9_1_1000000.gtf",
    package = "bambu")
bambuAnnotations <- prepareAnnotations(gtf.file)
se <- bambu(reads = test.bam, annotations = bambuAnnotations,
    genome = fa.file)

## -----------------------------------------------------------------------------
library(bambu)
library(NanoporeRNASeq)
data("SGNexSamples")
SGNexSamples
library(ExperimentHub)
NanoporeData <- query(ExperimentHub(), c("NanoporeRNA", "GRCh38","BAM"))
bamFiles <- Rsamtools::BamFileList(NanoporeData[["EH3808"]],
    NanoporeData[["EH3809"]],NanoporeData[["EH3810"]], NanoporeData[["EH3811"]],
    NanoporeData[["EH3812"]], NanoporeData[["EH3813"]])

## -----------------------------------------------------------------------------
# get path to fasta file
genomeSequenceData <- query(ExperimentHub(), c("NanoporeRNA", "GRCh38","FASTA"))
genomeSequence <- genomeSequenceData[["EH7260"]]

## -----------------------------------------------------------------------------
library(BSgenome.Hsapiens.NCBI.GRCh38)
genomeSequenceBsgenome <- BSgenome.Hsapiens.NCBI.GRCh38

## -----------------------------------------------------------------------------
gtf.file <- system.file("extdata", "Homo_sapiens.GRCh38.91_chr9_1_1000000.gtf",
    package = "bambu")
annotation <- prepareAnnotations(gtf.file)

## -----------------------------------------------------------------------------
txdb <- system.file("extdata", "Homo_sapiens.GRCh38.91_chr9_1_1000000.gtf",
    package = "bambu")
annotation <- prepareAnnotations(txdb)

## -----------------------------------------------------------------------------
data("HsChr22BambuAnnotation")
HsChr22BambuAnnotation

## ---- results = "hide"--------------------------------------------------------
se <- bambu(reads = bamFiles, annotations = HsChr22BambuAnnotation,
    genome = genomeSequence, ncore = 1)

## -----------------------------------------------------------------------------
se

## -----------------------------------------------------------------------------
colData(se)$condition <- as.factor(SGNexSamples$cellLine)

## ---- results = "hide"--------------------------------------------------------
seUnextended <- bambu(reads = bamFiles,
    annotations = HsChr22BambuAnnotation,
    genome = genomeSequence, discovery = FALSE)

## -----------------------------------------------------------------------------
seUnextended

## ---- fig.width = 8, fig.height = 6-------------------------------------------
library(ggplot2)
plotBambu(se, type = "heatmap")

## ---- fig.width = 8, fig.height = 6-------------------------------------------
plotBambu(se, type = "pca")

## ---- fig.width = 8, fig.height = 10------------------------------------------
plotBambu(se, type = "annotation", gene_id = "ENSG00000099968")

## -----------------------------------------------------------------------------
seGene <- transcriptToGeneExpression(se)
seGene

## ---- fig.width = 8, fig.height = 6-------------------------------------------
colData(seGene)$groupVar <- SGNexSamples$cellLine
plotBambu(seGene, type = "heatmap")

## -----------------------------------------------------------------------------
save.dir <- tempdir()
writeBambuOutput(se, path = save.dir, prefix = "NanoporeRNASeq_")

## -----------------------------------------------------------------------------
save.file <- tempfile(fileext = ".gtf")
writeToGTF(rowRanges(se), file = save.file)

## -----------------------------------------------------------------------------
library(DESeq2)
dds <- DESeqDataSetFromMatrix(round(assays(seGene)$counts),
                                    colData = colData(se),
                                    design = ~ condition)
dds.deseq <- DESeq(dds)
deGeneRes <- DESeq2::results(dds.deseq, independentFiltering = FALSE)
head(deGeneRes[order(deGeneRes$padj),])

## -----------------------------------------------------------------------------
summary(deGeneRes)

## ---- fig.width = 8, fig.height = 6-------------------------------------------
library(apeglm)
resLFC <- lfcShrink(dds.deseq, coef = "condition_MCF7_vs_K562", type = "apeglm")
plotMA(resLFC, ylim = c(-3,3))

## -----------------------------------------------------------------------------
library(DEXSeq)
dxd <- DEXSeqDataSet(countData = round(assays(se)$counts), 
sampleData = as.data.frame(colData(se)), 
design = ~sample + exon + condition:exon,
featureID = rowData(se)$TXNAME, 
groupID = rowData(se)$GENEID)
dxr <- DEXSeq(dxd)
head(dxr)

## ----fig.width = 8, fig.height = 6--------------------------------------------
plotMA(dxr, cex = 0.8 )

## ---- eval = FALSE------------------------------------------------------------
#  se <- bambu(reads = bamFiles, rcOutDir = "./bambu/",
#      annotations = annotaiton,
#      genome = genomeSequence)

## ---- eval = FALSE------------------------------------------------------------
#  bambu(reads, annotations, genome,
#      opt.discovery = list(min.readCount = 5))

## ---- eval = FALSE------------------------------------------------------------
#  bambu(reads, annotations, genome,
#      opt.discovery = list(min.sampleNumber = 5))

## ---- eval = FALSE------------------------------------------------------------
#  bambu(reads, annotations, genome,
#      opt.discovery = list(min.readFractionByGene = 0.1))

## ---- eval = FALSE------------------------------------------------------------
#  bambu(reads, annotations, genome, opt.em = list(bias = FALSE))

## ---- eval = FALSE------------------------------------------------------------
#  bambu(reads, annotations, genome, ncore = 8)

## -----------------------------------------------------------------------------
sessionInfo()

