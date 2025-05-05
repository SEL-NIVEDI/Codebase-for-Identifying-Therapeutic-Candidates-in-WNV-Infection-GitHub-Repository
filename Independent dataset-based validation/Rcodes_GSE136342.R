install.packages("tidyverse")
BiocManager::install("airway")
BiocManager::install("DESeq2")
library(DESeq2)
library(tidyverse)
library(airway)
library(ggplot2)

#Step1: Preparing count Data ------

setwd("E:/VARSHA RAMESH/Transcriptomics/West Nile Virus/GSE136342")

#read in counts data

counts.dat <- read.delim("GSE136342_allnormalized_deseq_counts.txt", sep = "\t")
rownames(counts.dat) <- counts.dat[, 1]
counts.dat$GeneID <- NULL
counts.dat <- counts.dat[, -c(21:35)]
head(counts.dat)

#read in sample info

colData <- read.csv("GSE136342_phenodata.csv")
rownames(colData) <- colData[, 1]

#make sure the rownames in colData matches to column names in counts_data
all(colnames(counts.dat)%in% rownames(colData))
#are they in the same order?
all(colnames(counts.dat) == rownames(colData))

#Step2: Construct a DESeqDataSet object -------

counts.dat <- na.omit(counts.dat)

counts.dat <- round(counts.dat)

dds <- DESeqDataSetFromMatrix(countData = counts.dat,
                              colData = colData,
                              design = ~ State)

#pre-filtering: removing rows with low gene counts
#keeping rows that have at least 10 reads total
#dds <- dds[rowSums(counts(dds) >= 5)]
nrow(dds)

#keep <- rowSums(counts(dds)) >=10
#dds <- dds[keep,]

#set the factor level. (when reference 'ref=' level is not mentioned, it will automatically chose reference according to the alphabetical order) )
dds$State <- relevel(dds$State, ref = "Mock")

#NOTE: collapse technical replicates (if any). Never collapse biological replicates

#Step3: Run DESeq ---------------------

dds <- DESeq(dds)
res <- results(dds)

#Explore results ----------

summary(res) #currently using padj <0.1 as threshold
dim(res)

filtered <- subset(res, pvalue <0.05)
write.csv(filtered, "DE results_new.csv")

up <- subset(filtered, log2FoldChange > 1)
up1 <- as.data.frame(up)
write.csv(up1, "up-regulated_new.csv")

down <- subset(filtered, log2FoldChange < -1)
down1 <- as.data.frame(down)
write.csv(down1, 'down-regulated_new.csv')
###########################

#####UNIQUE AND COMMON DEGs############

setwd("E:/VARSHA RAMESH/Transcriptomics/West Nile Virus/WGCNA_GSE43190/Results_new")
microarray_data <- read.csv("Genes of magenta module(positive correlated).csv")

setwd("E:/VARSHA RAMESH/Transcriptomics/West Nile Virus/GSE136342")
RNAseq_data <- read.csv("up-regulated_new.csv")

microarray_DEGs <- microarray_data$genes.magenta.module
RNAseq_DEGs <- RNAseq_data$X

##unque genes of microarray dataset (severe disease)
microarray_unique_genes <- setdiff(microarray_DEGs, RNAseq_DEGs)
SevereDisease_genes <- as.data.frame(microarray_unique_genes)

##unique genes of RNA-seq dataset (after 24hrs)
RNAseq_unique_genes <- setdiff(RNAseq_DEGs, microarray_DEGs)

#Common DEGs
Common_genes <- intersect(microarray_DEGs, RNAseq_DEGs)
list(Common_genes)
