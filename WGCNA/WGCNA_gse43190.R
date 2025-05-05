# Only run the following commands once to install WGCNA and flashClust on your computer

install.packages("BiocManager")
BiocManager::install("WGCNA")
install.packages("flashClust") 
BiocManager::install("hclust")

# Load WGCNA and flashClust libraries every time you open R
library(WGCNA)
library(flashClust)
#library(hclust)
cor = WGCNA::cor

#Set your current working directory (where all your files are)
#setwd("E:/VARSHA RAMESH/Transcriptomics/West Nile Virus/GSE46681")
setwd("E:/VARSHA RAMESH/Transcriptomics/West Nile Virus/GSE43190")

# Uploading data into R and formatting it for WGCNA --------------

# This creates an object called "datExpr" that contains the normalized counts file output from DESeq2
datExpr = read.csv("normalized expression data.csv")
# "head" the file to preview it
head(datExpr) # You see that genes are listed in a column named "X" and samples are in columns


#######Extract Gene names corresponding to Illumina ID#########
library(illuminaHumanv4.db)

ls("package:illuminaHumanv4.db")
Ref_id <- datExpr$X
Symbols <- unlist(mget(Ref_id, illuminaHumanv4SYMBOL, ifnotfound = NA))
rma <- cbind(Symbols, datExpr)

#remove the unnecessary columns
rma$X <- NULL
rownames(rma) <- NULL

#remove NA values and duplicates
rma <- na.omit(rma)
rma <- rma[!duplicated(rma$Symbols),]

rownames(rma) <- rma$Symbols
rma$Symbols <- NULL

write.csv(rma, "Expression Data_GSE43190.csv")

#final normalized gene expression dataset
ExprsDat <- t(rma)[ ,1:10000] #transpose the matrix (row to culums & colums to rows)


# Manipulate file so it matches the format WGCNA needs 
#row.names(datExpr) = datExpr$X
#datExpr$X = NULL
#datExpr = as.data.frame(t(datExpr))[, 1:5000] # now samples are rows and genes are columns
#dim(datExpr) 


# Run this to check if there are gene outliers
gsg = goodSamplesGenes(rma, verbose = 3)
gsg$allOK 

#If the last statement returns TRUE, all genes have passed the cuts. If not, we remove the offending genes and samples from the data with the following:

#if (!gsg$allOK)
#	{if (sum(!gsg$goodGenes)>0)
#		printFlush(paste("Removing genes:", paste(names(datExpr)[!gsg$goodGenes], collapse= ", ")));
#		if (sum(!gsg$goodSamples)>0)
#			printFlush(paste("Removing samples:", paste(rownames(datExpr)[!gsg$goodSamples], collapse=", ")))
#		datExpr= datExpr[gsg$goodSamples, gsg$goodGenes]
#		}

#Create an object called "datTraits" that contains your trait data
setwd("E:/VARSHA RAMESH/Transcriptomics/West Nile Virus/GSE43190")

Traits.dat = read.csv("GSE43190_traits.csv")
head(Traits.dat)

#form a data frame analogous to expression data that will hold the clinical traits.
rownames(Traits.dat) = Traits.dat$Samples
Traits.dat$Samples <- NULL
Traits <- Traits.dat
Traits.dat$Condition <- NULL
table(rownames(Traits.dat)==rownames(ExprsDat)) #should return TRUE if datasets align correctly, otherwise your names are out of order

# You have finished uploading and formatting expression and trait data
# Expression data is in datExpr, corresponding traits are datTraits

setwd("E:/VARSHA RAMESH/Transcriptomics/West Nile Virus/WGCNA_GSE43190/WGCNA UNSIGNED")
save(ExprsDat, Traits.dat, file="SamplesAndTraits.RData")
load("SamplesAndTraits.RData")

# Cluster samples by expression ----------------------------------------------------------------
##Setting threshold for outlier detection

A = adjacency(t(ExprsDat),type="unsigned") # this calculates the whole network connectivity
k = as.numeric(apply(A,2,sum))-1 # standardized connectivity
Z.k = scale(k)
thresholdZ.k = -2.5 # often -2.5
outlierColor = ifelse(Z.k<thresholdZ.k,"red","black")
sampleTree = flashClust(as.dist(1-A), method = "average")
Traits.dat_numeric <- as.numeric(as.factor(Traits$Traits)) ## The plot shows asymptomatic & severe disease condition as traits.
# Convert traits to a color representation where red indicates high values
traitColors <- data.frame(numbers2colors(Traits.dat_numeric, signed = FALSE))
dimnames(traitColors)[[2]] = paste(names(Traits.dat))
datColors <- data.frame(outlier = outlierColor,traitColors)

plotDendroAndColors(sampleTree,groupLabels=names(datColors),
                    colors=datColors,main="Sample Outlier Detection")

# Day "0" outliers have been identified. You could exclude these samples with the code below, but this scientists had a biological reason to NOT exclude these samples. It's up to you. Justify whatever decision you make.

# Remove outlying samples 
remove.samples = Z.k<thresholdZ.k | is.na(Z.k)
datExprOut = ExprsDat[!remove.samples,]
Traits.datOut = Traits.dat[!remove.samples,]
Traits = Traits[!remove.samples,]
save(datExprOut, Traits.datOut, file="SamplesAndTraits_OutliersRemoved.RData")


# Choose a soft threshold power- USE A SUPERCOMPUTER IRL ------------------------------------
powers = c(c(1:10), seq(from =10, to=20, by=2)) #choosing a set of soft-thresholding powers
sft = pickSoftThreshold(datExprOut, powerVector=powers, verbose =5, networkType="unsigned") #call network topology analysis function
sft.data <- sft$fitIndices

#sizeGrWindow(9,5)
#par(mfrow= c(2,1))
#cex1= 0.9
#plot(sft$fitIndices[,1], -sign(sft$fitIndices[,3])*sft$fitIndices[,2], xlab= "Soft Threshold (power)", ylab="Scale Free Topology Model Fit, signed R^2", type= "n", main= paste("Scale independence")) 
#text(sft$fitIndices[,1], -sign(sft$fitIndices[,3])*sft$fitIndices[,2], labels=powers, cex=cex1, col="red")
#abline(h=0.90, col="red")
#plot(sft$fitIndices[,1], sft$fitIndices[,5], xlab= "Soft Threshold (power)", ylab="Mean Connectivity", type="n", main = paste("Mean connectivity"))
#text(sft$fitIndices[,1], sft$fitIndices[,5], labels=powers, cex=cex1, col="red")
library(ggplot2)
##### OR #####
a1 <- ggplot(sft.data, aes(Power, SFT.R.sq, label = Power)) +
  geom_point() +
  geom_text(nudge_y = 0.1) + 
  geom_hline(yintercept = 0.9, color = 'red') +
  labs(x = 'Soft Threshhold (Power)', y = 'Sacle free topology model fit, signed R^2', title = 'Scale independence') +
  theme_classic() +
  theme(plot.title = element_text(hjust = 0.5))

a2 <- ggplot(sft.data, aes(Power, mean.k., label = Power)) +
  geom_point() +
  geom_text(nudge_y = 0.1) + 
  labs(x = 'Soft Threshhold (Power)', y = 'Mean Connectivity',  title = ('Mean Connectivity')) +
  theme_classic() +
  theme(plot.title = element_text(hjust = 0.5))

grid.arrange(a1, a2, ncol = 2)
##############

#from this plot, we would choose a power of 18 because it's the lowest power for which the scale free topology index reaches 0.90

enableWGCNAThreads()
sft$powerEstimate
softPower = 6
adjacency = adjacency(datExprOut, power = softPower, type = "unsigned") #specify network type

# Construct Networks- USE A SUPERCOMPUTER IRL -----------------------------

#translate the adjacency into topological overlap matrix and calculate the corresponding dissimilarity:
TOM = TOMsimilarity(adjacency, TOMType="unsigned") # specify network type
dissTOM = 1-TOM


# Generate Modules --------------------------------------------------------

# Generate a clustered gene tree
geneTree = flashClust(as.dist(dissTOM), method="average")

plot(geneTree, xlab="", sub="", main= "Gene Clustering on TOM-based dissimilarity", labels= FALSE, hang=0.04)

#This sets the minimum number of genes to cluster into a module
minModuleSize = 30 
dynamicMods = cutreeDynamic(dendro= geneTree, distM= dissTOM, deepSplit=2, pamRespectsDendro= FALSE, minClusterSize = minModuleSize)

dynamicColors= labels2colors(dynamicMods)
MEList= moduleEigengenes(datExprOut, colors= dynamicColors,softPower = 3)
MEs= MEList$eigengenes
MEDiss= 1-cor(MEs)
METree= flashClust(as.dist(MEDiss), method= "average")

save(dynamicMods, MEList, MEs, MEDiss, METree, file= "Network_allSamples_signed_RLDfiltered.RData")

#plots tree showing how the eigengenes cluster together
plot(METree, main= "Clustering of module eigengenes", xlab= "", sub= "")

#set a threhold for merging modules. In this example we are not merging so MEDissThres=0.0
MEDissThres = 0.25
merge = mergeCloseModules(datExprOut, dynamicColors, cutHeight= MEDissThres, verbose =3)
mergedColors = merge$colors
mergedMEs = merge$newMEs


#plot dendrogram with module colors below it
plotDendroAndColors(geneTree, cbind(dynamicColors, mergedColors), 
                    c("Dynamic Tree Cut ", "Merged dynamic"), 
                    dendroLabels= FALSE, 
                    hang=0.03, 
                    addGuide= TRUE, 
                    guideHang=0.05)

moduleColors = mergedColors
colorOrder = c("grey", standardColors(50))
moduleLabels = match(moduleColors, colorOrder)-1
MEs = mergedMEs

save(MEs, moduleLabels, moduleColors, geneTree, file= "Network_allSamples_signed_nomerge_RLDfiltered.RData")


##############-start-########################
####################

#to check
datExprOut_datframe <- as.data.frame(datExprOut)
cor = WGCNA::cor

mergingThresh = 0.25
net = blockwiseModules(datExprOut_datframe, corType = "pearson", maxBlockSize = 5000, 
                       networkType = "unsigned", power = 6, 
                       minModuleSize = 30, mergeCutHeight = 0.25, 
                       numericLabels = FALSE)

softPower = 6  # Based on the analysis

# Construct the network and detect modules
cor <- stats::cor
net = blockwiseModules(datExprOut,
                       power = softPower,
                       TOMType = "unsigned",
                       minModuleSize = 30,
                       reassignThreshold = 0,
                       mergeCutHeight = 0.25,
                       numericLabels = FALSE,
                       corType = "pearson",
                       verbose = 3)


moduleLabelsAutomatic = net$colors

# Convert labels to colors for plotting
moduleColorsAutomatic = labels2colors(moduleLabelsAutomatic)

# A data frame with module eigengenes can be obtained as follows
MEsAutomatic = net$MEs

blocknumber = 1

net_unmerged = net$unmergedColors
net_merged = net$colors


datcolors = as.data.frame(cbind(net_unmerged, net_merged))[net$blockGenes[[blocknumber]], ]

dev.off()
plotDendroAndColors(net$dendrograms[[blocknumber]], 
                    cbind(net_unmerged, net_merged)[net$blockGenes[[blocknumber]], ], 
                    c("unmerged", "merged"), 
                    dendroLabels = FALSE, hang = 0.03, 
                    addGuide = TRUE, guideHang = 0.05)


###----trial------###
# previous analysis
moduleLabelsBlockwise = matchLabels(net$colors, moduleLabelsAutomatic)
# Convert labels to colors for plotting
moduleColorsBlockwise = labels2colors(moduleLabelsBlockwise)

# measure agreement with single block automatic procedure
mean(moduleLabelsBlockwise == moduleLabelsAutomatic)

blockNumber = 1
# Plot the dendrogram for the chosen block
plotDendroAndColors(net$dendrograms[[blockNumber]], moduleColorsBlockwise[net$blockGenes[[blockNumber]]], 
                    "Module colors", main = paste("Dendrogram and module colors in block", blockNumber), 
                    dendroLabels = FALSE, hang = 0.03, addGuide = TRUE, guideHang = 0.05)

###----trial end-----###


###################
###############-end-#####################



library(CorLevelPlot)
library(dplyr)
library(tidyverse)
library(gridExtra)


# 6A. Correlate traits ---------------------------------------------------------------------

traits.corr <- Traits %>% 
  mutate(disease_status_bin = ifelse(grepl('Severe', Condition), 1, 0)) %>%
  dplyr::select(3)

#binarize categorical variables
Traits$Traits <- factor(Traits$Traits, 
                        levels = c("Asymptomatic", "Encephalitis", "Meningitis", "Encephalitis and Meningitis"))
severity.out <- binarizeCategoricalColumns(Traits$Traits,
                                           includePairwise = FALSE,
                                           includeLevelVsAll = TRUE,
                                           minCount = 1)
traits.new <- cbind(traits.corr, severity.out)


#change the colnames of the traits.new data
colnames(traits.new)[colnames(traits.new) == "disease_status_bin"] <- "Overall_severity_status"
colnames(traits.new)[colnames(traits.new) == "data.Encephalitis.vs.all"] <- "Encephalitis"
colnames(traits.new)[colnames(traits.new) == "data.Meningitis.vs.all"] <- "Meningitis"
colnames(traits.new)[colnames(traits.new) == "data.Encephalitis and Meningitis.vs.all"] <- "Encephalitis and Meningitis"



#Define number of genes and samples
nGenes = ncol(datExprOut)
nSamples = nrow(datExprOut)

#Recalculate MEs with color labels
MEs0 = moduleEigengenes(datExprOut, net_merged)$eigengenes
MEs = orderMEs(MEs0)
summary(MEs)
moduleTraitCor = cor(MEs, traits.new, use= "p")
moduleTraitPvalue = corPvalueStudent(moduleTraitCor, nSamples)

setwd("E:/VARSHA RAMESH/Transcriptomics/West Nile Virus/WGCNA_GSE43190/WGCNA UNSIGNED/Plots")

#Print correlation heatmap between modules and traits
textMatrix= paste(signif(moduleTraitCor, 2), "(", 
                  signif(moduleTraitPvalue, 1), ")", sep= "")
dim(textMatrix)= dim(moduleTraitCor)
sizeGrWindow(9,5)

pdf("heatmap_plot_1.pdf", width = 10, height = 10)
tiff("heatmap_plot.tiff", width = 2420, height = 2610, res = 300) 
png("heatmap_plot.png", width = 2420, height = 2610, res = 300)

par(mar = c(10, 9, 7, 9))
#display the correlation values with a heat map plot
labeledHeatmap(Matrix= moduleTraitCor, 
               xLabels= names(traits.new), 
               yLabels= names(MEs), 
               ySymbols= names(MEs), 
               colorLabels= FALSE, 
               colors= greenWhiteRed(50), 
               textMatrix= textMatrix, 
               setStdMargins= FALSE, 
               cex.text= 0.6,
               cex.lab.x = 0.8,
               cex.lab.y = 0.8,
               zlim= c(-1,1),
               las = 1,
               main= paste("Module-trait relationships"))

dev.off()
#########################

#To generate Eigengene adjacency heatmap

png("Eigengene adjacency heatmap.png", width = 2420, height = 2610, res = 300)
tiff("Eigengene adjacency heatmap.tiff", width = 2420, height = 2610, res = 300)
par(mar = c(10,9,7,9))

plotEigengeneNetworks(MEs, 
                      "Eigengene adjacency heatmap", 
                      marDendro = c(0, 4, 1, 2), 
                      marHeatmap = c(3, 4, 1, 2), 
                      cex.lab = 0.8, 
                      xLabelsAngle = 90)
dev.off()

##--------------Module membership trial------------------##

#Defining the variable Peso10dias containing the column Peso10dias of datTrait

Overall_sev_Stat = as.data.frame(traits.new$Overall_severity_status)
#Encephalitis = as.data.frame(traits.new$Encephalitis)

#names (colors) of the modules
modNames = substring(names(MEs), 3)

geneModuleMembership = as.data.frame(cor(datExprOut, MEs, use = "p"))
MMPvalue = as.data.frame(corPvalueStudent(as.matrix(geneModuleMembership), nSamples))

names(geneModuleMembership) = paste("MM", modNames, sep="")
names(MMPvalue) = paste("p.MM", modNames, sep="")

geneTraitSignificance = as.data.frame(cor(datExprOut, Overall_sev_Stat, use = "p"))
##geneTraitSignificance = as.data.frame(cor(datExprOut, Encephalitis, use = "p"))

GSPvalue = as.data.frame(corPvalueStudent(as.matrix(geneTraitSignificance), nSamples))

names(geneTraitSignificance) = paste("GS.", names(Overall_sev_Stat), sep="")
#names(geneTraitSignificance) = paste("GS.", names(Encephalitis), sep="")

names(GSPvalue) = paste("p.GS.", names(Overall_sev_Stat), sep="")
#names(GSPvalue) = paste("p.GS.", names(Encephalitis), sep="")


module = c("salmon")  #########################putting the color below the plot
column = match(module, modNames)
moduleGenes = net_merged==module

sizeGrWindow(7, 7)
par(mfrow = c(1,1))
verboseScatterplot(abs(geneModuleMembership[moduleGenes, column]),
                   abs(geneTraitSignificance[moduleGenes, 1]),
                   xlab = paste("Module Membership in", module, "module"),
                   ylab = "Gene significance for Overall_Severity_Status",
                   main = paste("Module membership vs. gene significance\n"),
                   cex.main = 1.2, cex.lab = 1.2, cex.axis = 1.2, col = module)

#Identifying most important genes for one determined characteristic inside of the cluster
geneInfo0 = data.frame(EST = colnames(datExprOut),
                       moduleColor = net_merged,
                       geneTraitSignificance,
                       GSPvalue)

salmon_genes = geneInfo0$moduleColor == "salmon"
important.salmon.Genes <- geneInfo0[salmon_genes, ]
write.csv(important.salmon.Genes, "Important SALMON genes.csv")

brown_genes = geneInfo0$moduleColor == "brown"
important.brown.Genes = geneInfo0[brown_genes, ]
write.csv(important.brown.Genes, "Important BROWN genes.csv")
#turquoise_genes_meningt = geneInfo0$moduleColor == "turquoise"
#important.turquoise.Genes = geneInfo0[turquoise_genes_meningt, ]

magenta_genes = geneInfo0$moduleColor == "magenta"
important.magenta.Genes = geneInfo0[magenta_genes, ]
write.csv(important.magenta.Genes, "Important MAGENTA genes.csv")

modOrder = order(-abs(cor(MEs, Overall_sev_Stat, use = "p")))
for (mod in 1:ncol(geneModuleMembership))
{
  oldNames = names(geneInfo0)
  geneInfo0 = data.frame(geneInfo0, geneModuleMembership[, modOrder[mod]], 
                         MMPvalue[, modOrder[mod]]);
  names(geneInfo0) = c(oldNames, paste("MM.", modNames[modOrder[mod]], sep=""),
                       paste("p.MM.", modNames[modOrder[mod]], sep=""))
}
geneOrder = order(geneInfo0$moduleColor, -abs(geneInfo0$GS.traits.new.Overall_severity_status))
geneInfo = geneInfo0[geneOrder, ]
#geneOrder_meningt = order(geneInfo0$moduleColor, -abs(geneInfo0$GS.traits.new.Meningitis))
#geneInfo_meningt = geneInfo0[geneOrder_meningt, ]
#write.csv(geneInfo_meningt, "GeneInfo_1000(Meningitis).csv")


#if you want to write the information in a csv file, just uncomment line below
setwd("E:/VARSHA RAMESH/Transcriptomics/West Nile Virus/WGCNA_GSE43190/WGCNA UNSIGNED/Plots")

write.csv(geneInfo, file = "geneInfo_10000_unsigned.csv")   #all genes corresponding to each module saved

Magenta_08 <- subset(geneInfo, geneInfo[,"MM.magenta"] > 0.7)
head(Magenta_08)
write.csv(Magenta_08, "Hub genes of MAGENTA module_07(16).csv")

Turquoise_08 <- subset(geneInfo, geneInfo[, "MM.turquoise"] > 0.7)
head(Turquoise_08)
write.csv(Turquoise_08, "Hub genes of TURQUOISE module_07(102).csv")

Blue_08 <- subset(geneInfo, geneInfo[, "MM.blue"] > 0.7)
head(Blue_08)
write.csv(Blue_08, "Hub genes of BLUE module_07(15).csv")



###############---------------END------------------#####################
