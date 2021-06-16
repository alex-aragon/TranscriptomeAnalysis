setwd('/Users/alejandroaragon/Desktop')

library(edgeR)
library(limma)
library(tximport)
library(readr)
library(tidyverse)
library(viridis)
library(ggplot2)
library(pheatmap)

###GETTING THE DATA###

CriSAMPLES <- read.table(file.path("RNASeqAnalysis/5.Quantification/0.kallistoQuantification", "kallisto_CriSamples.txt"), header = TRUE)
CriFiles <- file.path("RNASeqAnalysis/5.Quantification/0.kallistoQuantification/1.CeratopterisSamples_kallistoQuantification/", CriSAMPLES$dir, "abundance.h5")
names(CriFiles) <- paste0("SAMPLE", 1:9)
all(file.exists(CriFiles))
tx2gene <- read_csv(file.path("RNASeqAnalysis/5.Quantification", "CriTranscriptome_tx2gene.csv"))

CriTXITable <- tximport(CriFiles, type = "kallisto", tx2gene = tx2gene)


CriCounts <- CriTXITable$counts
CriLength <- CriTXITable$length #effective length for normalization

CriLength <- CriLength/exp(rowMeans(log(CriLength)))
normCriCounts <- CriCounts/CriLength
eff.lib <- calcNormFactors(normCriCounts) * colSums(normCriCounts)

CriLength <- sweep(CriLength, 2, eff.lib, "*")
CriLength <- log(CriLength)





###DATA PRE-PROCESSING###

##Filtering low expressing genes##
#Without normalization

CriDGE <- DGEList(CriCounts)
group <- as.factor(c("Leaves" ,"Leaves" , "Leaves" , "Root", "Root", "Root", "RootTip", "RootTip", "RootTip")) 
CriDGE$samples$group <- group

CriGenes_cpm <- cpm(CriDGE)
CriGenes_lcpm <- cpm(CriDGE, log=TRUE)

L <- mean(CriDGE$samples$lib.size)*1e-6
M <- median(CriDGE$samples$lib.size)*1e-6

expressedCriGenes_list <- filterByExpr(CriDGE, group=group)
expressedCriGenes <- CriDGE[expressedCriGenes_list, , keep.lib.sizes=FALSE] 

table(rowSums(CriDGE$counts==0)==9)
dim(expressedCriGenes)

lcpm.cutoff <- log2(10/M + 2/L)
 
 

 

#Normalization with TMM
TMMCriGenes <- calcNormFactors(expressedCriGenes, method = "TMM")
TMMCriGenes$samples$norm.factors
TMMCriGenes_lcpm <- cpm(TMMCriGenes, log=TRUE)





#Normalization with gene effective lengths (TPM)
normCriDGE <- scaleOffset(CriDGE, CriLength)
group <- as.factor(c("Leaves" ,"Leaves" , "Leaves" , "Root", "Root", "Root", "RootTip", "RootTip", "RootTip")) 
normCriDGE$samples$group <- group

norm_expressedCriGenes_list <- filterByExpr(normCriDGE, group=group)
norm_expressedCriGenes <- normCriDGE[norm_expressedCriGenes_list, , keep.lib.sizes=FALSE] 

normCriGenes_cpm <- cpm(norm_expressedCriGenes, offset = norm_expressedCriGenes$offset, log=FALSE)
normCriGenes_lcpm <- cpm(norm_expressedCriGenes, offset = norm_expressedCriGenes$offset, log=TRUE)





##Density plots of before and after filtering low expressed genes##
nsamples <- ncol(CriDGE)
col <- viridis_pal()(9)

tiff('densityPlot_Raw&Filtered_ExpressedCriGenes.tiff', units="in", width=10, height=5, res=600, compression = 'lzw')

par(mfrow=c(1,2))

plot(density(CriGenes_lcpm[,1]), col=col[1], lwd=2, xlim=c(-5,9), ylim=c(0,4), las=2, main="", xlab="")
title(main="A. Raw Data", xlab="Log-CPM")
abline(v=lcpm.cutoff, lty=3) 
for (i in 2:nsamples){
  den <- density(CriGenes_lcpm[,i])
  lines(den$x, den$y, col=col[i], lwd=2)
}
legend("topright", colnames(CriDGE), text.col=col, bty="n", y.intersp = 0.65)


expressedCri_lcpm <- cpm(expressedCriGenes, log=TRUE)
plot(density(expressedCri_lcpm[,1]), col=col[1], lwd=2, xlim=c(-5,9), ylim=c(0,0.5), las=2, main="", xlab="") 
title(main="B. Filtered Data", xlab="Log-CPM")
abline(v=lcpm.cutoff, lty=3)
for (i in 2:nsamples){
  den <- density(expressedCri_lcpm[,i])
  lines(den$x, den$y, col=col[i], lwd=2)
}
dev.off()





##Effects of the normalization
tiff('boxPlot_NormalizationCriGenes.tiff', units="in", width=10, height=5, res=600, compression = 'lzw')

par(mfrow=c(1,3))

boxplot(expressedCri_lcpm, las=2, col=col, ylim=c(-5,16), main="")
title(main="C. Unnormalised Data", ylab="Log-CPM")

boxplot(TMMCriGenes_lcpm, las=2, col=col, ylim=c(-5,16), main="")
title(main="D. Normalized Data TMM", ylab="Log-CPM")

boxplot(normCriGenes_lcpm, las=2, col=col, ylim=c(-5,16), main="")
title(main="E. Normalized Data Offset", ylab="Log-CPM")

dev.off()





##MDS to cluster samples##

col.group <- group
levels(col.group) <- c("#35b779", "#482878", "#fde725") 
col.group <- as.character(col.group)

shape.group <- group
levels(shape.group) <- (c(21, 22, 25))
shape.group <- as.numeric(as.character(shape.group))

samplenames <- c("S1","S2","S3","S4","S5","S6","S7","S8","S9")

par(mfrow=c(1,1))

tiff('MDS_WholeCriTranscriptome.tiff', units="in", width=5, height=5, res=600, compression = 'lzw')
MDS_WholeCriTranscriptome <- plotMDS(CriGenes_lcpm, top=7727, pch=shape.group, cex=4, col=col.group, bg=col.group, xlim=c(-4,4), ylim=c(-4,4), xlab="", ylab="")
text(MDS_WholeCriTranscriptome, labels=samplenames, cex= 1)
title(main="Complete Dataset", xlab="logFC2 DIM1", ylab="logFC2 DIM2", font.main=2, cex.main=2, font.lab=2, cex.lab=1.3)
legend("topright", levels(group), text.col=c("#35b779", "#482878", "#000000"), text.font=2, y.intersp = 1, 
       bty="n", pch=c(21, 22, 25), col=c("#35b779", "#482878", "#fde725"), pt.bg=c("#35b779", "#482878", "#fde725"), pt.cex = 2)
dev.off()



tiff('MDS_TMMCriTranscriptome.tiff', units="in", width=5, height=5, res=600, compression = 'lzw')
MDS_TMMCriTranscriptome <- plotMDS(TMMCriGenes_lcpm, top=1495, pch=shape.group, cex=4, col=col.group, bg=col.group, xlim=c(-4,4), ylim=c(-4,4), xlab="", ylab="")
text(MDS_TMMCriTranscriptome, labels=samplenames, cex= 1)
title(main="TMM Normalized Dataset", xlab="logFC2 DIM1", ylab="logFC2 DIM2", font.main=2, cex.main=2, font.lab=2, cex.lab=1.3)
legend("topright", levels(group), text.col=c("#35b779", "#482878", "#000000"), text.font=2, y.intersp = 1, 
       bty="n", pch=c(21, 22, 25), col=c("#35b779", "#482878", "#fde725"), pt.bg=c("#35b779", "#482878", "#fde725"), pt.cex = 2)
dev.off()



tiff('MDS_TPMCriTranscriptome.tiff', units="in", width=5, height=5, res=600, compression = 'lzw')
MDS_NormCriTranscriptome <- plotMDS(normCriGenes_lcpm, top=1495, pch=shape.group, cex=4, col=col.group, bg=col.group, xlim=c(-4,4), ylim=c(-4,4), xlab="", ylab="")
text(MDS_NormCriTranscriptome, labels=samplenames, cex= 1)
title(main="TPM Normalized Dataset", xlab="logFC2 DIM1", ylab="logFC2 DIM2", font.main=2, cex.main=2, font.lab=2, cex.lab=1.3)
legend("topright", levels(group), text.col=c("#35b779", "#482878", "#000000"), text.font=2, y.intersp = 1, 
       bty="n", pch=c(21, 22, 25), col=c("#35b779", "#482878", "#fde725"), pt.bg=c("#35b779", "#482878", "#fde725"), pt.cex = 2)
dev.off()





##Heatmap & dendogram of sample distances##
pseudocounts_allDGECounts <- log2(CriDGE$counts + 1)
allSAMPLEsDists <- dist(t(pseudocounts_allDGECounts))
allSAMPLEsDists_Matrix <- as.matrix(allSAMPLEsDists)


rownames(allSAMPLEsDists_Matrix) <- paste0("SAMPLE", 1:9, sep = " - ", CriDGE$samples$group)
colnames(allSAMPLEsDists_Matrix) <- NULL

heatColor <- viridis_pal(option = "A", direction = 1)(256)
par(mfrow=c(1,1))

tiff('HeatmapEuclideanDistances_Whole-Dataset.tiff', units="in", width=7, height=5, res=600, compression = 'lzw')
pheatmap(allSAMPLEsDists_Matrix,
         clustering_distance_rows = allSAMPLEsDists,
         clustering_distance_cols = allSAMPLEsDists,
         col = heatColor, 
         main ="Heatmap of Euclidean distances between samples\nWhole dataset")
dev.off()



pseudocounts_expDGECounts <- log2(norm_expressedCriGenes$counts + 1)
expSAMPLEsDists <- dist(t(pseudocounts_expDGECounts))
expSAMPLEsDists_Matrix <- as.matrix(expSAMPLEsDists)

rownames(expSAMPLEsDists_Matrix) <- paste0("SAMPLE", 1:9, sep = " - ", CriDGE$samples$group)
colnames(expSAMPLEsDists_Matrix) <- NULL

tiff('HeatmapEuclideanDistances_Normalized-Dataset.tiff', units="in", width=7, height=5, res=600, compression = 'lzw')
pheatmap(expSAMPLEsDists_Matrix,
         clustering_distance_rows = expSAMPLEsDists,
         clustering_distance_cols = expSAMPLEsDists,
         col = heatColor, main ="Heatmap of Euclidean distances between samples\nFiltered & normalized dataset")
dev.off()
