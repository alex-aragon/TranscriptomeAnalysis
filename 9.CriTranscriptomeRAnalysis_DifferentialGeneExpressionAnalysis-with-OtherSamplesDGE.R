setwd('/Users/alejandroaragon/Desktop')

library(edgeR)
library(limma)
library(tximport)
library(readr)
library(tidyverse)
library(dplyr)
library(colorspace)

library(viridis)
library(ggplot2)
#library(pheatmap)
library(trinotateR)
library(topGO)
library(UpSetR)


############## DIFFERENTIAL GENE EXPRESSION ANALYSIS with OTHER SAMPLES ##############
CriSAMPLES <- read.table(file.path("RNASeqAnalysis/5.Quantification/0.kallistoQuantification", "kallisto_OtherCriSamples.txt"), header = TRUE)
CriFiles <- file.path("RNASeqAnalysis/5.Quantification/0.kallistoQuantification/", CriSAMPLES$dir, "abundance.h5")
names(CriFiles) <- paste0("SAMPLE", 1:24)
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



group <- as.factor(c("Leaves" ,"Leaves" , "Leaves" , 
                     "Root", "Root", "Root", 
                     "RootTip", "RootTip", "RootTip",
                     "LeafTip", "LeafTip", "LeafTip",
                     "Root2", "Root2", "Root2",
                     "RootTip2", "RootTip2", "RootTip2",
                     "Shoot", "Shoot", "Shoot",
                     "Leaf2", "Leaf2", "Leaf2"))
                     #"Control", "Control",
                     #"24D1day", "24D1day",
                     #"24D2day", "24D2day")) 

CriDGE <- DGEList(CriCounts)
CriDGE$samples$group <- group



##Filtering low expressing genes##
#Without normalization
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
group <- as.factor(c("Leaves" ,"Leaves" , "Leaves" , 
                     "Root", "Root", "Root", 
                     "RootTip", "RootTip", "RootTip",
                     "LeafTip", "LeafTip", "LeafTip",
                     "Root2", "Root2", "Root2",
                     "RootTip2", "RootTip2", "RootTip2",
                     "Shoot", "Shoot", "Shoot",
                     "Leaves2", "Leaves2", "Leaves2"))
                     #"Control", "Control",
                     #"24D1day", "24D1day",
                     #"24D2day", "24D2day")) 
normCriDGE$samples$group <- group

norm_expressedCriGenes_list <- filterByExpr(normCriDGE, group=group)
norm_expressedCriGenes <- normCriDGE[norm_expressedCriGenes_list, , keep.lib.sizes=FALSE] 

normCriGenes_cpm <- cpm(norm_expressedCriGenes, offset = norm_expressedCriGenes$offset, log=FALSE)
normCriGenes_lcpm <- cpm(norm_expressedCriGenes, offset = norm_expressedCriGenes$offset, log=TRUE)


##Density plots of before and after filtering low expressed genes##
nsamples <- ncol(CriDGE)
col <- viridis_pal()(24)

tiff('densityPlot_Raw&Filtered_ExpressedCriGenes.tiff', units="in", width=10, height=5, res=600, compression = 'lzw')

par(mfrow=c(1,2))

plot(density(CriGenes_lcpm[,1]), col=col[1], lwd=2, xlim=c(-5,9), ylim=c(0,4), las=2, main="", xlab="")
title(main="A. Raw Data", xlab="Log-CPM")
abline(v=lcpm.cutoff, lty=3) 
for (i in 2:nsamples){
  den <- density(CriGenes_lcpm[,i])
  lines(den$x, den$y, col=col[i], lwd=2)
}
legend("topright", colnames(CriDGE), text.col=col, bty="n", cex = 0.7, y.intersp = 0.65)

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
tiff('boxPlot_NormalizationCriGenes.tiff', units="cm", width=40, height=10, res=600, compression = 'lzw')

par(mfrow=c(1,3))

boxplot(expressedCri_lcpm, las=2, col=col, ylim=c(-5,16), main="")
title(main="C. Unnormalised Data", ylab="Log-CPM")

boxplot(TMMCriGenes_lcpm, las=2, col=col, ylim=c(-5,16), main="")
title(main="D. Normalized Data TMM", ylab="Log-CPM")

boxplot(normCriGenes_lcpm, las=2, col=col, ylim=c(-5,16), main="")
title(main="E. Normalized Data Offset", ylab="Log-CPM")

dev.off()




col.group <- group
levels(col.group) <- c("#73d055", "#29af7f", "#3cbb75",
                       "#440154", "#482677",
                       "#fde725","#dce319", "#b8de29") 
                        #"#287d8e","#33638d", "#404788",
col.group <- as.character(col.group)

shape.group <- group
levels(shape.group) <- (c(24, 21, 21, 22, 22, 25, 25, 23))
shape.group <- as.numeric(as.character(shape.group))

samplenames <- c("S1","S2","S3","S4","S5","S6","S7","S8","S9","S10",
                 "S11","S12","S13","S14","S15","S16","S17","S18","S19","S20",
                 "S21","S22","S23","S24")

par(mfrow=c(1,1))

tiff('MDS_WholeCriTranscriptome.tiff', units="in", width=5, height=5, res=600, compression = 'lzw')
MDS_WholeCriTranscriptome <- plotMDS(CriGenes_lcpm, top=3500, pch=shape.group, cex=4, col=col.group, bg=col.group, xlim=c(-4,4), ylim=c(-4,4), xlab="", ylab="")
text(MDS_WholeCriTranscriptome, labels=samplenames, cex= 1)
title(main="Complete Dataset", xlab="logFC2 DIM1", ylab="logFC2 DIM2", font.main=2, cex.main=2, font.lab=2, cex.lab=1.3)
legend("topright", levels(group), 
       text.col=c("#73d055", "#29af7f", "#3cbb75", "#440154", "#482677","#000000","#000000", "#b8de29"), 
       text.font=2, y.intersp = 1, 
       bty="n", pch=c(24, 21, 21, 22, 22, 25, 25, 23), 
       col=c("#73d055", "#29af7f", "#3cbb75", "#440154", "#482677","#fde725","#dce319", "#b8de29"), 
       pt.bg=c("#73d055", "#29af7f", "#3cbb75", "#440154", "#482677","#fde725","#dce319", "#b8de29"), 
       pt.cex = 0.6, cex = 0.6)
dev.off()




tiff('MDS_TPMCriTranscriptome.tiff', units="in", width=5, height=5, res=600, compression = 'lzw')
MDS_NormCriTranscriptome <- plotMDS(normCriGenes_lcpm, top=1750, pch=shape.group, cex=4, col=col.group, bg=col.group, xlim=c(-4,4), ylim=c(-4,4), xlab="", ylab="")
text(MDS_NormCriTranscriptome, labels=samplenames, cex= 1)
title(main="TPM Normalized Dataset", xlab="logFC2 DIM1", ylab="logFC2 DIM2", font.main=2, cex.main=2, font.lab=2, cex.lab=1.3)
legend("bottomright", levels(group), 
       text.col=c("#73d055", "#29af7f", "#3cbb75", "#440154", "#482677","#000000","#000000", "#b8de29"), 
       text.font=2, y.intersp = 1, 
       bty="n", pch=c(24, 21, 21, 22, 22, 25, 25, 23), 
       col=c("#73d055", "#29af7f", "#3cbb75", "#440154", "#482677","#fde725","#dce319", "#b8de29"), 
       pt.bg=c("#73d055", "#29af7f", "#3cbb75", "#440154", "#482677","#fde725","#dce319", "#b8de29"), 
       pt.cex = 0.6, cex = 0.6)
dev.off()




batch_lcpm <- removeBatchEffect(CriGenes_lcpm, batch = c(1,1,1,1,1,1,1,1,1,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2))
normbatch_lcpm <- removeBatchEffect(normCriGenes_lcpm, batch = c(1,1,1,1,1,1,1,1,1,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2))

tiff('MDS_BATCHCriTranscriptome.tiff', units="in", width=5, height=5, res=600, compression = 'lzw')
MDS_NormCriTranscriptome <- plotMDS(batch_lcpm, top=1750, pch=shape.group, cex=4, col=col.group, bg=col.group, xlim=c(-5,5), ylim=c(-4,4), xlab="", ylab="")
text(MDS_NormCriTranscriptome, labels=samplenames, cex= 1)
title(main="Removing Batch Effect", xlab="logFC2 DIM1", ylab="logFC2 DIM2", font.main=2, cex.main=2, font.lab=2, cex.lab=1.3)
legend("bottomright", levels(group), 
       text.col=c("#73d055", "#29af7f", "#3cbb75", "#440154", "#482677","#000000","#000000", "#b8de29"), 
       text.font=2, y.intersp = 1, 
       bty="n", pch=c(24, 21, 21, 22, 22, 25, 25, 23), 
       col=c("#73d055", "#29af7f", "#3cbb75", "#440154", "#482677","#fde725","#dce319", "#b8de29"), 
       pt.bg=c("#73d055", "#29af7f", "#3cbb75", "#440154", "#482677","#fde725","#dce319", "#b8de29"), 
       pt.cex = 0.6, cex = 0.6)
dev.off()


tiff('MDS_BATCHNormalizedCriTranscriptome.tiff', units="in", width=5, height=5, res=600, compression = 'lzw')
MDS_NormCriTranscriptome <- plotMDS(normbatch_lcpm, top=1750, pch=shape.group, cex=4, col=col.group, bg=col.group, xlim=c(-5,5), ylim=c(-4,4), xlab="", ylab="")
text(MDS_NormCriTranscriptome, labels=samplenames, cex= 1)
title(main="Removing Batch Effect", xlab="logFC2 DIM1", ylab="logFC2 DIM2", font.main=2, cex.main=2, font.lab=2, cex.lab=1.3)
legend("bottomright", levels(group), 
       text.col=c("#73d055", "#29af7f", "#3cbb75", "#440154", "#482677","#000000","#000000", "#b8de29"), 
       text.font=2, y.intersp = 1, 
       bty="n", pch=c(24, 21, 21, 22, 22, 25, 25, 23), 
       col=c("#73d055", "#29af7f", "#3cbb75", "#440154", "#482677","#fde725","#dce319", "#b8de29"), 
       pt.bg=c("#73d055", "#29af7f", "#3cbb75", "#440154", "#482677","#fde725","#dce319", "#b8de29"), 
       pt.cex = 0.6, cex = 0.6)
dev.off()

##Heatmap & dendogram of sample distances##
pseudocounts_allDGECounts <- log2(CriDGE$counts + 1)
allSAMPLEsDists <- dist(t(pseudocounts_allDGECounts))
allSAMPLEsDists_Matrix <- as.matrix(allSAMPLEsDists)


rownames(allSAMPLEsDists_Matrix) <- paste0("SAMPLE", 1:24, sep = " - ", CriDGE$samples$group)
colnames(allSAMPLEsDists_Matrix) <- NULL

heatColor <- viridis_pal(option = "A", direction = 1)(256)
par(mfrow=c(1,1))

tiff('HeatmapEuclideanDistances_Whole-Dataset.tiff', units="in", width=7, height=5, res=600, compression = 'lzw')
pheatmap(allSAMPLEsDists_Matrix,
         clustering_distance_rows = allSAMPLEsDists,
         clustering_distance_cols = allSAMPLEsDists,
         col = heatColor, 
         main = "Heatmap of Euclidean distances between samples\nWhole dataset")
dev.off()
dev.off()



pseudocounts_expDGECounts <- log2(norm_expressedCriGenes$counts + 1)
expSAMPLEsDists <- dist(t(pseudocounts_expDGECounts))
expSAMPLEsDists_Matrix <- as.matrix(expSAMPLEsDists)

rownames(expSAMPLEsDists_Matrix) <- paste0("SAMPLE", 1:24, sep = " - ", CriDGE$samples$group)
colnames(expSAMPLEsDists_Matrix) <- NULL

tiff('HeatmapEuclideanDistances_Normalized-Dataset.tiff', units="in", width=7, height=5, res=600, compression = 'lzw')
pheatmap(expSAMPLEsDists_Matrix,
         clustering_distance_rows = expSAMPLEsDists,
         clustering_distance_cols = expSAMPLEsDists,
         col = heatColor, main ="Heatmap of Euclidean distances between samples\nFiltered & normalized dataset")
dev.off()
dev.off()










########### ANALYZING THE DATASET ###########

group <- as.factor(c("Leaves" ,"Leaves" , "Leaves" , 
                     "Root", "Root", "Root", 
                     "RootTip", "RootTip", "RootTip",
                     "LeafTip", "LeafTip", "LeafTip",
                     "Root2", "Root2", "Root2",
                     "RootTip2", "RootTip2", "RootTip2",
                     "Shoot", "Shoot", "Shoot",
                     "Leaves2", "Leaves2", "Leaves2"))
CriDGE <- DGEList(CriCounts)
CriDGE$samples$group <- group

#CriDGE$samples$batch <- batch
#batch <- c(0,0,0,0,0,0,0,0,0,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1) #,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2))

CriTdesign <- model.matrix(~0 + group)
colnames(CriTdesign) <- gsub("group", "", colnames(CriTdesign))

CriTdesign

contr.matrix <- makeContrasts(
  TheRootTip1 = RootTip,
  TheRootTip2 = RootTip2,
  TheShoot = Shoot,
  TheLeafTip = LeafTip,
  TheLeaves1 = Leaves,
  TheLeaves2 = Leaves2,
  TheRoot1 = Root,
  TheRoot2 = Root2,
  ARootTip = RootTip - RootTip2,
  ARootTip2 = RootTip2 - RootTip,
  ARoot = Root - Root2,
  ALeaves = Leaves - Leaves2,
  RootTip_vs_LeafTip = "((RootTip+RootTip2)/2) - LeafTip",
  RootTip_vs_Shoot = "((RootTip+RootTip2)/2) - Shoot",
  LeafTip_vs_Shoot = LeafTip - Shoot,
  RootTip_vs_Leaves = RootTip - Leaves,
  OnlyShoot = "Shoot - ((((RootTip+RootTip2)/2)+LeafTip)/2)",
  levels = colnames(CriTdesign))
contr.matrix

voom_normCriGenes <- voom(norm_expressedCriGenes, CriTdesign, plot = TRUE)
vExCriTfit <- lmFit(voom_normCriGenes, CriTdesign)
eBvExCriTfit <- eBayes(vExCriTfit)

resultsCriDGE <- decideTests(eBvExCriTfit, adjust.method = "BH", p.value = 0.05)
resultsCriDGE <- resultsCriDGE %>%
  dplyr::filter(logFC!=0)
write.fit(eBvExCriTfit, resultsCriDGE, file="Ceratopteris-richardii_Transcriptome_ResultsDGE_OtherSamples.txt")

versus_vExCriTfit <- contrasts.fit(vExCriTfit, contrasts=contr.matrix) 
versus_eBvExCriTfit <- eBayes(versus_vExCriTfit)
contrast_resultsCriDGE <- topTable(versus_eBvExCriTfit, coef=NULL, n=Inf, p = 0.05, adjust.method = "BH")

tiff('RootTip_VennDiagram.tiff', units="in", width=6.7, height=5, res=600, compression = 'lzw')
vennDiagram(resultsCriDGE[,6:7], 
            include=c("up","down"), 
            mar = c(0,4,0,4), 
            cex=1,
            lwd = 2,
            names = c("Aragón-Raygoza","Yu et al, 2020"),
            counts.col = c("#df488dff","#95d840ff"),
            circle.col = c("#481567ff","#fde725ff"),
            main = "\n\n\n\n\n\nUp & down regulated transcripts\nbetween Root Tip datasets")
dev.off()

vennCounts(resultsCriDGE[,6:7])
TheRootTip1_ResultsCriDGE <- topTable(versus_eBvExCriTfit, coef=1, n=Inf, p = 0.05, adjust.method = "BH")
TheRootTip2_ResultsCriDGE <- topTable(versus_eBvExCriTfit, coef=2, n=Inf, p = 0.05, adjust.method = "BH")
TheShoot_ResultsCriDGE <- topTable(versus_eBvExCriTfit, coef=3, n=Inf, p = 0.05, adjust.method = "BH")
TheLeafTip_ResultsCriDGE <- topTable(versus_eBvExCriTfit, coef=4, n=Inf, p = 0.05, adjust.method = "BH")
TheLeaves1_ResultsCriDGE <- topTable(versus_eBvExCriTfit, coef=5, n=Inf, p = 0.05, adjust.method = "BH")
TheLeaves2_ResultsCriDGE <- topTable(versus_eBvExCriTfit, coef=6, n=Inf, p = 0.05, adjust.method = "BH")
TheRoot1_ResultsCriDGE <- topTable(versus_eBvExCriTfit, coef=7, n=Inf, p = 0.05, adjust.method = "BH")
TheRoot2_ResultsCriDGE <- topTable(versus_eBvExCriTfit, coef=8, n=Inf, p = 0.05, adjust.method = "BH")



TheRootTip1_ResultsCriDGE <- TheRootTip1_ResultsCriDGE %>% 
  rownames_to_column(var = "ID")
TheRootTip2_ResultsCriDGE <- TheRootTip2_ResultsCriDGE %>% 
  rownames_to_column(var = "ID")
TheRoot1_ResultsCriDGE <- TheRoot1_ResultsCriDGE %>% 
  rownames_to_column(var = "ID")
TheRoot2_ResultsCriDGE <- TheRoot2_ResultsCriDGE %>% 
  rownames_to_column(var = "ID")
TheLeaves1_ResultsCriDGE <- TheLeaves1_ResultsCriDGE %>% 
  rownames_to_column(var = "ID")
TheLeaves2_ResultsCriDGE <- TheLeaves2_ResultsCriDGE %>% 
  rownames_to_column(var = "ID")
TheLeafTip_ResultsCriDGE <- TheLeafTip_ResultsCriDGE %>% 
  rownames_to_column(var = "ID")
TheShoot_ResultsCriDGE <- TheShoot_ResultsCriDGE %>% 
  rownames_to_column(var = "ID")




########### Upset Plots ###########
UpList <- list ("Root Tip" = TheRootTip1_ResultsCriDGE$ID, 
                #"Leaves and Shoot" = TheLeaves1_ResultsCriDGE$ID,
                "Root Tip 2" = TheRootTip2_ResultsCriDGE$ID,
                #"Differentiated Root 2" = TheRoot2_ResultsCriDGE$ID,
                "Leaf Tip" = TheLeafTip_ResultsCriDGE$ID,
                #"Leaves" = TheLeaves2_ResultsCriDGE$ID,
                "Shoot" = TheShoot_ResultsCriDGE$ID)
                #"Differentiated Root" = TheRoot1_ResultsCriDGE$ID)
tmp <- UpSetR::fromList(UpList)
tiff('CriDGE_UpSetPlot.tiff', units="in", width=10, height=5, res=600, compression = 'lzw')
UpSetR::upset(UpSetR::fromList(UpList),
              nsets = 8,
              nintersects = 200,
              keep.order = F,
              order.by = c("degree"), 
              decreasing = c(F),
              group.by = "degree",
              mainbar.y.label = "Number of Transcripts\n(Intersection)",
              sets.x.label = "Total of Transcripts",
              point.size = 3,
              text.scale = 2,
              queries = list(
                list(
                  query = intersects,
                  params = list("Root Tip","Root Tip 2","Leaf Tip", "Shoot"), 
                  color = "#482677ff", 
                  active = T),
                list(
                  query = intersects,
                  params = list("Root Tip","Root Tip 2","Leaf Tip"), 
                  color = "#238a8dff", 
                  active = T),
                list(
                  query = intersects,
                  params = list("Root Tip","Root Tip 2", "Shoot"), 
                  color = "#238a8dff", 
                  active = T),
                list(
                  query = intersects,
                  params = list("Root Tip","Leaf Tip", "Shoot"), 
                  color = "#238a8dff", 
                  active = T),
                list(
                  query = intersects,
                  params = list("Root Tip 2","Leaf Tip", "Shoot"), 
                  color = "#238a8dff", 
                  active = T),
                list(
                  query = intersects,
                  params = list("Root Tip","Shoot"), 
                  color = "#73d055ff", 
                  active = T),
                list(
                  query = intersects,
                  params = list("Root Tip","Root Tip 2"), 
                  color = "#73d055ff", 
                  active = T),
                list(
                  query = intersects,
                  params = list("Root Tip","Leaf Tip"), 
                  color = "#73d055ff", 
                  active = T),
                list(
                  query = intersects,
                  params = list("Root Tip 2","Shoot"), 
                  color = "#73d055ff", 
                  active = T),
                list(
                  query = intersects,
                  params = list("Root Tip 2","Leaf Tip"), 
                  color = "#73d055ff", 
                  active = T),
                list(
                  query = intersects,
                  params = list("Shoot","Leaf Tip"), 
                  color = "#73d055ff", 
                  active = T),
                list(
                  query = intersects,
                  params = list("Root Tip"), 
                  color = "#b8de29ff", 
                  active = T),
                list(
                  query = intersects,
                  params = list("Shoot"), 
                  color = "#b8de29ff", 
                  active = T),
                list(
                  query = intersects,
                  params = list("Root Tip 2"), 
                  color = "#b8de29ff", 
                  active = T),
                list(
                  query = intersects,
                  params = list("Leaf Tip"), 
                  color = "#b8de29ff", 
                  active = T)
              )
)
dev.off()





UpList <- list ("Root Tip" = TheRootTip1_ResultsCriDGE$ID, 
                "Root Tip 2" = TheRootTip2_ResultsCriDGE$ID,
                "Differentiated Root" = TheRoot2_ResultsCriDGE$ID,
                "Differentiated Root 2" = TheRoot1_ResultsCriDGE$ID)

tiff('CriDGE_UpSetPlot2.tiff', units="in", width=10, height=5, res=600, compression = 'lzw')
UpSetR::upset(UpSetR::fromList(UpList),
              nsets = 8,
              nintersects = 200,
              keep.order = F,
              order.by = c("degree"), 
              decreasing = c(F),
              group.by = "degree",
              mainbar.y.label = "Number of Transcripts\n(Intersection)",
              sets.x.label = "Total of Transcripts",
              point.size = 3,
              queries = list(
                list(
                  query = intersects,
                  params = list("Root Tip","Root Tip 2","Differentiated Root", "Differentiated Root 2"), 
                  color = "#482677ff", 
                  active = T),
                list(
                  query = intersects,
                  params = list("Root Tip","Root Tip 2","Differentiated Root"), 
                  color = "#238a8dff", 
                  active = T),
                list(
                  query = intersects,
                  params = list("Root Tip","Root Tip 2", "Differentiated Root 2"), 
                  color = "#238a8dff", 
                  active = T),
                list(
                  query = intersects,
                  params = list("Root Tip","Differentiated Root", "Differentiated Root 2"), 
                  color = "#238a8dff", 
                  active = T),
                list(
                  query = intersects,
                  params = list("Root Tip 2","Differentiated Root", "Differentiated Root 2"), 
                  color = "#238a8dff", 
                  active = T),
                list(
                  query = intersects,
                  params = list("Root Tip","Differentiated Root"), 
                  color = "#73d055ff", 
                  active = T),
                list(
                  query = intersects,
                  params = list("Root Tip","Root Tip 2"), 
                  color = "#73d055ff", 
                  active = T),
                list(
                  query = intersects,
                  params = list("Root Tip","Differentiated Root 2"), 
                  color = "#73d055ff", 
                  active = T),
                list(
                  query = intersects,
                  params = list("Root Tip 2","Differentiated Root"), 
                  color = "#73d055ff", 
                  active = T),
                list(
                  query = intersects,
                  params = list("Root Tip 2","Differentiated Root 2"), 
                  color = "#73d055ff", 
                  active = T),
                list(
                  query = intersects,
                  params = list("Differentiated Root","Differentiated Root 2"), 
                  color = "#73d055ff", 
                  active = T),
                list(
                  query = intersects,
                  params = list("Root Tip"), 
                  color = "#b8de29ff", 
                  active = T),
                list(
                  query = intersects,
                  params = list("Differentiated Root"), 
                  color = "#b8de29ff", 
                  active = T),
                list(
                  query = intersects,
                  params = list("Root Tip 2"), 
                  color = "#b8de29ff", 
                  active = T),
                list(
                  query = intersects,
                  params = list("Differentiated Root 2"), 
                  color = "#b8de29ff", 
                  active = T)
              )
)
dev.off()







UpList <- list ("Leaves" = TheLeaves1_ResultsCriDGE$ID,
                "Leaves 2" = TheLeaves2_ResultsCriDGE$ID,
                "Leaf Tip" = TheLeafTip_ResultsCriDGE$ID,
                "Shoot" = TheShoot_ResultsCriDGE$ID)
#"Differentiated Root" = TheRoot1_ResultsCriDGE$ID)

tiff('CriDGE_UpSetPlot3.tiff', units="in", width=10, height=5, res=600, compression = 'lzw')
UpSetR::upset(UpSetR::fromList(UpList),
              nsets = 8,
              nintersects = 200,
              keep.order = F,
              order.by = c("degree"), 
              decreasing = c(F),
              group.by = "degree",
              mainbar.y.label = "Number of Transcripts\n(Intersection)",
              sets.x.label = "Total of Transcripts",
              point.size = 3,
              queries = list(
                list(
                  query = intersects,
                  params = list("Leaves","Leaves 2","Leaf Tip", "Shoot"), 
                  color = "#482677ff", 
                  active = T),
                list(
                  query = intersects,
                  params = list("Leaves","Leaves 2","Leaf Tip"), 
                  color = "#238a8dff", 
                  active = T),
                list(
                  query = intersects,
                  params = list("Leaves","Leaves 2","Shoot"), 
                  color = "#238a8dff", 
                  active = T),
                list(
                  query = intersects,
                  params = list("Leaves","Leaf Tip", "Shoot"), 
                  color = "#238a8dff", 
                  active = T),
                list(
                  query = intersects,
                  params = list("Leaves 2","Leaf Tip", "Shoot"), 
                  color = "#238a8dff", 
                  active = T),
                list(
                  query = intersects,
                  params = list("Leaves","Leaves 2"), 
                  color = "#73d055ff", 
                  active = T),
                list(
                  query = intersects,
                  params = list("Leaves","Leaf Tip"), 
                  color = "#73d055ff", 
                  active = T),
                list(
                  query = intersects,
                  params = list("Leaves","Shoot"), 
                  color = "#73d055ff", 
                  active = T),
                list(
                  query = intersects,
                  params = list("Leaves 2","Leaf Tip"), 
                  color = "#73d055ff", 
                  active = T),
                list(
                  query = intersects,
                  params = list("Leaves 2","Shoot"), 
                  color = "#73d055ff", 
                  active = T),
                list(
                  query = intersects,
                  params = list("Leaf Tip", "Shoot"), 
                  color = "#73d055ff", 
                  active = T),
                list(
                  query = intersects,
                  params = list("Leaves"), 
                  color = "#b8de29ff", 
                  active = T),
                list(
                  query = intersects,
                  params = list("Shoot"), 
                  color = "#b8de29ff", 
                  active = T),
                list(
                  query = intersects,
                  params = list("Leaves 2"), 
                  color = "#b8de29ff", 
                  active = T),
                list(
                  query = intersects,
                  params = list("Leaf Tip"), 
                  color = "#b8de29ff", 
                  active = T)
              )
)
dev.off()








############# Gene Ontology Analysis ############
ARootTip_ResultsCriDGE <- topTable(versus_eBvExCriTfit, coef=9, n=Inf, p = 0.05, adjust.method = "BH")
ARootTip2_ResultsCriDGE <- topTable(versus_eBvExCriTfit, coef=10, n=Inf, p = 0.05, adjust.method = "BH")
ARoot_ResultsCriDGE <- topTable(versus_eBvExCriTfit, coef=11, n=Inf, p = 0.05, adjust.method = "BH")
ALeaf_ResultsCriDGE <- topTable(versus_eBvExCriTfit, coef=12, n=Inf, p = 0.05, adjust.method = "BH")
RootTip_vs_LeafTip_ResultsCriDGE <- topTable(versus_eBvExCriTfit, coef=13, n=Inf, p = 0.05, adjust.method = "BH")
RootTip_vs_Shoot_ResultsCriDGE <- topTable(versus_eBvExCriTfit, coef=14, n=Inf, p = 0.05, adjust.method = "BH")
LeafTip_vs_Shoot_ResultsCriDGE <- topTable(versus_eBvExCriTfit, coef=15, n=Inf, p = 0.05, adjust.method = "BH")
RootTip_vs_Leaves_ResultsCriDGE <- topTable(versus_eBvExCriTfit, coef=16, n=Inf, p = 0.05, adjust.method = "BH")
Shoot_ResultsCriDGE <- topTable(versus_eBvExCriTfit, coef=17, n=Inf, p = 0.05, adjust.method = "BH")

library(trinotateR)
library(topGO)

TrinotateAnnotation <- read_trinotate("1.CeratopterisTranscriptomeAnnotation-v2.xls")
TrinotateAnnotation <- TrinotateAnnotation %>%
  distinct(gene_id, .keep_all = TRUE)
summary_trinotate(TrinotateAnnotation)

GOx_CriAnno <- split_GO(TrinotateAnnotation, hit = "gene_ontology_BLASTP")
summary_GO(GOx_CriAnno)

BPGOx_CriAnno <- GOx_CriAnno %>%
  filter(ontology=="biological_process") %>%
  dplyr::select(gene,go) %>%
  dplyr::rename("GO:BP"="go") 
BPGOx_CriAnno_list <- split(BPGOx_CriAnno$`GO:BP`, BPGOx_CriAnno$gene)

########Shared Meristematic Genes #########
TheShoot_ResultsCriDGE <- TheShoot_ResultsCriDGE %>% 
  rownames_to_column(var = "ID")

TheShoot_ResultsCriDGE$Presence <- ifelse(TheShoot_ResultsCriDGE$ID%in%TheRootTip1_ResultsCriDGE$ID, as.character("Shared"), 
                                          ifelse(TheShoot_ResultsCriDGE$ID%in%TheRootTip1_ResultsCriDGE$ID, as.character("Shared"), 
                                                 ifelse(TheShoot_ResultsCriDGE$ID%in%TheLeafTip_ResultsCriDGE$ID, as.character("Shared"), as.character("Unique"))))

Shoot_SharedExp <- TheShoot_ResultsCriDGE %>%
  filter(Presence=="Shared") %>%
  dplyr::rename("gene_id" = "ID") 

Shoot_SharedExp <- Shoot_SharedExp %>%
  left_join(., TrinotateAnnotation, by = "gene_id")

Shoot_SharedExp <- Shoot_SharedExp %>%
  dplyr::select(gene_id, logFC, AveExpr, P.Value, adj.P.Val, Presence, Crichardii_Transcriptome_676_v2_BLASTX, Arabidopsis_Araport11_BLASTP, Azolla_FernBase_BLASTP, Salvinia_FernBase_BLASTP, gene_ontology_BLASTX)

Meristem_SharedUpGOs <- Shoot_SharedExp %>%
  filter(logFC>=2) %>%
  dplyr::rename("CriID" = "gene_id") 

MeristemUpSharedforGOs <- factor(as.integer(TrinotateAnnotation$gene_id %in% Meristem_SharedUpGOs$CriID))
names(MeristemUpSharedforGOs) <- TrinotateAnnotation$gene_id
head(MeristemUpSharedforGOs)

BP_GOdata <- new("topGOdata", ontology = "BP", allGenes = MeristemUpSharedforGOs,
                 annot = annFUN.gene2GO, gene2GO = BPGOx_CriAnno_list)

BP_Meristem_statTest <- runTest(BP_GOdata, algorithm = "elim", statistic = "fisher" )

BPUp_Meristem_resultTopGO <- GenTable(BP_GOdata, raw.P.value = BP_Meristem_statTest,
                                        orderBy = "raw.P.value", 
                                        topNodes = 30,
                                        numChar=1000)

BPUp_Meristem_resultTopGO <- cbind(BPUp_Meristem_resultTopGO, GeneRatio = BPUp_Meristem_resultTopGO$Significant/(as.numeric(nrow(Meristem_SharedUpGOs))))
BPUp_Meristem_resultTopGO$Name <- paste(BPUp_Meristem_resultTopGO$GO.ID, "-", BPUp_Meristem_resultTopGO$Term)

plot_Meristem_BP_SharedExpressedGOs <- ggplot(BPUp_Meristem_resultTopGO, aes(y=reorder(Name, GeneRatio), x = GeneRatio, size=Significant, color=-log10(as.numeric(raw.P.value)))) +
  geom_point(alpha=1) +
  scale_colour_viridis(option="D", limits = c(0,30)) +
  scale_size_continuous(limits = c(1, 1500), range = c(1,8)) +
  theme_minimal() +
  theme(panel.grid.major = element_line(size = 0.3, linetype = "dotted", colour = "gray", lineend = "butt")) + 
  theme(axis.text.x = element_text(angle = 45, hjust = 1, size = 8), axis.text.y = element_text(vjust = 1, size = 8)) +
  labs(size = "Gene Counts", color = "-log10(FDR)") +
  theme(plot.title = element_text(face = "bold")) +
  ylab("Gene Onthology\nCellular Component") +
  xlab("\nGene Ratio") +
  xlim(0,0.08) +
  theme(axis.title.x = element_text(face="bold"), axis.title.y = element_text(face="bold"), )

tiff('MeristemSharedUpGOs.tiff', units="cm", width=18, height=22, res=600, compression = 'lzw')
plot_Meristem_BP_SharedExpressedGOs
dev.off()





######## Unique Shoot Genes #########
Shoot_UniqueExp <- TheShoot_ResultsCriDGE %>%
  filter(Presence=="Unique") %>%
  dplyr::rename("gene_id" = "ID") 

Shoot_UniqueExp <- Shoot_UniqueExp %>%
  left_join(., TrinotateAnnotation, by = "gene_id")

Shoot_UniqueExp <- Shoot_UniqueExp %>%
  dplyr::select(gene_id, logFC, AveExpr, P.Value, adj.P.Val, Presence, Crichardii_Transcriptome_676_v2_BLASTX, Arabidopsis_Araport11_BLASTP, Azolla_FernBase_BLASTP, Salvinia_FernBase_BLASTP, gene_ontology_BLASTX)

Shoot_UniqueExp <- Shoot_UniqueExp %>%
  filter(logFC>=2) %>%
  dplyr::rename("CriID" = "gene_id") 

ShootUpUniqueforGOs <- factor(as.integer(TrinotateAnnotation$gene_id %in% Shoot_UniqueExp$CriID))
names(ShootUpUniqueforGOs) <- TrinotateAnnotation$gene_id
head(ShootUpUniqueforGOs)

BP_GOdata <- new("topGOdata", ontology = "BP", allGenes = ShootUpUniqueforGOs,
                 annot = annFUN.gene2GO, gene2GO = BPGOx_CriAnno_list)

BP_UniqueShoot_statTest <- runTest(BP_GOdata, algorithm = "elim", statistic = "fisher" )

BPUp_UniqueShoot_resultTopGO <- GenTable(BP_GOdata, raw.P.value = BP_UniqueShoot_statTest,
                                      orderBy = "raw.P.value", 
                                      topNodes = 15,
                                      numChar=1000)

BPUp_UniqueShoot_resultTopGO <- cbind (BPUp_UniqueShoot_resultTopGO, Dataset = "Shoot")
BPUp_UniqueShoot_resultTopGO <- cbind(BPUp_UniqueShoot_resultTopGO, GeneRatio = BPUp_UniqueShoot_resultTopGO$Significant/(as.numeric(nrow(Shoot_UniqueExp))))
BPUp_UniqueShoot_resultTopGO$Name <- paste(BPUp_UniqueShoot_resultTopGO$GO.ID, "-", BPUp_UniqueShoot_resultTopGO$Term)





######## Unique Leaf Tip Genes #########
TheLeafTip_ResultsCriDGE <- TheLeafTip_ResultsCriDGE %>% 
  rownames_to_column(var = "ID")

TheLeafTip_ResultsCriDGE$Presence <- ifelse(TheLeafTip_ResultsCriDGE$ID%in%TheRootTip1_ResultsCriDGE$ID, as.character("Shared"), 
                                          ifelse(TheLeafTip_ResultsCriDGE$ID%in%TheRootTip1_ResultsCriDGE$ID, as.character("Shared"), 
                                                 ifelse(TheLeafTip_ResultsCriDGE$ID%in%TheShoot_ResultsCriDGE$ID, as.character("Shared"), as.character("Unique"))))

LeafTip_UniqueExp <- TheLeafTip_ResultsCriDGE %>%
  filter(Presence=="Unique") %>%
  dplyr::rename("gene_id" = "ID") 

LeafTip_UniqueExp <- LeafTip_UniqueExp %>%
  left_join(., TrinotateAnnotation, by = "gene_id")

LeafTip_UniqueExp <- LeafTip_UniqueExp %>%
  dplyr::select(gene_id, logFC, AveExpr, P.Value, adj.P.Val, Presence, Crichardii_Transcriptome_676_v2_BLASTX, Arabidopsis_Araport11_BLASTP, Azolla_FernBase_BLASTP, Salvinia_FernBase_BLASTP, gene_ontology_BLASTX)

LeafTip_UniqueUpGOs <- LeafTip_UniqueExp %>%
  filter(logFC>=2) %>%
  dplyr::rename("CriID" = "gene_id") 


LeafTipUpUniqueforGOs <- factor(as.integer(TrinotateAnnotation$gene_id %in% LeafTip_UniqueUpGOs$CriID))
names(LeafTipUpUniqueforGOs) <- TrinotateAnnotation$gene_id
head(LeafTipUpUniqueforGOs)

BP_GOdata <- new("topGOdata", ontology = "BP", allGenes = LeafTipUpUniqueforGOs,
                 annot = annFUN.gene2GO, gene2GO = BPGOx_CriAnno_list)

BP_UniqueLeaf_statTest <- runTest(BP_GOdata, algorithm = "elim", statistic = "fisher" )

BPUp_UniqueLeaf_resultTopGO <- GenTable(BP_GOdata, raw.P.value = BP_UniqueLeaf_statTest,
                                         orderBy = "raw.P.value", 
                                         topNodes = 15,
                                         numChar=1000)

BPUp_UniqueLeaf_resultTopGO <- cbind (BPUp_UniqueLeaf_resultTopGO, Dataset = "Leaf")
BPUp_UniqueLeaf_resultTopGO <- cbind(BPUp_UniqueLeaf_resultTopGO, GeneRatio = BPUp_UniqueLeaf_resultTopGO$Significant/(as.numeric(nrow(LeafTip_UniqueUpGOs))))
BPUp_UniqueLeaf_resultTopGO$Name <- paste(BPUp_UniqueLeaf_resultTopGO$GO.ID, "-", BPUp_UniqueLeaf_resultTopGO$Term)





BP_UpperMeristem_resultTopGO <- full_join(BPUp_UniqueShoot_resultTopGO,BPUp_UniqueLeaf_resultTopGO)
BP_UpperMeristem_resultTopGO$Name <- factor(BP_UpperMeristem_resultTopGO$Name, levels = BP_UpperMeristem_resultTopGO$Name[order(BP_UpperMeristem_resultTopGO$Dataset)])

plot_UpperMeristem_BP_UniqueUpGOs <- ggplot(BP_UpperMeristem_resultTopGO, aes(x = Dataset, y = Name, size=Significant, color=-log10(as.numeric(raw.P.value)))) +
  geom_point(alpha=1) +
  scale_colour_viridis(option="D", limits = c(0,30)) +
  scale_size_continuous(limits = c(1, 10), range = c(1,8), labels = c(2,4,6,8,10)) +
  theme_minimal() +
  theme(panel.grid.major = element_line(size = 0.3, linetype = "dotted", colour = "gray", lineend = "butt")) + 
  theme(axis.text.y = element_text(angle =0, hjust = 1, size = 8), axis.text.x = element_text(angle =90, hjust = 1)) +
  labs(size = "Gene Counts", color = "-log10(FDR)") +
  theme(plot.title = element_text(face = "bold")) +
  ylab("Gene Onthology\nBiological Processes") +
  xlab("\nDataset") +
  theme(axis.title.x = element_text(face="bold"), axis.title.y = element_text(face="bold"), )

tiff('UpperMeristem_BP_UniqueUpGOs.tiff', units="cm", width=15, height=20, res=600, compression = 'lzw')
plot_UpperMeristem_BP_UniqueUpGOs
dev.off()





######## Differentialed Expressed Between Root Tip Datasets #######
RootTip_GOs <- ARootTip_ResultsCriDGE %>%
  filter(logFC>=2) %>%
  rownames_to_column(var = "CriID")

RootTip_GOs_list <- factor(as.integer(TrinotateAnnotation$gene_id %in% RootTip_GOs$CriID))
names(RootTip_GOs_list) <- TrinotateAnnotation$gene_id
head(RootTip_GOs_list)

BP_GOdata <- new("topGOdata", ontology = "BP", allGenes = RootTip_GOs_list,
                 annot = annFUN.gene2GO, gene2GO = BPGOx_CriAnno_list)

BP_RootTip_statTest <- runTest(BP_GOdata, algorithm = "elim", statistic = "fisher" )

BPUp_DERootTip1_resultTopGO <- GenTable(BP_GOdata, raw.P.value = BP_RootTip_statTest,
                                     orderBy = "raw.P.value", 
                                     topNodes = 13,
                                     numChar=1000)

BPUp_DERootTip1_resultTopGO <- BPUp_DERootTip1_resultTopGO[-c(2),]
BPUp_DERootTip1_resultTopGO <- cbind (BPUp_DERootTip1_resultTopGO, Dataset = "Aragón-Raygoza")
BPUp_DERootTip1_resultTopGO <- cbind(BPUp_DERootTip1_resultTopGO, GeneRatio = BPUp_DERootTip1_resultTopGO$Significant/(as.numeric(nrow(RootTip_GOs))))





RootTip2_GOs <- ARootTip2_ResultsCriDGE %>%
  filter(logFC>=2) %>%
  rownames_to_column(var = "CriID")

RootTip2_GOs_list <- factor(as.integer(TrinotateAnnotation$gene_id %in% RootTip2_GOs$CriID))
names(RootTip2_GOs_list) <- TrinotateAnnotation$gene_id
head(RootTip2_GOs_list)

BP_GOdata <- new("topGOdata", ontology = "BP", allGenes = RootTip2_GOs_list,
                 annot = annFUN.gene2GO, gene2GO = BPGOx_CriAnno_list)

BP_RootTip_statTest <- runTest(BP_GOdata, algorithm = "elim", statistic = "fisher" )

BPUp_DERootTip2_resultTopGO <- GenTable(BP_GOdata, raw.P.value = BP_RootTip_statTest,
                                     orderBy = "raw.P.value", 
                                     topNodes = 12,
                                     numChar=1000)

BPUp_DERootTip2_resultTopGO <- cbind(BPUp_DERootTip2_resultTopGO, GeneRatio = BPUp_DERootTip2_resultTopGO$Significant/(as.numeric(nrow(RootTip2_GOs))))
BPUp_DERootTip2_resultTopGO <- cbind (BPUp_DERootTip2_resultTopGO, Dataset = "Yu et al, 2020")





BP_RootTipBoth_resultTopGO <- full_join(BPUp_DERootTip1_resultTopGO,BPUp_DERootTip2_resultTopGO)
BP_RootTipBoth_resultTopGO$Term <- gsub(" involved in multidimensional cell growth", "", BP_RootTipBoth_resultTopGO$Term)
BP_RootTipBoth_resultTopGO$Term <- gsub("biosynthetic process", "biosynthesis", BP_RootTipBoth_resultTopGO$Term)
BP_RootTipBoth_resultTopGO$Name <- paste(BP_RootTipBoth_resultTopGO$GO.ID, "-", BP_RootTipBoth_resultTopGO$Term)
BP_RootTipBoth_resultTopGO$Name <- factor(BP_RootTipBoth_resultTopGO$Name, levels = BP_RootTipBoth_resultTopGO$Name[order(BP_RootTipBoth_resultTopGO$Dataset)])

plot_RootTipBoth_BP_UniqueDEUpGOs <- ggplot(BP_RootTipBoth_resultTopGO, aes(x = Dataset, y = Name, size=Significant, color=-log10(as.numeric(raw.P.value)))) +
  geom_point(alpha=1) +
  scale_colour_viridis(option="D", limits = c(0,30)) +
  scale_size_continuous(limits = c(1, 80), range = c(1,8)) +
  theme_minimal() +
  theme(panel.grid.major = element_line(size = 0.3, linetype = "dotted", colour = "gray", lineend = "butt")) + 
  theme(axis.text.y = element_text(angle =0, hjust = 1, size = 8), axis.text.x = element_text(angle =90, hjust = 1)) +
  labs(size = "Gene Counts", color = "-log10(FDR)") +
  theme(plot.title = element_text(face = "bold")) +
  ylab("Gene Onthology\nBiological Processes") +
  xlab("\nDataset") +
  theme(axis.title.x = element_text(face="bold"), axis.title.y = element_text(face="bold"), )

tiff('UniqueDEUpGOs_RootTipBothDataSet.tiff', units="cm", width=15, height=18, res=600, compression = 'lzw')
plot_RootTipBoth_BP_UniqueDEUpGOs
dev.off()





########Shared Root Tip Genes #########
TheRootTip1_ResultsCriDGE <- TheRootTip1_ResultsCriDGE %>% 
  rownames_to_column(var = "ID")
TheRootTip2_ResultsCriDGE <- TheRootTip2_ResultsCriDGE %>% 
  rownames_to_column(var = "ID")

TheRootTip1_ResultsCriDGE$Presence <- ifelse(TheRootTip1_ResultsCriDGE$ID%in%TheRootTip2_ResultsCriDGE$ID, as.character("Shared"), as.character("Unique"))
TheRootTip2_ResultsCriDGE$Presence <- ifelse(TheRootTip2_ResultsCriDGE$ID%in%TheRootTip1_ResultsCriDGE$ID, as.character("Shared"), as.character("Unique"))

RootTip_SharedExp <- TheRootTip1_ResultsCriDGE %>%
  filter(Presence=="Shared") %>%
  dplyr::rename("gene_id" = "ID") 

RootTip_SharedExp <- RootTip_SharedExp %>%
  left_join(., TrinotateAnnotation, by = "gene_id")

RootTip_SharedExp <- RootTip_SharedExp %>%
  dplyr::select(gene_id, logFC, AveExpr, P.Value, adj.P.Val, Presence, Crichardii_Transcriptome_676_v2_BLASTX, Arabidopsis_Araport11_BLASTP, Azolla_FernBase_BLASTP, Salvinia_FernBase_BLASTP, gene_ontology_BLASTX)

RAM_SharedUpGOs <- RootTip_SharedExp %>%
  filter(logFC>2) %>%
  dplyr::rename("CriID" = "gene_id") %>%
  na.omit(Crichardii_Transcriptome_676_v2_BLASTX)


RAMUpSharedforGOs <- factor(as.integer(TrinotateAnnotation$gene_id %in% RAM_SharedUpGOs$CriID))
names(RAMUpSharedforGOs) <- TrinotateAnnotation$gene_id
head(RAMUpSharedforGOs)


BP_GOdata <- new("topGOdata", ontology = "BP", allGenes = RAMUpSharedforGOs,
                 annot = annFUN.gene2GO, gene2GO = BPGOx_CriAnno_list)

BP_RootTipShared_statTest <- runTest(BP_GOdata, algorithm = "elim", statistic = "fisher" )

BPUp_RootTipShared_resultTopGO <- GenTable(BP_GOdata, raw.P.value = BP_RootTipShared_statTest,
                                           orderBy = "raw.P.value", 
                                           topNodes = 28,
                                           numChar=1000)

BPUp_RootTipShared_resultTopGO <- BPUp_RootTipShared_resultTopGO[-c(8:10,13,17:19,24),]
BPUp_RootTipShared_resultTopGO$raw.P.value <- gsub("< 1e-30", as.numeric("1.0e-30"), BPUp_RootTipShared_resultTopGO$raw.P.value)  
BPUp_RootTipShared_resultTopGO$Term <- gsub(" ending in seed dormancy", "", BPUp_RootTipShared_resultTopGO$Term)
BPUp_RootTipShared_resultTopGO$Term <- gsub("leaf ", "", BPUp_RootTipShared_resultTopGO$Term)
BPUp_RootTipShared_resultTopGO$Term <- gsub("vegetative to reproductive ", "", BPUp_RootTipShared_resultTopGO$Term)
BPUp_RootTipShared_resultTopGO$Term <- gsub(", DNA-templated", "", BPUp_RootTipShared_resultTopGO$Term)
BPUp_RootTipShared_resultTopGO$Term <- gsub(" mediated", "", BPUp_RootTipShared_resultTopGO$Term)

BPUp_RootTipShared_resultTopGO <- cbind(BPUp_RootTipShared_resultTopGO, GeneRatio = BPUp_RootTipShared_resultTopGO$Significant/(as.numeric(nrow(RAM_SharedUpGOs))))
BPUp_RootTipShared_resultTopGO$Name <- paste(BPUp_RootTipShared_resultTopGO$GO.ID, "-", BPUp_RootTipShared_resultTopGO$Term)

plot_RootTip_BP_SharedUpGOs <- ggplot(BPUp_RootTipShared_resultTopGO, aes(y=reorder(Name, GeneRatio), x = GeneRatio, size=Significant, color=-log10(as.numeric(raw.P.value)))) +
  geom_point(alpha=1) +
  scale_colour_viridis(option="D", limits = c(0,30)) +
  scale_size_continuous(limits = c(1, 800), range = c(1,8)) +
  theme_minimal() +
  theme(panel.grid.major = element_line(size = 0.3, linetype = "dotted", colour = "gray", lineend = "butt")) + 
  theme(axis.text.x = element_text(angle = 45, hjust = 1, size = 8), axis.text.y = element_text(vjust = 1, size = 10)) +
  labs(title = "GOs from shared & upregulated transcripts\nbetween the Root Tip dataset", size = "Gene Counts", color = "-log10(FDR)") +
  theme(plot.title = element_text(face = "bold")) +
  ylab("Gene Onthology\nBiological Processes") +
  xlab("\nGene Ratio") +
  xlim(0,0.1) +
  theme(axis.title.x = element_text(face="bold"), axis.title.y = element_text(face="bold"), )

tiff('RootTip_BP_SharedUpGOs.tiff', units="cm", width=20, height=18, res=600, compression = 'lzw')
plot_RootTip_BP_SharedUpGOs
dev.off()








######## Unique Root Tip Genes #########

######## Aragón-Raygoza Dataset ######## 
RootTip_Unique1Exp <- TheRootTip1_ResultsCriDGE %>%
  filter(Presence=="Unique") %>%
  dplyr::rename("gene_id" = "ID") 

RootTip_Unique1Exp <- RootTip_Unique1Exp %>%
  left_join(., TrinotateAnnotation, by = "gene_id")

RootTip_Unique1Exp <- RootTip_Unique1Exp %>%
  dplyr::select(gene_id, logFC, AveExpr, P.Value, adj.P.Val, Presence, Crichardii_Transcriptome_676_v2_BLASTX, Arabidopsis_Araport11_BLASTP, Azolla_FernBase_BLASTP, Salvinia_FernBase_BLASTP, gene_ontology_BLASTX)

RootTip1_UniqueUpGOs <- RootTip_Unique1Exp %>%
  filter(logFC>0) %>%
  dplyr::rename("CriID" = "gene_id") %>%
  na.omit(Crichardii_Transcriptome_676_v2_BLASTX)


RootTip1UniqueforGOs <- factor(as.integer(TrinotateAnnotation$gene_id %in% RootTip1_UniqueUpGOs$CriID))
names(RootTip1UniqueforGOs) <- TrinotateAnnotation$gene_id
head(RootTip1UniqueforGOs)


BP_GOdata <- new("topGOdata", ontology = "BP", allGenes = RootTip1UniqueforGOs,
                 annot = annFUN.gene2GO, gene2GO = BPGOx_CriAnno_list)

BP_RootTip1Unique_statTest <- runTest(BP_GOdata, algorithm = "elim", statistic = "fisher" )

BP_RootTip1Unique_resultTopGO <- GenTable(BP_GOdata, raw.P.value = BP_RootTip1Unique_statTest,
                                           orderBy = "raw.P.value", 
                                           topNodes = 20,
                                           numChar=1000)

BP_RootTip1Unique_resultTopGO <- BP_RootTip1Unique_resultTopGO[-c(6,9:17),]
BP_RootTip1Unique_resultTopGO$raw.P.value <- gsub("< 1e-30", as.numeric("1.0e-30"), BP_RootTip1Unique_resultTopGO$raw.P.value)  
BP_RootTip1Unique_resultTopGO$Term <- gsub(" family amino acid", "", BP_RootTip1Unique_resultTopGO$Term)
BP_RootTip1Unique_resultTopGO$Term <- gsub("secondary ", "", BP_RootTip1Unique_resultTopGO$Term)
BP_RootTip1Unique_resultTopGO$Term <- gsub("-activated", "", BP_RootTip1Unique_resultTopGO$Term)
BP_RootTip1Unique_resultTopGO <- cbind (BP_RootTip1Unique_resultTopGO, Dataset = "Aragón-Raygoza")
BP_RootTip1Unique_resultTopGO <- cbind(BP_RootTip1Unique_resultTopGO, GeneRatio = BP_RootTip1Unique_resultTopGO$Significant/(as.numeric(nrow(RootTip1_UniqueUpGOs))))
BP_RootTip1Unique_resultTopGO$Name <- paste(BP_RootTip1Unique_resultTopGO$GO.ID, "-", BP_RootTip1Unique_resultTopGO$Term)








######## Yu Paper Dataset ######## 
RootTip_Unique2Exp <- TheRootTip2_ResultsCriDGE %>%
  filter(Presence=="Unique") %>%
  dplyr::rename("gene_id" = "ID") 

RootTip_Uniqu2Exp <- RootTip_Unique2Exp %>%
  left_join(., TrinotateAnnotation, by = "gene_id")

RootTip_Unique2Exp <- RootTip_Unique2Exp %>%
  dplyr::select(gene_id, logFC, AveExpr, P.Value, adj.P.Val, Presence, Crichardii_Transcriptome_676_v2_BLASTX, Arabidopsis_Araport11_BLASTP, Azolla_FernBase_BLASTP, Salvinia_FernBase_BLASTP, gene_ontology_BLASTX)

RootTip2_UniqueUpGOs <- RootTip_Unique2Exp %>%
  filter(logFC>0) %>%
  dplyr::rename("CriID" = "gene_id") %>%
  na.omit(Crichardii_Transcriptome_676_v2_BLASTX)


RootTip2UniqueforGOs <- factor(as.integer(TrinotateAnnotation$gene_id %in% RootTip2_UniqueUpGOs$CriID))
names(RootTip2UniqueforGOs) <- TrinotateAnnotation$gene_id
head(RootTip2UniqueforGOs)


BP_GOdata <- new("topGOdata", ontology = "BP", allGenes = RootTip2UniqueforGOs,
                 annot = annFUN.gene2GO, gene2GO = BPGOx_CriAnno_list)

BP_RootTip2Unique_statTest <- runTest(BP_GOdata, algorithm = "elim", statistic = "fisher" )

BP_RootTip2Unique_resultTopGO <- GenTable(BP_GOdata, raw.P.value = BP_RootTip2Unique_statTest,
                                          orderBy = "raw.P.value", 
                                          topNodes = 10,
                                          numChar=1000)

BP_RootTip2Unique_resultTopGO$Term <- gsub(" biosynthetic process.*", " biosynthesis", BP_RootTip2Unique_resultTopGO$Term)
BP_RootTip2Unique_resultTopGO$Term <- gsub(" biosynthetic ", "", BP_RootTip2Unique_resultTopGO$Term)
BP_RootTip2Unique_resultTopGO$Term <- gsub("secondary ", "", BP_RootTip2Unique_resultTopGO$Term)
BP_RootTip2Unique_resultTopGO$Term <- gsub(" involved in", "for", BP_RootTip2Unique_resultTopGO$Term)
BP_RootTip2Unique_resultTopGO$Term <- gsub(", heterochronic", "", BP_RootTip2Unique_resultTopGO$Term)

BP_RootTip2Unique_resultTopGO <- cbind (BP_RootTip2Unique_resultTopGO, Dataset = "Yu et al., 2020")

BP_RootTip2Unique_resultTopGO <- cbind(BP_RootTip2Unique_resultTopGO, GeneRatio = BP_RootTip2Unique_resultTopGO$Significant/(as.numeric(nrow(RootTip2_UniqueUpGOs))))
BP_RootTip2Unique_resultTopGO$Name <- paste(BP_RootTip2Unique_resultTopGO$GO.ID, "-", BP_RootTip2Unique_resultTopGO$Term)



##########Joining both datasets tables ###########
BP_RootTipBoth_resultTopGO <- full_join(BP_RootTip1Unique_resultTopGO,BP_RootTip2Unique_resultTopGO)
BP_RootTipBoth_resultTopGO$Name <- factor(BP_RootTipBoth_resultTopGO$Name, levels = BP_RootTipBoth_resultTopGO$Name[order(BP_RootTipBoth_resultTopGO$Dataset)])


plot_RootTipBoth_BP_UniqueUpGOs <- ggplot(BP_RootTipBoth_resultTopGO, aes(x = Dataset, y = Name, size=Significant, color=-log10(as.numeric(raw.P.value)))) +
  geom_point(alpha=1) +
  scale_colour_viridis(option="D", limits = c(0,30)) +
  scale_size_continuous(limits = c(1, 20), range = c(1,8)) +
  theme_minimal() +
  theme(panel.grid.major = element_line(size = 0.3, linetype = "dotted", colour = "gray", lineend = "butt")) + 
  theme(axis.text.y = element_text(angle =0, hjust = 1, size = 10), axis.text.x = element_text(angle =90, hjust = 1)) +
  labs(size = "Gene Counts", color = "-log10(FDR)") +
  theme(plot.title = element_text(face = "bold")) +
  ylab("Gene Onthology\nBiological Processes") +
  xlab("\nDataset") +
  theme(axis.title.x = element_text(face="bold"), axis.title.y = element_text(face="bold"), )

tiff('UniqueUpGOs_RootTipBothDataSet.tiff', units="cm", width=15, height=18, res=600, compression = 'lzw')
plot_RootTipBoth_BP_UniqueUpGOs
dev.off()









######### GOs for Arabidopsis Annotation #########
GOs_CritoAth <- read.table(file.path("RNASeqAnalysis/6.DGE-Analysis/", "GOs_Cri-to-Ath.txt"), header = FALSE)
AtGOs <- read.table(file.path("RNASeqAnalysis/6.DGE-Analysis/", "AtGOs.txt"), header = FALSE)

BP_AtGOs <- AtGOs %>%
  filter(V3 == "P") %>%
  dplyr::rename(AT_Number = V1)

GOs_CritoAth <- GOs_CritoAth %>%
  dplyr::rename(CriID = V1, AT_Number = V2)

tmp <- left_join(GOs_CritoAth, BP_AtGOs, by = "AT_Number")

tmp <- tmp %>%
  group_by(CriID) %>%
  dplyr::rename("GO:BP" = V2) %>%
  dplyr::select(CriID,"GO:BP")

BPGOx_CriAnno_list <- split(tmp$`GO:BP`, tmp$CriID)


######## Aragón-Raygoza Dataset ######## 
BP_GOdata <- new("topGOdata", ontology = "BP", allGenes = RootTip1UniqueforGOs,
                 annot = annFUN.gene2GO, gene2GO = BPGOx_CriAnno_list)

BP_RootTip1Unique_statTest <- runTest(BP_GOdata, algorithm = "elim", statistic = "fisher" )

BP_RootTip1Unique_resultTopGO <- GenTable(BP_GOdata, raw.P.value = BP_RootTip1Unique_statTest,
                                          orderBy = "raw.P.value", 
                                          topNodes = 20,
                                          numChar=1000)





######## Yu Paper Dataset ######## 
BP_GOdata <- new("topGOdata", ontology = "BP", allGenes = RootTip2UniqueforGOs,
                 annot = annFUN.gene2GO, gene2GO = BPGOx_CriAnno_list)

BP_RootTip2Unique_statTest <- runTest(BP_GOdata, algorithm = "elim", statistic = "fisher" )

BP_RootTip2Unique_resultTopGO <- GenTable(BP_GOdata, raw.P.value = BP_RootTip2Unique_statTest,
                                          orderBy = "raw.P.value", 
                                          topNodes = 20,
                                          numChar=1000)




##########Joining both datasets tables ###########
BP_RootTipBoth_resultTopGO <- full_join(BP_RootTip1Unique_resultTopGO,BP_RootTip2Unique_resultTopGO)
BP_RootTipBoth_resultTopGO$Name <- factor(BP_RootTipBoth_resultTopGO$Name, levels = BP_RootTipBoth_resultTopGO$Name[order(BP_RootTipBoth_resultTopGO$Dataset)])


plot_RootTipBoth_BP_UniqueUpGOs <- ggplot(BP_RootTipBoth_resultTopGO, aes(x = Dataset, y = Name, size=Significant, color=-log10(as.numeric(raw.P.value)))) +
  geom_point(alpha=1) +
  scale_colour_viridis(option="D", limits = c(0,30)) +
  scale_size_continuous(limits = c(1, 20), range = c(1,8)) +
  theme_minimal() +
  theme(panel.grid.major = element_line(size = 0.3, linetype = "dotted", colour = "gray", lineend = "butt")) + 
  theme(axis.text.y = element_text(angle =0, hjust = 1, size = 8), axis.text.x = element_text(angle =90, hjust = 1)) +
  labs(size = "Gene Counts", color = "-log10(FDR)") +
  theme(plot.title = element_text(face = "bold")) +
  ylab("Gene Onthology\nBiological Processes") +
  xlab("\nDataset") +
  theme(axis.title.x = element_text(face="bold"), axis.title.y = element_text(face="bold"), )

tiff('UniqueUpGOs_RootTipBothDataSet.tiff', units="cm", width=15, height=18, res=600, compression = 'lzw')
plot_RootTipBoth_BP_UniqueUpGOs
dev.off()








sum(is.na(TrinotateAnnotation$Crichardii_Transcriptome_676_v2_BLASTX))
sum(is.na(RootTip_SharedExp$Crichardii_Transcriptome_676_v2_BLASTX))
