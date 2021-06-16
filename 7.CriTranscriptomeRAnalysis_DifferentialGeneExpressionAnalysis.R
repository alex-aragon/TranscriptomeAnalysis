setwd('/Users/alejandroaragon/Desktop')

library(edgeR)
library(limma)
library(tximport)
library(readr)
library(tidyverse)
library(viridis)
library(ggplot2)
library(pheatmap)
library(dplyr)
library(trinotateR)
library(topGO)
library(colorspace)


############## DIFFERENTIAL GENE EXPRESSION ANALYSIS ##############


########### PREPARING THE DATASET ###########

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

group <- as.factor(c("Leaves" ,"Leaves" , "Leaves" , "Root", "Root", "Root", "RootTip", "RootTip", "RootTip")) 
CriDGE <- DGEList(CriCounts)
CriDGE$samples$group <- group

normCriDGE <- scaleOffset(CriDGE, CriLength)
normCriDGE$samples$group <- group
norm_expressedCriGenes_list <- filterByExpr(normCriDGE, group=group)
norm_expressedCriGenes <- normCriDGE[norm_expressedCriGenes_list, , keep.lib.sizes=FALSE] 





########### ANALYZING THE DATASET ###########

CriTdesign <- model.matrix(~0+group)
colnames(CriTdesign) <- gsub("group", "", colnames(CriTdesign))
CriTdesign

contr.matrix <- makeContrasts(
  TheRootTip = RootTip,
  TheRoot = Root,
  TheLeaves = Leaves,
  RootTip_vs_Leaves = RootTip - Leaves, 
  RootTip_vs_Root = RootTip - Root, 
  DiffRoot_vs_Leaves = Root - Leaves,
  RAM = "RootTip - ((Root+Leaves)/2)",
  DiffRoot = "Root - ((RootTip+Leaves)/2)",
  OnlyLeaves = "Leaves - ((Root+RootTip)/2)",
  levels = colnames(CriTdesign))
contr.matrix

voom_normCriGenes <- voom(norm_expressedCriGenes, CriTdesign, plot = TRUE)
vExCriTfit <- lmFit(voom_normCriGenes, CriTdesign)
eBvExCriTfit <- eBayes(vExCriTfit)
plotSA(eBvExCriTfit)

resultsCriDGE <- decideTests(eBvExCriTfit, adjust.method = "BH", p.value = 0.05)
write.fit(eBvExCriTfit, resultsCriDGE, file="Ceratopteris-richardii_Transcriptome_ResultsDGE.txt")



versus_vExCriTfit <- contrasts.fit(vExCriTfit, contrasts=contr.matrix) 
versus_eBvExCriTfit <- eBayes(versus_vExCriTfit)


contrast_resultsCriDGE <- topTable(versus_eBvExCriTfit, coef=NULL, n=Inf, p = 0.05, adjust.method = "BH")

TheLeaves_ResultsCriDGE <- topTable(versus_eBvExCriTfit, coef=3, n=Inf, p = 0.05, adjust.method = "BH")
RTvsL_ResultsCriDGE <- topTable(versus_eBvExCriTfit, coef=4, n=Inf, p = 0.05, adjust.method = "BH")
RTvsR_ResultsCriDGE <- topTable(versus_eBvExCriTfit, coef=5, n=Inf, p = 0.05, adjust.method = "BH")
RvsL_ResultsCriDGE <- topTable(versus_eBvExCriTfit, coef=6, n=Inf, p = 0.05, adjust.method = "BH")
RAM_ResultsCriDGE <- topTable(versus_eBvExCriTfit, coef=7, n=Inf, p = 0.05, adjust.method = "BH")
DiffRoot_ResultsCriDGE <- topTable(versus_eBvExCriTfit, coef=8, n=Inf, p = 0.05, adjust.method = "BH")
Leaf_ResultsCriDGE <- topTable(versus_eBvExCriTfit, coef=9, n=Inf, p = 0.05, adjust.method = "BH")



RootTip_ResultsCriDGE <- topTable(eBvExCriTfit, coef=1, n=Inf, p = 0.05, adjust.method = "BH")
Roots_ResultsCriDGE <- topTable(eBvExCriTfit, coef=2, n=Inf, p = 0.05, adjust.method = "BH")
Leaves_ResultsCriDGE <- topTable(eBvExCriTfit, coef=3, n=Inf, p = 0.05, adjust.method = "BH")


RootTip_ResultsCriDGE <- RootTip_ResultsCriDGE %>% 
  rownames_to_column(var = "ID")
Roots_ResultsCriDGE <- Roots_ResultsCriDGE %>% 
  rownames_to_column(var = "ID")
Leaves_ResultsCriDGE <- Leaves_ResultsCriDGE %>% 
  rownames_to_column(var = "ID")

UpList <- list ("Root Tip" = RootTip_ResultsCriDGE$ID, "Differentiated Root" = Roots_ResultsCriDGE$ID, "Leaves and Shoot" = Leaves_ResultsCriDGE$ID)

tiff('CriDGE_UpSetPlot.tiff', units="in", width=7, height=5, res=600, compression = 'lzw')
UpSetR::upset(UpSetR::fromList(UpList),
              mainbar.y.label = "Number of Transcripts\n(Intersection)",
              sets.x.label = "Total of Transcripts",
              point.size = 3,
              queries = list(
                list(
                  query = intersects,
                  params = list("Root Tip","Differentiated Root","Leaves and Shoot"), 
                  color = "#481567ff", 
                  active = T),
                list(
                  query = intersects,
                  params = list("Root Tip","Leaves and Shoot"), 
                  color = "#1f968bff", 
                  active = T),
                list(
                  query = intersects,
                  params = list("Root Tip","Differentiated Root"), 
                  color = "#1f968bff", 
                  active = T),
                list(
                  query = intersects,
                  params = list("Differentiated Root","Leaves and Shoot"), 
                  color = "#1f968bff", 
                  active = T),
                list(
                  query = intersects,
                  params = list("Root Tip"), 
                  color = "#b8de29ff", 
                  active = T),
                list(
                  query = intersects,
                  params = list("Leaves and Shoot"), 
                  color = "#b8de29ff", 
                  active = T),
                list(
                  query = intersects,
                  params = list("Differentiated Root"), 
                  color = "#b8de29ff", 
                  active = T)
              )
)
dev.off()

############## GENE ONTOLOGY ANALYSIS ##############

#Whole transcriptome annotation
TrinotateAnnotation <- read_trinotate("1.CeratopterisTranscriptomeAnnotation-v2.xls")
TrinotateAnnotation <- TrinotateAnnotation %>%
  distinct(gene_id, .keep_all = TRUE)
summary_trinotate(TrinotateAnnotation)


#Obtaining gene ontology from Trinotate file
GOx_CriAnno <- split_GO(TrinotateAnnotation, hit = "gene_ontology_BLASTP")
summary_GO(GOx_CriAnno)

#Obtaining gene ontology for only biological process
BPGOx_CriAnno <- GOx_CriAnno %>%
  filter(ontology=="biological_process") %>%
  dplyr::select(gene,go) %>%
  dplyr::rename("GO:BP"="go") 
BPGOx_CriAnno_list <- split(BPGOx_CriAnno$`GO:BP`, BPGOx_CriAnno$gene)

#Obtaining gene ontology for only molecular function
MFGOx_CriAnno <- GOx_CriAnno %>%
  filter(ontology=="molecular_function") %>%
  dplyr::select(gene,go) %>%
  dplyr::rename("GO:MF"="go") 
MFGOx_CriAnno_list <- split(MFGOx_CriAnno$`GO:MF`, MFGOx_CriAnno$gene)

#Obtaining gene ontology for only cellular component
CCGOx_CriAnno <- GOx_CriAnno %>%
  filter(ontology=="cellular_component") %>%
  dplyr::select(gene,go) %>%
  dplyr::rename("GO:CC"="go") 
CCGOx_CriAnno_list <- split(CCGOx_CriAnno$`GO:CC`, CCGOx_CriAnno$gene)

########### BIOLOGICAL PROCCESSES ###########



########### BP-GOs Upregulated in the RAM ###########
#Running analysis for Biological Processes GOs of the root tip expressed genes with TopGO
RAM_UpGOs <- RAM_ResultsCriDGE %>%
  filter(logFC>=2) %>%
  rownames_to_column(var = "CriID")

RAMforUpGOs <- factor(as.integer(TrinotateAnnotation$gene_id %in% RAM_UpGOs$CriID))
names(RAMforUpGOs) <- TrinotateAnnotation$gene_id
head(RAMforUpGOs)

BP_GOdata <- new("topGOdata", ontology = "BP", allGenes = RAMforUpGOs,
              annot = annFUN.gene2GO, gene2GO = BPGOx_CriAnno_list)

BP_RootTip_statTest <- runTest(BP_GOdata, algorithm = "elim", statistic = "fisher" )

BPUp_RootTip_resultTopGO <- GenTable(BP_GOdata, raw.P.value = BP_RootTip_statTest,
                                orderBy = "raw.P.value", 
                                topNodes = 30,
                                numChar=1000)

BPUp_RootTip_resultTopGO <- cbind(BPUp_RootTip_resultTopGO, GeneRatio = BPUp_RootTip_resultTopGO$Significant/(as.numeric(nrow(RAM_UpGOs))))
BPUp_RootTip_resultTopGO$Name <- paste(BPUp_RootTip_resultTopGO$GO.ID, "-", BPUp_RootTip_resultTopGO$Term)

plot_RootTip1_BP_UpGOs <- ggplot(BPUp_RootTip_resultTopGO, aes(y=reorder(Name, GeneRatio), x = GeneRatio, size=Significant, color=-log10(as.numeric(raw.P.value)))) +
  geom_point(alpha=1) +
  scale_colour_viridis(option="D", limits = c(0,25.5)) +
  scale_size_continuous(limits = c(1, 200), range = c(1,8)) +
  theme_minimal() +
  theme(panel.grid.major = element_line(size = 0.3, linetype = "dotted", colour = "gray", lineend = "butt")) + 
  theme(axis.text.x = element_text(angle = 45, hjust = 1, size = 8), axis.text.y = element_text(vjust = 1, size = 8)) +
  labs(title = "\nGOs Upregulated in the Root Tip", size = "Gene Counts", color = "-log10(FDR)") +
  theme(plot.title = element_text(face = "bold")) +
  ylab("Gene Onthology\nBiological Processes") +
  xlab("\nGene Ratio") +
  xlim(0,0.05) +
  theme(axis.title.x = element_text(face="bold"), axis.title.y = element_text(face="bold"), )




BPUp_RootTip_resultTopGO <- GenTable(BP_GOdata, raw.P.value = BP_RootTip_statTest,
                                     orderBy = "raw.P.value", 
                                     topNodes = 16,
                                     numChar=1000)

BPUp_RootTip_resultTopGO <- BPUp_RootTip_resultTopGO[-c(11),]
BPUp_RootTip_resultTopGO$Term <- gsub("structural ", "", BPUp_RootTip_resultTopGO$Term)
BPUp_RootTip_resultTopGO$Term <- gsub("regulation of cyclin-dependent ", "", BPUp_RootTip_resultTopGO$Term)
BPUp_RootTip_resultTopGO$Term <- gsub("via break-induced replication", "", BPUp_RootTip_resultTopGO$Term)
BPUp_RootTip_resultTopGO$Term <- gsub("involved in nuclear cell cycle DNA replication", "", BPUp_RootTip_resultTopGO$Term)
BPUp_RootTip_resultTopGO <- cbind(BPUp_RootTip_resultTopGO, GeneRatio = BPUp_RootTip_resultTopGO$Significant/(as.numeric(nrow(RAM_UpGOs))))
BPUp_RootTip_resultTopGO$Name <- paste(BPUp_RootTip_resultTopGO$GO.ID, "-", BPUp_RootTip_resultTopGO$Term)

plot_RAM_BP_UpGOs_Main <- ggplot(BPUp_RootTip_resultTopGO, aes(y=reorder(Name, GeneRatio), x = GeneRatio, size=Significant, color=-log10(as.numeric(raw.P.value)))) +
  geom_point(alpha=1) +
  scale_colour_viridis(option="D", limits = c(0,25.5)) +
  scale_size_continuous(limits = c(1, 200), range = c(1,8)) +
  theme_minimal() +
  theme(panel.grid.major = element_line(size = 0.3, linetype = "dotted", colour = "gray", lineend = "butt")) + 
  theme(axis.text.x = element_text(angle = 45, hjust = 1, size = 8), axis.text.y = element_text(vjust = 1, size = 10)) +
  labs(title = "GOs Upregulated in the Root Tip", size = "Gene Counts", color = "-log10(FDR)") +
  theme(plot.title = element_text(face = "bold")) +
  ylab("Gene Onthology\nBiological Processes") +
  xlab("\nGene Ratio") +
  xlim(0,0.05) +
  theme(axis.title.x = element_text(face="bold"), axis.title.y = element_text(face="bold"), )

tiff('CriDGE_BPUpGOs-RootTip.tiff', units="cm", width=22, height=15, res=600, compression = 'lzw')
plot_RAM_BP_UpGOs_Main
dev.off()

tiff('CriDGE_BPUpGOs-RootTip_Complete.tiff', units="cm", width=22, height=18, res=600, compression = 'lzw')
plot_RAM_BP_UpGOs_Complete
dev.off()





########### BP-GOs Downregulated in the RAM ###########
RAM_DownGOs <- RAM_ResultsCriDGE %>%
  filter(logFC<=-2) %>%
  rownames_to_column(var = "CriID") 

RAMforDownGOs <- factor(as.integer(TrinotateAnnotation$gene_id %in% RAM_DownGOs$CriID))
names(RAMforDownGOs) <- TrinotateAnnotation$gene_id
head(RAMforDownGOs)

BP_GOdata <- new("topGOdata", ontology = "BP", allGenes = RAMforDownGOs,
                 annot = annFUN.gene2GO, gene2GO = BPGOx_CriAnno_list)

BP_RootTip_statTest <- runTest(BP_GOdata, algorithm = "elim", statistic = "fisher" )

BPDown_RootTip_resultTopGO <- GenTable(BP_GOdata, raw.P.value = BP_RootTip_statTest,
                                     orderBy = "raw.P.value", 
                                     topNodes = 31,
                                     numChar=1000)

BPDown_RootTip_resultTopGO <- BPDown_RootTip_resultTopGO[-c(9),]
BPDown_RootTip_resultTopGO <- cbind(BPDown_RootTip_resultTopGO, GeneRatio = BPDown_RootTip_resultTopGO$Significant/(as.numeric(nrow(RAM_DownGOs))))
BPDown_RootTip_resultTopGO$Name <- paste(BPDown_RootTip_resultTopGO$GO.ID, "-", BPDown_RootTip_resultTopGO$Term)

plot_RAM_BP_DownGOs <- ggplot(BPDown_RootTip_resultTopGO, aes(y=reorder(Name, GeneRatio), x = GeneRatio, size=Significant, color=-log10(as.numeric(raw.P.value)))) +
  geom_point(alpha=1) +
  scale_colour_viridis(option="D", limits = c(0,25.5)) +
  scale_size_continuous(limits = c(1, 200), range = c(1,8)) +
  theme_minimal() +
  theme(panel.grid.major = element_line(size = 0.3, linetype = "dotted", colour = "gray", lineend = "butt")) + 
  theme(axis.text.x = element_text(angle = 45, hjust = 1, size = 8), axis.text.y = element_text(vjust = 1, size = 8)) +
  labs(title = "GOs Downregulated in the Root Tip", size = "Gene Counts", color = "-log10(FDR)") +
  theme(plot.title = element_text(face = "bold")) +
  ylab("Gene Onthology\nBiological Processes") +
  xlab("\nGene Ratio") +
  xlim(0,0.056) +
  theme(axis.title.x = element_text(face="bold"), axis.title.y = element_text(face="bold"), )

tiff('CriDGE_BPDownGOs-RootTip.tiff', units="cm", width=18, height=21, res=600, compression = 'lzw')
plot_RAM_BP_DownGOs
dev.off()





########### BP-GOs Upregulated in the Differentiated Root ###########
#Running analysis for Biological Processes GOs of the differentiated root expressed genes with TopGO

DiffRoot_UpGOs <- DiffRoot_ResultsCriDGE %>%
  filter(logFC>=2) %>%
  rownames_to_column(var = "CriID")

DiffRootforUpGOs <- factor(as.integer(TrinotateAnnotation$gene_id %in% DiffRoot_UpGOs$CriID))
names(DiffRootforUpGOs) <- TrinotateAnnotation$gene_id
head(DiffRootforUpGOs)

BP_GOdata <- new("topGOdata", ontology = "BP", allGenes = DiffRootforUpGOs,
                 annot = annFUN.gene2GO, gene2GO = BPGOx_CriAnno_list)

BP_DiffRoot_statTest <- runTest(BP_GOdata, algorithm = "elim", statistic = "fisher" )

BPUp_DiffRoot_resultTopGO <- GenTable(BP_GOdata, raw.P.value = BP_DiffRoot_statTest,
                                     orderBy = "raw.P.value", 
                                     topNodes = 30,
                                     numChar=1000)
 
BPUp_DiffRoot_resultTopGO <- cbind(BPUp_DiffRoot_resultTopGO, GeneRatio = BPUp_DiffRoot_resultTopGO$Significant/(as.numeric(nrow(DiffRoot_UpGOs))))
BPUp_DiffRoot_resultTopGO$Name <- paste(BPUp_DiffRoot_resultTopGO$GO.ID, "-", BPUp_DiffRoot_resultTopGO$Term)

plot_DiffRoot_BP_UpGOs_Complete <- ggplot(BPUp_DiffRoot_resultTopGO, aes(y=reorder(Name, GeneRatio), x = GeneRatio, size=Significant, color=-log10(as.numeric(raw.P.value)))) +
  geom_point(alpha=1) +
  scale_colour_viridis(option="D", limits = c(0,25.5)) +
  scale_size_continuous(limits = c(1, 200), range = c(1,8)) +
  theme_minimal() +
  theme(panel.grid.major = element_line(size = 0.3, linetype = "dotted", colour = "gray", lineend = "butt")) + 
  theme(axis.text.x = element_text(angle = 45, hjust = 1, size = 8), axis.text.y = element_text(vjust = 1, size = 8)) +
  labs(title = "GOs Upregulated in the Differentiated Root", size = "Gene Counts", color = "-log10(FDR)") +
  theme(plot.title = element_text(face = "bold")) +
  ylab("Gene Onthology\nBiological Processes") +
  xlab("\nGene Ratio") +
  xlim(0,0.05) +
  theme(axis.title.x = element_text(face="bold"), axis.title.y = element_text(face="bold"), )



BPUp_DiffRoot_resultTopGO <- GenTable(BP_GOdata, raw.P.value = BP_DiffRoot_statTest,
                                      orderBy = "raw.P.value", 
                                      topNodes = 15,
                                      numChar=1000)

BPUp_DiffRoot_resultTopGO$Term <- gsub(", incompatible interaction", "", BPUp_DiffRoot_resultTopGO$Term)
BPUp_DiffRoot_resultTopGO <- cbind(BPUp_DiffRoot_resultTopGO, GeneRatio = BPUp_DiffRoot_resultTopGO$Significant/(as.numeric(nrow(DiffRoot_UpGOs))))
BPUp_DiffRoot_resultTopGO$Name <- paste(BPUp_DiffRoot_resultTopGO$GO.ID, "-", BPUp_DiffRoot_resultTopGO$Term)

plot_DiffRoot_BP_UpGOs_Main <- ggplot(BPUp_DiffRoot_resultTopGO, aes(y=reorder(Name, GeneRatio), x = GeneRatio, size=Significant, color=-log10(as.numeric(raw.P.value)))) +
  geom_point(alpha=1) +
  scale_colour_viridis(option="D", limits = c(0,25.5)) +
  scale_size_continuous(limits = c(1, 200), range = c(1,8)) +
  theme_minimal() +
  theme(panel.grid.major = element_line(size = 0.3, linetype = "dotted", colour = "gray", lineend = "butt")) + 
  theme(axis.text.x = element_text(angle = 45, hjust = 1, size = 8), axis.text.y = element_text(vjust = 1, size = 10)) +
  labs(title = "GOs Upregulated in the Differentiated Root", size = "Gene Counts", color = "-log10(FDR)") +
  theme(plot.title = element_text(face = "bold")) +
  ylab("Gene Onthology\nBiological Processes") +
  xlab("\nGene Ratio") +
  xlim(0,0.05) +
  theme(axis.title.x = element_text(face="bold"), axis.title.y = element_text(face="bold"), )

tiff('CriDGE_BPUpGOs-DiffRoot.tiff', units="cm", width=22, height=15, res=600, compression = 'lzw')
plot_DiffRoot_BP_UpGOs_Main
dev.off()

tiff('CriDGE_BPUpGOs-DiffRoot_Complete.tiff', units="cm", width=22, height=18, res=600, compression = 'lzw')
plot_DiffRoot_BP_UpGOs_Complete
dev.off()





########### BP-GOs Downregulated in the Differentiated Root ###########
DiffRoot_DownGOs <- DiffRoot_ResultsCriDGE %>%
  filter(logFC<=-2) %>%
  rownames_to_column(var = "CriID") 

DiffRootforDownGOs <- factor(as.integer(TrinotateAnnotation$gene_id %in% DiffRoot_DownGOs$CriID))
names(DiffRootforDownGOs) <- TrinotateAnnotation$gene_id
head(DiffRootforDownGOs)

BP_GOdata <- new("topGOdata", ontology = "BP", allGenes = DiffRootforDownGOs,
                 annot = annFUN.gene2GO, gene2GO = BPGOx_CriAnno_list)

BP_DiffRoot_statTest <- runTest(BP_GOdata, algorithm = "elim", statistic = "fisher" )

BPDown_DiffRoot_resultTopGO <- GenTable(BP_GOdata, raw.P.value = BP_DiffRoot_statTest,
                                       orderBy = "raw.P.value", 
                                       topNodes = 30,
                                       numChar=1000)

BPDown_DiffRoot_resultTopGO <- cbind(BPDown_DiffRoot_resultTopGO, GeneRatio = BPDown_DiffRoot_resultTopGO$Significant/(as.numeric(nrow(DiffRoot_DownGOs))))
BPDown_DiffRoot_resultTopGO$Name <- paste(BPDown_DiffRoot_resultTopGO$GO.ID, "-", BPDown_DiffRoot_resultTopGO$Term)

plot_DiffRoot_BP_DownGOs <- ggplot(BPDown_DiffRoot_resultTopGO, aes(y=reorder(Name, GeneRatio), x = GeneRatio, size=Significant, color=-log10(as.numeric(raw.P.value)))) +
  geom_point(alpha=1) +
  scale_colour_viridis(option="D", limits = c(0,25.5)) +
  scale_size_continuous(limits = c(1, 200), range = c(1,8)) +
  theme_minimal() +
  theme(panel.grid.major = element_line(size = 0.3, linetype = "dotted", colour = "gray", lineend = "butt")) + 
  theme(axis.text.x = element_text(angle = 45, hjust = 1, size = 8), axis.text.y = element_text(vjust = 1, size = 8)) +
  labs(title = "GOs Downregulated in\nthe Differentiated Root", size = "Gene Counts", color = "-log10(FDR)") +
  theme(plot.title = element_text(face = "bold")) +
  ylab("Gene Onthology\nBiological Processes") +
  xlab("\nGene Ratio") +
  xlim(0,0.056) +
  theme(axis.title.x = element_text(face="bold"), axis.title.y = element_text(face="bold"), )



tiff('CriDGE_BPDownGOs-DiffRoot.tiff', units="cm", width=21, height=20, res=600, compression = 'lzw')
plot_DiffRoot_BP_DownGOs
dev.off()





########### BP-GOs Upregulated in the Leaves and Shoot ###########
#Running analysis for Biological Processes GOs of the leaves expressed genes with TopGO
Leaves_UpGOs <- Leaf_ResultsCriDGE %>%
  filter(logFC>=2) %>%
  rownames_to_column(var = "CriID")

LeavesforUpGOs <- factor(as.integer(TrinotateAnnotation$gene_id %in% Leaves_UpGOs$CriID))
names(LeavesforUpGOs) <- TrinotateAnnotation$gene_id
head(LeavesforUpGOs)

BP_GOdata <- new("topGOdata", ontology = "BP", allGenes = LeavesforUpGOs,
                 annot = annFUN.gene2GO, gene2GO = BPGOx_CriAnno_list)

BP_Leaves_statTest <- runTest(BP_GOdata, algorithm = "elim", statistic = "fisher" )

BPUp_Leaves_resultTopGO <- GenTable(BP_GOdata, raw.P.value = BP_Leaves_statTest,
                                      orderBy = "raw.P.value", 
                                      topNodes = 30,
                                      numChar=1000)

BPUp_Leaves_resultTopGO <- cbind(BPUp_Leaves_resultTopGO, GeneRatio = BPUp_Leaves_resultTopGO$Significant/(as.numeric(nrow(Leaves_UpGOs))))
BPUp_Leaves_resultTopGO$Name <- paste(BPUp_Leaves_resultTopGO$GO.ID, "-", BPUp_Leaves_resultTopGO$Term)

plot_Leaves_BP_UpGOs_Complete <- ggplot(BPUp_Leaves_resultTopGO, aes(y=reorder(Name, GeneRatio), x = GeneRatio, size=Significant, color=-log10(as.numeric(raw.P.value)))) +
  geom_point(alpha=1) +
  scale_colour_viridis(option="D", limits = c(0,25.5)) +
  scale_size_continuous(limits = c(1, 200), range = c(1,8)) +
  theme_minimal() +
  theme(panel.grid.major = element_line(size = 0.3, linetype = "dotted", colour = "gray", lineend = "butt")) + 
  theme(axis.text.x = element_text(angle = 45, hjust = 1, size = 8), axis.text.y = element_text(vjust = 1, size = 8)) +
  labs(title = "GOs Upregulated in the Leaves & Shoot", size = "Gene Counts", color = "-log10(FDR)") +
  theme(plot.title = element_text(face = "bold")) +
  ylab("Gene Onthology\nBiological Processes") +
  xlab("\nGene Ratio") +
  xlim(0,0.05) +
  theme(axis.title.x = element_text(face="bold"), axis.title.y = element_text(face="bold"), )



BPUp_Leaves_resultTopGO <- GenTable(BP_GOdata, raw.P.value = BP_Leaves_statTest,
                                    orderBy = "raw.P.value", 
                                    topNodes = 16,
                                    numChar=1000)

BPUp_Leaves_resultTopGO <- BPUp_Leaves_resultTopGO[-c(4),]
BPUp_Leaves_resultTopGO$Term <- gsub("coupled ", "", BPUp_Leaves_resultTopGO$Term)
BPUp_Leaves_resultTopGO$Term <- gsub("metabolic ", "", BPUp_Leaves_resultTopGO$Term)
BPUp_Leaves_resultTopGO$Term <- gsub("carbon ", "", BPUp_Leaves_resultTopGO$Term)

BPUp_Leaves_resultTopGO <- cbind(BPUp_Leaves_resultTopGO, GeneRatio = BPUp_Leaves_resultTopGO$Significant/(as.numeric(nrow(Leaves_UpGOs))))
BPUp_Leaves_resultTopGO$Name <- paste(BPUp_Leaves_resultTopGO$GO.ID, "-", BPUp_Leaves_resultTopGO$Term)

plot_Leaves_BP_UpGOs_Main <- ggplot(BPUp_Leaves_resultTopGO, aes(y=reorder(Name, GeneRatio), x = GeneRatio, size=Significant, color=-log10(as.numeric(raw.P.value)))) +
  geom_point(alpha=1) +
  scale_colour_viridis(option="D", limits = c(0,25.5)) +
  scale_size_continuous(limits = c(1, 200), range = c(1,8)) +
  theme_minimal() +
  theme(panel.grid.major = element_line(size = 0.3, linetype = "dotted", colour = "gray", lineend = "butt")) + 
  theme(axis.text.x = element_text(angle = 45, hjust = 1, size = 8), axis.text.y = element_text(vjust = 1, size = 10)) +
  labs(title = "GOs Upregulated in the Leaves & Shoot", size = "Gene Counts", color = "-log10(FDR)") +
  theme(plot.title = element_text(face = "bold")) +
  ylab("Gene Onthology\nBiological Processes") +
  xlab("\nGene Ratio") +
  xlim(0,0.05) +
  theme(axis.title.x = element_text(face="bold"), axis.title.y = element_text(face="bold"), )



tiff('CriDGE_BPUpGOs-Leaves.tiff', units="cm", width=22, height=15, res=600, compression = 'lzw')
plot_Leaves_BP_UpGOs_Main
dev.off()

tiff('CriDGE_BPUpGOs-Leaves_Complete.tiff', units="cm", width=21, height=18, res=600, compression = 'lzw')
plot_Leaves_BP_UpGOs_Complete
dev.off()





########### BP-GOs Downregulated in the Leaves and Shoot ###########
Leaves_DownGOs <- Leaf_ResultsCriDGE %>%
  filter(logFC<=-2) %>%
  rownames_to_column(var = "CriID") 

LeavesforDownGOs <- factor(as.integer(TrinotateAnnotation$gene_id %in% Leaves_DownGOs$CriID))
names(LeavesforDownGOs) <- TrinotateAnnotation$gene_id
head(LeavesforDownGOs)

BP_GOdata <- new("topGOdata", ontology = "BP", allGenes = LeavesforDownGOs,
                 annot = annFUN.gene2GO, gene2GO = BPGOx_CriAnno_list)

BP_Leaves_statTest <- runTest(BP_GOdata, algorithm = "elim", statistic = "fisher" )

BPDown_Leaves_resultTopGO <- GenTable(BP_GOdata, raw.P.value = BP_Leaves_statTest,
                                        orderBy = "raw.P.value", 
                                        topNodes = 30,
                                        numChar=1000)

BPDown_Leaves_resultTopGO <- cbind(BPDown_Leaves_resultTopGO, GeneRatio = BPDown_Leaves_resultTopGO$Significant/(as.numeric(nrow(Leaves_DownGOs))))
BPDown_Leaves_resultTopGO$Name <- paste(BPDown_Leaves_resultTopGO$GO.ID, "-", BPDown_Leaves_resultTopGO$Term)

plot_Leaves_BP_DownGOs <- ggplot(BPDown_Leaves_resultTopGO, aes(y=reorder(Term, GeneRatio), x = GeneRatio, size=Significant, color=-log10(as.numeric(raw.P.value)))) +
  geom_point(alpha=1) +
  scale_colour_viridis(option="D", limits = c(0,25.5)) +
  scale_size_continuous(limits = c(1, 200), range = c(1,8)) +
  theme_minimal() +
  theme(panel.grid.major = element_line(size = 0.3, linetype = "dotted", colour = "gray", lineend = "butt")) + 
  theme(axis.text.x = element_text(angle = 45, hjust = 1, size = 8), axis.text.y = element_text(vjust = 1, size = 8)) +
  labs(title = "GOs Downregulated in the Leaves & Shoot", size = "Gene Counts", color = "-log10(FDR)") +
  theme(plot.title = element_text(face = "bold")) +
  ylab("Gene Onthology\nBiological Processes") +
  xlab("\nGene Ratio") +
  xlim(0,0.056) +
  theme(axis.title.x = element_text(face="bold"), axis.title.y = element_text(face="bold"), )

tiff('CriDGE_BPDownGOs-Leaves.tiff', units="cm", width=18.5, height=21, res=600, compression = 'lzw')
plot_Leaves_BP_DownGOs
dev.off()



########### MOLECULAR FUNCTION ###########


########### MF-GOs Upregulated in the RAM ###########
MF_GOdata <- new("topGOdata", ontology = "MF", allGenes = RAMforUpGOs,
                 annot = annFUN.gene2GO, gene2GO = MFGOx_CriAnno_list)

MF_RootTip_statTest <- runTest(MF_GOdata, algorithm = "elim", statistic = "fisher" )

MFUp_RootTip_resultTopGO <- GenTable(MF_GOdata, raw.P.value = MF_RootTip_statTest,
                                     orderBy = "raw.P.value", 
                                     topNodes = 30,
                                     numChar=1000)

MFUp_RootTip_resultTopGO <- cbind(MFUp_RootTip_resultTopGO, GeneRatio = MFUp_RootTip_resultTopGO$Significant/(as.numeric(nrow(RAM_UpGOs))))
MFUp_RootTip_resultTopGO$Name <- paste(MFUp_RootTip_resultTopGO$GO.ID, "-", MFUp_RootTip_resultTopGO$Term)

plot_RAM_MF_UpGOs <- ggplot(MFUp_RootTip_resultTopGO, aes(y=reorder(Name, GeneRatio), x = GeneRatio, size=Significant, color=-log10(as.numeric(raw.P.value)))) +
  geom_point(alpha=1) +
  scale_colour_viridis(option="D", limits = c(0,30)) +
  scale_size_continuous(limits = c(1, 220), range = c(1,8)) +
  theme_minimal() +
  theme(panel.grid.major = element_line(size = 0.3, linetype = "dotted", colour = "gray", lineend = "butt")) + 
  theme(axis.text.x = element_text(angle = 45, hjust = 1, size = 8), axis.text.y = element_text(vjust = 1, size = 8)) +
  labs(title = "GOs Upregulated in the Root Tip", size = "Gene Counts", color = "-log10(FDR)") +
  theme(plot.title = element_text(face = "bold")) +
  ylab("Gene Onthology\nMolecular Function") +
  xlab("\nGene Ratio") +
  xlim(0,0.1) +
  theme(axis.title.x = element_text(face="bold"), axis.title.y = element_text(face="bold"), )





########### MF-GOs Downregulated in the RAM ###########
MF_GOdata <- new("topGOdata", ontology = "MF", allGenes = RAMforDownGOs,
                 annot = annFUN.gene2GO, gene2GO = MFGOx_CriAnno_list)

MF_RootTip_statTest <- runTest(MF_GOdata, algorithm = "elim", statistic = "fisher" )

MFDown_RootTip_resultTopGO <- GenTable(MF_GOdata, raw.P.value = MF_RootTip_statTest,
                                     orderBy = "raw.P.value", 
                                     topNodes = 30,
                                     numChar=1000)

MFDown_RootTip_resultTopGO <- cbind(MFDown_RootTip_resultTopGO, GeneRatio = MFDown_RootTip_resultTopGO$Significant/(as.numeric(nrow(RAM_DownGOs))))
MFDown_RootTip_resultTopGO$Term <- gsub("oxidoreductase activity.*", "oxidoreductase activity", MFDown_RootTip_resultTopGO$Term)
MFDown_RootTip_resultTopGO$Term <- gsub("protein serine/threonine ", "", MFDown_RootTip_resultTopGO$Term)
MFDown_RootTip_resultTopGO$Term <- gsub("-1,5-bisphosphat", "", MFDown_RootTip_resultTopGO$Term)
MFDown_RootTip_resultTopGO$Name <- paste(MFDown_RootTip_resultTopGO$GO.ID, "-", MFDown_RootTip_resultTopGO$Term)


plot_RAM_MF_DownGOs <- ggplot(MFDown_RootTip_resultTopGO, aes(y=reorder(Name, GeneRatio), x = GeneRatio, size=Significant, color=-log10(as.numeric(raw.P.value)))) +
  geom_point(alpha=1) +
  scale_colour_viridis(option="D", limits = c(0,30)) +
  scale_size_continuous(limits = c(1, 220), range = c(1,8)) +
  theme_minimal() +
  theme(panel.grid.major = element_line(size = 0.3, linetype = "dotted", colour = "gray", lineend = "butt")) + 
  theme(axis.text.x = element_text(angle = 45, hjust = 1, size = 8), axis.text.y = element_text(vjust = 1, size = 8)) +
  labs(title = "GOs Downregulated in the Root Tip", size = "Gene Counts", color = "-log10(FDR)") +
  theme(plot.title = element_text(face = "bold")) +
  ylab("Gene Onthology\nMolecular Function") +
  xlab("\nGene Ratio") +
  xlim(0,0.1) +
  theme(axis.title.x = element_text(face="bold"), axis.title.y = element_text(face="bold"), )





########### MF-GOs Upregulated in the Differentiated Root ###########
MF_GOdata <- new("topGOdata", ontology = "MF", allGenes = DiffRootforUpGOs,
                 annot = annFUN.gene2GO, gene2GO = MFGOx_CriAnno_list)

MF_DiffRoot_statTest <- runTest(MF_GOdata, algorithm = "elim", statistic = "fisher" )

MFUp_DiffRoot_resultTopGO <- GenTable(MF_GOdata, raw.P.value = MF_DiffRoot_statTest,
                                     orderBy = "raw.P.value", 
                                     topNodes = 30,
                                     numChar=1000)

MFUp_DiffRoot_resultTopGO <- cbind(MFUp_DiffRoot_resultTopGO, GeneRatio = MFUp_DiffRoot_resultTopGO$Significant/(as.numeric(nrow(DiffRoot_UpGOs))))
MFUp_DiffRoot_resultTopGO$Term <- gsub("oxidoreductase activity.*", "oxidoreductase activity", MFUp_DiffRoot_resultTopGO$Term)
MFUp_DiffRoot_resultTopGO$Term <- gsub("protein serine/threonine ", "", MFUp_DiffRoot_resultTopGO$Term)
MFUp_DiffRoot_resultTopGO$Name <- paste(MFUp_DiffRoot_resultTopGO$GO.ID, "-", MFUp_DiffRoot_resultTopGO$Term)

plot_DiffRoot_MF_UpGOs <- ggplot(MFUp_DiffRoot_resultTopGO, aes(y=reorder(Name, GeneRatio), x = GeneRatio, size=Significant, color=-log10(as.numeric(raw.P.value)))) +
  geom_point(alpha=1) +
  scale_colour_viridis(option="D", limits = c(0,30)) +
  scale_size_continuous(limits = c(1, 220), range = c(1,8)) +
  theme_minimal() +
  theme(panel.grid.major = element_line(size = 0.3, linetype = "dotted", colour = "gray", lineend = "butt")) + 
  theme(axis.text.x = element_text(angle = 45, hjust = 1, size = 8), axis.text.y = element_text(vjust = 1, size = 8)) +
  labs(title = "GOs Upregulated in\nthe Differentiated Root", size = "Gene Counts", color = "-log10(FDR)") +
  theme(plot.title = element_text(face = "bold")) +
  ylab("Gene Onthology\nMolecular Function") +
  xlab("\nGene Ratio") +
  xlim(0,0.1) +
  theme(axis.title.x = element_text(face="bold"), axis.title.y = element_text(face="bold"), )





########### MF-GOs Downregulated in the Differentiated Root ###########
MF_GOdata <- new("topGOdata", ontology = "MF", allGenes = DiffRootforDownGOs,
                 annot = annFUN.gene2GO, gene2GO = MFGOx_CriAnno_list)

MF_DiffRoot_statTest <- runTest(MF_GOdata, algorithm = "elim", statistic = "fisher" )

MFDown_DiffRoot_resultTopGO <- GenTable(MF_GOdata, raw.P.value = MF_DiffRoot_statTest,
                                       orderBy = "raw.P.value", 
                                       topNodes = 30,
                                       numChar=1000)

MFDown_DiffRoot_resultTopGO <- cbind(MFDown_DiffRoot_resultTopGO, GeneRatio = MFDown_DiffRoot_resultTopGO$Significant/(as.numeric(nrow(DiffRoot_DownGOs))))
MFDown_DiffRoot_resultTopGO$Term <- gsub("oxidoreductase activity.*", "oxidoreductase activity", MFDown_DiffRoot_resultTopGO$Term)
MFDown_DiffRoot_resultTopGO$Term <- gsub(", plus-end-directed", "", MFDown_DiffRoot_resultTopGO$Term)
MFDown_DiffRoot_resultTopGO$Term <- gsub(", thioredoxin disulfide as acceptor", "", MFDown_DiffRoot_resultTopGO$Term)
MFDown_DiffRoot_resultTopGO$Term <- gsub("C16:0-DCA-CoA ", "", MFDown_DiffRoot_resultTopGO$Term)
MFDown_DiffRoot_resultTopGO$Term <- gsub("transferring electrons within the cyclic electron transport pathway of ", "", MFDown_DiffRoot_resultTopGO$Term)
MFDown_DiffRoot_resultTopGO$Term <- gsub("protein serine/threonine", "", MFDown_DiffRoot_resultTopGO$Term)
MFDown_DiffRoot_resultTopGO$Name <- paste(MFDown_DiffRoot_resultTopGO$GO.ID, "-", MFDown_DiffRoot_resultTopGO$Term)

plot_DiffRoot_MF_DownGOs <- ggplot(MFDown_DiffRoot_resultTopGO, aes(y=reorder(Name, GeneRatio), x = GeneRatio, size=Significant, color=-log10(as.numeric(raw.P.value)))) +
  geom_point(alpha=1) +
  scale_colour_viridis(option="D", limits = c(0,30)) +
  scale_size_continuous(limits = c(1, 220), range = c(1,8)) +
  theme_minimal() +
  theme(panel.grid.major = element_line(size = 0.3, linetype = "dotted", colour = "gray", lineend = "butt")) + 
  theme(axis.text.x = element_text(angle = 45, hjust = 1, size = 8), axis.text.y = element_text(vjust = 1, size = 8)) +
  labs(title = "GOs Downregulated in\nthe Differentiated Root", size = "Gene Counts", color = "-log10(FDR)") +
  theme(plot.title = element_text(face = "bold")) +
  ylab("Gene Onthology\nMolecular Function") +
  xlab("\nGene Ratio") +
  xlim(0,0.1) +
  theme(axis.title.x = element_text(face="bold"), axis.title.y = element_text(face="bold"), )





########### MF-GOs Upregulated in the Leaves and Shoot ###########
MF_GOdata <- new("topGOdata", ontology = "MF", allGenes = LeavesforUpGOs,
                 annot = annFUN.gene2GO, gene2GO = MFGOx_CriAnno_list)

MF_Leaves_statTest <- runTest(MF_GOdata, algorithm = "elim", statistic = "fisher" )

MFUp_Leaves_resultTopGO <- GenTable(MF_GOdata, raw.P.value = MF_Leaves_statTest,
                                      orderBy = "raw.P.value", 
                                      topNodes = 30,
                                      numChar=1000)

MFUp_Leaves_resultTopGO <- cbind(MFUp_Leaves_resultTopGO, GeneRatio = MFUp_Leaves_resultTopGO$Significant/(as.numeric(nrow(Leaves_UpGOs))))
MFUp_Leaves_resultTopGO$raw.P.value <- gsub("< 1e-30", as.numeric("1.0e-30"), MFUp_Leaves_resultTopGO$raw.P.value)  
MFUp_Leaves_resultTopGO$Term <- gsub("oxidoreductase activity.*", "oxidoreductase activity", MFUp_Leaves_resultTopGO$Term)
MFUp_Leaves_resultTopGO$Term <- gsub("electron transporter.*", "electron transporter, photosynthesis activity", MFUp_Leaves_resultTopGO$Term)
MFUp_Leaves_resultTopGO$Term <- gsub(", rotational mechanism", "", MFUp_Leaves_resultTopGO$Term)
MFUp_Leaves_resultTopGO$Term <- gsub("-long-chain-", "-", MFUp_Leaves_resultTopGO$Term)
MFUp_Leaves_resultTopGO$Term <- gsub("-1,5-bisphosphate", "", MFUp_Leaves_resultTopGO$Term)
MFUp_Leaves_resultTopGO$Term <- gsub("glyceraldehyde-3-phosphate.*", "NADP+ phosphorylating activity", MFUp_Leaves_resultTopGO$Term)
MFUp_Leaves_resultTopGO$Name <- paste(MFUp_Leaves_resultTopGO$GO.ID, "-", MFUp_Leaves_resultTopGO$Term)

plot_Leaves_MF_UpGOs <- ggplot(MFUp_Leaves_resultTopGO, aes(y=reorder(Name, GeneRatio), x = GeneRatio, size=Significant, color=-log10(as.numeric(raw.P.value)))) +
  geom_point(alpha=1) +
  scale_colour_viridis(option="D", limits = c(0,30)) +
  scale_size_continuous(limits = c(1, 220), range = c(1,8)) +
  theme_minimal() +
  theme(panel.grid.major = element_line(size = 0.3, linetype = "dotted", colour = "gray", lineend = "butt")) + 
  theme(axis.text.x = element_text(angle = 45, hjust = 1, size = 8), axis.text.y = element_text(vjust = 1, size = 8)) +
  labs(title = "GOs Upregulated in\nthe Leaves & Shoot", size = "Gene Counts", color = "-log10(FDR)") +
  theme(plot.title = element_text(face = "bold")) +
  ylab("Gene Onthology\nMolecular Function") +
  xlab("\nGene Ratio") +
  xlim(0,0.1) +
  theme(axis.title.x = element_text(face="bold"), axis.title.y = element_text(face="bold"), )





########### MF-GOs Downregulated in the Leaves and Shoot ###########
MF_GOdata <- new("topGOdata", ontology = "MF", allGenes = LeavesforDownGOs,
                 annot = annFUN.gene2GO, gene2GO = MFGOx_CriAnno_list)

MF_Leaves_statTest <- runTest(MF_GOdata, algorithm = "elim", statistic = "fisher" )

MFDown_Leaves_resultTopGO <- GenTable(MF_GOdata, raw.P.value = MF_Leaves_statTest,
                                        orderBy = "raw.P.value", 
                                        topNodes = 30,
                                        numChar=1000)

MFDown_Leaves_resultTopGO <- cbind(MFDown_Leaves_resultTopGO, GeneRatio = MFDown_Leaves_resultTopGO$Significant/(as.numeric(nrow(Leaves_DownGOs))))
MFDown_Leaves_resultTopGO$Term <- gsub("oxidoreductase activity.*", "oxidoreductase activity", MFDown_Leaves_resultTopGO$Term)
MFDown_Leaves_resultTopGO$Term <- gsub("polygalacturonate 4-alpha-", "", MFDown_Leaves_resultTopGO$Term)
MFDown_Leaves_resultTopGO$Term <- gsub(", RNA polymerase II-specific", "", MFDown_Leaves_resultTopGO$Term)
MFDown_Leaves_resultTopGO$Name <- paste(MFDown_Leaves_resultTopGO$GO.ID, "-", MFDown_Leaves_resultTopGO$Term)

plot_Leaves_MF_DownGOs <- ggplot(MFDown_Leaves_resultTopGO, aes(y=reorder(Name, GeneRatio), x = GeneRatio, size=Significant, color=-log10(as.numeric(raw.P.value)))) +
  geom_point(alpha=1) +
  scale_colour_viridis(option="D", limits = c(0,30)) +
  scale_size_continuous(limits = c(1, 220), range = c(1,8)) +
  theme_minimal() +
  theme(panel.grid.major = element_line(size = 0.3, linetype = "dotted", colour = "gray", lineend = "butt")) + 
  theme(axis.text.x = element_text(angle = 45, hjust = 1, size = 8), axis.text.y = element_text(vjust = 1, size = 8)) +
  labs(title = "GOs Downregulated in\nthe Leaves & Shoot", size = "Gene Counts", color = "-log10(FDR)") +
  theme(plot.title = element_text(face = "bold")) +
  ylab("Gene Onthology\nMolecular Function") +
  xlab("\nGene Ratio") +
  xlim(0,0.1) +
  theme(axis.title.x = element_text(face="bold"), axis.title.y = element_text(face="bold"), )





########### CELLULAR COMPONENT ###########
########### CC-GOs Upregulated in the RAM ###########

CC_GOdata <- new("topGOdata", ontology = "CC", allGenes = RAMforUpGOs,
                 annot = annFUN.gene2GO, gene2GO = CCGOx_CriAnno_list)

CC_RootTip_statTest <- runTest(CC_GOdata, algorithm = "elim", statistic = "fisher" )

CCUp_RootTip_resultTopGO <- GenTable(CC_GOdata, raw.P.value = CC_RootTip_statTest,
                                     orderBy = "raw.P.value", 
                                     topNodes = 30,
                                     numChar=1000)

CCUp_RootTip_resultTopGO <- cbind(CCUp_RootTip_resultTopGO, GeneRatio = CCUp_RootTip_resultTopGO$Significant/(as.numeric(nrow(RAM_UpGOs))))
CCUp_RootTip_resultTopGO$raw.P.value <- gsub("< 1e-30", as.numeric("1.0e-30"), CCUp_RootTip_resultTopGO$raw.P.value)  
CCUp_RootTip_resultTopGO$Name <- paste(CCUp_RootTip_resultTopGO$GO.ID, "-", CCUp_RootTip_resultTopGO$Term)

plot_RAM_CC_UpGOs <- ggplot(CCUp_RootTip_resultTopGO, aes(y=reorder(Name, GeneRatio), x = GeneRatio, size=Significant, color=-log10(as.numeric(raw.P.value)))) +
  geom_point(alpha=1) +
  scale_colour_viridis(option="D", limits = c(0,30)) +
  scale_size_continuous(limits = c(1, 500), range = c(1,8)) +
  theme_minimal() +
  theme(panel.grid.major = element_line(size = 0.3, linetype = "dotted", colour = "gray", lineend = "butt")) + 
  theme(axis.text.x = element_text(angle = 45, hjust = 1, size = 8), axis.text.y = element_text(vjust = 1, size = 8)) +
  labs(title = "GOs Upregulated in the Root Tip", size = "Gene Counts", color = "-log10(FDR)") +
  theme(plot.title = element_text(face = "bold")) +
  ylab("Gene Onthology\nCellular Component") +
  xlab("\nGene Ratio") +
  xlim(0,0.15) +
  theme(axis.title.x = element_text(face="bold"), axis.title.y = element_text(face="bold"), )





########### CC-GOs Downregulated in the RAM ###########
CC_GOdata <- new("topGOdata", ontology = "CC", allGenes = RAMforDownGOs,
                 annot = annFUN.gene2GO, gene2GO = CCGOx_CriAnno_list)

CC_RootTip_statTest <- runTest(CC_GOdata, algorithm = "elim", statistic = "fisher" )

CCDown_RootTip_resultTopGO <- GenTable(CC_GOdata, raw.P.value = CC_RootTip_statTest,
                                       orderBy = "raw.P.value", 
                                       topNodes = 30,
                                       numChar=1000)

CCDown_RootTip_resultTopGO <- cbind(CCDown_RootTip_resultTopGO, GeneRatio = CCDown_RootTip_resultTopGO$Significant/(as.numeric(nrow(RAM_DownGOs))))
CCDown_RootTip_resultTopGO$Name <- paste(CCDown_RootTip_resultTopGO$GO.ID, "-", CCDown_RootTip_resultTopGO$Term)

plot_RAM_CC_DownGOs <- ggplot(CCDown_RootTip_resultTopGO, aes(y=reorder(Name, GeneRatio), x = GeneRatio, size=Significant, color=-log10(as.numeric(raw.P.value)))) +
  geom_point(alpha=1) +
  scale_colour_viridis(option="D", limits = c(0,30)) +
  scale_size_continuous(limits = c(1, 500), range = c(1,8)) +
  theme_minimal() +
  theme(panel.grid.major = element_line(size = 0.3, linetype = "dotted", colour = "gray", lineend = "butt")) + 
  theme(axis.text.x = element_text(angle = 45, hjust = 1, size = 8), axis.text.y = element_text(vjust = 1, size = 8)) +
  labs(title = "GOs Downregulated in the Root Tip", size = "Gene Counts", color = "-log10(FDR)") +
  theme(plot.title = element_text(face = "bold")) +
  ylab("Gene Onthology\nCellular Component") +
  xlab("\nGene Ratio") +
  xlim(0,0.15) +
  theme(axis.title.x = element_text(face="bold"), axis.title.y = element_text(face="bold"), )





########### CC-GOs Upregulated in the Differentiated Root ###########
CC_GOdata <- new("topGOdata", ontology = "CC", allGenes = DiffRootforUpGOs,
                 annot = annFUN.gene2GO, gene2GO = CCGOx_CriAnno_list)

CC_DiffRoot_statTest <- runTest(CC_GOdata, algorithm = "elim", statistic = "fisher" )

CCUp_DiffRoot_resultTopGO <- GenTable(CC_GOdata, raw.P.value = CC_DiffRoot_statTest,
                                      orderBy = "raw.P.value", 
                                      topNodes = 15,
                                      numChar=1000)

CCUp_DiffRoot_resultTopGO <- CCUp_DiffRoot_resultTopGO[-c(11,12,13),]
CCUp_DiffRoot_resultTopGO <- cbind(CCUp_DiffRoot_resultTopGO, GeneRatio = CCUp_DiffRoot_resultTopGO$Significant/(as.numeric(nrow(DiffRoot_UpGOs))))
CCUp_DiffRoot_resultTopGO$Name <- paste(CCUp_DiffRoot_resultTopGO$GO.ID, "-", CCUp_DiffRoot_resultTopGO$Term)

plot_DiffRoot_CC_UpGOs <- ggplot(CCUp_DiffRoot_resultTopGO, aes(y=reorder(Name, GeneRatio), x = GeneRatio, size=Significant, color=-log10(as.numeric(raw.P.value)))) +
  geom_point(alpha=1) +
  scale_colour_viridis(option="D", limits = c(0,30)) +
  scale_size_continuous(limits = c(1, 500), range = c(1,8)) +
  theme_minimal() +
  theme(panel.grid.major = element_line(size = 0.3, linetype = "dotted", colour = "gray", lineend = "butt")) + 
  theme(axis.text.x = element_text(angle = 45, hjust = 1, size = 8), axis.text.y = element_text(vjust = 1, size = 8)) +
  labs(title = "GOs Upregulated in\nthe Differentiated Root", size = "Gene Counts", color = "-log10(FDR)") +
  theme(plot.title = element_text(face = "bold")) +
  ylab("Gene Onthology\nCellular Component") +
  xlab("\nGene Ratio") +
  xlim(0,0.15) +
  theme(axis.title.x = element_text(face="bold"), axis.title.y = element_text(face="bold"), )





########### CC-GOs Downregulated in the Differentiated Root ###########
CC_GOdata <- new("topGOdata", ontology = "CC", allGenes = DiffRootforDownGOs,
                 annot = annFUN.gene2GO, gene2GO = CCGOx_CriAnno_list)

CC_DiffRoot_statTest <- runTest(CC_GOdata, algorithm = "elim", statistic = "fisher" )

CCDown_DiffRoot_resultTopGO <- GenTable(CC_GOdata, raw.P.value = CC_DiffRoot_statTest,
                                        orderBy = "raw.P.value", 
                                        topNodes = 30,
                                        numChar=1000)

CCDown_DiffRoot_resultTopGO <- cbind(CCDown_DiffRoot_resultTopGO, GeneRatio = CCDown_DiffRoot_resultTopGO$Significant/(as.numeric(nrow(DiffRoot_DownGOs))))
CCDown_DiffRoot_resultTopGO$Name <- paste(CCDown_DiffRoot_resultTopGO$GO.ID, "-", CCDown_DiffRoot_resultTopGO$Term)

plot_DiffRoot_CC_DownGOs <- ggplot(CCDown_DiffRoot_resultTopGO, aes(y=reorder(Name, GeneRatio), x = GeneRatio, size=Significant, color=-log10(as.numeric(raw.P.value)))) +
  geom_point(alpha=1) +
  scale_colour_viridis(option="D", limits = c(0,30)) +
  scale_size_continuous(limits = c(1, 220), range = c(1,8)) +
  theme_minimal() +
  theme(panel.grid.major = element_line(size = 0.3, linetype = "dotted", colour = "gray", lineend = "butt")) + 
  theme(axis.text.x = element_text(angle = 45, hjust = 1, size = 8), axis.text.y = element_text(vjust = 1, size = 8)) +
  labs(title = "GOs Downregulated in\nthe Differentiated Root", size = "Gene Counts", color = "-log10(FDR)") +
  theme(plot.title = element_text(face = "bold")) +
  ylab("Gene Onthology\nCellular Component") +
  xlab("\nGene Ratio") +
  xlim(0,0.15) +
  theme(axis.title.x = element_text(face="bold"), axis.title.y = element_text(face="bold"), )





########### CC-GOs Upregulated in the Leaves and Shoot ###########
CC_GOdata <- new("topGOdata", ontology = "CC", allGenes = LeavesforUpGOs,
                 annot = annFUN.gene2GO, gene2GO = CCGOx_CriAnno_list)

CC_Leaves_statTest <- runTest(CC_GOdata, algorithm = "elim", statistic = "fisher" )

CCUp_Leaves_resultTopGO <- GenTable(CC_GOdata, raw.P.value = CC_Leaves_statTest,
                                    orderBy = "raw.P.value", 
                                    topNodes = 30,
                                    numChar=1000)

CCUp_Leaves_resultTopGO <- cbind(CCUp_Leaves_resultTopGO, GeneRatio = CCUp_Leaves_resultTopGO$Significant/(as.numeric(nrow(Leaves_UpGOs))))
CCUp_Leaves_resultTopGO$raw.P.value <- gsub("< 1e-30", as.numeric("1.0e-30"), CCUp_Leaves_resultTopGO$raw.P.value)  
CCUp_Leaves_resultTopGO$Term <- gsub(", catalytic core.*", "", CCUp_Leaves_resultTopGO$Term)
CCUp_Leaves_resultTopGO$Term <- gsub(" dehydrogenase complex.*", " dehydrogenase complex", CCUp_Leaves_resultTopGO$Term)
CCUp_Leaves_resultTopGO$Term <- gsub("protein complex", "", CCUp_Leaves_resultTopGO$Term)
CCUp_Leaves_resultTopGO$Name <- paste(CCUp_Leaves_resultTopGO$GO.ID, "-", CCUp_Leaves_resultTopGO$Term)

plot_Leaves_CC_UpGOs <- ggplot(CCUp_Leaves_resultTopGO, aes(y=reorder(Name, GeneRatio), x = GeneRatio, size=Significant, color=-log10(as.numeric(raw.P.value)))) +
  geom_point(alpha=1) +
  scale_colour_viridis(option="D", limits = c(0,30)) +
  scale_size_continuous(limits = c(1, 500), range = c(1,8)) +
  theme_minimal() +
  theme(panel.grid.major = element_line(size = 0.3, linetype = "dotted", colour = "gray", lineend = "butt")) + 
  theme(axis.text.x = element_text(angle = 45, hjust = 1, size = 8), axis.text.y = element_text(vjust = 1, size = 8)) +
  labs(title = "GOs Upregulated in\nthe Leaves & Shoot", size = "Gene Counts", color = "-log10(FDR)") +
  theme(plot.title = element_text(face = "bold")) +
  ylab("Gene Onthology\nCellular Component") +
  xlab("\nGene Ratio") +
  xlim(0,0.15) +
  theme(axis.title.x = element_text(face="bold"), axis.title.y = element_text(face="bold"), )





########### CC-GOs Downregulated in the Leaves and Shoot ###########
CC_GOdata <- new("topGOdata", ontology = "CC", allGenes = LeavesforDownGOs,
                 annot = annFUN.gene2GO, gene2GO = CCGOx_CriAnno_list)

CC_Leaves_statTest <- runTest(CC_GOdata, algorithm = "elim", statistic = "fisher" )

CCDown_Leaves_resultTopGO <- GenTable(CC_GOdata, raw.P.value = CC_Leaves_statTest,
                                      orderBy = "raw.P.value", 
                                      topNodes = 29,
                                      numChar=1000)

CCDown_Leaves_resultTopGO <- cbind(CCDown_Leaves_resultTopGO, GeneRatio = CCDown_Leaves_resultTopGO$Significant/(as.numeric(nrow(Leaves_DownGOs))))
CCDown_Leaves_resultTopGO$Name <- paste(CCDown_Leaves_resultTopGO$GO.ID, "-", CCDown_Leaves_resultTopGO$Term)

plot_Leaves_CC_DownGOs <- ggplot(CCDown_Leaves_resultTopGO, aes(y=reorder(Name, GeneRatio), x = GeneRatio, size=Significant, color=-log10(as.numeric(raw.P.value)))) +
  geom_point(alpha=1) +
  scale_colour_viridis(option="D", limits = c(0,30)) +
  scale_size_continuous(limits = c(1, 500), range = c(1,8)) +
  theme_minimal() +
  theme(panel.grid.major = element_line(size = 0.3, linetype = "dotted", colour = "gray", lineend = "butt")) + 
  theme(axis.text.x = element_text(angle = 45, hjust = 1, size = 8), axis.text.y = element_text(vjust = 1, size = 8)) +
  labs(title = "GOs Downregulated in\nthe Leaves & Shoot", size = "Gene Counts", color = "-log10(FDR)") +
  theme(plot.title = element_text(face = "bold")) +
  ylab("Gene Onthology\nCellular Component") +
  xlab("\nGene Ratio") +
  xlim(0,0.15) +
  theme(axis.title.x = element_text(face="bold"), axis.title.y = element_text(face="bold"), )







########### ANALYZING CERATOPTERIS TRANSCRIPTION FACTORS ########### 

GOs_CritoAth <- read.table(file.path("RNASeqAnalysis/6.DGE-Analysis/", "GOs_Cri-to-Ath.txt"), header = FALSE)
AthTFs <- read.table(file.path("RNASeqAnalysis/6.DGE-Analysis/", "Ath_TF_list.txt"), header = TRUE)

GOs_CritoAth <- GOs_CritoAth  %>% 
  dplyr::rename(AT_Number = V2,
                CriID = V1)

AthTFs <- AthTFs %>%
  dplyr::rename(AT_Number = Gene_ID)

CriTFs <- left_join(GOs_CritoAth, AthTFs, by = "AT_Number")

CriTFs <- CriTFs %>%
  dplyr::select(CriID,Family) %>%
  dplyr::distinct() %>%
  na.omit() %>%
  arrange(Family)





RAM_ResultsCriDGE <- topTable(versus_eBvExCriTfit, coef=7, n=Inf, p = 0.05, adjust.method = "BH")
DiffRoot_ResultsCriDGE <- topTable(versus_eBvExCriTfit, coef=8, n=Inf, p = 0.05, adjust.method = "BH")
Leaf_ResultsCriDGE <- topTable(versus_eBvExCriTfit, coef=9, n=Inf, p = 0.05, adjust.method = "BH")




############# TF Family Analysis in the RAM ############# 
RAM_ResultsCriDGE <- RAM_ResultsCriDGE %>%
  rownames_to_column(var = "CriID")
  
CriTFs_RAM <- left_join(RAM_ResultsCriDGE, CriTFs, by = "CriID") %>%
  na.omit(Family) %>%
  arrange(Family,logFC)

CriTFs_RAM <- cbind(CriTFs_RAM, Expression = "NA")
CriTFs_RAM$Expression <- fifelse((as.numeric(CriTFs_RAM$logFC))>0, "Upregulated", "Downregulated")

#Checking TFs with FC > or < to 2
tmp0 <- CriTFs_RAM[ which( CriTFs_RAM$logFC >= 2 | CriTFs_RAM$logFC <= -2) , ]

tmp0 <- tmp0 %>%
group_by(Family,Expression) %>% 
  summarise_each(funs(mean)) 

#Checking TFs with FC > or < to 0
#tmp0 <- CriTFs_RAM %>%
  #group_by(Family,Expression) %>% 
  #summarise_each(funs(mean)) 
  

tmp1 <- as.data.frame((table(CriTFs_RAM$Family,CriTFs_RAM$Expression))) %>%
  dplyr::rename(Family = Var1,
                Expression = Var2)

tmp <- left_join(tmp0, tmp1, by = c("Family", "Expression"))

plot_RAM_TFs <- ggplot(tmp, aes(x=Family, y=Expression, size=Freq, color=logFC)) +
  geom_point(alpha=1) +
  scale_colour_continuous_diverging(palette = "Purple-Green", limits = c(-6,6)) +
  theme_minimal() +
  theme(panel.grid.major = element_line(size = 0.3, linetype = "dotted", colour = "gray", lineend = "butt")) + 
  theme(axis.text.x = element_text(angle = 90, hjust = 1)) +
  theme(plot.margin = unit(c(6,1,6,1),"cm")) +
  labs(title = "Transcription Factors Expressed in the Root Tip", size = "TFs Count", color = "Average logFC") +
  theme(plot.title = element_text(face = "bold")) +
  xlab("\nTranscription Factor Family") +
  theme(axis.title.x = element_text(face="bold"), axis.title.y = element_text(face="bold"))
 
tiff('CriDGE_CriTFs-RAM.tiff', units="in", width=10, height=7, res=600, compression = 'lzw')
plot_RAM_TFs
dev.off()




#Checking TFs with FC between 1 to 1.99 (+ and -)
tmp0 <- CriTFs_RAM[ which( CriTFs_RAM$logFC >= 1 & CriTFs_RAM$logFC < 2 | CriTFs_RAM$logFC <= -1 & CriTFs_RAM$logFC > -2) , ]

tmp0 <- tmp0 %>%
  group_by(Family,Expression) %>% 
  summarise_each(funs(mean))  

tmp1 <- as.data.frame((table(CriTFs_RAM$Family,CriTFs_RAM$Expression))) %>%
  dplyr::rename(Family = Var1,
                Expression = Var2)

tmp <- left_join(tmp0, tmp1, by = c("Family", "Expression"))

plot_RAM_TFsFC1 <- ggplot(tmp, aes(x=Family, y=Expression, size=Freq, color=logFC)) +
  geom_point(alpha=1) +
  scale_colour_continuous_diverging(palette = "Purple-Green", limits = c(-6,6)) +
  theme_minimal() +
  theme(panel.grid.major = element_line(size = 0.3, linetype = "dotted", colour = "gray", lineend = "butt")) + 
  theme(axis.text.x = element_text(angle = 90, hjust = 1)) +
  theme(plot.margin = unit(c(6,1,6,1),"cm")) +
  labs(title = "Transcription Factors Expressed in the Root Tip", size = "TFs Count", color = "Average logFC") +
  theme(plot.title = element_text(face = "bold")) +
  xlab("\nTranscription Factor Family") +
  theme(axis.title.x = element_text(face="bold"), axis.title.y = element_text(face="bold"))

tiff('CriDGE_CriTFs_FC1-RAM.tiff', units="in", width=12, height=7, res=600, compression = 'lzw')
plot_RAM_TFsFC1
dev.off()





############# TF Family Analysis in the Differentiated Root ############# 
DiffRoot_ResultsCriDGE <- DiffRoot_ResultsCriDGE %>%
  rownames_to_column(var = "CriID")

CriTFs_DiffRoot <- left_join(DiffRoot_ResultsCriDGE, CriTFs, by = "CriID") %>%
  na.omit(Family) %>%
  arrange(Family,logFC)

CriTFs_DiffRoot <- cbind(CriTFs_DiffRoot, Expression = "NA")
CriTFs_DiffRoot$Expression <- fifelse((as.numeric(CriTFs_DiffRoot$logFC))>0, "Upregulated", "Downregulated")

#Checking TFs with FC > or < to 0
#tmp0 <- CriTFs_DiffRoot %>%
  #group_by(Family,Expression) %>% 
  #summarise_each(funs(mean)) 

#Checking TFs with FC > or < to 2
tmp0 <- CriTFs_DiffRoot[ which( CriTFs_DiffRoot$logFC >= 2 | CriTFs_DiffRoot$logFC <= -2) , ]

tmp0 <- tmp0 %>%
  group_by(Family,Expression) %>% 
  summarise_each(funs(mean)) 


tmp1 <- as.data.frame((table(CriTFs_DiffRoot$Family,CriTFs_DiffRoot$Expression))) %>%
  dplyr::rename(Family = Var1,
                Expression = Var2)

tmp <- left_join(tmp0, tmp1, by = c("Family", "Expression"))

plot_DiffRoot_TFs <- ggplot(tmp, aes(x=Family, y=Expression, size=Freq, color=logFC)) +
  geom_point(alpha=1) +
  scale_colour_continuous_diverging(palette = "Purple-Green", limits = c(-6,6)) +
  theme_minimal() +
  theme(panel.grid.major = element_line(size = 0.3, linetype = "dotted", colour = "gray", lineend = "butt")) + 
  theme(axis.text.x = element_text(angle = 90, hjust = 1)) +
  theme(plot.margin = unit(c(6,1,6,1),"cm")) +
  labs(title = "Transcription Factors Expressed in the Differentiated Root", size = "TFs Count", color = "Average logFC") +
  theme(plot.title = element_text(face = "bold")) +
  xlab("\nTranscription Factor Family") +
  theme(axis.title.x = element_text(face="bold"), axis.title.y = element_text(face="bold"))

tiff('CriDGE_CriTFs-DiffRoot.tiff', units="in", width=10, height=7, res=600, compression = 'lzw')
plot_DiffRoot_TFs
dev.off()



#Checking TFs with FC between 1 to 1.99 (+ and -)
tmp0 <- CriTFs_DiffRoot[ which( CriTFs_DiffRoot$logFC >= 1 & CriTFs_DiffRoot$logFC < 2 | CriTFs_DiffRoot$logFC <= -1 & CriTFs_DiffRoot$logFC > -2) , ]

tmp0 <- tmp0 %>%
  group_by(Family,Expression) %>% 
  summarise_each(funs(mean)) 


tmp1 <- as.data.frame((table(CriTFs_DiffRoot$Family,CriTFs_DiffRoot$Expression))) %>%
  dplyr::rename(Family = Var1,
                Expression = Var2)

tmp <- left_join(tmp0, tmp1, by = c("Family", "Expression"))

plot_DiffRoot_TFsFC1 <- ggplot(tmp, aes(x=Family, y=Expression, size=Freq, color=logFC)) +
  geom_point(alpha=1) +
  scale_colour_continuous_diverging(palette = "Purple-Green", limits = c(-6,6)) +
  theme_minimal() +
  theme(panel.grid.major = element_line(size = 0.3, linetype = "dotted", colour = "gray", lineend = "butt")) + 
  theme(axis.text.x = element_text(angle = 90, hjust = 1)) +
  theme(plot.margin = unit(c(6,1,6,1),"cm")) +
  labs(title = "Transcription Factors Expressed in the Differentiated Root", size = "TFs Count", color = "Average logFC") +
  theme(plot.title = element_text(face = "bold")) +
  xlab("\nTranscription Factor Family") +
  theme(axis.title.x = element_text(face="bold"), axis.title.y = element_text(face="bold"))

tiff('CriDGE_CriTFs_FC1-DiffRoot.tiff', units="in", width=12, height=7, res=600, compression = 'lzw')
plot_DiffRoot_TFsFC1
dev.off()





############# TF Family Analysis in the Leaves and Shoot ############# 
Leaf_ResultsCriDGE <- Leaf_ResultsCriDGE %>%
  rownames_to_column(var = "CriID")

CriTFs_Leaf <- left_join(Leaf_ResultsCriDGE, CriTFs, by = "CriID") %>%
  na.omit(Family) %>%
  arrange(Family,logFC)

CriTFs_Leaf <- cbind(CriTFs_Leaf, Expression = "NA")
CriTFs_Leaf$Expression <- fifelse((as.numeric(CriTFs_Leaf$logFC))>0, "Upregulated", "Downregulated")

#Checking TFs with FC > or < to 0
#tmp0 <- CriTFs_Leaf %>%
  #group_by(Family,Expression) %>% 
  #summarise_each(funs(mean)) 

#Checking TFs with FC > or < to 2
tmp0 <- CriTFs_Leaf[ which( CriTFs_Leaf$logFC >= 2 | CriTFs_Leaf$logFC <= -2) , ]

tmp0 <- tmp0 %>%
  group_by(Family,Expression) %>% 
  summarise_each(funs(mean)) 


tmp1 <- as.data.frame((table(CriTFs_Leaf$Family,CriTFs_Leaf$Expression))) %>%
  dplyr::rename(Family = Var1,
                Expression = Var2)

tmp <- left_join(tmp0, tmp1, by = c("Family", "Expression"))


plot_Leaf_TFs <- ggplot(tmp, aes(x=Family, y=Expression, size=Freq, color=logFC)) +
  geom_point(alpha=1) +
  scale_colour_continuous_diverging(palette = "Purple-Green", limits = c(-6,6)) +
  theme_minimal() +
  theme(panel.grid.major = element_line(size = 0.3, linetype = "dotted", colour = "gray", lineend = "butt")) + 
  theme(axis.text.x = element_text(angle = 90, hjust = 1)) +
  theme(plot.margin = unit(c(6,1,6,1),"cm")) +
  labs(title = "Transcription Factors Expressed in the Leaves and Shoot", size = "TFs Count", color = "Average logFC") +
  theme(plot.title = element_text(face = "bold")) +
  xlab("\nTranscription Factor Family") +
  theme(axis.title.x = element_text(face="bold"), axis.title.y = element_text(face="bold"))

tiff('CriDGE_CriTFs-Leaves.tiff', units="in", width=10, height=7, res=600, compression = 'lzw')
plot_Leaf_TFs
dev.off()



#Checking TFs with FC between 1 to 1.99 (+ and -)
tmp0 <- CriTFs_Leaf[ which( CriTFs_Leaf$logFC >= 1 & CriTFs_Leaf$logFC < 2 | CriTFs_Leaf$logFC <= -1 & CriTFs_Leaf$logFC > 2) , ]

tmp0 <- tmp0 %>%
  group_by(Family,Expression) %>% 
  summarise_each(funs(mean)) 


tmp1 <- as.data.frame((table(CriTFs_Leaf$Family,CriTFs_Leaf$Expression))) %>%
  dplyr::rename(Family = Var1,
                Expression = Var2)

tmp <- left_join(tmp0, tmp1, by = c("Family", "Expression"))


plot_Leaf_TFsFC1 <- ggplot(tmp, aes(x=Family, y=Expression, size=Freq, color=logFC)) +
  geom_point(alpha=1) +
  scale_colour_continuous_diverging(palette = "Purple-Green", limits = c(-6,6)) +
  theme_minimal() +
  theme(panel.grid.major = element_line(size = 0.3, linetype = "dotted", colour = "gray", lineend = "butt")) + 
  theme(axis.text.x = element_text(angle = 90, hjust = 1)) +
  theme(plot.margin = unit(c(6,1,6,1),"cm")) +
  labs(title = "Transcription Factors Expressed in the Leaves and Shoot", size = "TFs Count", color = "Average logFC") +
  theme(plot.title = element_text(face = "bold")) +
  xlab("\nTranscription Factor Family") +
  theme(axis.title.x = element_text(face="bold"), axis.title.y = element_text(face="bold"))

tiff('CriDGE_CriTFs_FC1-Leaves.tiff', units="in", width=12, height=7, res=600, compression = 'lzw')
plot_Leaf_TFsFC1
dev.off()






#################  ANALYZING UNIQUE GENES IN EACH TISSUE #################
TrinotateAnnotation <- read_trinotate("1.CeratopterisTranscriptomeAnnotation-v2.xls")
TrinotateAnnotation <- TrinotateAnnotation %>%
  distinct(gene_id, .keep_all = TRUE)
summary_trinotate(TrinotateAnnotation)

TrinotateAnnotation$Crichardii_Transcriptome_676_v2_BLASTX <- gsub("[[:punct:]]Ceric.*$", "", TrinotateAnnotation$Crichardii_Transcriptome_676_v2_BLASTX)

CriGenome_GOs <- read.table(file.path("RNASeqAnalysis/6.DGE-Analysis/", "CeratopterisGenome_GOs.txt"), header = FALSE, sep = "\t")

CriGenome_GOs <- CriGenome_GOs %>%
  dplyr::rename(Crichardii_Transcriptome_676_v2_BLASTX = V1,
                gene_ontology = V2)

CriGenomeAnnotation <- left_join(TrinotateAnnotation, CriGenome_GOs, by = "Crichardii_Transcriptome_676_v2_BLASTX") %>%
  dplyr::select(gene_id,transcript_id,prot_id,Crichardii_Transcriptome_676_v2_BLASTX,gene_ontology)
CriGenomeAnnotation$gene_ontology <- gsub(" ", "`", CriGenomeAnnotation$gene_ontology)
CriGenomeAnnotation <- dplyr::rename(CriGenomeAnnotation, gene_ontology_BLASTX = gene_ontology)

GOx_CriGenoAnno <- split_GO(CriGenomeAnnotation, hit = "gene_ontology_BLASTX")
summary_GO(GOx_CriGenoAnno)

CriGenoAnno_list <- split(GOx_CriGenoAnno$go, BPGOx_CriAnno$gene)




######### Unique Genes in the Root Tip ##########


RootTip_ResultsCriDGE$Presence <- ifelse(RootTip_ResultsCriDGE$ID%in%Roots_ResultsCriDGE$ID, as.character("Shared"),
                                   ifelse(RootTip_ResultsCriDGE$ID%in%Leaves_ResultsCriDGE$ID, as.character("Shared"), as.character("Unique")))


RootTip_UniqueExp <- RootTip_ResultsCriDGE %>%
  filter(Presence=="Unique") %>%
  dplyr::rename("gene_id" = "ID") 

RootTip_UniqueExp <- RootTip_UniqueExp %>%
  left_join(., CriGenomeAnnotation, by = "gene_id")

RootTip_UniqueExp <- RootTip_UniqueExp %>%
  dplyr::select(gene_id, logFC, AveExpr, P.Value, adj.P.Val, Crichardii_Transcriptome_676_v2_BLASTX, gene_ontology_BLASTX)



GOs_CritoAth <- read.table(file.path("RNASeqAnalysis/6.DGE-Analysis/", "GOs_Cri-to-Ath.txt"), header = FALSE)
AtGOs <- read.table(file.path("RNASeqAnalysis/6.DGE-Analysis/", "AtGOs.txt"), header = FALSE)

BP_AtGOs <- AtGOs %>%
  filter(V3 == "P") %>%
  dplyr::rename(AT_Number = V1)

MF_AtGOs <- AtGOs %>%
  filter(V3 == "F") %>%
  dplyr::rename(AT_Number = V1)

GOs_CritoAth <- GOs_CritoAth %>%
  dplyr::rename(CriID = V1, AT_Number = V2)

tmp <- left_join(GOs_CritoAth, BP_AtGOs, by = "AT_Number")

tmp <- tmp %>%
  group_by(CriID) %>%
  dplyr::rename("GO:BP" = V2) %>%
  dplyr::select(CriID,"GO:BP")

BPGOx_CriAnno_list <- split(tmp$`GO:BP`, tmp$CriID)

RAM_UniqueUpGOs <- RootTip_UniqueExp %>%
  filter(logFC>0) %>%
  dplyr::rename("CriID" = "gene_id") 

RAMforUniqueUpGOs <- factor(as.integer(TrinotateAnnotation$gene_id %in% RAM_UniqueUpGOs$CriID))
names(RAMforUniqueUpGOs) <- TrinotateAnnotation$gene_id
head(RAMforUniqueUpGOs)

BP_GOdata <- new("topGOdata", ontology = "BP", allGenes = RAMforUniqueUpGOs,
                 annot = annFUN.gene2GO, gene2GO = BPGOx_CriAnno_list)

BP_RootTipUnique_statTest <- runTest(BP_GOdata, algorithm = "elim", statistic = "fisher" )

BPUp_RootTipUnique_resultTopGO <- GenTable(BP_GOdata, raw.P.value = BP_RootTipUnique_statTest,
                                     orderBy = "raw.P.value", 
                                     topNodes = 15,
                                     numChar=1000)

BPUp_RootTipUnique_resultTopGO$Term <- gsub("positive regulation of ", "", BPUp_RootTipUnique_resultTopGO$Term)
BPUp_RootTipUnique_resultTopGO <- cbind(BPUp_RootTipUnique_resultTopGO, GeneRatio = BPUp_RootTipUnique_resultTopGO$Significant/(as.numeric(nrow(RAM_UpGOs))))
BPUp_RootTipUnique_resultTopGO$Name <- paste(BPUp_RootTipUnique_resultTopGO$GO.ID, "-", BPUp_RootTipUnique_resultTopGO$Term)

plot_RAM_BP_UniqueUpGOs <- ggplot(BPUp_RootTipUnique_resultTopGO, aes(y=reorder(Name, GeneRatio), x = GeneRatio, size=Significant, color=-log10(as.numeric(raw.P.value)))) +
  geom_point(alpha=1) +
  scale_colour_viridis(option="D", limits = c(0,25.5)) +
  scale_size_continuous(limits = c(1, 200), range = c(1,8)) +
  theme_minimal() +
  theme(panel.grid.major = element_line(size = 0.3, linetype = "dotted", colour = "gray", lineend = "butt")) + 
  theme(axis.text.x = element_text(angle = 45, hjust = 1, size = 8), axis.text.y = element_text(vjust = 1, size = 8)) +
  labs(title = "GOs from Unique & Upregulated\nTranscripts in the Root Tip", size = "Gene Counts", color = "-log10(FDR)") +
  theme(plot.title = element_text(face = "bold")) +
  ylab("Gene Onthology\nBiological Processes") +
  xlab("\nGene Ratio") +
  xlim(0,0.06) +
  theme(axis.title.x = element_text(face="bold"), axis.title.y = element_text(face="bold"), )

tiff('CriDGE_UniqueBPUpGOs-RootTip.tiff', units="cm", width=18, height=18, res=600, compression = 'lzw')
plot_RAM_BP_UniqueUpGOs
dev.off()





######### Unique Genes in the Differentiated Root ##########
Roots_ResultsCriDGE$Presence <- ifelse(Roots_ResultsCriDGE$ID%in%RootTip_ResultsCriDGE$ID, as.character("Shared"),
                                         ifelse(Roots_ResultsCriDGE$ID%in%Leaves_ResultsCriDGE$ID, as.character("Shared"), as.character("Unique")))

DiffRoot_UniqueExp <- Roots_ResultsCriDGE %>%
  filter(Presence=="Unique") %>%
  dplyr::rename("gene_id" = "ID") 

DiffRoot_UniqueExp <- DiffRoot_UniqueExp %>%
  left_join(., CriGenomeAnnotation, by = "gene_id")

DiffRoot_UniqueExp <- DiffRoot_UniqueExp %>%
  dplyr::select(gene_id, logFC, AveExpr, P.Value, adj.P.Val, Crichardii_Transcriptome_676_v2_BLASTX, gene_ontology_BLASTX)

DiffRoot_UniqueUpGOs <- DiffRoot_UniqueExp %>%
  filter(logFC>0) %>%
  dplyr::rename("CriID" = "gene_id") 

DiffRootforUniqueUpGOs <- factor(as.integer(TrinotateAnnotation$gene_id %in% DiffRoot_UniqueUpGOs$CriID))
names(DiffRootforUniqueUpGOs) <- TrinotateAnnotation$gene_id
head(DiffRootforUniqueUpGOs)

BP_GOdata <- new("topGOdata", ontology = "BP", allGenes = DiffRootforUniqueUpGOs,
                 annot = annFUN.gene2GO, gene2GO = BPGOx_CriAnno_list)

BP_DiffRootUnique_statTest <- runTest(BP_GOdata, algorithm = "elim", statistic = "fisher" )

BPUp_DiffRootUnique_resultTopGO <- GenTable(BP_GOdata, raw.P.value = BP_DiffRootUnique_statTest,
                                     orderBy = "raw.P.value", 
                                     topNodes = 15,
                                     numChar=1000)

BPUp_DiffRootUnique_resultTopGO$Term <- gsub("jasmonic acid and", "JA &", BPUp_DiffRootUnique_resultTopGO$Term)
BPUp_DiffRootUnique_resultTopGO <- cbind(BPUp_DiffRootUnique_resultTopGO, GeneRatio = BPUp_DiffRootUnique_resultTopGO$Significant/(as.numeric(nrow(DiffRoot_UniqueUpGOs))))
BPUp_DiffRootUnique_resultTopGO$Name <- paste(BPUp_DiffRootUnique_resultTopGO$GO.ID, "-", BPUp_DiffRootUnique_resultTopGO$Term)

plot_DiffRoot_BP_UniqueUpGOs <- ggplot(BPUp_DiffRootUnique_resultTopGO, aes(y=reorder(Name, GeneRatio), x = GeneRatio, size=Significant, color=-log10(as.numeric(raw.P.value)))) +
  geom_point(alpha=1) +
  scale_colour_viridis(option="D", limits = c(0,25.5)) +
  scale_size_continuous(limits = c(1, 200), range = c(1,8)) +
  theme_minimal() +
  theme(panel.grid.major = element_line(size = 0.3, linetype = "dotted", colour = "gray", lineend = "butt")) + 
  theme(axis.text.x = element_text(angle = 45, hjust = 1, size = 8), axis.text.y = element_text(vjust = 1, size = 8)) +
  labs(title = "GOs from Unique & Upregulated\nTranscripts in the Differentiated Root", size = "Gene Counts", color = "-log10(FDR)") +
  theme(plot.title = element_text(face = "bold")) +
  ylab("Gene Onthology\nBiological Processes") +
  xlab("\nGene Ratio") +
  xlim(0,0.06) +
  theme(axis.title.x = element_text(face="bold"), axis.title.y = element_text(face="bold"), )

tiff('CriDGE_UniqueBPUpGOs-DiffRoot.tiff', units="cm", width=18, height=18, res=600, compression = 'lzw')
plot_DiffRoot_BP_UniqueUpGOs
dev.off()





######### Unique Genes in the Leaves and Shoot ##########
Leaves_ResultsCriDGE$Presence <- ifelse(Leaves_ResultsCriDGE$ID%in%RootTip_ResultsCriDGE$ID, as.character("Shared"),
                                       ifelse(Leaves_ResultsCriDGE$ID%in%Roots_ResultsCriDGE$ID, as.character("Shared"), as.character("Unique")))

Leaves_UniqueExp <- Leaves_ResultsCriDGE %>%
  filter(Presence=="Unique") %>%
  dplyr::rename("gene_id" = "ID") 

Leaves_UniqueExp <- Leaves_UniqueExp %>%
  left_join(., CriGenomeAnnotation, by = "gene_id")

Leaves_UniqueExp <- Leaves_UniqueExp %>%
  dplyr::select(gene_id, logFC, AveExpr, P.Value, adj.P.Val, Crichardii_Transcriptome_676_v2_BLASTX, gene_ontology_BLASTX)

Leaves_UniqueUpGOs <- Leaves_UniqueExp %>%
  filter(logFC>0) %>%
  dplyr::rename("CriID" = "gene_id") 

LeavesforUniqueUpGOs <- factor(as.integer(TrinotateAnnotation$gene_id %in% Leaves_UniqueUpGOs$CriID))
names(LeavesforUniqueUpGOs) <- TrinotateAnnotation$gene_id
head(LeavesforUniqueUpGOs)

BP_GOdata <- new("topGOdata", ontology = "BP", allGenes = LeavesforUniqueUpGOs,
                 annot = annFUN.gene2GO, gene2GO = BPGOx_CriAnno_list)

BP_LeavesUnique_statTest <- runTest(BP_GOdata, algorithm = "elim", statistic = "fisher" )

BPUp_LeavesUnique_resultTopGO <- GenTable(BP_GOdata, raw.P.value = BP_LeavesUnique_statTest,
                                      orderBy = "raw.P.value", 
                                      topNodes = 15,
                                      numChar=1000)

BPUp_LeavesUnique_resultTopGO$Term <- gsub(" in multidimensional cell growth", "", BPUp_LeavesUnique_resultTopGO$Term)
BPUp_LeavesUnique_resultTopGO <- cbind(BPUp_LeavesUnique_resultTopGO, GeneRatio = BPUp_LeavesUnique_resultTopGO$Significant/(as.numeric(nrow(Leaves_UniqueUpGOs))))
BPUp_LeavesUnique_resultTopGO$Name <- paste(BPUp_LeavesUnique_resultTopGO$GO.ID, "-", BPUp_LeavesUnique_resultTopGO$Term)

plot_Leaves_BP_UniqueUpGOs <- ggplot(BPUp_LeavesUnique_resultTopGO, aes(y=reorder(Name, GeneRatio), x = GeneRatio, size=Significant, color=-log10(as.numeric(raw.P.value)))) +
  geom_point(alpha=1) +
  scale_colour_viridis(option="D", limits = c(0,25.5)) +
  scale_size_continuous(limits = c(1, 200), range = c(1,8)) +
  theme_minimal() +
  theme(panel.grid.major = element_line(size = 0.3, linetype = "dotted", colour = "gray", lineend = "butt")) + 
  theme(axis.text.x = element_text(angle = 45, hjust = 1, size = 8), axis.text.y = element_text(vjust = 1, size = 8)) +
  labs(title = "GOs from Unique & Upregulated\nTranscripts in the Leaves & Shoot", size = "Gene Counts", color = "-log10(FDR)") +
  theme(plot.title = element_text(face = "bold")) +
  ylab("Gene Onthology\nBiological Processes") +
  xlab("\nGene Ratio") +
  xlim(0,0.06) +
  theme(axis.title.x = element_text(face="bold"), axis.title.y = element_text(face="bold"), )

tiff('CriDGE_UniqueBPUpGOs-Leaves.tiff', units="cm", width=18, height=18, res=600, compression = 'lzw')
plot_Leaves_BP_UniqueUpGOs
dev.off()




