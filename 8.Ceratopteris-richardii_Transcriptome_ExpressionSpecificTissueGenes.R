setwd('/Users/alejandroaragon/Desktop')

library(edgeR)
library(limma)
library(tximport)
library(readr)
library(tidyverse)
library(viridis)
library(ggplot2)
library(dplyr)
library(colorspace)


########### DIFFERENTIAL GENE EXPRESSION ANALYSIS ########### 



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

versus_vExCriTfit <- contrasts.fit(vExCriTfit, contrasts=contr.matrix) 
versus_eBvExCriTfit <- eBayes(versus_vExCriTfit)

RAM_ResultsCriDGE <- topTable(versus_eBvExCriTfit, coef=7, n=Inf, p = 0.05, adjust.method = "BH")
DiffRoot_ResultsCriDGE <- topTable(versus_eBvExCriTfit, coef=8, n=Inf, p = 0.05, adjust.method = "BH")
Leaf_ResultsCriDGE <- topTable(versus_eBvExCriTfit, coef=9, n=Inf, p = 0.05, adjust.method = "BH")





###################################
#### CHECKING KNOWN ROOT GENES ####
###################################

######### THE ENDODERMIS ##########

Endodermis <- c("Crichardii-Hnn_ARJA-v1.0_transcript23688",
                "Crichardii-Hnn_ARJA-v1.0_transcript116363",
                "Crichardii-Hnn_ARJA-v1.0_transcript30597",
                "Crichardii-Hnn_ARJA-v1.0_transcript83572",
                "Crichardii-Hnn_ARJA-v1.0_transcript130227",
                "Crichardii-Hnn_ARJA-v1.0_transcript44715",
                "Crichardii-Hnn_ARJA-v1.0_transcript42446",
                "Crichardii-Hnn_ARJA-v1.0_transcript77150",
                "Crichardii-Hnn_ARJA-v1.0_transcript76580",
                "Crichardii-Hnn_ARJA-v1.0_transcript68371",
                "Crichardii-Hnn_ARJA-v1.0_transcript68372",
                "Crichardii-Hnn_ARJA-v1.0_transcript48787",
                "Crichardii-Hnn_ARJA-v1.0_transcript34969",
                "Crichardii-Hnn_ARJA-v1.0_transcript85335",
                "Crichardii-Hnn_ARJA-v1.0_transcript45117")

RAM_ResultsCriDGE <- RAM_ResultsCriDGE %>%
  rownames_to_column(var = "CriID")
CriEndodermisRAM <- subset(RAM_ResultsCriDGE, CriID %in% Endodermis)
CriEndodermisRAM <- cbind(CriEndodermisRAM, Tissue = "Root Tip")

DiffRoot_ResultsCriDGE <- DiffRoot_ResultsCriDGE %>%
  rownames_to_column(var = "CriID")
CriEndodermisDiffRoot <- subset(DiffRoot_ResultsCriDGE, CriID %in% Endodermis)
CriEndodermisDiffRoot <- cbind(CriEndodermisDiffRoot, Tissue = "Differentiated Root")

Leaf_ResultsCriDGE <- Leaf_ResultsCriDGE %>%
  rownames_to_column(var = "CriID")
CriEndodermisLeaves <- subset(Leaf_ResultsCriDGE, CriID %in% Endodermis)
CriEndodermisLeaves <- cbind(CriEndodermisLeaves, Tissue = "Leaves and Shoot")

CriEndodermis <- full_join(CriEndodermisRAM,CriEndodermisDiffRoot) %>%
  full_join(. ,CriEndodermisLeaves)

SCR <- c("Crichardii-Hnn_ARJA-v1.0_transcript23688")
SHR <- c("Crichardii-Hnn_ARJA-v1.0_transcript116363")
SCL28 <- c("Crichardii-Hnn_ARJA-v1.0_transcript45117")
RBR <- c("Crichardii-Hnn_ARJA-v1.0_transcript30597",
         "Crichardii-Hnn_ARJA-v1.0_transcript83572")

CriEndodermis$GeneName <- ifelse(CriEndodermis$CriID%in%SCR, "SCR",
                                 ifelse(CriEndodermis$CriID%in%SHR, "SHR", 
                                        ifelse(CriEndodermis$CriID%in%RBR, "RBR",
                                                      ifelse(CriEndodermis$CriID%in%SCL28, "SCL28", "CYCD"))))
CriEndodermis$Name <- paste(CriEndodermis$GeneName, "-", CriEndodermis$CriID)
CriEndodermis$Name <- gsub("Crichardii-Hnn_ARJA-v1.0_t", "T", CriEndodermis$Name)

tiff('CriDGE_Endodermis.tiff', units="cm", width=25, height=13, res=600, compression = 'lzw')
ggplot(CriEndodermis, aes(x=reorder(Name, GeneName), y = Tissue, color=logFC, size=-log10(as.numeric(adj.P.Val)))) +
  geom_point(alpha=1) +
  scale_colour_continuous_diverging(palette = "Purple-Green", limits = c(-7,7)) +
  scale_size_continuous(range = c(1,10)) +
  theme_minimal() +
  theme(panel.grid.major = element_line(size = 0.3, linetype = "dotted", colour = "gray", lineend = "butt")) + 
  theme(axis.text.x = element_text(angle = 60, hjust = 1, size = 8), axis.text.y = element_text(vjust = 1, size = 8)) +
  labs(title = "\nExpression of Ground Tissue-Specifying Orthologs", color = "Fold Change", size = "-log10(FDR)") +
  theme(plot.title = element_text(face = "bold")) +
  xlab("\nGenes") +
  theme(axis.title.x = element_text(face="bold"), axis.title.y = element_text(face="bold"), )
dev.off()




########### THE ROOT CAP ############
RootCap <- c("Crichardii-Hnn_ARJA-v1.0_transcript125422",
             "Crichardii-Hnn_ARJA-v1.0_transcript125421",
             "Crichardii-Hnn_ARJA-v1.0_transcript34113",
             "Crichardii-Hnn_ARJA-v1.0_transcript48177",
             "Crichardii-Hnn_ARJA-v1.0_transcript107875",
             "Crichardii-Hnn_ARJA-v1.0_transcript68964",
             "Crichardii-Hnn_ARJA-v1.0_transcript24375",
             "Crichardii-Hnn_ARJA-v1.0_transcript28691",
             "Crichardii-Hnn_ARJA-v1.0_transcript69831",
             "Crichardii-Hnn_ARJA-v1.0_transcript98223",
             "Crichardii-Hnn_ARJA-v1.0_transcript37933",
             "Crichardii-Hnn_ARJA-v1.0_transcript62331",
             "Crichardii-Hnn_ARJA-v1.0_transcript75969",
             "Crichardii-Hnn_ARJA-v1.0_transcript99946",
             "Crichardii-Hnn_ARJA-v1.0_transcript99950",
             "Crichardii-Hnn_ARJA-v1.0_transcript99947",
             "Crichardii-Hnn_ARJA-v1.0_transcript89165",
             "Crichardii-Hnn_ARJA-v1.0_transcript48864",
             "Crichardii-Hnn_ARJA-v1.0_transcript50074",
             "Crichardii-Hnn_ARJA-v1.0_transcript38379")
         
CUC <- c("Crichardii-Hnn_ARJA-v1.0_transcript84216",
         "Crichardii-Hnn_ARJA-v1.0_transcript42129",
         "Crichardii-Hnn_ARJA-v1.0_transcript129801",
         "Crichardii-Hnn_ARJA-v1.0_transcript99311")

SMB <- c("Crichardii-Hnn_ARJA-v1.0_transcript125422",
         "Crichardii-Hnn_ARJA-v1.0_transcript125421",
         "Crichardii-Hnn_ARJA-v1.0_transcript34113",
         "Crichardii-Hnn_ARJA-v1.0_transcript48177")

CR4 <- c("Crichardii-Hnn_ARJA-v1.0_transcript107875",
         "Crichardii-Hnn_ARJA-v1.0_transcript68964",
         "Crichardii-Hnn_ARJA-v1.0_transcript24375",
         "Crichardii-Hnn_ARJA-v1.0_transcript28691")

ARF10 <- c("Crichardii-Hnn_ARJA-v1.0_transcript69831",
           "Crichardii-Hnn_ARJA-v1.0_transcript98223",
           "Crichardii-Hnn_ARJA-v1.0_transcript37933",
           "Crichardii-Hnn_ARJA-v1.0_transcript62331",
           "Crichardii-Hnn_ARJA-v1.0_transcript75969")

PIN <- c("Crichardii-Hnn_ARJA-v1.0_transcript99946",
          "Crichardii-Hnn_ARJA-v1.0_transcript99950",
          "Crichardii-Hnn_ARJA-v1.0_transcript99947",
          "Crichardii-Hnn_ARJA-v1.0_transcript89165",
          "Crichardii-Hnn_ARJA-v1.0_transcript48864",
          "Crichardii-Hnn_ARJA-v1.0_transcript50074",
          "Crichardii-Hnn_ARJA-v1.0_transcript38379")

CriRootCapRAM <- subset(RAM_ResultsCriDGE, CriID %in% RootCap)
CriRootCapRAM <- cbind(CriRootCapRAM, Tissue = "Root Tip")

CriRootCapDiffRoot <- subset(DiffRoot_ResultsCriDGE, CriID %in% RootCap)
CriRootCapDiffRoot <- cbind(CriRootCapDiffRoot, Tissue = "Differentiated Root")

CriRootCapLeaves <- subset(Leaf_ResultsCriDGE, CriID %in% RootCap)
CriRootCapLeaves <- cbind(CriRootCapLeaves, Tissue = "Leaves and Shoot")

CriRootCap <- full_join(CriRootCapRAM,CriRootCapDiffRoot) %>%
  full_join(. ,CriRootCapLeaves)

CriRootCap$GeneName <- ifelse(CriRootCap$CriID%in%SMB, "SMB",
                              ifelse(CriRootCap$CriID%in%CUC, "CUC",
                                     ifelse(CriRootCap$CriID%in%FEZ, "FEZ",
                                            ifelse(CriRootCap$CriID%in%CR4, "CR4", 
                                                  ifelse(CriRootCap$CriID%in%ARF10, "ARF10/16","PIN")))))

CriRootCap$Name <- paste(CriRootCap$GeneName, "-", CriRootCap$CriID)
CriRootCap$Name <- gsub("Crichardii-Hnn_ARJA-v1.0_t", "T", CriRootCap$Name)


tiff('CriDGE_RootCap.tiff', units="cm", width=25, height=12, res=600, compression = 'lzw')
ggplot(CriRootCap, aes(x=reorder(Name, GeneName), y = Tissue, color=logFC, size=-log10(as.numeric(adj.P.Val)))) +
  geom_point(alpha=1) +
  scale_colour_continuous_diverging(palette = "Purple-Green", limits = c(-7,7)) +
  scale_size_continuous(range = c(1,10)) +
  theme_minimal() +
  theme(panel.grid.major = element_line(size = 0.3, linetype = "dotted", colour = "gray", lineend = "butt")) + 
  theme(axis.text.x = element_text(angle = 60, hjust = 1, size = 8), axis.text.y = element_text(vjust = 1, size = 8)) +
  labs(title = "\nExpression of Root Cap-Specifying Orthologs", color = "\nFold Change", size = "-log10(FDR)") +
  theme(plot.title = element_text(face = "bold")) +
  xlab("\nGenes") +
  theme(axis.title.x = element_text(face="bold"), axis.title.y = element_text(face="bold"), )
dev.off()




######### QUIESCENT CENTER ##########
QuiescentCenter <- c("Crichardii-Hnn_ARJA-v1.0_transcript69634",
                     "Crichardii-Hnn_ARJA-v1.0_transcript110487",
                     "Crichardii-Hnn_ARJA-v1.0_transcript110489",
                     "Crichardii-Hnn_ARJA-v1.0_transcript9747",
                     "Crichardii-Hnn_ARJA-v1.0_transcript124117",
                     "Crichardii-Hnn_ARJA-v1.0_transcript43706",
                     "Crichardii-Hnn_ARJA-v1.0_transcript69630",
                     "Crichardii-Hnn_ARJA-v1.0_transcript69632",
                     "Crichardii-Hnn_ARJA-v1.0_transcript69631",
                     "Crichardii-Hnn_ARJA-v1.0_transcript39955",
                     "Crichardii-Hnn_ARJA-v1.0_transcript134494",
                     "Crichardii-Hnn_ARJA-v1.0_transcript140100",
                     "Crichardii-Hnn_ARJA-v1.0_transcript86378",
                     "Crichardii-Hnn_ARJA-v1.0_transcript45206",
                     "Crichardii-Hnn_ARJA-v1.0_transcript86941")

PLT <- c("Crichardii-Hnn_ARJA-v1.0_transcript110487",
         "Crichardii-Hnn_ARJA-v1.0_transcript110489",
         "Crichardii-Hnn_ARJA-v1.0_transcript9747",
         "Crichardii-Hnn_ARJA-v1.0_transcript124117",
         "Crichardii-Hnn_ARJA-v1.0_transcript43706")

OtherPLT <- c("Crichardii-Hnn_ARJA-v1.0_transcript69630",
         "Crichardii-Hnn_ARJA-v1.0_transcript69632",
         "Crichardii-Hnn_ARJA-v1.0_transcript69631",
         "Crichardii-Hnn_ARJA-v1.0_transcript69634")

ANT <- c("Crichardii-Hnn_ARJA-v1.0_transcript124117",
         "Crichardii-Hnn_ARJA-v1.0_transcript43706")

QuiescentCenter <- c("Crichardii-Hnn_ARJA-v1.0_transcript39955",
                     "Crichardii-Hnn_ARJA-v1.0_transcript134494",
                     "Crichardii-Hnn_ARJA-v1.0_transcript140100",
                     "Crichardii-Hnn_ARJA-v1.0_transcript86378",
                     "Crichardii-Hnn_ARJA-v1.0_transcript45206",
                     "Crichardii-Hnn_ARJA-v1.0_transcript86941")

WOX13 <- c("Crichardii-Hnn_ARJA-v1.0_transcript39955",
         "Crichardii-Hnn_ARJA-v1.0_transcript134494",
         "Crichardii-Hnn_ARJA-v1.0_transcript140100")

WOX9  <- c("Crichardii-Hnn_ARJA-v1.0_transcript86378",
         "Crichardii-Hnn_ARJA-v1.0_transcript45206")

WUS <- c("Crichardii-Hnn_ARJA-v1.0_transcript86941")

CriQC_RAM <- subset(RAM_ResultsCriDGE, CriID %in% QuiescentCenter)
CriQC_RAM <- cbind(CriQC_RAM, Tissue = "Root Tip")

CriQC_DiffRoot <- subset(DiffRoot_ResultsCriDGE, CriID %in% QuiescentCenter)
CriQC_DiffRoot <- cbind(CriQC_DiffRoot, Tissue = "Differentiated Root")

CriQC_Leaves <- subset(Leaf_ResultsCriDGE, CriID %in% QuiescentCenter)
CriQC_Leaves <- cbind(CriQC_Leaves, Tissue = "Leaves and Shoot")

CriQC <- full_join(CriQC_RAM,CriQC_DiffRoot) %>%
  full_join(. ,CriQC_Leaves)

CriQC$GeneName <- ifelse(CriQC$CriID%in%PLT, "PLT",
                         ifelse(CriQC$CriID%in%ANT, "ANT",
                                ifelse(CriQC$CriID%in%OtherPLT, "Other-PLTs",
                                       ifelse(CriQC$CriID%in%WOX13, "WOX13",
                                              ifelse(CriQC$CriID%in%WOX9, "WOX9","WUS")))))
                
CriQC$GeneName <- ifelse(CriQC$CriID%in%WOX13, "WOX13",
                         ifelse(CriQC$CriID%in%WOX9, "WOX9","WUS"))

CriQC$Name <- paste(CriQC$GeneName, "-", CriQC$CriID)
CriQC$Name <- gsub("Crichardii-Hnn_ARJA-v1.0_t", "T", CriQC$Name)


tiff('CriDGE_QuiescentCenter.tiff', units="cm", width=15, height=13, res=600, compression = 'lzw')
ggplot(CriQC, aes(x=reorder(Name, GeneName), y = Tissue, color=logFC, size=-log10(as.numeric(adj.P.Val)))) +
  geom_point(alpha=1) +
  scale_colour_continuous_diverging(palette = "Purple-Green", limits = c(-7,7)) +
  scale_size_continuous(range = c(1,10)) +
  theme_minimal() +
  theme(panel.grid.major = element_line(size = 0.3, linetype = "dotted", colour = "gray", lineend = "butt")) + 
  theme(axis.text.x = element_text(angle = 60, hjust = 1, size = 8), axis.text.y = element_text(vjust = 1, size = 8)) +
  labs(title = "\nExpression of Stem Cell Maintenance Orthologs", color = "\nFold Change", size = "-log10(FDR)") +
  theme(plot.title = element_text(face = "bold")) +
  xlab("\nGenes") +
  theme(axis.title.x = element_text(face="bold"), axis.title.y = element_text(face="bold"), )
dev.off()




########### VASCULATURE ############
Vasculature <- c("Crichardii-Hnn_ARJA-v1.0_transcript108746",
                 "Crichardii-Hnn_ARJA-v1.0_transcript93770",
                 "Crichardii-Hnn_ARJA-v1.0_transcript48366",
                 "Crichardii-Hnn_ARJA-v1.0_transcript48502",
                 "Crichardii-Hnn_ARJA-v1.0_transcript79660",
                 "Crichardii-Hnn_ARJA-v1.0_transcript45096",
                 "Crichardii-Hnn_ARJA-v1.0_transcript42081",
                 "Crichardii-Hnn_ARJA-v1.0_transcript67498",
                 "Crichardii-Hnn_ARJA-v1.0_transcript68466",
                 "Crichardii-Hnn_ARJA-v1.0_transcript68469",
                 "Crichardii-Hnn_ARJA-v1.0_transcript85733",
                 "Crichardii-Hnn_ARJA-v1.0_transcript49412",
                 "Crichardii-Hnn_ARJA-v1.0_transcript35649",
                 "Crichardii-Hnn_ARJA-v1.0_transcript83583",
                 "Crichardii-Hnn_ARJA-v1.0_transcript68467",
                 "Crichardii-Hnn_ARJA-v1.0_transcript68468",
                 "Crichardii-Hnn_ARJA-v1.0_transcript62177",
                 "Crichardii-Hnn_ARJA-v1.0_transcript68465",
                 "Crichardii-Hnn_ARJA-v1.0_transcript88651",
                 "Crichardii-Hnn_ARJA-v1.0_transcript88652",
                 "Crichardii-Hnn_ARJA-v1.0_transcript97546",
                 "Crichardii-Hnn_ARJA-v1.0_transcript108746",
                 "Crichardii-Hnn_ARJA-v1.0_transcript93770",
                 "Crichardii-Hnn_ARJA-v1.0_transcript48366",
                 "Crichardii-Hnn_ARJA-v1.0_transcript48502",
                 "Crichardii-Hnn_ARJA-v1.0_transcript79660",
                 "Crichardii-Hnn_ARJA-v1.0_transcript45096",
                 "Crichardii-Hnn_ARJA-v1.0_transcript42081",
                 "Crichardii-Hnn_ARJA-v1.0_transcript67498",
                 "Crichardii-Hnn_ARJA-v1.0_transcript87948",
                 "Crichardii-Hnn_ARJA-v1.0_transcript85120",
                 "Crichardii-Hnn_ARJA-v1.0_transcript47519",
                 "Crichardii-Hnn_ARJA-v1.0_transcript46682",
                 "Crichardii-Hnn_ARJA-v1.0_transcript47161",
                 "Crichardii-Hnn_ARJA-v1.0_transcript36244",
                 "Crichardii-Hnn_ARJA-v1.0_transcript18743",
                 "Crichardii-Hnn_ARJA-v1.0_transcript30936",
                 "Crichardii-Hnn_ARJA-v1.0_transcript30483",
                 "Crichardii-Hnn_ARJA-v1.0_transcript76958",
                 "Crichardii-Hnn_ARJA-v1.0_transcript41977",
                 "Crichardii-Hnn_ARJA-v1.0_transcript49563",
                 "Crichardii-Hnn_ARJA-v1.0_transcript56914",
                 "Crichardii-Hnn_ARJA-v1.0_transcript117870",
                 "Crichardii-Hnn_ARJA-v1.0_transcript58539",
                 "Crichardii-Hnn_ARJA-v1.0_transcript45770",
                 "Crichardii-Hnn_ARJA-v1.0_transcript42025",
                 "Crichardii-Hnn_ARJA-v1.0_transcript79024",
                 "Crichardii-Hnn_ARJA-v1.0_transcript91098",
                 "Crichardii-Hnn_ARJA-v1.0_transcript127925",
                 "Crichardii-Hnn_ARJA-v1.0_transcript24867")

KAN <- c("Crichardii-Hnn_ARJA-v1.0_transcript108746",
         "Crichardii-Hnn_ARJA-v1.0_transcript93770",
         "Crichardii-Hnn_ARJA-v1.0_transcript48366",
         "Crichardii-Hnn_ARJA-v1.0_transcript48502",
         "Crichardii-Hnn_ARJA-v1.0_transcript79660",
         "Crichardii-Hnn_ARJA-v1.0_transcript45096",
         "Crichardii-Hnn_ARJA-v1.0_transcript42081",
         "Crichardii-Hnn_ARJA-v1.0_transcript67498")

PHB <- c("Crichardii-Hnn_ARJA-v1.0_transcript68466",
         "Crichardii-Hnn_ARJA-v1.0_transcript68469",
         "Crichardii-Hnn_ARJA-v1.0_transcript85733",
         "Crichardii-Hnn_ARJA-v1.0_transcript49412",
         "Crichardii-Hnn_ARJA-v1.0_transcript35649",
         "Crichardii-Hnn_ARJA-v1.0_transcript83583",
         "Crichardii-Hnn_ARJA-v1.0_transcript68467",
         "Crichardii-Hnn_ARJA-v1.0_transcript68468",
         "Crichardii-Hnn_ARJA-v1.0_transcript62177",
         "Crichardii-Hnn_ARJA-v1.0_transcript68465")

WOL <- c("Crichardii-Hnn_ARJA-v1.0_transcript88651",
         "Crichardii-Hnn_ARJA-v1.0_transcript88652",
         "Crichardii-Hnn_ARJA-v1.0_transcript97546")

TMO5 <- c("Crichardii-Hnn_ARJA-v1.0_transcript87948",
          "Crichardii-Hnn_ARJA-v1.0_transcript85120",
          "Crichardii-Hnn_ARJA-v1.0_transcript47519")

LHW <- c("Crichardii-Hnn_ARJA-v1.0_transcript46682",
         "Crichardii-Hnn_ARJA-v1.0_transcript47161",
         "Crichardii-Hnn_ARJA-v1.0_transcript36244")

LBD <- c("Crichardii-Hnn_ARJA-v1.0_transcript18743",
         "Crichardii-Hnn_ARJA-v1.0_transcript30936",
         "Crichardii-Hnn_ARJA-v1.0_transcript30483",
         "Crichardii-Hnn_ARJA-v1.0_transcript76958",
         "Crichardii-Hnn_ARJA-v1.0_transcript41977",
         "Crichardii-Hnn_ARJA-v1.0_transcript49563",
         "Crichardii-Hnn_ARJA-v1.0_transcript56914",
         "Crichardii-Hnn_ARJA-v1.0_transcript117870",
         "Crichardii-Hnn_ARJA-v1.0_transcript58539",
         "Crichardii-Hnn_ARJA-v1.0_transcript45770",
         "Crichardii-Hnn_ARJA-v1.0_transcript42025",
         "Crichardii-Hnn_ARJA-v1.0_transcript79024",
         "Crichardii-Hnn_ARJA-v1.0_transcript91098")

AHL <- c("Crichardii-Hnn_ARJA-v1.0_transcript127925",
         "Crichardii-Hnn_ARJA-v1.0_transcript24867")

CriVasculatureRAM <- subset(RAM_ResultsCriDGE, CriID %in% Vasculature)
CriVasculatureRAM <- cbind(CriVasculatureRAM, Tissue = "Root Tip")

CriVasculatureDiffRoot <- subset(DiffRoot_ResultsCriDGE, CriID %in% Vasculature)
CriVasculatureDiffRoot <- cbind(CriVasculatureDiffRoot, Tissue = "Differentiated Root")

CriVasculatureLeaves <- subset(Leaf_ResultsCriDGE, CriID %in% Vasculature)
CriVasculatureLeaves <- cbind(CriVasculatureLeaves, Tissue = "Leaves and Shoot")

CriVasculature <- full_join(CriVasculatureRAM,CriVasculatureDiffRoot) %>%
  full_join(. ,CriVasculatureLeaves)

CriVasculature$GeneName <- ifelse(CriVasculature$CriID%in%KAN, "KAN",
                              ifelse(CriVasculature$CriID%in%PHB, "C3-HD-Zip", 
                                     ifelse(CriVasculature$CriID%in%WOL, "WOL", 
                                                   ifelse(CriVasculature$CriID%in%TMO5, "TMO5",
                                                          ifelse(CriVasculature$CriID%in%LHW, "LHW",
                                                                 ifelse(CriVasculature$CriID%in%LBD, "LBD","AHL"))))))

CriVasculature$Name <- paste(CriVasculature$GeneName, "-", CriVasculature$CriID)


tiff('CriDGE_Vasculature.tiff', units="cm", width=35, height=12, res=600, compression = 'lzw')
ggplot(CriVasculature, aes(x=reorder(Name, GeneName), y = Tissue, color=logFC, size=-log10(as.numeric(adj.P.Val)))) +
  geom_point(alpha=1) +
  scale_colour_continuous_diverging(palette = "Purple-Green", limits = c(-7,7)) +
  scale_size_continuous(range = c(1,10)) +
  theme_minimal() +
  theme(panel.grid.major = element_line(size = 0.3, linetype = "dotted", colour = "gray", lineend = "butt")) + 
  theme(axis.text.x = element_text(angle = 60, hjust = 1, size = 8), axis.text.y = element_text(vjust = 1, size = 8)) +
  labs(title = "\nVasculature", color = "\nFold Change", size = "-log10(FDR)") +
  theme(plot.title = element_text(face = "bold")) +
  xlab("\nGenes") +
  theme(axis.title.x = element_text(face="bold"), axis.title.y = element_text(face="bold"), )
dev.off()





MYB <- c("Crichardii-Hnn_ARJA-v1.0_transcript38695",
         "Crichardii-Hnn_ARJA-v1.0_transcript38590",
         "Crichardii-Hnn_ARJA-v1.0_transcript35406",
         "Crichardii-Hnn_ARJA-v1.0_transcript76904",
         "Crichardii-Hnn_ARJA-v1.0_transcript34930",
         "Crichardii-Hnn_ARJA-v1.0_transcript146723",
         "Crichardii-Hnn_ARJA-v1.0_transcript77323",
         "Crichardii-Hnn_ARJA-v1.0_transcript26407",
         "Crichardii-Hnn_ARJA-v1.0_transcript127489",
         "Crichardii-Hnn_ARJA-v1.0_transcript41752",
         "Crichardii-Hnn_ARJA-v1.0_transcript24626",
         "Crichardii-Hnn_ARJA-v1.0_transcript115742",
         "Crichardii-Hnn_ARJA-v1.0_transcript85333",
         "Crichardii-Hnn_ARJA-v1.0_transcript22853",
         "Crichardii-Hnn_ARJA-v1.0_transcript45098",
         "Crichardii-Hnn_ARJA-v1.0_transcript141397",
         "Crichardii-Hnn_ARJA-v1.0_transcript84710",
         "Crichardii-Hnn_ARJA-v1.0_transcript46590",
         "Crichardii-Hnn_ARJA-v1.0_transcript27931",
         "Crichardii-Hnn_ARJA-v1.0_transcript115436",
         "Crichardii-Hnn_ARJA-v1.0_transcript154428",
         "Crichardii-Hnn_ARJA-v1.0_transcript46423",
         "Crichardii-Hnn_ARJA-v1.0_transcript84328",
         "Crichardii-Hnn_ARJA-v1.0_transcript23357",
         "Crichardii-Hnn_ARJA-v1.0_transcript146414",
         "Crichardii-Hnn_ARJA-v1.0_transcript136126",
         "Crichardii-Hnn_ARJA-v1.0_transcript21699")

CriVasculatureRAM <- subset(RAM_ResultsCriDGE, CriID %in% MYB)
CriVasculatureRAM <- cbind(CriVasculatureRAM, Tissue = "Root Tip")

CriVasculatureDiffRoot <- subset(DiffRoot_ResultsCriDGE, CriID %in% MYB)
CriVasculatureDiffRoot <- cbind(CriVasculatureDiffRoot, Tissue = "Differentiated Root")

CriVasculatureLeaves <- subset(Leaf_ResultsCriDGE, CriID %in% MYB)
CriVasculatureLeaves <- cbind(CriVasculatureLeaves, Tissue = "Leaves and Shoot")

CriVasculature <- full_join(CriVasculatureRAM,CriVasculatureDiffRoot) %>%
  full_join(. ,CriVasculatureLeaves)

CriVasculature$GeneName <- ifelse(CriVasculature$CriID%in%MYB, "MYB", NA)
CriVasculature$Name <- paste(CriVasculature$GeneName, "-", CriVasculature$CriID)


tiff('CriDGE_Vasculature2.tiff', units="cm", width=25, height=12, res=600, compression = 'lzw')
ggplot(CriVasculature, aes(x=reorder(Name, GeneName), y = Tissue, color=logFC, size=-log10(as.numeric(adj.P.Val)))) +
  geom_point(alpha=1) +
  scale_colour_continuous_diverging(palette = "Purple-Green", limits = c(-7,7)) +
  scale_size_continuous(range = c(1,10)) +
  theme_minimal() +
  theme(panel.grid.major = element_line(size = 0.3, linetype = "dotted", colour = "gray", lineend = "butt")) + 
  theme(axis.text.x = element_text(angle = 60, hjust = 1, size = 8), axis.text.y = element_text(vjust = 1, size = 8)) +
  labs(title = "\nVasculature", color = "\nFold Change", size = "-log10(FDR)") +
  theme(plot.title = element_text(face = "bold")) +
  xlab("\nGenes") +
  theme(axis.title.x = element_text(face="bold"), axis.title.y = element_text(face="bold"), )
dev.off()











FEZ <- c("Crichardii-Hnn_ARJA-v1.0_transcript109535",
         "Crichardii-Hnn_ARJA-v1.0_transcript49026",
         "Crichardii-Hnn_ARJA-v1.0_transcript63659",
         "Crichardii-Hnn_ARJA-v1.0_transcript96676",
         "Crichardii-Hnn_ARJA-v1.0_transcript51240",
         "Crichardii-Hnn_ARJA-v1.0_transcript94432",
         "Crichardii-Hnn_ARJA-v1.0_transcript51293",
         "Crichardii-Hnn_ARJA-v1.0_transcript51292",
         "Crichardii-Hnn_ARJA-v1.0_transcript42662",
         "Crichardii-Hnn_ARJA-v1.0_transcript89820",
         "Crichardii-Hnn_ARJA-v1.0_transcript57589",
         "Crichardii-Hnn_ARJA-v1.0_transcript41992")