library(readxl)
library(data.table)
library(qusage)
library(stringr)
library(lmerTest)
library(dplyr)
library(tidyr)

###################### Set up ######################

# Import reads
#reads <- as.data.frame(readxl::read_xlsx("D:/COVID_TISSUE_PROJECT/PlacentaPaper/GeoMX/Q3_norm_Filter_1_25_08.xlsx", sheet = "TargetCountMatrix"))
reads <- as.data.frame(readxl::read_xlsx("D:/COVID_TISSUE_PROJECT/PlacentaPaper/GeoMX/TG-AllTargets-Q3_NoFilter.xlsx", sheet = "TargetCountMatrix"))
rownames(reads) <- reads[,1]
reads <- reads[,2:ncol(reads)]

#nucleusCounts <- as.data.frame(readxl::read_xlsx("D:/COVID_TISSUE_PROJECT/PlacentaPaper/GeoMX/Q3_norm_Filter_1_25_08.xlsx", sheet = "SegmentProperties"))
nucleusCounts <- as.data.frame(readxl::read_xlsx("D:/COVID_TISSUE_PROJECT/PlacentaPaper/GeoMX/TG-AllTargets-Q3_NoFilter.xlsx", sheet = "SegmentProperties"))
#nucleusCounts <- nucleusCounts[,c("AOINucleiCount","AOISurfaceArea")]
surfaceArea <- nucleusCounts$AOISurfaceArea


# Import annotation file
annotationFile <- as.data.frame(readxl::read_xlsx("D:/COVID_TISSUE_PROJECT/PlacentaPaper/GeoMX/NS_annotationFile.xlsx", sheet = "SegmentProperties"))

#nucleusCounts <- nucleusCounts[!(annotationFile$Comment=="GeoMx excluded"),]
annotationFile <- annotationFile[!(annotationFile$Comment=="GeoMx excluded"),]

# Ensure annotation file and read file are referring to the same samples
setcolorder(reads, paste0(annotationFile$`Scan name`," | ",annotationFile$`ROI (label)`," | ",annotationFile$`Segment (Name/ Label)`))
identical(colnames(reads), paste0(annotationFile$`Scan name`," | ",annotationFile$`ROI (label)`," | ",annotationFile$`Segment (Name/ Label)`))

### To do ###
# - Boxplots of gene sets per region
# - DEG (multi-effects model using what was used in synovial sarcoma paper)
# - GO/GSEA
# - PCA/LDA
# - SpatialDecon

### QC

# Remove rows/columns from reads and annotation files that Matt has flagged
reads <- reads[,!(annotationFile$Comment!="NA")]
nucleusCounts <- nucleusCounts[!(annotationFile$Comment!="NA"),]
annotationFile <- annotationFile[!(annotationFile$Comment!="NA"),]

# Make new virus column

annotationFile$VirusPositive <- annotationFile$`Nucleocapsid IHC`
annotationFile$VirusPositive[annotationFile$`Nucleocapsid IHC` %in% c("Strong")] <- "Positive"
annotationFile$VirusPositive[annotationFile$`Nucleocapsid IHC` %in% c("Negative","Weak")] <- "Negative"


###################### Average area per compartment ######################

annotationFile$surfaceArea <- nucleusCounts$AOISurfaceArea
annotationFile$AOINucleiCount <- nucleusCounts$AOINucleiCount

tapply(annotationFile$surfaceArea, annotationFile$Compartment, mean)

tapply(annotationFile$AOINucleiCount, annotationFile$Compartment, mean)

table(annotationFile$DiseaseState2[annotationFile$Compartment == "Decidua"])
table(annotationFile$DiseaseState2[annotationFile$Compartment == "Villous stroma"])
table(annotationFile$DiseaseState2[annotationFile$Compartment == "Trophoblast"])
table(annotationFile$DiseaseState2[annotationFile$Compartment == "Macrophage"])

table(annotationFile$Compartment[annotationFile$DiseaseState2 == "COVID+ Placenta"])
table(annotationFile$Compartment[annotationFile$DiseaseState2 == "COVID+ Mother CHI"])
table(annotationFile$Compartment[annotationFile$DiseaseState2 == "COVID+ Mother CV"])
table(annotationFile$Compartment[annotationFile$DiseaseState2 == "CHI"])
table(annotationFile$Compartment[annotationFile$DiseaseState2 == "VUE"])
table(annotationFile$Compartment[annotationFile$DiseaseState2 == "Normal"])


###################### Analysis ######################

## Heatmap

#a <- str_detect(rownames(reads),"^IFIT")
#genes <- rownames(reads)[a == TRUE]

genes <- read.table("D:/COVID_TISSUE_PROJECT/QuantSeq/PlacentaQuantSeq/Macrophage_chemotaxis.txt", header = TRUE, sep = '\t', row.names=NULL)
genes <- genes$MGI.Gene.Marker.ID
genes <- toupper(genes)
genes <- as.data.frame(genes)
genes <- unique(genes[,1])

restrictionFactors <- c("IFITM3","ZBP1","MYD88","IFITM2","STAT2","NRN1","BST2","RETREG1","JADE2","IFIT1","MSR1","FNDC4",
                        "NT5C3A","ISG20","TMEM268","IFIT3","ETV6","TAGAP","TENT5C","TENT5A","SPATS2L","MAX","CCND3","RAB27A",
                        "UPP2","ST3GAL4","CLEC4D","GBP3","B4GALT5","FZD5","ARNTL","TRIM21","NAPA","CRP","APOL2","ERLIN1","CASP7","DDX60")

genes <- read.table(file = "D:/COVID_TISSUE_PROJECT/QuantSeq/PlacentaQuantSeq/allograftRejection.txt", sep = "\t")
genes <- genes[,1]

genes <- restrictionFactors

ROI <- c("Trophoblast")

readsHM <- reads[,annotationFile[,"Compartment"] %in% ROI]
annoHM <- annotationFile[annotationFile[,"Compartment"] %in% ROI,]

logTransformedData <- log2(readsHM+1)
scaledData <- as.data.frame(t(apply(na.omit(logTransformedData[genes,]), MARGIN = 1, scale)))

#scaledData <- na.omit(logTransformedData[genes,])

upperLimit <- 3
lowerLimit <- -3
scaledData[scaledData > upperLimit] <- upperLimit
scaledData[scaledData < lowerLimit] <- lowerLimit

ha = HeatmapAnnotation(
  Disease = annoHM$DiseaseState,
  Compartment = annoHM$Compartment,
  Virus = annoHM$`Nucleocapsid IHC`
)
split = annoHM$DiseaseState

ComplexHeatmap::Heatmap(as.matrix(scaledData), col = colorRampPalette(rev(brewer.pal(11, "RdBu")))(25), border = FALSE, rect_gp = gpar(col = "grey80", lwd = 1), name = "Correlation Coefficient",
                        column_title = "", column_title_gp = gpar(fontsize = 15, fontface = "bold"), clustering_distance_rows = "euclidean", clustering_distance_columns = "euclidean",
                        cluster_columns = FALSE, column_split = split,
                        top_annotation = ha,
                        cluster_rows = TRUE,
                        column_names_gp = grid::gpar(fontsize = 9),
                        row_names_gp = grid::gpar(fontsize = 9),
                        heatmap_legend_param = list(
                          title = "Normalised Counts", title_gp=gpar(fontsize=11, fontface="bold"),
                          legend_height = unit(6, "cm"), title_position = "leftcenter-rot"
                        )
)


### Boxplots of genes

regions <- unique(annotationFile$Compartment)
whichRegion <- 1
print(regions[whichRegion])

readsTrophoblast <- reads[,annotationFile[,"Compartment"] == regions[whichRegion]]
annoTrophoblast <- annotationFile[annotationFile[,"Compartment"] == regions[whichRegion],]

# SARS-CoV-2 entry genes
genes <- c("ACE2","BSG","NRP1","TMPRSS2")
# Collagen
b <- str_detect(rownames(reads),"^COL")
genes <- rownames(reads)[b == TRUE][1:10]
# B cells
genes <- c("CD19","MS4A1")
# T cells
genes <- c("CD3D","CD4","CD8A")
# Macrophage
genes <- c("CD68","CD163","MRC1","IFNG","IFNA2")
# Immune
genes <- c("PDCD1","CD274","TIGIT","HAVCR2","LAG3","CXCL10","CXCR3")


genes <- "GAPDH"

a <- as.data.frame(t(readsTrophoblast[genes,])) %>%
  gather(Channel, A, c(1:length(genes)))

a <- cbind(a,cbind(annoTrophoblast$DiseaseState2, annoTrophoblast$`Nucleocapsid IHC`))
colnames(a) <- c("Gene","Expression","Disease","Virus")

a$Virus[a$Virus %in% c("Strong","Weak")] <- "Positive"
a$Virus[a$Virus %in% c("NA")] <- "Negative"

#a$Disease[a$Disease == "COVID+ Placenta"] <- paste0(a$Disease[a$Disease == "COVID+ Placenta"],' - ', a$Virus[a$Disease == "COVID+ Placenta"] )

a$Disease <- factor(a$Disease, levels = c("COVID+ Placenta", "COVID+ Mother CHI", "COVID+ Mother CV","CHI","VUE","Normal"))


#a$Disease <- factor(a$Disease, levels = c("COVID+ Placenta - Positive", "COVID+ Placenta - Negative", "COVID+ Mother CHI", "COVID+ Mother CV","CHI","VUE","Normal"))


my_comparisons <- list( c("COVID+ Placenta", "COVID+ Mother CHI"), c("COVID+ Placenta", "COVID+ Mother CV"), c("COVID+ Placenta", "VUE"), c("COVID+ Placenta", "CHI"), c("COVID+ Placenta", "Normal") )
#my_comparisons_macrophage <- list( c("COVID+ Placenta", "COVID+ Mother"), c("COVID+ Placenta", "CHI"))

#my_comparisons <- list( c("COVID+ Placenta - Positive", "COVID+ Placenta - Negative"), c("COVID+ Placenta - Positive", "COVID+ Mother CHI"), c("COVID+ Placenta - Positive", "COVID+ Mother CV"))


#stat.test <- a %>%
#  wilcox_test(Expression ~ Disease) %>%
#  add_xy_position()


ggboxplot(a, x = "Disease", y = "Expression", ncol = 4, 
          palette = "nejm", fill = "Disease", color = "black", size = 1.05, scales = "free") +
  xlab("") + ylab("Normalised Expression") + theme_bw() + theme(axis.text.x = element_text(color = "black", size = 11),
                                                                     axis.text.y = element_text(color = "black", size = 11),
                                                                     axis.title.x = element_text(color = "black", size = 12),
                                                                     axis.title.y = element_text(color = "black", size = 13),
                                                                     legend.position = "none",
                                                                     axis.ticks = element_line(size = 1),
                                                                     axis.line = element_line(size = 1),
                                                                     panel.grid.major = element_blank(),
                                                                     panel.grid.minor = element_blank(),
                                                                     strip.background = element_blank(), strip.text = element_text(color = "black", size = 12)) +
  stat_compare_means(comparisons = my_comparisons, label = "p.signif", hide.ns = TRUE) + rotate_x_text(45)
  #stat_pvalue_manual(stat.test, hide.ns = TRUE)


ggboxplot(a, x = "Disease", y = "Expression",ncol = 6, palette = c("#D43F3A","#EEA236","#46B8DA","#F8DAAF","#AEE0EF","#CCCCCC"), fill = "Disease", color = "black", size = 1.05, scales = "free",
          bxp.errorbar = TRUE, bxp.errorbar.width = 0) +
  xlab("") + ylab(paste0(genes," Expression")) + ggdist::theme_ggdist() +
  geom_jitter(
    width = 0,
    color = "black",
    size = 2
  ) + theme(panel.grid.major = element_blank(),
            panel.grid.minor = element_blank(),
            legend.position = "none",
            axis.text.x = element_text(size = 18, colour = "black"),
            axis.text.y = element_text(size = 18, colour = "black"),
            axis.title.y = element_text(size = 20, colour = "black"),
            axis.title.x = element_text(size = 20, colour = "black"),
            axis.ticks.length=unit(.25, "cm"),
            axis.line.x = element_line(color = "black", linewidth = rel(0.5)),
            axis.line.y = element_line(color = "black", linewidth = rel(0.5)),
            axis.ticks = element_line(color="black"))  + rotate_x_text(60) +
  stat_compare_means(comparisons = my_comparisons, label = "p.signif", tip.length = 0)






### Combined boxplots of gene sets

genes <- read.table("D:/COVID_TISSUE_PROJECT/QuantSeq/PlacentaQuantSeq/Macrophage_activation.txt", header = TRUE, sep = '\t', row.names=NULL)
genes <- genes$MGI.Gene.Marker.ID
genes <- toupper(genes)
genes <- as.data.frame(genes)
genes <- unique(genes[,1])

#genes <- read.table(file = "D:/COVID_TISSUE_PROJECT/QuantSeq/PlacentaQuantSeq/Macrophage_activation.txt", sep = "\t")
genes <- read.table(file = "D:/COVID_TISSUE_PROJECT/QuantSeq/PlacentaQuantSeq/InterferonAlphaGeneSet.txt", sep = "\t")
genes <- genes[,1]

genes <- as.character(read.table(file = 'D:/COVID_TISSUE_PROJECT/QuantSeq/PlacentaQuantSeq/HSIAO_HOUSEKEEPING_GENES.tsv', sep = ','))
genes <- rownames(reads)[sample(1:nrow(reads),500, replace = F)]

genes <- read.table(file = "D:/COVID_TISSUE_PROJECT/QuantSeq/PlacentaQuantSeq/UV_RESPONSE.txt", sep = "\t", fill = TRUE)
genes <- genes[,1]

genes <- sample(rownames(reads),500)

#genes <- rownames(reads)[sample(1:nrow(reads),5, replace = FALSE)]
#genes <- restrictionFactors
#genes <- "ACE2"

ROI <- c("Macrophage")

readsBP <- reads[,annotationFile[,"Compartment"] %in% ROI]
annoBP <- annotationFile[annotationFile[,"Compartment"] %in% ROI,]

#readsBP <- reads
#annoBP <- annotationFile

annoBP$DS <- annoBP$DiseaseState2

annoBP$DS[annoBP$`Nucleocapsid IHC` %in% c("Strong")] <- "COVID+ Placenta, SARS-CoV-2+"

a <- colMeans(na.omit(readsBP[genes,]))

a <- cbind(a,annoBP$DS,annoBP$Number, annoBP$Compartment)
colnames(a) <- c("Expression","Disease","Sample","Compartment")
a <- as.data.frame(a)
a$Expression <- as.numeric(a$Expression)

#a <- a[a$Disease %in% c("COVID+ Placenta, SARS-CoV-2+","COVID+ Placenta","COVID+ Mother CHI","CHI"),]
#a <- a[a$Disease %in% c("COVID+ Placenta","COVID+ Mother CHI","COVID+ Mother CV"),]
#a <- a[a$Disease %in% c("COVID+ Placenta, SARS-CoV-2+","COVID+ Placenta","COVID+ Mother CHI","COVID+ Mother CV","CHI","VUE","Normal"),]


a$Disease <- factor(a$Disease, levels = c("COVID+ Placenta","COVID+ Placenta, SARS-CoV-2+","COVID+ Mother CHI","COVID+ Mother CV","CHI","VUE","Normal"))
#a$Disease <- factor(a$Disease, levels = c("COVID+ Placenta","COVID+ Placenta, SARS-CoV-2+","COVID+ Mother CHI","CHI"))

my_comparisons <- list( c("COVID+ Placenta", "CHI"),c("COVID+ Placenta, SARS-CoV-2+", "COVID+ Placenta"), c("COVID+ Placenta, SARS-CoV-2+", "COVID+ Mother CHI"),c("COVID+ Placenta, SARS-CoV-2+", "CHI"))
#my_comparisons <- list( c("COVID+ Placenta", "CHI"),c("COVID+ Placenta", "COVID+ Mother"))
#my_comparisons <- list( c("COVID+ Placenta", "CHI"),c("COVID+ Placenta", "COVID+ Mother"), c("COVID+ Placenta", "VUE"), c("COVID+ Placenta", "Normal") )
#my_comparisons <- list(c("COVID+ Placenta, SARS-CoV-2+", "COVID+ Placenta"),c("COVID+ Placenta, SARS-CoV-2+", "COVID+ Mother"))
#my_comparisons <- list(c("Macrophage", "Decidua"),c("Macrophage", "Villous stroma"),c("Macrophage", "Trophoblast"))
#my_comparisons <- list( c("COVID+ Placenta, SARS-CoV-2+", "COVID+ Placenta"), c("COVID+ Placenta, SARS-CoV-2+", "COVID+ Mother CHI"), c("COVID+ Placenta, SARS-CoV-2+", "CHI"))


#my_comparisons <- list(c("COVID+ Placenta", "COVID+ Placenta, SARS-CoV-2+"))
#my_comparisons <- list(c("COVID+ Placenta", "COVID+ Mother CHI"))


ggboxplot(a, x = "Disease", y = "Expression",ncol = 6, palette = "npg", fill = "Disease", color = "black", size = 1.1, scales = "free", # , palette = "locuszoom"
          bxp.errorbar = TRUE, bxp.errorbar.width = 0.3) +
  xlab("") + ylab("Gene Set Expression") + theme_bw() + theme(axis.text.x = element_text(color = "black", size = 11),
                                                              panel.grid.major = element_blank(),
                                                              panel.grid.minor = element_blank(),
                                                                     axis.text.y = element_text(color = "black", size = 11),
                                                                     axis.title.x = element_text(color = "black", size = 12),
                                                                     axis.title.y = element_text(color = "black", size = 13),
                                                                     legend.position = "none",
                                                                     axis.ticks = element_line(size = 1),
                                                                     axis.line = element_line(size = 1),
                                                                     strip.background = element_blank(), strip.text = element_text(color = "black", size = 12)) + rotate_x_text(60) +
  stat_compare_means(comparisons = my_comparisons) #, label = "p.signif")

genes[100]

## For correlation matrix

for(i in 1:length(unique(a$Sample))){

  print(paste0(unique(a$Sample)[i]," - ",mean(a$Expression[a$Sample == unique(a$Sample)[i]])))

}





#my_comparisons

a <- lmerTest::lmer(formula = Expression ~ Disease + (1|Sample),
                    data    = a)



### Statistical models - homemade

diseases <- c("COVID+ Placenta","COVID+ Mother")
randomAge <- sample(20:40, 177, replace = T)
statResults <- as.data.frame(phonTools::zeros(nrow(reads),4))
colnames(statResults) <- c("Gene","p-value","logChange","adjustedP")
statResults$Gene <- rownames(reads)

for(i in 1:nrow(reads)){
  dataToTest <- as.data.frame(cbind(t(reads[i,]),annotationFile$DiseaseState,randomAge, annotationFile$Number))
  dataToTest[,1] <- as.numeric(dataToTest[,1])
  dataToTest[,3] <- as.numeric(dataToTest[,3])
  colnames(dataToTest) <- c("gene","disease","age","sample")
  dataToTest <- dataToTest[dataToTest$disease %in% diseases,]

  a <- lmerTest::lmer(formula = gene ~ disease + age + (1|sample),
                      data    = dataToTest)
  statResults$logChange[i] <- log2(mean(dataToTest[dataToTest$disease == diseases[1],'gene'])/mean(dataToTest[dataToTest$disease == diseases[2],'gene']) )
  statResults$`p-value`[i] <- parameters::p_value(a)
}

statResults$adjustedP <- p.adjust(statResults$`p-value`, method = "bonferroni")


colnames(dataToTest) <- c("gene","disease","age","sample")
dataToTest <- dataToTest[dataToTest$disease %in% diseases,]

a <- lmerTest::lmer(formula = gene ~ disease + age + (1|sample),
               data    = dataToTest)

log2(mean(dataToTest[dataToTest$disease == diseases[1],'gene'])/mean(dataToTest[dataToTest$disease == diseases[2],'gene']) )

summary(a)
anova(a)
b <- parameters::p_value(a)

### Statistical models - glmmSeq

readsTrophoblast <- reads[,annotationFile[,"Compartment"] == regions[whichRegion]]
annoTrophoblast <- annotationFile[annotationFile[,"Compartment"] == regions[whichRegion],]

library(glmmSeq)

randomAge <- sample(20:40, 72, replace = T)
annoTrophoblast <- cbind(annoTrophoblast,randomAge)

# Filter out the two disease states of interest

diseases <- c("COVID+ Placenta","Normal")

readsTrophoblastModel <- readsTrophoblast[,annoTrophoblast$DiseaseState %in% diseases]
annoTrophoblastModel <- annoTrophoblast[annoTrophoblast$DiseaseState %in% diseases,]

annoTrophoblastModel$VirusPositive2 <- annoTrophoblastModel$`Nucleocapsid IHC`
annoTrophoblastModel$VirusPositive2[annoTrophoblastModel$`Nucleocapsid IHC` %in% c("Strong")] <- TRUE
annoTrophoblastModel$VirusPositive2[annoTrophoblastModel$`Nucleocapsid IHC` %in% c("Negative","Weak")] <- FALSE

# Calculate dispersion

disp <- apply(readsTrophoblastModel, 1, function(x){
  (var(x, na.rm=TRUE)-mean(x, na.rm=TRUE))/(mean(x, na.rm=TRUE)**2)
})

head(disp)

# Run model

statResults <- glmmSeq::glmmSeq(~ VirusPositive  + (1|Number),
                                countdata = readsTrophoblastModel,
                                metadata = annoTrophoblastModel,
                                id = "Number",
                                cores = 7,
                                dispersion = disp,
                                removeDuplicatedMeasures = FALSE,
                                removeSingles=FALSE,
                                progress = TRUE)


stats = data.frame(statResults@stats)


foldChange <- zeros(nrow(stats),1)

for(i in 1:nrow(stats)){
  foldChange[i,1] = log2(mean(as.numeric(readsTrophoblastModel[rownames(stats)[i],annoTrophoblastModel$VirusPositive == TRUE]))/mean(as.numeric(readsTrophoblastModel[rownames(stats)[i],annoTrophoblastModel$VirusPositive == FALSE])) )
}

# add to stats
# save stats as csv

#glmmSeq::fcPlot(glmmResult = statResults,
#       x1Label = "DiseaseState",
#       x2Label = "randomAge",
#       x2Values = c("COVID+ Placenta", "COVID+ Mother"),
#       pCutoff = 0.05,
#       useAdjusted = FALSE,
#       plotCutoff = 1,
#       graphics = "plotly")

stats <- cbind(stats,foldChange)

write.table(stats, file = "D:/COVID_TISSUE_PROJECT/PlacentaPaper/GeoMX/DEG/glmmSeq_COVID_placenta_virus_high_vs_virus_low_trophoblast.csv", sep = ",", quote = FALSE)

keyvals.colour <- as.character(zeros(nrow(stats)) * NA)

keyvals.colour[stats$foldChange > 0] <- "red"
keyvals.colour[stats$foldChange < 0] <- "deepskyblue"
keyvals.colour[stats$P_VirusPositive > 0.05] <- "grey80"
keyvals.colour[is.na(keyvals.colour)] <- 'grey80'

names(keyvals.colour) <- rownames(stats)

EnhancedVolcano(stats,
                lab = rownames(stats),
                x = 'foldChange',
                y = 'P_VirusPositive',
                cutoffLineType = 'blank',
                #selectLab = c('Nucleocapsid','Spike',"ORF1ab"),
                cutoffLineWidth = 0.8,
                colCustom = keyvals.colour,
                pointSize = 2.0,
                pCutoff = 0.05,
                FCcutoff = 0.5,
                labFace = 'bold',
                #boxedLabels = TRUE,
                #drawConnectors = TRUE,
                widthConnectors = 1.0,
                colConnectors = 'black',
                title = "",
                subtitle = "",
                caption = "",
                labSize = 4.0,
                colAlpha = 0.6,
                gridlines.major = FALSE,
                gridlines.minor = FALSE) + theme(legend.position = "none") + theme(
                  axis.title.x = element_text(size = 12),
                  axis.text.x = element_text(size = 12),
                  axis.title.y = element_text(size = 12),
                  axis.text.y = element_text(size = 12)) + xlab("Log2FC") # + ylim(0, max(abs(stats$foldChange)) + 0.5)


# DEseq 2

dds <- DESeqDataSetFromMatrix(countData = round(readsTrophoblastModel),
                              colData = annoTrophoblastModel,
                              design = ~ DiseaseState + Gestational_Age)

DESeq2::sizeFactors(dds) <- 1

dds <- DESeq(dds)
res <- results(dds, contrast=c("DiseaseState","COVID+ Placenta","Normal"))

write.table(res, file = "D:/COVID_TISSUE_PROJECT/PlacentaPaper/GeoMX/DEG/deSeq_COVID_placenta_decidua_vs_Normal_GestationalAgeCoVariate.csv", sep = ",", quote = FALSE)

keyvals.colour <- as.character(zeros(nrow(res)) * NA)
keyvals.colour[res$log2FoldChange > 0] <- "red"
keyvals.colour[res$log2FoldChange < 0] <- "deepskyblue"
keyvals.colour[res$pvalue > 0.05] <- "grey80"
keyvals.colour[is.na(keyvals.colour)] <- 'grey80'

names(keyvals.colour) <- rownames(res)

EnhancedVolcano(res,
                lab = rownames(res),
                x = 'log2FoldChange',
                y = 'pvalue',
                cutoffLineType = 'blank',
                selectLab = c('MMP11','TIMP2',"IL6ST"),
                cutoffLineWidth = 0.8,
                colCustom = keyvals.colour,
                pointSize = 2.0,
                pCutoff = 0.05,
                FCcutoff = 1.2,
                labFace = 'bold',
                #boxedLabels = TRUE,
                #drawConnectors = TRUE,
                widthConnectors = 1.0,
                colConnectors = 'black',
                title = "",
                subtitle = "",
                caption = "",
                labSize = 4.0,
                colAlpha = 0.6,
                gridlines.major = FALSE,
                gridlines.minor = FALSE) + theme(legend.position = "none") + theme(
                  axis.title.x = element_text(size = 12),
                  axis.text.x = element_text(size = 12),
                  axis.title.y = element_text(size = 12),
                  axis.text.y = element_text(size = 12)) + xlab("Log2FC") # + ylim(0, max(abs(res$log2FoldChange)) + 0.5)


## Spatial deconvolution ##

library(SpatialDecon)

per.observation.mean.neg = reads["NegProbe", ]
bg = sweep(reads * 0, 2, as.numeric(per.observation.mean.neg), "+")
dim(bg)

#placentaSC = read.csv("D:/COVID_TISSUE_PROJECT/PlacentaPaper/GeoMX/SingleCellData/CellProfileLibrary-master/FetalMaternal_Fetal_Placenta_10x.csv")
placentaSC = read.csv("D:/COVID_TISSUE_PROJECT/PlacentaPaper/GeoMX/SingleCellData/CellProfileLibrary-master/FetalMaternal_Adult_Decidua_10x.csv")
placentaSC = as.data.frame(placentaSC)
rownames(placentaSC) <- placentaSC[,1]
rownamesPlacenta <- rownames(placentaSC)
placentaSC <- placentaSC[,c(2:ncol(placentaSC))]
colnamesPlacenta <- colnames(placentaSC)
placentaSC = lapply(placentaSC, as.numeric)
placentaSC <- matrix(unlist(placentaSC), ncol = length(placentaSC))
rownames(placentaSC) <- rownamesPlacenta
colnames(placentaSC) <- colnamesPlacenta


res = spatialdecon(norm = as.matrix(reads),
                   bg = bg,
                   X = placentaSC,
                   align_genes = TRUE) #,
                   #cell_counts = annotationFile$nuclei)
str(res)

heatmapData <- as.data.frame(res$beta)

my_comparisons <- list( c("COVID+ Placenta", "COVID+ Mother"), c("COVID+ Placenta", "VUE"), c("COVID+ Placenta", "CHI"), c("COVID+ Placenta", "Normal") )

region <- "Macrophage"
cell <- "decidual.macrophage.type.1"

regionCellData <- heatmapData[,annotationFile$Compartment == region]

regionCellDataToPlot <- as.data.frame(cbind(t(regionCellData[cell,]),annotationFile[annotationFile$Compartment == region, c("DiseaseState")]))
colnames(regionCellDataToPlot) <- c("CellType","Disease")
regionCellDataToPlot$CellType <- as.numeric(regionCellDataToPlot$CellType)

ggboxplot(regionCellDataToPlot, x = "Disease", y = "CellType",ncol = 6, palette = "locuszoom", fill = "Disease", color = "black", size = 1.05, scales = "free",
          bxp.errorbar = TRUE, bxp.errorbar.width = 0.3, title = region) +
  xlab("") + ylab(paste0(cell," abundance")) + theme_classic() + theme(axis.text.x = element_text(color = "black", size = 11),
                                                                   axis.text.y = element_text(color = "black", size = 11),
                                                                   axis.title.x = element_text(color = "black", size = 12),
                                                                   axis.title.y = element_text(color = "black", size = 13),
                                                                   legend.position = "none",
                                                                   axis.ticks = element_line(size = 1),
                                                                   axis.line = element_line(size = 1),
                                                                   strip.background = element_blank(), strip.text = element_text(color = "black", size = 12)) + rotate_x_text(60) +
  stat_compare_means(comparisons = my_comparisons) #, label = "p.signif")




annoHM <- annotationFile[annotationFile$Compartment == region,]

ha = HeatmapAnnotation(
  Disease = annoHM$DiseaseState,
  Compartment = annoHM$Compartment,
  Virus = annoHM$`Nucleocapsid IHC`
)
split = annoHM$DiseaseState

ComplexHeatmap::Heatmap(as.matrix(regionCellData), col = colorRampPalette(rev(brewer.pal(11, "RdBu")))(25), border = FALSE, rect_gp = gpar(col = "grey80", lwd = 1), name = "Correlation Coefficient",
                        column_title = "", column_title_gp = gpar(fontsize = 15, fontface = "bold"), clustering_distance_rows = "euclidean", clustering_distance_columns = "euclidean",
                        cluster_columns = FALSE, column_split = split,
                        top_annotation = ha,
                        cluster_rows = TRUE,
                        column_names_gp = grid::gpar(fontsize = 9),
                        row_names_gp = grid::gpar(fontsize = 9),
                        heatmap_legend_param = list(
                          title = "Cell estimates", title_gp=gpar(fontsize=11, fontface="bold"),
                          legend_height = unit(6, "cm"), title_position = "leftcenter-rot"
                        )
)


heatmap(res$beta, cexCol = 0.5, cexRow = 0.7, margins = c(10,7))


rdecon = reverseDecon(norm = as.matrix(reads),
                      beta = res$beta)

heatmap(pmax(pmin(rdecon$resids, 2), -2))


plot(rdecon$cors, rdecon$resid.sd, col = 0)
showgenes = c("CXCL14", "LYZ", "NKG7")
text(rdecon$cors[setdiff(names(rdecon$cors), showgenes)],
     rdecon$resid.sd[setdiff(names(rdecon$cors), showgenes)],
     setdiff(names(rdecon$cors), showgenes), cex = 0.5)
text(rdecon$cors[showgenes], rdecon$resid.sd[showgenes],
     showgenes, cex = 0.75, col = 2)

###### Spatial deconvolution - immune #########

resImmune = spatialdecon(norm = as.matrix(reads),
                   bg = bg,
                   X = AA,
                   align_genes = TRUE,
                   cell_counts = nucleusCounts)

heatmap(resImmune$beta, cexCol = 0.5, cexRow = 0.7, margins = c(10,7))

heatmapData <- as.data.frame(resImmune$prop_of_all)

my_comparisons <- list( c("COVID+ Placenta", "COVID+ Mother"), c("COVID+ Placenta", "VUE"), c("COVID+ Placenta", "CHI"), c("COVID+ Placenta", "Normal") )
#my_comparisons <- list( c("COVID+ Placenta", "COVID+ Mother"), c("COVID+ Placenta", "CHI") )

regionList <- c("Macrophage","Villous stroma","Trophoblast","Decidua")

region <- regionList[4]
cell <- "B.naive"

regionCellData <- heatmapData[,annotationFile$Compartment == region]
colnames(regionCellData) <- annotationFile[annotationFile$Compartment == region, c("DiseaseState")]

regionCellDataToPlot <- as.data.frame(cbind(t(regionCellData[cell,]),annotationFile[annotationFile$Compartment == region, c("DiseaseState")]))
colnames(regionCellDataToPlot) <- c("CellType","Disease")
regionCellDataToPlot$CellType <- as.numeric(regionCellDataToPlot$CellType)

ggboxplot(regionCellDataToPlot, x = "Disease", y = "CellType",ncol = 6, palette = "locuszoom", fill = "Disease", color = "black", size = 1.05, scales = "free",
          bxp.errorbar = TRUE, bxp.errorbar.width = 0.3, title = region) +
  xlab("") + ylab(paste0(cell," abundance")) + theme_classic() + theme(axis.text.x = element_text(color = "black", size = 11),
                                                                       axis.text.y = element_text(color = "black", size = 11),
                                                                       axis.title.x = element_text(color = "black", size = 12),
                                                                       axis.title.y = element_text(color = "black", size = 13),
                                                                       legend.position = "none",
                                                                       axis.ticks = element_line(size = 1),
                                                                       axis.line = element_line(size = 1),
                                                                       strip.background = element_blank(), strip.text = element_text(color = "black", size = 12)) + rotate_x_text(60) +
  stat_compare_means(comparisons = my_comparisons) #, label = "p.signif")

# Stacked

regionCellDataStacked <- data.frame(type=rownames(regionCellData),setNames(stack(regionCellData),c("cell","disease")))
regionCellDataStacked$disease <- as.character(regionCellDataStacked$disease)
disease <- "VUE"
regionCellDataStacked$disease[regionCellDataStacked$disease %in% grep("VUE", as.character(regionCellDataStacked$disease), value=TRUE)] <- disease

ggboxplot(regionCellDataStacked, x = "disease", y = "cell",ncol = 6, palette = "locuszoom", fill = "disease", color = "black", facet.by = "type", size = 1.05, scales = "free",
          bxp.errorbar = TRUE, bxp.errorbar.width = 0.3, title = region) +
  xlab("") + ylab(paste0("Abundance")) + theme_classic() + theme(axis.text.x = element_text(color = "black", size = 11),
                                                                       axis.text.y = element_text(color = "black", size = 11),
                                                                       axis.title.x = element_text(color = "black", size = 12),
                                                                       axis.title.y = element_text(color = "black", size = 13),
                                                                       legend.position = "none",
                                                                       axis.ticks = element_line(size = 1),
                                                                       axis.line = element_line(size = 1),
                                                                       strip.background = element_blank(), strip.text = element_text(color = "black", size = 12)) + rotate_x_text(60) +
  stat_compare_means(comparisons = my_comparisons) #, label = "p.signif")

library(reshape)
melt(regionCellData, id=c("id","time"))

# Heatmap

annoHM <- annotationFile[annotationFile$Compartment == region,]

ha = HeatmapAnnotation(
  Disease = annoHM$DiseaseState,
  Compartment = annoHM$Compartment,
  Virus = annoHM$`Nucleocapsid IHC`
)
split = annoHM$DiseaseState

ComplexHeatmap::Heatmap(as.matrix(regionCellData), col = colorRampPalette(rev(brewer.pal(11, "RdBu")))(25), border = FALSE, rect_gp = gpar(col = "grey80", lwd = 1), name = "Correlation Coefficient",
                        column_title = "", column_title_gp = gpar(fontsize = 15, fontface = "bold"), clustering_distance_rows = "euclidean", clustering_distance_columns = "euclidean",
                        cluster_columns = FALSE, column_split = split,
                        top_annotation = ha,
                        cluster_rows = TRUE,
                        column_names_gp = grid::gpar(fontsize = 9),
                        row_names_gp = grid::gpar(fontsize = 9),
                        heatmap_legend_param = list(
                          title = "Cell estimates", title_gp=gpar(fontsize=11, fontface="bold"),
                          legend_height = unit(6, "cm"), title_position = "leftcenter-rot"
                        )
)

### Plotting deconvolution

# Stacked barchart

region <- c("Macrophage","Villous stroma","Trophoblast","Decidua")
i <- 4

TIL_barplot(resImmune$prop_of_all[,annotationFile$Compartment == region[i]],cex.names = 0.5, main = region[i]) #, draw_legend = TRUE)


barplot(resImmune$prop_of_all[,annotationFile$Compartment == region[i]], col = cellcols)
legend("topright", rownames(resImmune$prop_of_all[,annotationFile$Compartment == region[i]]), fill = cellcols, inset=c(-0.2,0))





# 2D projection
pc = prcomp(reads)$x[, c(1, 2)] # change this to highest variable genes like regular quantseq

par(mar = c(5,5,1,1))
layout(mat = (matrix(c(1, 2), 1)), widths = c(6, 2))
florets(x = pc[, 1], y = pc[, 2],
        b = resImmune$beta, cex = 2,
        xlab = "PC1", ylab = "PC2")
par(mar = c(0,0,0,0))
frame()
legend("center", fill = cellcols[rownames(resImmune$beta)],
       legend = rownames(resImmune$beta), cex = 0.7)


#### Reverse deconvolution ####

rdecon = reverseDecon(norm = as.matrix(reads),
                      beta = resImmune$beta)


heatmap(pmax(pmin(rdecon$resids, 2), -2))

# Base plot
plot(rdecon$cors, rdecon$resid.sd, col = 0, xlab="Correlation between predicted and observed expression", ylab="SD of residuals from predicted expression")
showgenes = c("CXCL14", "LYZ", "NKG7")
text(rdecon$cors[setdiff(names(rdecon$cors), showgenes)],
     rdecon$resid.sd[setdiff(names(rdecon$cors), showgenes)],
     setdiff(names(rdecon$cors), showgenes), cex = 0.5)
text(rdecon$cors[showgenes], rdecon$resid.sd[showgenes],
     showgenes, cex = 0.75, col = 2)



# Nicer plot
plot(rdecon$cors, rdecon$resid.sd, pch = 19, col = alpha("#0000FF", 0.2), xlab="Correlation between predicted and observed expression", ylab="SD of residuals from predicted expression",
     lwd=1.5)
abline(v = mean(na.omit(rdecon$cors)), col="black", lwd=1.8, lty=2)
abline(h = mean(na.omit(rdecon$resid.sd)), col="black", lwd=1.8, lty=2)
showgenes = c("HLA-C", "TIMP1", "COL3A1")
text(rdecon$cors[showgenes], rdecon$resid.sd[showgenes],
     showgenes, cex = 1, col = "black")



#### PCA ####

#pca <- prcomp(t(log2(reads)))

ntop = 500
rv <- rowVars(as.matrix(reads))

rv <- rowVars(as.matrix(reads[,annotationFile$Compartment == compartment & annotationFile$DiseaseState2 %in% c("COVID+ Placenta","COVID+ Mother CHI")]))

select <- order(rv, decreasing=TRUE)[seq_len(min(ntop, length(rv)))]

pca <- prcomp(t(log2(reads[select,annotationFile$Compartment == compartment & annotationFile$DiseaseState2 %in% c("COVID+ Placenta","COVID+ Mother CHI")])))

#pca <- glmpca(reads[select,annotationFile$Compartment == compartment & annotationFile$DiseaseState2 %in% c("COVID+ Placenta","COVID+ Mother CHI")],2)


pcaData <- as.data.frame(pca$x[,1:2])
#pcaData <- as.data.frame(pca$factors[,1])
pcaData <- cbind(pcaData,annotationFile$DiseaseState2[annotationFile$Compartment == compartment & annotationFile$DiseaseState2 %in% c("COVID+ Placenta","COVID+ Mother CHI")],
                 clusters = clusterTest2_covidCHITrophoblasts$communities,
                 patient = annotationFile$Number[annotationFile$Compartment == compartment & annotationFile$DiseaseState2 %in% c("COVID+ Placenta","COVID+ Mother CHI")],
                 virus = annotationFile$`Nucleocapsid IHC`[annotationFile$Compartment == compartment & annotationFile$DiseaseState2 %in% c("COVID+ Placenta","COVID+ Mother CHI")])
colnames(pcaData) <- c("PC1","PC2","Disease","Clusters","Patient","Virus")

pcaData$Clusters <- as.factor(pcaData$Clusters)
pcaData$Patient <- as.factor(pcaData$Patient)
pcaData$Virus <- as.factor(pcaData$Virus)

pcaData$Virus[pcaData$Virus == "NA"] <- "Negative"

eigs <- pca$sdev^2

b <- rbind(
  SD = sqrt(eigs),
  Proportion = eigs/sum(eigs),
  Cumulative = cumsum(eigs)/sum(eigs))

lineWidth <- 1
transparency <- 0.5
pointLines <- 1.5

pmain <- ggplot(pcaData, aes(x = PC1, y = PC2))+
  geom_point(aes(fill = Clusters), shape=21, color="black", size = 4, stroke=pointLines) + theme_bw() + theme(legend.title = element_blank(),
                                         panel.grid.major = element_blank(),
                                         panel.grid.minor = element_blank()) +
  xlab(paste0("PC 1 (",format(b[2,1]*100, digits = 3),"%)")) + ylab(paste0("PC 2 (",format(b[2,2]*100, digits = 3),"%)")) + scale_fill_lancet()

xdens <- axis_canvas(pmain, axis = "x")+
  geom_density(data = pcaData, aes(x = PC1, fill = Clusters),
               alpha = transparency, size = lineWidth) + scale_fill_lancet()

ydens <- axis_canvas(pmain, axis = "y", coord_flip = TRUE)+
  geom_density(data = pcaData, aes(x = PC2, fill = Clusters),
               alpha = transparency, size = lineWidth) +
  coord_flip() + scale_fill_lancet()

p1 <- insert_xaxis_grob(pmain, xdens, grid::unit(.2, "null"), position = "top")
p2<- insert_yaxis_grob(p1, ydens, grid::unit(.2, "null"), position = "right")
ggdraw(p2)

# Loadings
loadings <- as.data.frame(pca$rotation)

loadingsOrderedPC1 <- loadings[order(loadings$PC1),]

loadingsOrderedPC2 <- loadings[order(loadings$PC2),]

PC1genes <- cbind(rownames(loadingsOrderedPC1), loadingsOrderedPC1[,1])

PC2genes <- cbind(rownames(loadingsOrderedPC2), loadingsOrderedPC2[,2])

write.table(PC1genes, file = "D:/COVID_TISSUE_PROJECT/PlacentaPaper/GeoMX/PCA_Loadings/PC1_2000_genes_TrophoblastRegion_NS.txt", sep = "\t", quote = FALSE, row.names = FALSE)

write.table(PC2genes, file = "D:/COVID_TISSUE_PROJECT/PlacentaPaper/GeoMX/PCA_Loadings/PC2_2000_genes_TrophoblastRegion_NS.txt", sep = "\t", quote = FALSE, row.names = FALSE)

###### tSNE

compartment <- "Trophoblast"

clusterTest <- FastPG::fastCluster(t(reads[select,annotationFile$Compartment == compartment]))

clusterTest2_covidCHI <- FastPG::fastCluster(t(reads[select,annotationFile$Compartment == compartment & annotationFile$DiseaseState2 %in% c("COVID+ Placenta","COVID+ Mother CHI")]), k = 10)

#clusterTest2_covidCHI$communities[clusterTest2_covidCHI$communities == 1] <- 0
clusterTest2_covidCHI$communities[clusterTest2_covidCHI$communities == 3] <- 2

annotationFile$TrophoblastClustersCOVIDCHI <- "NA"
annotationFile$TrophoblastClustersCOVIDCHI[annotationFile$Compartment == compartment & annotationFile$DiseaseState2 %in% c("COVID+ Placenta","COVID+ Mother CHI")] <- clusterTest2_covidCHI$communities

# annotationFile[annotationFile$Compartment == compartment & annotationFile$DiseaseState2 %in% c("COVID+ Placenta","COVID+ Mother CHI"),]


saveRDS(clusterTest2_covidCHI,file = "D:/COVID_TISSUE_PROJECT/PlacentaPaper/RDS_Objects/macrophageRegionClusteringNS_COVID_CHI.rds")
saveRDS(clusterTest2_covidCHI,file = "D:/COVID_TISSUE_PROJECT/PlacentaPaper/RDS_Objects/trophoblastRegionClusteringNS_COVID_CHI.rds")

clusterTest2_covidCHITrophoblasts <- readRDS("D:/COVID_TISSUE_PROJECT/PlacentaPaper/RDS_Objects/trophoblastRegionClusteringNS_COVID_CHI.rds")

trophoblastClusteringMetadataImportant <- annotationFile[annotationFile$Compartment == "Trophoblast" & annotationFile$DiseaseState2 %in% c("COVID+ Placenta","COVID+ Mother CHI"),c("Scan name","ROI (label)","DiseaseState2","TrophoblastClustersCOVIDCHI")]
write.table(trophoblastClusteringMetadataImportant, file = "D:/COVID_TISSUE_PROJECT/PlacentaPaper/RDS_Objects/trophoblastClusterRegionsForAnnotation.csv", sep = ",")

trophoblastClusteringMetadataImportant <- read.table(file = "D:/COVID_TISSUE_PROJECT/PlacentaPaper/RDS_Objects/trophoblastClusterRegionsForAnnotation.csv", sep = ",")
macrophageMetadataImportant <- read.table(file = "D:/COVID_TISSUE_PROJECT/PlacentaPaper/RDS_Objects/macrophageClusterRegionsForAnnotation_Matt.csv", sep = ",")
macrophageMetadataImportant <- macrophageMetadataImportant[,c(1:3)]
colnames(macrophageMetadataImportant) <- macrophageMetadataImportant[1,]
macrophageMetadataImportant <- macrophageMetadataImportant[c(2:nrow(macrophageMetadataImportant)),]
#####
# Matching trophoblast regions to macrophage regions

trophoblastClusteringMetadataImportant$scanROI <- paste0(trophoblastClusteringMetadataImportant$Scan.name,'-',trophoblastClusteringMetadataImportant$ROI..label.)
macrophageMetadataImportant$scanROI <- paste0(macrophageMetadataImportant$`Scan name`,'-',macrophageMetadataImportant$`ROI (label)`)

trophoblastClusteringMetadataImportant <- trophoblastClusteringMetadataImportant[trophoblastClusteringMetadataImportant$scanROI %in% macrophageMetadataImportant$scanROI, ]
macrophageMetadataImportant <- macrophageMetadataImportant[macrophageMetadataImportant$scanROI %in% trophoblastClusteringMetadataImportant$scanROI, ]

cor(as.numeric(macrophageMetadataImportant$MacrophageClustersCOVIDCHI), as.numeric(trophoblastClusteringMetadataImportant$TrophoblastClustersCOVIDCHI))


trophoblastClusteringMetadataImportant$TrophoblastClustersCOVIDCHI[macrophageMetadataImportant$MacrophageClustersCOVIDCHI == 0]
trophoblastClusteringMetadataImportant$TrophoblastClustersCOVIDCHI[macrophageMetadataImportant$MacrophageClustersCOVIDCHI == 1]
trophoblastClusteringMetadataImportant$TrophoblastClustersCOVIDCHI[macrophageMetadataImportant$MacrophageClustersCOVIDCHI == 2]

#write.table(macrophageMetadataImportant, file = "D:/COVID_TISSUE_PROJECT/PlacentaPaper/RDS_Objects/macrophageClusterRegionsForAnnotation_Matt.csv", sep = ",")

#annotationFile$cluster <- clusterTest$communities

rtsne_out <- Rtsne::Rtsne(t(reads[select,annotationFile$Compartment == compartment]), verbose = TRUE, perplexity = 10, num_threads = 14)

rtsne_plot <- data.frame(x = rtsne_out$Y[,1],
                         y = rtsne_out$Y[,2],
                         Compartment = annotationFile$Compartment[annotationFile$Compartment == compartment],
                         COVID_positive = annotationFile$`Nucleocapsid IHC`[annotationFile$Compartment == compartment],
                         Disease = annotationFile$DiseaseState[annotationFile$Compartment == compartment],
                         clusters = clusterTest$communities,
                         Patient = annotationFile$Number[annotationFile$Compartment == compartment])

rtsne_plot$Compartment <- as.factor(rtsne_plot$Compartment)
rtsne_plot$COVID_positive <- as.factor(rtsne_plot$COVID_positive)
rtsne_plot$Disease <- as.factor(rtsne_plot$Disease)
rtsne_plot$clusters <- as.factor(rtsne_plot$clusters)
rtsne_plot$Patient <- as.factor(rtsne_plot$Patient)


colorIn <- rtsne_plot$Disease
colorLabel <- "Region"

  ggplot2::ggplot(rtsne_plot) + ggplot2::geom_point(aes(x=x, y=y, color=colorIn), size = 4) +
    guides(colour = guide_legend(override.aes = list(size=6))) +
    theme_bw() + labs(y= "tSNE 2", x = "tSNE 1")

##### Matched virus high vs.low #####

sampleIDs <- c(185,96,46,392,5458)
whichSample <- 5

regions <- unique(annotationFile$Compartment)
whichRegion <- 1

readsSingleSample <- reads[,annotationFile[,"Number"] == sampleIDs[whichSample]]
annoSingleSample <- annotationFile[annotationFile[,"Number"] == sampleIDs[whichSample],]

readsSingleSampleTrophoblast <- readsSingleSample[,annoSingleSample[,"Compartment"] == regions[whichRegion]]
annoSingleSampleTrophoblast <- annoSingleSample[annoSingleSample[,"Compartment"] == regions[whichRegion],]


dds <- DESeqDataSetFromMatrix(countData = round(readsSingleSampleTrophoblast),
                                colData = annoSingleSampleTrophoblast,
                                design = ~ VirusPositive)

DESeq2::sizeFactors(dds) <- 1

dds <- DESeq(dds)
res <- results(dds, contrast=c("VirusPositive","Positive","Negative"))

keyvals.colour <- as.character(zeros(nrow(res)) * NA)
keyvals.colour[res$log2FoldChange > 0] <- "red"
keyvals.colour[res$log2FoldChange < 0] <- "deepskyblue"
keyvals.colour[res$pvalue > 0.05] <- "grey80"
keyvals.colour[is.na(keyvals.colour)] <- 'grey80'

names(keyvals.colour) <- rownames(res)

EnhancedVolcano(res,
                lab = rownames(res),
                x = 'log2FoldChange',
                y = 'pvalue',
                cutoffLineType = 'blank',
                #selectLab = c('Nucleocapsid','Spike',"ORF1ab"),
                cutoffLineWidth = 0.8,
                colCustom = keyvals.colour,
                pointSize = 2.0,
                pCutoff = 0.05,
                FCcutoff = 1.2,
                labFace = 'bold',
                #boxedLabels = TRUE,
                #drawConnectors = TRUE,
                widthConnectors = 1.0,
                colConnectors = 'black',
                title = "",
                subtitle = "",
                caption = "",
                labSize = 4.0,
                colAlpha = 0.6,
                gridlines.major = FALSE,
                gridlines.minor = FALSE) + theme(legend.position = "none") + theme(
                  axis.title.x = element_text(size = 12),
                  axis.text.x = element_text(size = 12),
                  axis.title.y = element_text(size = 12),
                  axis.text.y = element_text(size = 12)) + xlab("Log2FC") # + ylim(0, max(abs(res$log2FoldChange)) + 0.5)


write.table(res, file = "D:/COVID_TISSUE_PROJECT/PlacentaPaper/GeoMX/DEG/deSeq_COVID_placenta_virus_high_vs_virus_low_trophoblast___Sample_5458.csv", sep = ",", quote = FALSE)

##### Trophoblast virus vs. normal #####

readsTrophoblast <- reads[,annotationFile[,"Compartment"] == regions[whichRegion]]
annoTrophoblast <- annotationFile[annotationFile[,"Compartment"] == regions[whichRegion],]

dds <- DESeqDataSetFromMatrix(countData = round(readsTrophoblast),
                              colData = annoTrophoblast,
                              design = ~ DiseaseState)

DESeq2::sizeFactors(dds) <- 1

dds <- DESeq(dds)
res <- results(dds, contrast=c("DiseaseState","COVID+ Placenta","CHI"))


keyvals.colour <- as.character(zeros(nrow(res)) * NA)
keyvals.colour[res$log2FoldChange > 0] <- "red"
keyvals.colour[res$log2FoldChange < 0] <- "deepskyblue"
keyvals.colour[res$pvalue > 0.05] <- "grey80"
keyvals.colour[is.na(keyvals.colour)] <- 'grey80'

names(keyvals.colour) <- rownames(res)

EnhancedVolcano(res,
                lab = rownames(res),
                x = 'log2FoldChange',
                y = 'pvalue',
                cutoffLineType = 'blank',
                #selectLab = c('Nucleocapsid','Spike',"ORF1ab"),
                cutoffLineWidth = 0.8,
                colCustom = keyvals.colour,
                pointSize = 2.0,
                pCutoff = 0.05,
                FCcutoff = 1.2,
                labFace = 'bold',
                #boxedLabels = TRUE,
                #drawConnectors = TRUE,
                widthConnectors = 1.0,
                colConnectors = 'black',
                title = "",
                subtitle = "",
                caption = "",
                labSize = 4.0,
                colAlpha = 0.6,
                gridlines.major = FALSE,
                gridlines.minor = FALSE) + theme(legend.position = "none") + theme(
                  axis.title.x = element_text(size = 12),
                  axis.text.x = element_text(size = 12),
                  axis.title.y = element_text(size = 12),
                  axis.text.y = element_text(size = 12)) + xlab("Log2FC") # + ylim(0, max(abs(res$log2FoldChange)) + 0.5)

write.table(res, file = "D:/COVID_TISSUE_PROJECT/PlacentaPaper/GeoMX/DEG/deSeq_COVID_placenta_vs_CHI_trophoblast.csv", sep = ",", quote = FALSE)

##### Which region has the highest interferons? #####

genes <- read.table("D:/COVID_TISSUE_PROJECT/QuantSeq/PlacentaQuantSeq/Macrophage_chemotaxis.txt", header = TRUE, sep = '\t', row.names=NULL)
genes <- genes$MGI.Gene.Marker.ID
genes <- toupper(genes)
genes <- as.data.frame(genes)
genes <- unique(genes[,1])


genes <- read.table(file = "D:/COVID_TISSUE_PROJECT/QuantSeq/PlacentaQuantSeq/InterferonAlphaGeneSet.txt", sep = "\t")
genes <- genes[,1]

b <- reads[genes,]
b <- colMeans(b)

c <- cbind(as.data.frame(b),annotationFile$Compartment)
colnames(c) <- c("GeneSet","Region")

#my_comparisons <- list( c("COVID+ Placenta", "COVID+ Mother"), c("COVID+ Placenta", "VUE"), c("COVID+ Placenta", "CHI"), c("COVID+ Placenta", "Normal") )

ggboxplot(c, x = "Region", y = "GeneSet",ncol = 6, palette = "npg", fill = "Region", color = "black", size = 1.05, scales = "free", title = "Interferon Alpha Signature") +
  xlab("") + ylab("Mean Geneset Expression") + theme_classic() + theme(axis.text.x = element_text(color = "black", size = 11),
                                                                     axis.text.y = element_text(color = "black", size = 11),
                                                                     axis.title.x = element_text(color = "black", size = 12),
                                                                     axis.title.y = element_text(color = "black", size = 13),
                                                                     legend.position = "none",
                                                                     axis.ticks = element_line(size = 1),
                                                                     axis.line = element_line(size = 1),
                                                                     strip.background = element_blank(), strip.text = element_text(color = "black", size = 12)) + rotate_x_text(45)

 # stat_compare_means(comparisons = my_comparisons_macrophage, label = "p.signif")


##### Quality control #####


keratins <- rownames(reads)

keratins <- keratins[str_detect(keratins, "KRT")]

geneOfInterest <- "KRT7"

b <- reads[geneOfInterest,]
b <- colMeans(b)
c <- as.data.frame(cbind(as.numeric(t(b)),annotationFile$Compartment))
colnames(c) <- c("Gene","Region")
c$Gene <- as.numeric(c$Gene)


#my_comparisons <- list( c("COVID+ Placenta", "COVID+ Mother"), c("COVID+ Placenta", "VUE"), c("COVID+ Placenta", "CHI"), c("COVID+ Placenta", "Normal") )

ggboxplot(c, x = "Region", y = "Gene",ncol = 6, palette = "npg", fill = "Region", color = "black", size = 1.05, scales = "free", title = "KRT7") +
  xlab("") + ylab("Normalised Expression") + theme_classic() + theme(axis.text.x = element_text(color = "black", size = 11),
                                                                       axis.text.y = element_text(color = "black", size = 11),
                                                                       axis.title.x = element_text(color = "black", size = 12),
                                                                       axis.title.y = element_text(color = "black", size = 13),
                                                                       legend.position = "none",
                                                                       axis.ticks = element_line(size = 1),
                                                                       axis.line = element_line(size = 1),
                                                                       strip.background = element_blank(), strip.text = element_text(color = "black", size = 12)) + rotate_x_text(45)





##### XXXX #####

