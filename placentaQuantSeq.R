library(DESeq2)
library(IHW)
library(biomaRt)
library(stringr)
library(apeglm)

### Load Data ###

dirName <- "F:/QuantSeq/COVID-19_Including_Placenta_and_Gary_Reynolds_Liver_April_2021/Placenta/counts/"

samples <- list.files(path = dirName)

for (j in 1:length(samples)){
  if (j == 1){
    aa <- read.table(paste0(dirName,samples[j]), sep = "\t")
    gene_names <- aa[,1]
    aa$V1 <- NULL
  } else {
    bb <- read.table(paste0(dirName,samples[j]), sep = "\t")
    bb$V1 <- NULL
    aa <- cbind(aa,bb)
  }
}

sampleNames <-  substr(samples,1,nchar(samples)-4)

rownames(aa) <- gene_names
colnames(aa) <- sampleNames

### Organise - Remove metrics, split human and viral genes, rename genes, re-merge viral and human, save ###

numViralGenes <- 12
n <- dim(aa)[1]
rawReads <- aa[1:(n-5),]
readMetrics <- aa[(n-5):n,]
virusReads <- rawReads[(nrow(rawReads)-(numViralGenes-1)):nrow(rawReads),]
humanReads <- rawReads[1:(nrow(rawReads)-(numViralGenes)),]

rownames(virusReads) <- c("ORF1ab","ORF1a","Spike","Nucleocapsid","ORF3a","Membrane Glycoprotein","ORF8","ORF7a","Envelope Protein","ORF6","ORF7b","ORF10")

bb <- convertEnsembleToHGNC(inputMatrix = humanReads)

humanReads <- bb
rm(bb)


totalReads <- rbind(humanReads,virusReads)

write.table(totalReads, file = "F:/QuantSeq/COVID-19_Including_Placenta_and_Gary_Reynolds_Liver_April_2021/Placenta/totalReads_preReplicateAddition.txt", sep = "\t", quote = FALSE)

## Merge technical replicates ##

uniqueSamples <-  unique(substr(sampleNames,1,6))

for(i in 1:length(uniqueSamples)){

  sampleCases <- grep(uniqueSamples[i], colnames(totalReads))

  sampleTemp <- totalReads[,sampleCases]

  sampleTemp <- rowSums(sampleTemp)

  if(i == 1){

    combinedReads <- as.data.frame(sampleTemp)
    colnames(combinedReads) <- uniqueSamples[i]

  } else {
    combinedReads <- cbind(combinedReads,sampleTemp)
    colnames(combinedReads)[i] <- uniqueSamples[i]

  }

}


write.table(combinedReads, file = "F:/QuantSeq/COVID-19_Including_Placenta_and_Gary_Reynolds_Liver_April_2021/Placenta/combinedReads.txt", sep = "\t", quote = FALSE)

combinedReads <- read.table(file = "F:/QuantSeq/COVID-19_Including_Placenta_and_Gary_Reynolds_Liver_April_2021/Placenta/combinedReads.txt", sep = "\t")

### Import into DESeq2, normalise ###

library(readxl)

metadata <- read_xlsx("F:/QuantSeq/COVID-19_Including_Placenta_and_Gary_Reynolds_Liver_April_2021/Placenta/metadataPlacenta.xlsx")

colnames(combinedReads) <- metadata$...1

combinedReads <- combinedReads[1:(nrow(combinedReads)-12),]

dds <- DESeqDataSetFromMatrix(countData = combinedReads,
                              colData = metadata,
                              design = ~ condition2)

dds <- estimateSizeFactors(dds)

sizeBG <- sizeFactors(dds)

normalized_counts <- counts(dds, normalized=TRUE)

write.table(normalized_counts, file = "F:/QuantSeq/COVID-19_Including_Placenta_and_Gary_Reynolds_Liver_April_2021/Placenta/combinedNormalisedReads.txt", sep = "\t", quote = FALSE)

normalized_counts <- read.table(file = "F:/QuantSeq/COVID-19_Including_Placenta_and_Gary_Reynolds_Liver_April_2021/Placenta/combinedNormalisedReads.txt", sep = "\t")
#rownames(normalized_counts) <- normalized_counts[,1]
# <- normalized_counts[,c(2:ncol(normalized_counts))]
colnames(normalized_counts) <- normalized_counts[1,]
normalized_counts <- normalized_counts[c(2:nrow(normalized_counts)),]
genes <- rownames(normalized_counts)
normalized_counts <- as.data.frame(normalized_counts)
normalized_counts <- apply(normalized_counts, 2, as.numeric)
rownames(normalized_counts) <- genes


### Compare viral reads, total viral reads, and breakdown ###

virusReadsNormed <- normalized_counts[(nrow(normalized_counts)-11):nrow(normalized_counts),]
#rownames(virusReadsNormed) <- virusReadsNormed[,1]
virusReadsNormed <- virusReadsNormed[,c(1:ncol(virusReadsNormed))]
virusReadsNormed <- as.data.frame(apply(virusReadsNormed, MARGIN = 2, as.numeric))

#virusReadsNormed <- virusReadsNormed[,c(1:15,17:21)]

viralReadsTotal <- colSums(virusReadsNormed)

viralReadsTotal <- as.data.frame(viralReadsTotal)

viralReadsTotal <- cbind(viralReadsTotal,metadata)

colnames(viralReadsTotal) <- c("Reads","Sample","Disease","Disease2","Disease3")

viralReadsTotal$Reads[16] <- 0

viralReadsTotal$Disease <- factor(viralReadsTotal$Disease, levels = c("COVID_Placenta", "COVID_Mother", "CHI","VUE","Normal"))
viralReadsTotal$Disease2 <- factor(viralReadsTotal$Disease2, levels = c("COVID_Placenta", "COVID_Mother_CHI","COVID_Mother_CV", "CHI","VUE","Normal"))
viralReadsTotal$Disease3 <- factor(viralReadsTotal$Disease3, levels = c("COVID_CHI", "COVID_CV", "CHI","VUE","Normal"))

my_comparisons <- list(c("COVID_Placenta", "COVID_Mother_CHI"), c("COVID_Placenta", "COVID_Mother_CV"))
my_comparisons <- list(c("COVID_Placenta", "COVID_Mother"))
my_comparisons <- list(c("COVID_CHI", "COVID_CV"))

ggboxplot(viralReadsTotal, x = "Disease3", y = "Reads", fill = "Disease3", size = 0.95, palette = "nejm") + theme_classic() + xlab("") + ylab("Total Viral Reads") +
  rotate_x_text(45) + theme(legend.position = "none") +
  theme(panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        axis.line = element_line(colour = 'black', size = 1),
        panel.background = element_rect(fill = "#faffff",
                                        colour = "black",
                                        size = 0.5, linetype = "solid"),
        axis.ticks = element_line(colour = "black", size = 1),
        axis.title.x= element_text(size = 16),
        axis.text.x= element_text(size = 14, colour = "black"),
        axis.ticks.x= element_blank(),
        axis.text.y = element_text(size = 14, colour = "black"),
        axis.title.y = element_text(size = 16)) +
  stat_compare_means(comparisons = my_comparisons) +
  scale_x_discrete(labels=c("Mother+CHI", "Mother+CV", "CHI", "VUE","Normal"))
my_comparisons <- list(c("COVID_Placenta", "COVID_Mother"))

ggboxplot(viralReadsTotal, x = "Disease", y = "Reads", fill = "Disease", size = 1, add = "jitter", palette = "nejm") + theme_classic() + xlab("") + ylab("Total Viral Reads") +
  rotate_x_text(45) + theme(legend.position = "none") +
  theme(
    axis.ticks = element_line(size = 1),
    axis.line = element_line(size = 1),
    axis.text.x = element_text(color = "black", size = 11),
    axis.text.y = element_text(color = "black", size = 11),
    axis.title.x = element_text(color = "black", size = 13),
    axis.title.y = element_text(color = "black", size = 13)
  ) +      # Add global p-value
  stat_compare_means(label = "p.signif", method = "wilcoxon", label.y = 1600, comparisons = my_comparisons)


# Virus heatmap

#virusReadsLogged <- log2(virusReadsNormed[,metadata$condition %in% c("COVID_Mother","COVID_Placenta")] + 1)
virusReadsLogged <- log2(virusReadsNormed + 1)

#colnames(virusReadsLogged) = c("P152","P153-1","P153-2","P154","P155-1","P155-2","P156-1","P156-2","P158-1","P158-2","P159-1","P159-2")

metadata$condition

rownames(virusReadsLogged) <- c("ORF1ab","ORF1a","Spike","Nucleocapsid","ORF3a","Membrane","ORF8","ORF7a","Envelope","ORF6","ORF7b","ORF10")
#colnames(virusReadsLogged) <- metadata$...1[c(1:15,17:21)]
colnames(virusReadsLogged) <- metadata$...1
#virusReadsLogged[virusReadsLogged > 3] <- 3


metadata$condition[metadata$condition == "COVID_Mother"] <- "Mother+Placenta-"
metadata$condition[metadata$condition == "COVID_Placenta"] <- "Mother+Placenta+"

#split = metadata$condition[metadata$condition %in% c("COVID_Mother","COVID_Placenta")]
split = metadata$condition2

virusReadsLogged[virusReadsLogged < 2] <- 0

virusReadsLogged[c(1,3,4),16] <- 0

#metadata$condition2 <- factor(metadata$condition2, levels = c("COVID_Placenta",
#                                                              "COVID_Mother_CHI",
#                                                              "COVID_Mother_CV",
#                                                              "CHI",
#                                                              "VUE",
#                                                              "Normal"))

metadata$condition3 <- factor(metadata$condition3, levels = c("COVID_CHI",
                                                              "COVID_CV",
                                                              "CHI",
                                                              "VUE",
                                                              "Normal"))

#split = metadata$condition[metadata$condition %in% c("COVID_Mother","COVID_Placenta")]
split <- metadata$condition3

#virusReadsNormed[metadata$condition == "Normal"]

ComplexHeatmap::Heatmap(as.matrix(virusReadsLogged), col = brewer.pal(9,"Blues"), border = FALSE, rect_gp = gpar(col = "grey80", lwd = 2), name = "log2(Normalized Counts + 1)",   #grey80
                        column_title = "", column_title_gp = gpar(fontsize = 12, fontface = "bold"), clustering_distance_rows = "euclidean", clustering_distance_columns = "euclidean",
                        cluster_columns = FALSE, column_split = split,
                        cluster_rows = TRUE,show_column_names = FALSE,
                        column_names_gp = grid::gpar(fontsize = 11),
                        row_dend_side = c("right"),
                        row_names_side = c("left"),
                        row_names_gp = grid::gpar(fontsize = 11),
                        heatmap_legend_param = list(
                          title = "log2(Normalized Counts + 1)", title_gp=gpar(fontsize=11, fontface="bold"),
                          legend_height = unit(6, "cm"), title_position = "leftcenter-rot"
                        )
)

ComplexHeatmap::Heatmap(as.matrix(virusReadsLogged), col = viridis(100), border = FALSE, rect_gp = gpar(col = "grey80", lwd = 2), name = "log2(Normalized Counts + 1)",
                        column_title = "", column_title_gp = gpar(fontsize = 12, fontface = "bold"), clustering_distance_rows = "euclidean", clustering_distance_columns = "euclidean",
                        cluster_columns = FALSE, column_split = split,
                        cluster_rows = TRUE,show_column_names = FALSE,
                        column_names_gp = grid::gpar(fontsize = 11),
                        row_names_gp = grid::gpar(fontsize = 11),
                        heatmap_legend_param = list(
                          title = "log2(Normalized Counts + 1)", title_gp=gpar(fontsize=11, fontface="bold"),
                          legend_height = unit(6, "cm"), title_position = "leftcenter-rot"
                        )
)



# map_signif_level = TRUE



virusReadsLong <- as.data.frame(virusReadsNormed) %>% gather(Genes, Reads, 1:21)

virusReadsLong <- cbind(virusReadsLong,rep(rownames(virusReadsNormed), times = 21))

colnames(virusReadsLong) <- c("Samples","Reads","Gene")

ggbarplot(virusReadsLong, x = "Gene", y = "Reads", fill = "Samples") + theme_bw() + xlab("Genes") + ylab("Viral Reads") + rotate_x_text(45)

### Differential expression testing ###

dds <- DESeq(dds)
res <- results(dds, contrast=c("condition2","COVID_Placenta","COVID_Mother_CHI"))
#res <- results(dds, contrast=c("condition2","CHI","VUE"))

write.csv(res,"D:/COVID_TISSUE_PROJECT/QuantSeq/PlacentaQuantSeq/DEG/CHI_vs_CV.csv")

res <- read.csv("D:/COVID_TISSUE_PROJECT/QuantSeq/PlacentaQuantSeq/DEG/COVID_Placenta_vs_COVID_Mother.csv")

#rownames(res) <- res[,1]
#res <- res[,c(2:ncol(res))]

# reduce number of dots that are not significant

res2 <- res
res <- res[!is.na(res$padj),]
notSignif <- res$padj > 0.05
perc.70 <- round(sum(notSignif) * 0.7)
button.5 <- which(notSignif == TRUE)
sampled.70 <- sample(button.5, perc.70)
res <- res[-sampled.70, ]


keyvals.colour <- as.character(zeros(nrow(res)) * NA)

keyvals.colour[res$log2FoldChange > 0] <- "red"
keyvals.colour[res$log2FoldChange < 0] <- "deepskyblue"
keyvals.colour[res$pvalue > 0.05] <- "grey80"
keyvals.colour[is.na(keyvals.colour)] <- 'grey80'

names(keyvals.colour)[keyvals.colour == 'red'] <- 'Early-stage'
names(keyvals.colour)[keyvals.colour == 'deepskyblue'] <- 'Late-Stage'
names(keyvals.colour)[keyvals.colour == 'grey80'] <- 'NS'

#selectLab = c('VCAM1','KCTD12','ADAM12',
#              'CXCL12','CACNB2','SPARCL1','DUSP1','SAMHD1','MAOA')

EnhancedVolcano(res,
                lab = rownames(res),
                x = 'log2FoldChange',
                y = 'padj',
                cutoffLineType = 'blank',
                cutoffLineWidth = 0.8,
                colCustom = keyvals.colour,
                selectLab = c('GBP1','IER3','IFI27','IL1RN',
                              'CXCL10','IFIT2','PLIN2','HLA-DRB6','MX2','MT-ND4L'),
                pointSize = 1.5,
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
                colAlpha = 0.7,
                gridlines.major = FALSE,
                gridlines.minor = FALSE) + theme_bw() + theme(legend.position = "none") + theme(
                  panel.grid.major = element_blank(),
                  panel.grid.minor = element_blank(),
                  axis.title.x = element_text(size = 12),
                  axis.text.x = element_text(size = 12),
                  axis.title.y = element_text(size = 12),
                  axis.text.y = element_text(size = 12)) + xlab("Log2FC") + ylim(0, max(-log10(res$padj)) + 0.5)




### Principal component analysis and loadings ###



res.pca <- prcomp(t(normalized_counts), scale = FALSE)

a <- as.data.frame(res.pca$x[,1:2])
a <- cbind(a,metadata$condition)
colnames(a) <- c("x","y","Disease")

eigs <- res.pca$sdev^2

b <- rbind(
  SD = sqrt(eigs),
  Proportion = eigs/sum(eigs),
  Cumulative = cumsum(eigs)/sum(eigs))

pmain <- ggplot(a, aes(x = x, y = y))+
  geom_point(aes(fill = Disease), shape=21, color="black", size = 4) +
  scale_fill_nejm() + theme_bw() + theme(legend.title = element_blank()) +
  xlab(paste0("PC 1 (",format(b[2,1]*100, digits = 3),"%)")) + ylab(paste0("PC 2 (",format(b[2,2]*100, digits = 3),"%)"))

xdens <- axis_canvas(pmain, axis = "x")+
  geom_density(data = a, aes(x = x, fill = Disease),
               alpha = 0.9, size = 0.2)+
  scale_fill_nejm()

ydens <- axis_canvas(pmain, axis = "y", coord_flip = TRUE)+
  geom_density(data = a, aes(x = y, fill = Disease),
               alpha = 0.9, size = 0.2)+
  coord_flip()+
  scale_fill_nejm()

p1 <- insert_xaxis_grob(pmain, xdens, grid::unit(.2, "null"), position = "top")
p2<- insert_yaxis_grob(p1, ydens, grid::unit(.2, "null"), position = "right")
ggdraw(p2)






vsd <- vst(dds, blind=FALSE)

pcaData <- plotPCA(vsd, intgroup=c("condition"), ntop = 1000, returnData = TRUE)

percentVar <- round(100 * attr(pcaData, "percentVar"))
ggplot(pcaData, aes(PC1, PC2, color=condition)) +
  geom_point(size=3) +
  xlab(paste0("PC1: ",percentVar[1],"% variance")) +
  ylab(paste0("PC2: ",percentVar[2],"% variance"))


pcaData$condition2 <- metadata$condition2

lineWidth <- 1
transparency <- 0.5
pointLines <- 1.5

qualpalette <- qualpalr::qualpal(n=6)

pmain <- ggplot(pcaData, aes(x = PC1, y = PC2))+
  geom_point(aes(fill = condition2), shape=21, color="black", size = 4, stroke=pointLines) + theme_bw() + theme(legend.title = element_blank(),
                                         panel.grid.major = element_blank(),
                                         panel.grid.minor = element_blank()) +
  xlab(paste0("PC 1 (",format(percentVar[1], digits = 3),"%)")) + ylab(paste0("PC 2 (",format(percentVar[2], digits = 3),"%)")) + scale_fill_manual(values = qualpalette$hex)

xdens <- axis_canvas(pmain, axis = "x")+
  geom_density(data = pcaData, aes(x = PC1, fill = condition),
               alpha = transparency, size = lineWidth)+
  scale_fill_manual(values = qualpalette$hex)

ydens <- axis_canvas(pmain, axis = "y", coord_flip = TRUE)+
  geom_density(data = pcaData, aes(x = PC2, fill = condition),
               alpha = transparency, size = lineWidth)+
  coord_flip() + scale_fill_manual(values = qualpalette$hex)

p1 <- insert_xaxis_grob(pmain, xdens, grid::unit(.2, "null"), position = "top")
p2<- insert_yaxis_grob(p1, ydens, grid::unit(.2, "null"), position = "right")
ggdraw(p2)




ntop = 1000

rv <- rowVars(assay(vsd))

select <- order(rv, decreasing=TRUE)[seq_len(min(ntop, length(rv)))]

# GLMPCA / Freeman-Tukey

# Freeman-Tukey
library(expandFunctions)
pca <- expandFunctions::freemanTukey(assay(vsd)[select,])
pcaFT <- as.data.frame(matrix(pca,ncol =21,byrow = T))
pca <- prcomp(t(pcaFT))
pcaDataFT <- as.data.frame(pca$x[,1:2])
pcaDataFT$condition <- metadata$condition2
pcaDataFT$patientID <- metadata$...1

# GLMPCA
pca <- glmpca::glmpca(assay(vsd)[select,],4)
pcaDataGLM <- pca$factors
pcaDataGLM$condition <- metadata$condition2
pcaDataGLM$patientID <- metadata$...1

pmain <- ggplot(pcaDataFT, aes(x = PC1, y = PC2))+
  geom_point(aes(fill = condition), shape=21, color="black", size = 4, stroke=pointLines) + theme_bw() + theme(legend.title = element_blank(),
                                                                                                                panel.grid.major = element_blank(),
                                                                                                                panel.grid.minor = element_blank()) + scale_fill_manual(values = qualpalette$hex)

xdens <- axis_canvas(pmain, axis = "x")+
  geom_density(data = pcaDataFT, aes(x = PC1, fill = condition),
               alpha = transparency, size = lineWidth)+
  scale_fill_manual(values = qualpalette$hex)

ydens <- axis_canvas(pmain, axis = "y", coord_flip = TRUE)+
  geom_density(data = pcaDataFT, aes(x = PC2, fill = condition),
               alpha = transparency, size = lineWidth)+
  coord_flip() + scale_fill_manual(values = qualpalette$hex)

p1 <- insert_xaxis_grob(pmain, xdens, grid::unit(.2, "null"), position = "top")
p2<- insert_yaxis_grob(p1, ydens, grid::unit(.2, "null"), position = "right")
ggdraw(p2)




# + xlab(paste0("PC 2 (",format(percentVar[1], digits = 3),"%)")) + ylab(paste0("PC 3 (",format(percentVar[2], digits = 3),"%)"))






pca <- prcomp(t(assay(vsd)[select,]))

loadings <- as.data.frame(pca$rotation)

loadingsOrderedPC1 <- loadings[order(loadings$PC1),]

loadingsOrderedPC2 <- loadings[order(loadings$PC2),]


PC1genes <- cbind(rownames(loadingsOrderedPC1), loadingsOrderedPC1[,1])

PC2genes <- cbind(rownames(loadingsOrderedPC2), loadingsOrderedPC2[,2])


write.table(PC1genes, file = "D:/COVID_TISSUE_PROJECT/QuantSeq/PlacentaQuantSeq/PC1genes_JustHuman.txt", sep = "\t", quote = FALSE, row.names = FALSE)

write.table(PC2genes, file = "D:/COVID_TISSUE_PROJECT/QuantSeq/PlacentaQuantSeq/PC2genes_JustHuman.txt", sep = "\t", quote = FALSE, row.names = FALSE)


write.csv(loadings,"D:/COVID_TISSUE_PROJECT/QuantSeq/PlacentaQuantSeq/QS_loadings.csv", row.names = TRUE)


loadingsOrdered <- loadings[order(loadings[,1], loadings[,2]),]






genesOfInterest <- c("COL1A1","COL1A2","COL3A1","COL4A1","COL6A1","COL6A3","COL15A1")

genesOfInterest <- c("IL17A","IL6","IL2","CXCL8","TNF","IFNA1","IFNG","IL1B","NRP1","ACE2")

normalized_counts <- as.data.frame(normalized_counts)

genesOfInterestCounts <- normalized_counts[genesOfInterest,]

genesOfInterestCountsLong <- as.data.frame(genesOfInterestCounts) %>% gather(Genes, Reads, 1:21)

genesOfInterestCountsLong <- cbind(genesOfInterestCountsLong,rep(genesOfInterest, times = 21))

genesOfInterestCountsLong <- cbind(genesOfInterestCountsLong,rep(metadata$condition, each = length(genesOfInterest)))

colnames(genesOfInterestCountsLong) <- c("Samples","Reads","Gene","Disease")

ggboxplot(genesOfInterestCountsLong, x = "Disease", y = "Reads", facet.by = "Gene", scales = "free", ncol = length(genesOfInterest), fill = "Disease", palette = "nejm") + theme_bw() + rotate_x_text(45)


# GO and gene lists

dds <- DESeq(dds)

  a <- bulkGO(DESeqObject = dds,comparisonsLeft = c("COVID_Placenta","COVID_Placenta","COVID_Placenta","COVID_Placenta","COVID_Mother","COVID_Mother","COVID_Mother"),
            comparisonsRight = c("COVID_Mother","CHI","VUE","Normal","CHI","VUE","Normal"),
            plotHeatmap = FALSE, directory = "D:/COVID_TISSUE_PROJECT/QuantSeq/PlacentaQuantSeq", reduceTerms = TRUE, saveGeneList = TRUE)

## virDirect

  # Load in virus counts

  dirName <- "F:/QuantSeq/COVID-19_Including_Placenta_and_Gary_Reynolds_Liver_April_2021/virDirect_Placenta/counts/"

  samples <- list.files(path = dirName)

  for (j in 1:length(samples)){
    if (j == 1){
      aa <- t(read.table(paste0(dirName,samples[j]), sep = "\t"))
      aa <- as.data.frame(aa[c(2:length(aa)),])
      #gene_names <- aa[,1]
      #aa$V1 <- NULL
    } else {
      bb <- t(read.table(paste0(dirName,samples[j]), sep = "\t"))
      bb <- as.data.frame(bb[c(2:length(bb)),])
      #bb$V1 <- NULL
      aa <- cbind(aa,bb)
    }
  }

  sampleNames <-  substr(samples,1,nchar(samples)-30)

  colnames(aa) <- sampleNames

  # Get virus names from SAM file

  c <- read.table(file = "F:/QuantSeq/COVID-19_Including_Placenta_and_Gary_Reynolds_Liver_April_2021/virDirect_Placenta/alignedVirus/129_20_placenta_L001_ds_merged/129_20_placenta_L001_ds_mergedAligned.out.sam", header = TRUE, fill = TRUE, row.names = NULL)

  c <- c[c(1:1894),]

  library(stringi)

  d <- stri_split_fixed(c$X.HD,"|")

  e <- c$X.HD

  for(i in 1:length(d)){

    e[i] <- d[[i]][5]

  }

  rownames(aa) = make.names(e, unique=TRUE)

  totalReads <- aa

  # make numeric

  for(i in 1:ncol(totalReads)){
  totalReads[,i] <- as.numeric(totalReads[,i])
  }


  # remove duplicates

  uniqueSamples <-  unique(substr(sampleNames,1,6))

  for(i in 1:length(uniqueSamples)){

    sampleCases <- grep(uniqueSamples[i], colnames(totalReads))

    sampleTemp <- totalReads[,sampleCases]

    sampleTemp <- rowSums(sampleTemp)

    if(i == 1){

      combinedReads <- as.data.frame(sampleTemp)
      colnames(combinedReads) <- uniqueSamples[i]

    } else {

      combinedReads <- cbind(combinedReads,sampleTemp)
      colnames(combinedReads)[i] <- uniqueSamples[i]

    }

  }

  # Normalise counts

  metadata <- read_xlsx("F:/QuantSeq/COVID-19_Including_Placenta_and_Gary_Reynolds_Liver_April_2021/Placenta/metadataPlacenta.xlsx")

  colnames(combinedReads) <- metadata$...1

  dds <- DESeqDataSetFromMatrix(countData = combinedReads,
                                colData = metadata,
                                design = ~ condition)

  dds <- estimateSizeFactors(dds)

  sizeFactors(dds)

  normalized_counts <- counts(dds, normalized=TRUE)

  dds <- DESeq(dds)
  res <- results(dds, contrast=c("condition","COVID_Placenta","VUE"))



  keyvals.colour <- as.character(zeros(nrow(res)) * NA)

  keyvals.colour[res$log2FoldChange > 0] <- "red"
  keyvals.colour[res$log2FoldChange < 0] <- "deepskyblue"
  keyvals.colour[res$pvalue > 0.05] <- "grey80"
  keyvals.colour[is.na(keyvals.colour)] <- 'grey80'

  names(keyvals.colour)[keyvals.colour == 'red'] <- 'Early-stage'
  names(keyvals.colour)[keyvals.colour == 'deepskyblue'] <- 'Late-Stage'
  names(keyvals.colour)[keyvals.colour == 'grey80'] <- 'NS'

  #selectLab = c('VCAM1','KCTD12','ADAM12',
  #              'CXCL12','CACNB2','SPARCL1','DUSP1','SAMHD1','MAOA')

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
                    axis.text.y = element_text(size = 12)) + xlab("Log2FC") + ylim(0, 3)



  i = 10

  top5 <- as.data.frame(normalized_counts[order(normalized_counts[,i], decreasing = TRUE)[1:5],i])
  top5 <- cbind(top5,rownames(top5))
  colnames(top5) <- c("Reads","Virus")


  ggbarplot(top5, x = "Virus", y = "Reads", fill = "Virus", palette = "jco") + theme(legend.position = "none") + rotate_x_text(90) + ylab("Reads") + xlab("")



  normalized_counts[rownames(normalized_counts) %like% "Human_herpesvirus_5", ]

  # Human_papillomavirus_type_126,

  #ebv <- cbind(normalized_counts["Human_herpesvirus_4_complete_wild_type_genome",], colnames(normalized_counts), metadata$condition)
  ebv <- cbind(normalized_counts["Human_papillomavirus_type_126",], colnames(normalized_counts), metadata$condition)
  colnames(ebv) <- c("Reads","Sample","Condition")
  ebv <- as.data.frame(ebv)
  ebv$Reads <- as.numeric(ebv$Reads)

  ggboxplot(ebv, x = "Condition", y = "Reads", fill = "Condition", palette = "jco", size = 1.3, alpha = 0.3, add = "jitter") + theme(legend.position = "none") +
    rotate_x_text(45) + ylab("Reads") + xlab("")


  library(data.table)


  hpv <- normalized_counts[rownames(normalized_counts) %like% "Human_herpesvirus_5", ]
  hpv <- colSums(hpv)
  hpv <- cbind(hpv, colnames(normalized_counts), metadata$condition)
  colnames(hpv) <- c("Reads","Sample","Condition")
  hpv <- as.data.frame(hpv)
  hpv$Reads <- as.numeric(hpv$Reads)

  ggboxplot(hpv, x = "Condition", y = "Reads", fill = "Condition", palette = "jco", size = 1.3, alpha = 0.3, add = "jitter") + theme(legend.position = "none") +
    rotate_x_text(45) + ylab("CMV Reads") + xlab("")


## Inflammatory and Interferon Heatmap ##

interferonAlpha <- read.table("D:/COVID_TISSUE_PROJECT/QuantSeq/PlacentaQuantSeq/Macrophage_chemotaxis.txt", header = TRUE, sep = '\t', row.names=NULL)
interferonAlpha <- interferonAlpha$MGI.Gene.Marker.ID
interferonAlpha <- toupper(interferonAlpha)
interferonAlpha <- as.data.frame(interferonAlpha)
interferonAlpha <- unique(interferonAlpha[,1])


#inflammatory <- read.table(file = "D:/COVID_TISSUE_PROJECT/QuantSeq/PlacentaQuantSeq/inflammatoryGeneSet.txt", sep = "\t")

interferonGamma <- read.table(file = "D:/COVID_TISSUE_PROJECT/QuantSeq/PlacentaQuantSeq/InterferonGammaGeneSet.txt", sep = "\t")

interferonAlpha <- read.table(file = "D:/COVID_TISSUE_PROJECT/QuantSeq/PlacentaQuantSeq/InterferonAlphaGeneSet.txt", sep = "\t")

interferonBeta <- as.character(read.table(file = "D:/COVID_TISSUE_PROJECT/QuantSeq/PlacentaQuantSeq/DER_IFN_BETA_RESPONSE_UP.tsv", sep = ","))

inflamTotal <- rbind(interferonAlpha, interferonGamma, interferonBeta)

#inflamTotal <- interferonAlpha


inflamTotal <- unique(inflamTotal)

inflamTotal <- inflamTotal[,1]

filteredInflam <- normalized_counts[inflamTotal,]


c <- 0

for(i in 1:nrow(filteredInflam)){

  if(sum(filteredInflam[i,]) == 0){
    c <- append(c, i)
  }


}

c <- c[2:length(c)]


filteredInflam <- filteredInflam[-c, ]

#filteredInflamZ <- log2(filteredInflam+1)
#filteredInflamZ <- log2(t(apply(filteredInflam, MARGIN = 1, scale))+1)
filteredInflamZ <- t(apply(filteredInflam, MARGIN = 1, scale))


filteredInflamZ[is.na(filteredInflamZ)] <- 0

filteredInflamZ[filteredInflamZ > 2] <- 2
filteredInflamZ[filteredInflamZ < -2] <- -2


#split = factor(metadata$condition, levels = c("COVID_Placenta","COVID_Mother","CHI","VUE","Normal"))
split = factor(metadata$condition2, levels = c("COVID_Placenta","COVID_Mother_CHI","COVID_Mother_CV","CHI","VUE","Normal"))

#column_ha = HeatmapAnnotation(Condition = metadata$condition)


ids <- c(1,4,5,15,16,29,31,32,34,35,36,39,41,42,55,81,82,105,109,110,111,112,113,136)

labels <-  rownames(filteredInflamZ)[ids]
ha = rowAnnotation(foo = anno_mark(at = ids, labels = labels))

ComplexHeatmap::Heatmap(as.matrix(filteredInflamZ), col=viridis(100), border = TRUE, name = "Scaled Expression", row_names_gp = gpar(fontsize = 8),
                        column_title = "", column_title_gp = gpar(fontsize = 15, fontface = "bold"), column_split = split, right_annotation = ha,
                        clustering_distance_rows = "euclidean", cluster_columns = FALSE, heatmap_legend_param = list(
                        title = "Normalised Expression", title_gp=gpar(fontsize=11, fontface="bold"),
                        legend_height = unit(6, "cm"), title_position = "leftcenter-rot"))



# , bottom_annotation = column_ha
# col=viridis(100)
# colorRampPalette(c("navy", "white", "firebrick3"))(50)

# rect_gp = gpar(col = "black", lwd = 1)

zscore <- function(inputData){
  meanb <- mean(inputData)

  stdb <- sd(inputData)

  zScore <- ((inputData-meanb)/stdb)

  return(zScore)

}

## GSEA prep ##


expressionFileGSEA <- function(data, fileOut){

  gseaExpressionFile <- cbind(rownames(data), rep(NA, each = length(rownames(data))), data)

  gseaExpressionFile <- rbind(c("NAME","description",colnames(data)),gseaExpressionFile)

  write.table(gseaExpressionFile, file = fileOut, sep = "\t", row.names = FALSE, col.names = FALSE, quote = FALSE)

}


expressionFileGSEA(data = normalized_counts, fileOut = "D:/COVID_TISSUE_PROJECT/QuantSeq/PlacentaQuantSeq/GSEA/Expression.gct.txt")


write.table(t(as.data.frame(metadata$condition)), file = "D:/COVID_TISSUE_PROJECT/QuantSeq/PlacentaQuantSeq/GSEA/phenotype.txt", sep = "\t", row.names = FALSE, col.names = FALSE, quote = FALSE)

### Combined boxplots of gene sets

genes <- read.table("D:/COVID_TISSUE_PROJECT/QuantSeq/PlacentaQuantSeq/Immunosuppression.txt", header = TRUE, sep = '\t', row.names=NULL)
genes <- genes$MGI.Gene.Marker.ID
genes <- toupper(genes)
genes <- as.data.frame(genes)
genes <- unique(genes[,1])

genes <- read.table(file = "D:/COVID_TISSUE_PROJECT/QuantSeq/PlacentaQuantSeq/InterferonAlphaGeneSet.txt", sep = "\t")
genes <- genes[,1]


genes <- as.character(read.table(file = 'D:/COVID_TISSUE_PROJECT/QuantSeq/PlacentaQuantSeq/HSIAO_HOUSEKEEPING_GENES.tsv', sep = ','))

genes <- restrictionFactors

a <- colMeans(na.omit(normalized_counts[genes,]))

a <- cbind(a,metadata$condition2, metadata$...1)
colnames(a) <- c("Expression","Disease","Patient")
a <- as.data.frame(a)
a$Expression <- as.numeric(a$Expression)

#a <- a[a$Disease %in% c("Mother+Placenta+","Mother+Placenta-CHI"),]

#a$Disease <- factor(a$Disease, levels = c("COVID_Placenta", "COVID_Mother", "CHI", "VUE", "Normal"))
a$Disease <- factor(a$Disease, levels = c("COVID_Placenta", "COVID_Mother_CHI", "COVID_Mother_CV", "CHI", "VUE", "Normal"))

#my_comparisons <- list( c("COVID_Placenta", "COVID_Mother"), c("COVID_Placenta", "CHI"), c("COVID_Placenta", "VUE"), c("COVID_Placenta", "Normal") )
my_comparisons <- list( c("COVID_Placenta", "COVID_Mother_CHI"),c("COVID_Placenta", "CHI"), c("COVID_Placenta", "Normal"))
#my_comparisons <- list( c("COVID_Placenta", "CHI"), c("COVID_Placenta", "Normal") )

ggboxplot(a, x = "Disease", y = "Expression",ncol = 6, palette = "locuszoom", fill = "Disease", color = "black", size = 1.05, scales = "free",
          bxp.errorbar = TRUE, bxp.errorbar.width = 0.3) +
  xlab("")   + ylab("Gene Set Expression") + theme_bw() + theme(axis.text.x = element_text(color = "black", size = 11),
                                                                   axis.text.y = element_text(color = "black", size = 11),
                                                                   axis.title.x = element_text(color = "black", size = 12),
                                                                   axis.title.y = element_text(color = "black", size = 13),
                                                                   legend.position = "none",
                                                                   axis.ticks = element_line(size = 1),
                                                                   axis.line = element_line(size = 1),
                                                                    panel.grid.major = element_blank(),
                                                                panel.grid.minor = element_blank(),
                                                                   strip.background = element_blank(), strip.text = element_text(color = "black", size = 12)) + rotate_x_text(60) +
  stat_compare_means(comparisons = my_comparisons, label = "p.signif")

# restriction factor vs. virus load

a <- cbind(a, viralReadsTotal)
colnames(a)[4] <- "VirusReads"
a <- a[a$Disease %in% c("COVID_Placenta"),]

ggscatter(a, x = "Expression", y = "VirusReads",
          fill = "red", color = "black", shape = 19, size = 3, # Points color, shape and size
          add = "reg.line",  # Add regressin line - reg.line
          add.params = list(color = "blue", fill = "lightgray"), # Customize reg. line
          conf.int = TRUE, # Add confidence interval
          cor.coef = TRUE, # Add correlation coefficient. see ?stat_cor
          cor.coeff.args = list(method = "pearson", label.x = 15, label.sep = "\n"),
          ylab = "Viral Reads", xlab = "Restriction Factors"
) + theme_bw() + theme(panel.grid.major = element_blank(),
                       panel.grid.minor = element_blank())

## Deconvolution ##

bgZero <- as.data.frame(phonTools::zeros(nrow(normalized_counts),ncol(normalized_counts))+0.1)
colnames(bgZero) <- colnames(normalized_counts)
rownames(bgZero) <- rownames(normalized_counts)

bgZero = sweep(normalized_counts * 0, 2, as.numeric(sizeBG), "+")

res = spatialdecon(norm = as.matrix(normalized_counts),
                   bg = bgZero,
                   X = placentaSC,
                   align_genes = TRUE)

colMain <- colorRampPalette(rev(brewer.pal(11, "RdBu")))(25)

heatmap(res$beta, cexCol = 0.5, cexRow = 0.7, col=colMain)

heatmapData <- as.data.frame(res$prop_of_all * 100)

cell <- "monocyte"

regionCellDataToPlot <- as.data.frame(cbind(t(heatmapData[cell,]),metadata$condition))
colnames(regionCellDataToPlot) <- c("CellType","Disease")
regionCellDataToPlot$CellType <- as.numeric(regionCellDataToPlot$CellType)

my_comparisons <- list( c("COVID_Placenta", "COVID_Mother"), c("COVID_Placenta", "VUE"), c("COVID_Placenta", "CHI"), c("COVID_Placenta", "Normal") )

ggboxplot(regionCellDataToPlot, x = "Disease", y = "CellType",ncol = 6, palette = "locuszoom", fill = "Disease", color = "black", size = 1.05, scales = "free",
          bxp.errorbar = TRUE, bxp.errorbar.width = 0.3) +
  xlab("") + ylab(paste0(cell," abundance [%]")) + theme_classic() + theme(axis.text.x = element_text(color = "black", size = 11),
                                                                       axis.text.y = element_text(color = "black", size = 11),
                                                                       axis.title.x = element_text(color = "black", size = 12),
                                                                       axis.title.y = element_text(color = "black", size = 13),
                                                                       legend.position = "none",
                                                                       axis.ticks = element_line(size = 1),
                                                                       axis.line = element_line(size = 1),
                                                                       strip.background = element_blank(), strip.text = element_text(color = "black", size = 12)) + rotate_x_text(60) +
  stat_compare_means(comparisons = my_comparisons) #, label = "p.signif")


## Deconvolution - immune ##

bgZero <- as.data.frame(phonTools::zeros(nrow(normalized_counts),ncol(normalized_counts))+0.2)
colnames(bgZero) <- colnames(normalized_counts)
rownames(bgZero) <- rownames(normalized_counts)

resImmune = spatialdecon(norm = as.matrix(normalized_counts),
                   bg = bgZero,
                   X = AA,
                   align_genes = TRUE)

colMain <- colorRampPalette(rev(brewer.pal(11, "RdBu")))(25)

heatmap(resImmune$beta, cexCol = 0.5, cexRow = 0.7, col=colMain)

heatmapData <- as.data.frame(resImmune$prop_of_all * 100)

cell <- "B.naive"

regionCellDataToPlot <- as.data.frame(cbind(t(heatmapData[cell,]),metadata$condition))
colnames(regionCellDataToPlot) <- c("CellType","Disease")
regionCellDataToPlot$CellType <- as.numeric(regionCellDataToPlot$CellType)

my_comparisons <- list( c("COVID_Placenta", "COVID_Mother"), c("COVID_Placenta", "VUE"), c("COVID_Placenta", "CHI"), c("COVID_Placenta", "Normal") )

ggboxplot(regionCellDataToPlot, x = "Disease", y = "CellType",ncol = 6, palette = "locuszoom", fill = "Disease", color = "black", size = 1.05, scales = "free",
          bxp.errorbar = TRUE, bxp.errorbar.width = 0.3) +
  xlab("") + ylab(paste0(cell," abundance [%]")) + theme_classic() + theme(axis.text.x = element_text(color = "black", size = 11),
                                                                           axis.text.y = element_text(color = "black", size = 11),
                                                                           axis.title.x = element_text(color = "black", size = 12),
                                                                           axis.title.y = element_text(color = "black", size = 13),
                                                                           legend.position = "none",
                                                                           axis.ticks = element_line(size = 1),
                                                                           axis.line = element_line(size = 1),
                                                                           strip.background = element_blank(), strip.text = element_text(color = "black", size = 12)) + rotate_x_text(60) +
  stat_compare_means(comparisons = my_comparisons) #, label = "p.signif")

## CIBERSORT

designDoc <- metadata$condition
names(designDoc) <- metadata$...1


mean_by_cluster <- placentaSC
ciberRes <- RNAMagnet::runCIBERSORT(as.matrix(normalized_counts), mean_by_cluster, designDoc, mc.cores = 4)



mean_by_cluster <- AA
ciberResImmune <- RNAMagnet::runCIBERSORT(as.matrix(normalized_counts), mean_by_cluster, designDoc, mc.cores = 4)



ggboxplot(ciberRes, x = "SampleClass", y = "Fraction", facet.by = "CellType", ncol = 6, palette = "locuszoom", fill = "SampleClass", color = "black", size = 1.05, scales = "free",
          bxp.errorbar = TRUE, bxp.errorbar.width = 0.3) +
  xlab("") + ylab("Cell Abundance [%]") + theme_classic() + theme(axis.text.x = element_text(color = "black", size = 11),
                                                                           axis.text.y = element_text(color = "black", size = 11),
                                                                           axis.title.x = element_text(color = "black", size = 12),
                                                                           axis.title.y = element_text(color = "black", size = 13),
                                                                           legend.position = "none",
                                                                           axis.ticks = element_line(size = 1),
                                                                           axis.line = element_line(size = 1),
                                                                           strip.background = element_blank(), strip.text = element_text(color = "black", size = 12)) + rotate_x_text(60) +
  stat_compare_means(comparisons = my_comparisons) #, label = "p.signif")

## New palette testing


NicheDataColors <-
  c(Erythroblasts = "#bc7c7c", Chondrocytes = "#a6c7f7", Osteoblasts = "#0061ff",
    `Fibro/Chondro p.` = "#70a5f9", `pro-B` = "#7b9696", `Arteriolar ECs` = "#b5a800",
    `B cell` = "#000000", `large pre-B.` = "#495959", `Sinusoidal ECs` = "#ffee00",
    Fibroblasts = "#70a5f9", `Endosteal fibro.` = "#264570", `Arteriolar fibro.` = "#567fba",
    `Stromal fibro.` = "#465f82", `small pre-B.` = "#323d3d", `Adipo-CAR` = "#ffb556",
    `Ng2+ MSCs` = "#ab51ff", Neutrophils = "#1f7700", `T cells` = "#915400",
    `NK cells` = "#846232", `Schwann cells` = "#ff00fa", `Osteo-CAR` = "#ff0000",
    `Dendritic cells` = "#44593c", Myofibroblasts = "#dddddd", Monocytes = "#8fff68",
    `Smooth muscle` = "#ff2068", `Ery prog.` = "#f9a7a7", `Mk prog.` = "#f9e0a7",
    `Ery/Mk prog.` = "#f9cda7", `Gran/Mono prog.` = "#e0f9a7", `Neutro prog.` = "#c6f9a7",
    `Mono prog.` = "#f4f9a7", LMPPs = "#a7f9e9", `Eo/Baso prog.` = "#a7b7f9",
    HSPC = "#c6f9a7")

library(scales)
show_col(sample(NicheDataColors, 7))
show_col(NicheDataColors)

## SARS-CoV-2 restriction factors

restrictionFactors <- c("IFITM3","ZBP1","MYD88","IFITM2","STAT2","NRN1","BST2","RETREG1","JADE2","IFIT1","MSR1","FNDC4",
                        "NT5C3A","ISG20","TMEM268","IFIT3","ETV6","TAGAP","TENT5C","TENT5A","SPATS2L","MAX","CCND3","RAB27A",
                        "UPP2","ST3GAL4","CLEC4D","GBP3","B4GALT5","FZD5","ARNTL","TRIM21","NAPA","CRP","APOL2","ERLIN1","CASP7","DDX60")

normalized_counts <- read.table(file = "F:/QuantSeq/COVID-19_Including_Placenta_and_Gary_Reynolds_Liver_April_2021/Placenta/combinedNormalisedReads.txt", sep = "\t")

#colnames(normalized_counts) <- normalized_counts[1,]
#normalized_counts <- normalized_counts[c(2:nrow(normalized_counts)),]
#rownames(normalized_counts) <- normalized_counts[,1]
#normalized_counts <- normalized_counts[,c(2:ncol(normalized_counts))]

#namesOfgenes <- rownames(normalized_counts)

#normalized_counts <- as.data.frame(apply(normalized_counts, MARGIN = 2, as.numeric))

#rownames(normalized_counts) <-namesOfgenes

readsLogged <- log2(normalized_counts[restrictionFactors,]+1)

namesOfgenes2 <- rownames(readsLogged)



readsLogged <- t(apply(readsLogged, MARGIN = 1, scale))

rownames(readsLogged) <- namesOfgenes2

#colnames(virusReadsLogged) = c("P152","P153-1","P153-2","P154","P155-1","P155-2","P156-1","P156-2","P158-1","P158-2","P159-1","P159-2")

#colnames(readsLogged) <- metadata$...1

#virusReadsLogged[virusReadsLogged > 3] <- 3

#split = metadata$condition[metadata$condition %in% c("COVID_Mother","COVID_Placenta")]
split = metadata$condition2


readsLogged[metadata$condition == "Normal"]

ComplexHeatmap::Heatmap(as.matrix(readsLogged), col = brewer.pal(10,"GnBu"), border = FALSE, rect_gp = gpar(col = "grey80", lwd = 2), name = "log2(Normalized Counts + 1)",
                        column_title = "", column_title_gp = gpar(fontsize = 12, fontface = "bold"), clustering_distance_rows = "euclidean", clustering_distance_columns = "euclidean",
                        cluster_columns = FALSE, column_split = split,
                        cluster_rows = TRUE,
                        column_names_gp = grid::gpar(fontsize = 11),
                        row_names_gp = grid::gpar(fontsize = 11),
                        heatmap_legend_param = list(
                          title = "log2(Normalized Counts + 1)", title_gp=gpar(fontsize=11, fontface="bold"),
                          legend_height = unit(6, "cm"), title_position = "leftcenter-rot"
                        )
)


### Plot GO results

GOresults <- read.table("D:/COVID_TISSUE_PROJECT/PlacentaPaper/LegitFigures/Ex4/GO_Biological_Process_2021_table.txt", sep = "\t")
colnames(GOresults) <- GOresults[1,]
GOresults <- GOresults[c(2:nrow(GOresults)),]


GOresults$Log10AdjP <- -log10(as.numeric(GOresults$`Adjusted P-value`))


GOresultsShort <- GOresults[c(1:10),]

ggbarplot(GOresultsShort, x = "Term", y = "Log10AdjP",
          fill = c("#D0595A"),               # change fill color by cyl
          color = "black",            # Set bar border colors to white
          sort.val = "asc",          # Sort the value in dscending order
          sort.by.groups = FALSE,     # Don't sort inside each group
          x.text.angle = 90           # Rotate vertically x axis texts
) + rotate() + theme_bw() + theme(panel.grid.major = element_blank(),
                                  panel.grid.minor = element_blank()) + xlab("")

## CXCL8, CXCL10 and IL8 bar plots


genes <- "CXCL8"

a <- colMeans(na.omit(normalized_counts[genes,]))

a <- cbind(a,metadata$condition2, metadata$...1)
colnames(a) <- c("Expression","Disease","Patient")
a <- as.data.frame(a)
a$Expression <- as.numeric(a$Expression)

#a <- a[a$Disease %in% c("Mother+Placenta+","Mother+Placenta-CHI"),]

#a$Disease <- factor(a$Disease, levels = c("COVID_Placenta", "COVID_Mother", "CHI", "VUE", "Normal"))
a$Disease <- factor(a$Disease, levels = c("COVID_Placenta", "COVID_Mother_CHI", "COVID_Mother_CV", "CHI", "VUE", "Normal"))

#my_comparisons <- list( c("COVID_Placenta", "COVID_Mother"), c("COVID_Placenta", "CHI"), c("COVID_Placenta", "VUE"), c("COVID_Placenta", "Normal") )

#my_comparisons <- list( c("COVID_Placenta", "COVID_Mother_CHI"), c("COVID_Placenta", "CHI"), c("COVID_Placenta", "Normal"))
my_comparisons <- list( c("COVID_Placenta", "CHI"), c("COVID_Placenta", "Normal"))

#my_comparisons <- list( c("COVID_Placenta", "CHI"), c("COVID_Placenta", "Normal") )

ggboxplot(a, x = "Disease", y = "Expression",ncol = 6, palette = "locuszoom", fill = "Disease", color = "black", size = 1.05, scales = "free",
          bxp.errorbar = TRUE, bxp.errorbar.width = 0.3) +
  xlab("")   + ylab("Gene Set Expression") + theme_bw() + theme(axis.text.x = element_text(color = "black", size = 11),
                                                                axis.text.y = element_text(color = "black", size = 11),
                                                                axis.title.x = element_text(color = "black", size = 12),
                                                                axis.title.y = element_text(color = "black", size = 13),
                                                                legend.position = "none",
                                                                axis.ticks = element_line(size = 1),
                                                                axis.line = element_line(size = 1),
                                                                panel.grid.major = element_blank(),
                                                                panel.grid.minor = element_blank(),
                                                                strip.background = element_blank(), strip.text = element_text(color = "black", size = 12)) + rotate_x_text(60) +
  stat_compare_means(comparisons = my_comparisons, label = "p.signif")










