# Load high resolution object
covidPlacentaHighRes <- readRDS("D:/COVID_TISSUE_PROJECT/PlacentaPaper/RDS_Objects/baseHighResClusteringMISSILeObject.rds")

# Add metadata
metadataInput <- readxl::read_xlsx(path = "D:/COVID_TISSUE_PROJECT/PlacentaPaper/metadataPlacenta.xlsx")
metadataInput <- metadataInput[c(1:40),]
covidPlacentaHighRes <- MISSILe::createFactoredMetadataItem(covidPlacentaHighRes, baseItem = "allRegions",
                                                            newItem = metadataInput$PatientID, newItemName = "PatientID")

covidPlacentaHighRes <- MISSILe::createFactoredMetadataItem(covidPlacentaHighRes, baseItem = "allRegions",
                                                            newItem = metadataInput$RegionAreaPixels, newItemName = "RegionArea")

covidPlacentaHighRes <- MISSILe::createFactoredMetadataItem(covidPlacentaHighRes, baseItem = "allRegions",
                                                            newItem = metadataInput$FibrinAreaPixels, newItemName = "FibrinArea")

covidPlacentaHighRes <- addMISSILeMetadata(MISSILeObject = covidPlacentaHighRes, metaData = covidPlacenta@meta.data$LowResClusteringIdentities, metaDataName = "LowResolutionClustering")



covidPlacentaHighRes@meta.data$Disease[covidPlacentaHighRes@meta.data$Disease == "Unknown"] <- "COVID_Placenta"

covidPlacentaHighRes@meta.data$LowResolutionClustering[covidPlacentaHighRes@meta.data$LowResolutionClustering == "Hofbauer cells"] <- "M1 Macrophage"

covidPlacentaHighRes@meta.data$LowResolutionClustering <- droplevels(covidPlacentaHighRes@meta.data$LowResolutionClustering)
covidPlacentaHighRes@meta.data$Disease <- droplevels(covidPlacentaHighRes@meta.data$Disease)

# High resolution clustering

colnames(covidPlacentaHighRes@Expression$MultiIHC@counts)

markers <- c(c(1:16),c(19:28),c(30:31))

covidPlacentaHighRes <- clusterMISSILe(MISSILeObject = covidPlacentaHighRes, markers = markers, includeMetaData = c("size"), numNeighbours = 11, expSet = "counts")

# Heatmap of clustering

levels(covidPlacentaHighRes@meta.data$PhenographClustering)

plotClusterIntensities(MISSILeObject = covidPlacentaHighRes, markers = markers, clusteringName = 'PhenographClustering', heatType = "hello", includeMetaData = c("size"),
                       enrichmentLimits = c(-4,4), expSet = "counts", ignoreClusters = c(0,19,4,5,35,73,74,75,76,77))


clusterMISSILe(MISSILeObject = covidPlacentaHighRes, markers = markers, includeMetaData = c("size"), numNeighbours = 11, expSet = "counts", ignoreClusters = c("Unclassified","Noise"),
               ignoreMetadata = "LowResolutionClustering")

Rphenograph_out <- FastPG::fastCluster(data = as.matrix(dataToCluster), k = 12, num_threads = 14, verbose = TRUE)

oldClustering <- as.character(covidPlacentaHighRes@meta.data$LowResolutionClustering)

oldClustering[oldClustering %in% c("Noise")] <- "76"
oldClustering[oldClustering %in% c("Unclassified")] <- "77"

oldClustering[!(oldClustering %in% c("76","77"))] <- as.character(Rphenograph_out$communities)

oldClustering <- as.numeric(oldClustering)

covidPlacentaHighRes <- addMISSILeMetadata(MISSILeObject = covidPlacentaHighRes, metaData = oldClustering, metaDataName = "PhenographClustering")




# Reclustering protocol

colnames(covidPlacentaHighRes@Expression$MultiIHC@counts)

levels(covidPlacentaHighRes@meta.data$LowResolutionClustering)

functionalMarkersCD4 <- c(1,5,8,10,16,18,21,22,23,24,25,27,28,30,31)
functionalMarkersCD8 <- c(1,5,8,10,18,21,22,24,25,27,30,31)
functionalMarkersEpi <- c(3,8,10,16,21,23,25,30)
functionalMarkersB <- c(1,8,10,16,21,23,25,27,30,31)
functionalMarkersNeu <- c(1,3,8,10,16,21,22,23,25,27,30)
functionalMarkersM1 <- c(8,10,19,21,25,27,30)


functionalMarkersM2 <- c(7,8,10,19,21,25,27,30)

functionalMarkers <- functionalMarkersM2

cellType <- "M2 Macrophage"

cellTypeCounts <- covidPlacentaHighRes@Expression$MultiIHC@counts[covidPlacentaHighRes@meta.data$LowResolutionClustering %in% cellType,]

cellTypeCounts <- cellTypeCounts[,functionalMarkers]

Rphenograph_out <- FastPG::fastCluster(data = as.matrix(cellTypeCounts), k = 20, num_threads = 14, verbose = TRUE)


# Plot cluster intensities

expData <- cellTypeCounts

clusterID <- Rphenograph_out$communities

unique.clusters <- unique(clusterID)

cluster.intensities <- as.data.frame(zeros(length(unique.clusters),ncol(expData)))

for(i in 1:length(unique.clusters)){
  single.clusters <- expData[clusterID == unique.clusters[i],1:ncol(expData)]
  cluster.intensities[i,] <- colMeans(single.clusters)
  rm(single.clusters)
}

colnames(cluster.intensities) <- colnames(expData)
rownames(cluster.intensities) <- unique.clusters

#if(!is.null(ignoreClusters)){
#  cluster.intensities <- cluster.intensities[-which(rownames(cluster.intensities) %in% ignoreClusters),]
#  clusterID <- clusterID[! clusterID %in% ignoreClusters]
#}

length(unique(Rphenograph_out$communities))

enrichmentLimits <- c(-2,2)

enrichment <- as.data.frame(apply(cluster.intensities, MARGIN = 2, scale))

rownames(enrichment) <- rownames(cluster.intensities)

enrichment[enrichment > enrichmentLimits[2]] <- enrichmentLimits[2]

enrichment[enrichment < enrichmentLimits[1]] <- enrichmentLimits[1]

enrichment <- enrichment[ order(as.numeric(row.names(enrichment))), ]

row_ha = HeatmapAnnotation(Cells_per_cluster = anno_barplot(as.numeric(table(clusterID)), bar_width = 0.5, gp = gpar(fill = "black")), annotation_label = "Cells per cluster")

cluster.intensities <- cluster.intensities[ order(as.numeric(row.names(cluster.intensities))), ]

cluster.intensities <- t(round(cluster.intensities))

ComplexHeatmap::Heatmap(t(enrichment), col = colorRampPalette(rev(brewer.pal(5, "RdBu")))(25), rect_gp = gpar(col = "black", lwd = 1),
                        column_order = order(as.numeric(rownames(enrichment))), top_annotation = row_ha,
                        heatmap_legend_param = list(
                          title = "Enrichment", title_gp=gpar(fontsize=11, fontface="bold"),
                          legend_height = unit(6, "cm"), title_position = "leftcenter-rot"
                        ),
                        cell_fun = function(j, i, x, y, w, h, col) { # add text to each grid
                          grid.text(cluster.intensities[i, j], x, y,
                                    gp = gpar(fontsize = 8))
                        }
)

# M1 Macrophage function

CD80clusters <- c(22,27,28) # > 7500
IP10clusters <- c(0,4,20,22,27) # > 6000
GPNMBclusters <- c(0,1,4,28) # > 7000
PD1clusters <- c(1,20) # > 1300
PDL1clusters <- c(0,9,10) # > 3000
Ki67clusters <- c(11) # > 2000

m1Metdata <- covidPlacentaHighRes@meta.data[covidPlacentaHighRes@meta.data$LowResolutionClustering %in% cellType,]

diseases <- levels(m1Metdata$Disease3)

macrophageFunction <- phonTools::zeros(6, length(diseases))

functionalMarkers <- c("CD80clusters","IP10clusters","GPNMBclusters","PD1clusters","PDL1clusters","Ki67clusters")

  for(j in 1:ncol(macrophageFunction)){

  macrophageFunction[2,j] <- nrow(m1Metdata[Rphenograph_out$communities %in% IP10clusters & m1Metdata$Disease3 == diseases[j],]) / nrow(covidPlacentaHighRes@meta.data[covidPlacentaHighRes@meta.data$Disease3 == diseases[j],])

  }

macrophageFunctionscaled <- macrophageFunction

macrophageFunctionscaled <- as.data.frame(macrophageFunctionscaled)

for(i in 1:nrow(macrophageFunction)){

  macrophageFunctionscaled[i,] <- scale(macrophageFunction[i,])

}


rownames(macrophageFunctionscaled) <- c("CD80","IP10/CXCL10","GPNMB","PD1","PDL1","Ki67")
colnames(macrophageFunctionscaled) <- diseases


macrophageFunctionscaled$markers <- rownames(macrophageFunctionscaled)

macrophageFunctionscaled <- macrophageFunctionscaled[match(c("GPNMB","PDL1","PD1","Ki67","CD80","IP10/CXCL10"), macrophageFunctionscaled$markers),]

macrophageFunctionscaled <- macrophageFunctionscaled[,c(1:6)]


library(dplyr)

#macrophageFunctionscaled <- as.data.frame(macrophageFunctionscaled)

#macrophageFunctionscaled %>%
#slice(3,5,4,6,1,2)


#macrophageFunctionscaled <- macrophageFunctionscaled[order(c("GPNMB","PDL1","PD1","Ki67","CD80","CXCL10")), ]

ComplexHeatmap::Heatmap(as.matrix(macrophageFunctionscaled),
                        col = colorRampPalette(c("blue","black","yellow"))(100),
                        rect_gp = gpar(col = "black", lwd = 2), cluster_rows = FALSE, cluster_columns = TRUE,
                        heatmap_legend_param = list(
                          title = "Likelihood", title_gp=gpar(fontsize=11, fontface="bold"),
                          legend_height = unit(6, "cm"), title_position = "leftcenter-rot"
                        ))

saveRDS(macrophageFunction, "D:/COVID_TISSUE_PROJECT/PlacentaPaper/RDS_Objects/macrophageFunctionM2_unscaled.rds")
saveRDS(macrophageFunctionscaled, "D:/COVID_TISSUE_PROJECT/PlacentaPaper/RDS_Objects/macrophageFunctionM2.rds")
saveRDS(macrophageFunctionscaled, "D:/COVID_TISSUE_PROJECT/PlacentaPaper/RDS_Objects/macrophageFunction.rds")

macrophageFunctionscaled <- readRDS("D:/COVID_TISSUE_PROJECT/PlacentaPaper/RDS_Objects/macrophageFunction.rds")

# column_order = order(as.numeric(rownames(enrichment)))

# Renaming clusters in overall character vector

clusterInput <- as.data.frame(readxl::read_xlsx("D:/COVID_TISSUE_PROJECT/PlacentaPaper/reclustering.xlsx"))

reclusterIDs <- as.factor(Rphenograph_out$communities)

levels(reclusterIDs) <- clusterInput[,17]

#epitheliumRename <- highResClusteringRename[highResClusteringRename %in% cellType]
#highResClusteringRename <- as.character(covidPlacentaHighRes@meta.data$LowResolutionClustering)

saveRDS(highResClusteringRename,"D:/COVID_TISSUE_PROJECT/PlacentaPaper/RDS_Objects/highResClusteringAll.rds")

#highResClusteringRenameBase <- as.character(covidPlacentaHighRes@meta.data$LowResolutionClustering)
#highResClusteringRename[highResClusteringRenameBase %in% "B cells"] <- "B cells"

highResClusteringRename[highResClusteringRename %in% "M1 Macrophages"] <- as.character(reclusterIDs)

unique(highResClusteringRename)





## Finding virus cells using higher res clustering

highResClusteringRename <- readRDS("D:/COVID_TISSUE_PROJECT/PlacentaPaper/RDS_Objects/highResClusteringAll.rds")

highResClusteringRename[highResClusteringRename == "SARS-CoV-2+ Neutrophils"] <- "M1 Macrophages"

sarsM1 <- highResClusteringRename[highResClusteringRename == "SARS-CoV-2+ Macrophages"]




highResClusteringRename[covidPlacentaHighRes@meta.data$LowResolutionClustering == "M1 Macrophage"] <- "M1 Macrophages"

highResClusteringRename[sarsM1 == "SARS-CoV-2+ Macrophages"] <- "SARS-CoV-2+ Macrophages"




tempHighRes <- highResClusteringRename[covidPlacentaHighRes@meta.data$LowResolutionClustering == "M1 Macrophage"]

tempHighRes[Rphenograph_out$communities %in% c(23)] <- "SARS-CoV-2+ Macrophages"

highResClusteringRename[covidPlacentaHighRes@meta.data$LowResolutionClustering == "M1 Macrophage"] <- tempHighRes




highResClusteringRename[covidPlacentaHighRes@meta.data$HighResolutionClustering == "SARS-CoV-2+ Macrophages"] <- "M1 Macrophages"





# Clean up epithelium that is sars-cov-2+ but is actually noise

epitheliumDisease <- covidPlacentaHighRes@meta.data$Disease[covidPlacentaHighRes@meta.data$LowResolutionClustering == "Epithelium"]

Rphenograph_out$communities[Rphenograph_out$communities == 2 & epitheliumDisease %in% c("COVID_Mother","Normal")] <- 1

# Save clustering and rename

unique(highResClusteringRename)

highResClusteringRename[highResClusteringRename == "CD39+ IP10+ Epithelium"] <- "CD39+ Epithelium"
highResClusteringRename[highResClusteringRename == "IP10+ PD1+ Epithelium"] <- "PD1+ Epithelium"
highResClusteringRename[highResClusteringRename == "IP10+ Epithelium"] <- "Epithelium"
highResClusteringRename[highResClusteringRename == "PD1+ GrB+ B cells"] <- "PD1+ B cells"

saveRDS(highResClusteringRename, "D:/COVID_TISSUE_PROJECT/PlacentaPaper/RDS_Objects/highResClusteringAll.rds")

# Intermediate clustering

bb <- unique(interResClusteringRename)
phenotypes <- grep('Epithelium', bb, value=TRUE)
phenotypes

interResClusteringRename <- highResClusteringRename

interResClusteringRename[interResClusteringRename == "GPNMB+ PD1+ CD39+ M2 Macrophages"] <- "GPNMB+ PD1+ M2 Macrophages"


covidPlacentaHighRes <- addMISSILeMetadata(MISSILeObject = covidPlacentaHighRes, metaData = interResClusteringRename, metaDataName = "InterResolutionClustering1")




## Abundances of virus cells - Fig. 2

covidPlacentaHighRes <- addMISSILeMetadata(MISSILeObject = covidPlacentaHighRes, metaData = highResClusteringRename, metaDataName = "HighResolutionClustering")

unique(covidPlacentaHighRes@meta.data$HighResolutionClustering)

bb <- unique(covidPlacentaHighRes@meta.data$InterResolutionClustering1)

phenotypes <- grep('SARS-CoV-2', bb, value=TRUE)

phenotypes <- c("SARS-CoV-2+ Macrophages", "SARS-CoV-2+ Epithelium")

my_comparisons <- list( c("COVID_Placenta", "COVID_Mother"), c("COVID_Placenta", "VUE"), c("COVID_Placenta", "CHI"), c("COVID_Placenta", "Normal") )

a <- plotCellAbundanceGrouped(MISSILeObject = covidPlacentaHighRes, Idents = "InterResolutionClustering2", ignoreClusters = "Noise",
                              cluster = phenotypes,
                              condition = "Disease3", normalisation = "totalCells")

a$Disease[a$Disease %in% c("COVID_Mother_CHI","COVID_Mother_CV")] <- "COVID_Mother_CHI_CV"
a$Disease[a$Disease %in% c("CHI","VUE","Normal")] <- "All controls"
a$Disease <- as.factor(a$Disease)
a$Disease <- factor(a$Disease, levels = c("COVID_Placenta","COVID_Mother_CHI_CV","All controls"))


bxp <- ggboxplot(
  a, x = "Disease", y = "Frequency", fill = "Disease", color = "black", palette = "npg", alpha = 0.4,
  scales = "free", facet.by = "Phenotype", ncol = 8, shape = 17, size = 1.2
)

bxp <- ggboxplot(
  a, x = "Disease", y = "Frequency", fill = "Phenotype", color = "black", palette = "npg", alpha = 0.4,
  scales = "free", ncol = 8, shape = 17, size = 1
)
#, add = "jitter"

bxp + xlab("") + ylab("Abundance [%]") + theme_bw() +
  theme(strip.text = element_text(colour = "black", face = "bold", size = 9)) +
  theme(strip.background = element_rect(fill=NA, colour = NA))  + theme(axis.text.x = element_text(color = "black", size = 9),
                                                                                               axis.text.y = element_text(color = "black", size = 9),
                                                                                               axis.title.x = element_text(color = "black", size = 11),
                                                                                               axis.title.y = element_text(color = "black", size = 12),
                                                                                               legend.position = "none",panel.grid.major = element_blank(),
                                                                                               panel.grid.minor = element_blank()) + rotate_x_text(45) +
  stat_compare_means(comparisons = my_comparisons, label = "p.signif", hide.ns = TRUE, vjust = 0.5) +
  geom_jitter(color = "black", width = 0, height = 0)


pdf("D:/COVID_TISSUE_PROJECT/PlacentaPaper/REVIEWER_COMMENTS/Updated_Figures/Fig2/virus_mIHC.pdf", width = 4, height = 5)

ggpubr::ggboxplot(a, x = "Disease", y = "Frequency",
                  fill = "Phenotype") + rotate_x_text(45) + scale_fill_manual(values = c("#706EC4","black","black","black","black")) + 
  xlab("") + ylab("Abundance [%]") + theme(legend.position = "none") +
  geom_jitter(
    aes(x = Disease, y = Frequency, fill = Phenotype),
    position = position_jitterdodge(jitter.width = 0, dodge.width = 0.75),
    color = "black",
    size = 2
  ) + ggdist::theme_ggdist() + theme(legend.position = "none",
                                                                                                 axis.text.x = element_text(size = 18, colour = "black"),
                                                                                                 axis.text.y = element_text(size = 18, colour = "black"),
                                                                                                 axis.title.y = element_text(size = 20, colour = "black"),
                                                                                                 axis.ticks.length=unit(.25, "cm"),
                                                                                                 axis.line.x = element_line(color = "black", linewidth = rel(0.5)),
                                                                                                 axis.line.y = element_line(color = "black", linewidth = rel(0.5)),
                                                                                                 axis.ticks = element_line(color="black"))

dev.off()

a_stats <- a[a$Phenotype == "SARS-CoV-2+ Macrophages",] 
a_stats <- a_stats[a_stats$Disease %in% c("COVID_Placenta","COVID_Mother_CHI_CV"),]

wilcox.test(Frequency ~ Disease,data = a_stats)

## ACE2

bb <- unique(covidPlacentaHighRes@meta.data$InterResolutionClustering1)

phenotypes <- grep('Epithelium', bb, value=TRUE)


ace2SARS <- as.data.frame(cbind(as.character(covidPlacentaHighRes@meta.data$InterResolutionClustering1[covidPlacentaHighRes@meta.data$InterResolutionClustering1 %in% phenotypes]),covidPlacentaHighRes@Expression$MultiIHC@counts[covidPlacentaHighRes@meta.data$InterResolutionClustering1 %in% phenotypes,"ACE2"]))
ace2Epi <- as.data.frame(cbind(as.character(covidPlacentaHighRes@meta.data$Disease[covidPlacentaHighRes@meta.data$InterResolutionClustering1 == "Epithelium"]),covidPlacentaHighRes@Expression$MultiIHC@counts[covidPlacentaHighRes@meta.data$InterResolutionClustering1 == "Epithelium","ACE2"]))
colnames(ace2Epi) <- c("Disease","ACE2")
ace2Epi$ACE2 <- as.numeric(ace2Epi$ACE2)

colnames(ace2SARS) <- c("Disease","ACE2")
ace2SARS$ACE2 <- as.numeric(ace2SARS$ACE2)

my_sum <- ace2Epi %>%
  group_by(Disease) %>%
  summarise(
    n=n(),
    mean=mean(ACE2),
    sd=sd(ACE2)
  ) %>%
  mutate( se=sd/sqrt(n))  %>%
  mutate( ic=se * qt((1-0.05)/2 + .5, n-1))


my_sum <- ace2SARS %>%
  group_by(Disease) %>%
  summarise(
    n=n(),
    mean=mean(ACE2),
    sd=sd(ACE2)
  ) %>%
  mutate( se=sd/sqrt(n))  %>%
  mutate( ic=se * qt((1-0.05)/2 + .5, n-1))

my_comparisons <- list( c("SARS-CoV-2+ Epithelium", "Ki67+ Epithelium"), c("SARS-CoV-2+ Epithelium", "Epithelium") )

ggplot(my_sum, aes(x=Disease, y=mean, fill=Disease)) +
  geom_bar( stat="identity", color="black", position=position_dodge(), size = 1.2)  +
  geom_errorbar(aes(ymin=mean-sd, ymax=mean+sd), width=.2, size = 1.2,
                position=position_dodge(.9)) + theme_bw() + ylab("ACE2 Expression") + xlab("")


## Fibrin

phenotypes <- c("Epithelium")

my_comparisons <- list( c("COVID_Placenta", "COVID_Mother"), c("COVID_Placenta", "VUE"), c("COVID_Placenta", "CHI"), c("COVID_Placenta", "Normal") )

a <- plotCellAbundanceGrouped(MISSILeObject = covidPlacentaHighRes, Idents = "LowResolutionClustering", ignoreClusters = c("Noise"), cluster = phenotypes, condition = "Disease", normalisation = "area")

fibrinPercent <- (metadataInput$FibrinAreaPixels / metadataInput$RegionAreaPixels)*100

metadataInput$Disease[metadataInput$Disease == "Unknown"] <- "COVID_Placenta"

fibrinDisease <- as.data.frame(cbind(metadataInput$Disease2,fibrinPercent,metadataInput$PatientID))
colnames(fibrinDisease) <- c("Disease","FibrinPercent","PatientID")
fibrinDisease$FibrinPercent <- as.numeric(fibrinDisease$FibrinPercent)

fibrinDisease <- fibrinDisease[fibrinDisease$Disease %in% c("COVID_Placenta","COVID_Mother_CV","COVID_Mother_CHI"),]

my_comparisons <- list( c("COVID_Placenta", "COVID_Mother_CHI"),c("COVID_Placenta", "COVID_Mother_CV"))

fibrinDisease$Disease <- factor(fibrinDisease$Disease, levels = c("COVID_Placenta","COVID_Mother_CHI","COVID_Mother_CV"))

ggboxplot(
  fibrinDisease, x = "Disease", y = "FibrinPercent", fill = "Disease", color = "black", palette = "locuszoom",
  shape = 17, size = 1
) +  scale_x_discrete(labels=c("CHI" = "CHI",
                               "COVID_Mother" = "CM", "COVID_Placenta" = "CP",
                               "Normal" = "Normal",
                               "VUE" = "VUE"))  + xlab("") + ylab("Fibrin Percent") + theme_bw() +
  theme(strip.text = element_text(colour = "black", face = "bold", size = 8)) +
  theme(strip.background = element_rect(fill=NA, colour = "black", linetype="solid"))  + theme(axis.text.x = element_text(color = "black", size = 9),
                                                                                               axis.text.y = element_text(color = "black", size = 9),
                                                                                               axis.title.x = element_text(color = "black", size = 11),
                                                                                               axis.title.y = element_text(color = "black", size = 12),
                                                                                               legend.position = "none",panel.grid.major = element_blank(),
                                                                                               panel.grid.minor = element_blank()) + rotate_x_text(45) +
  stat_compare_means(comparisons = my_comparisons)




fibrinDisease <- cbind(fibrinDisease,metadataInput$PatientID)


fibrinDisease$Epithelium <- phonTools::zeros(nrow(fibrinDisease))

fibrinDisease$Epithelium[fibrinDisease$Disease == "VUE"] <- a$Frequency[a$Disease == "VUE"]

ggscatter(fibrinDisease, x = "Epithelium", y = "FibrinPercent", color = "Disease",
          add = "reg.line") + ylim(0,50)

# Distances

bb <- unique(covidPlacentaHighRes@meta.data$InterResolutionClustering1)
phenotypes1 <- grep('M1 Macrophages', bb, value=TRUE)
phenotypes2 <- grep('Epithelium', bb, value=TRUE)

covidPlacentaHighRes@current.identity$ActiveIdents <- "InterResolutionClustering1"

allDistancesHigh <- minimalDistances(MISSILeObject = covidPlacentaHighRes, phenotypes = NULL, secondPhenotypes = c(phenotypes1,phenotypes2))

# Plot across diseases
id <- 111
specificDistances <- allDistancesHigh[[id]]
cellTypeIDs <- names(allDistancesHigh)[id]
specificDistances <- as.data.frame(specificDistances)
colnames(specificDistances) <- c("Meannn","sd")
specificDistances$Disease <- metadataInput$Disease
specificDistances$regions <- metadataInput$allRegions
specificDistances$Disease[specificDistances$Disease %in% c("VUE","CHI")] <- "Disease Controls"

my_comparisons <- list( c("COVID_Placenta", "COVID_Mother"), c("COVID_Placenta", "Disease Controls"), c("COVID_Placenta", "Normal") )

ggboxplot(specificDistances, x = "Disease", y = "Meannn", color =  "Disease", palette = "npg", size = 1.2) + stat_compare_means(comparisons = my_comparisons,label = "p.format") +
  theme_bw() + xlab("") + ylab(paste0("Distance between ", cellTypeIDs))

# Plot across cell interactions

names(allDistancesHigh)
ids <- c(7:11)
diseaseState <- "COVID_Placenta"
specificDistances <- allDistancesHigh[ids]
for(i in 1:length(specificDistances)){

  if(i == 1){
    cellTypeDistances <- as.data.frame(specificDistances[[i]][metadataInput$Disease == diseaseState,])
    cellTypeDistances$cellType <- rep(names(specificDistances)[i], nrow(cellTypeDistances))
  } else {
    cellTypeDistancesTemp <- as.data.frame(specificDistances[[i]][metadataInput$Disease == diseaseState,])
    cellTypeDistancesTemp$cellType <- rep(names(specificDistances)[i], nrow(cellTypeDistancesTemp))
    cellTypeDistances <- rbind(cellTypeDistances,cellTypeDistancesTemp)
    rm(cellTypeDistancesTemp)
  }

}

colnames(cellTypeDistances) <- c("Mean","sd","Celltype")

ggboxplot(cellTypeDistances, x = "Celltype", y = "Mean", color =  "Celltype", size = 1.2) +
  theme_bw() + xlab("") + ylab(paste0(diseaseState)) + rotate_x_text(45) + stat_compare_means(label = "p.format") + theme(legend.position = "none") +
  geom_hline(yintercept = mean(na.omit(cellTypeDistances$Mean)), linetype = 2)  # Add horizontal line at base mean
#  stat_compare_means(method = "anova", label.y = 15000)+        # Add global annova p-value
#  stat_compare_means(label = "p.signif", method = "t.test",
#                     ref.group = ".all.")


## Point pattern

covidPlacentaHighRes@current.identity$ActiveIdents <- "LowResolutionClustering"
covidPlacentaHighRes@current.identity$ActiveIdents <- "HighResolutionClustering"
covidPlacentaHighRes@current.identity$ActiveIdents <- "InterResolutionClustering1"

plotSpatialCoordinates(MISSILeObject = covidPlacentaHighRes, region = 39, palette = "Set1", pointSize = 2, crop = c(0.6,0.6,0.8,0.8))

## Interactions

covidPlacentaHighRes <- createDelaunayMISSILe(MISSILeObject = covidPlacentaHighRes, maxNeighbourDistance = 100, ignoreClusters = c("Unclassified","Noise"))

saveRDS(covidPlacentaHighRes,"D:/COVID_TISSUE_PROJECT/PlacentaPaper/RDS_Objects/highResClusteringMISSILeObject_withDelaunay.rds")

covidPlacentaHighRes <- prepInteractions(MISSILeObject = covidPlacentaHighRes, ignoreClusters = c("Unclassified","Noise"))

saveRDS(covidPlacentaHighRes,"D:/COVID_TISSUE_PROJECT/PlacentaPaper/RDS_Objects/highResClusteringMISSILeObject_withPrepInteractions.rds")

covidPlacentaHighRes <- combineInteractionData(MISSILeObject = covidPlacentaHighRes, ignoreClusters = c("Unclassified","Noise"), condition = "Disease3")

likRat <- likelihoodRatio(MISSILeObject = covidPlacentaHighRes, combined = FALSE)

relFre <- relativeFrequencies(MISSILeObject = covidPlacentaHighRes, combined = FALSE)

drops <- c("Unclassified","Noise")
id <- "COVID_Placenta"

whichInteraction <- likRat

ComplexHeatmap::Heatmap(whichInteraction[[id]][!(rownames(whichInteraction[[id]]) %in% drops) , !(colnames(whichInteraction[[id]]) %in% drops)], rect_gp = gpar(col = "grey", lwd = 2),
                        col = colorRamp2(c(0, 0.2, 0.4), c("deepskyblue3", "white", "red2")), cluster_rows = TRUE, cluster_columns = TRUE,
                        heatmap_legend_param = list(
                          title = "Likelihood", title_gp=gpar(fontsize=11, fontface="bold"),
                          legend_height = unit(6, "cm"), title_position = "leftcenter-rot"
                        ))

# Plot likelihood per disease

whichInteraction <- relFre

location <- c(2,11)
likelihoodTargeted <- zeros(length(whichInteraction))
for(i in 1:length(whichInteraction)){
  likelihoodTargeted[i] <- whichInteraction[[i]][location[1],location[2]]
}
interactionCells <- paste0(rownames(whichInteraction[[i]])[location[1]]," - ", colnames(whichInteraction[[i]])[location[2]])

#toPlot <- as.data.frame(cbind(likelihoodTargeted,names(whichInteraction)))
toPlot <- as.data.frame(cbind(likelihoodTargeted,metadataInput$Disease2))

colnames(toPlot) <- c("Likelihood","Condition")

toPlot$Likelihood <- as.numeric(toPlot$Likelihood)

toPlot$Condition <- factor(toPlot$Condition, levels = c("COVID_Placenta","COVID_Mother_CHI","COVID_Mother_CV","CHI","VUE", "Normal"))

my_comparisons <- list( c("COVID_Placenta", "COVID_Mother_CHI"), c("COVID_Placenta", "COVID_Mother_CV"))
my_comparisons <- list( c("COVID_Placenta", "COVID_Mother_CHI"), c("COVID_Placenta", "COVID_Mother_CV"), c("COVID_Placenta", "CHI"))

ggboxplot(toPlot, x = "Condition", y = "Likelihood", size = 1.2, fill = "Condition", color = "black", palette = "locuszoom" ) + xlab("") +
  ylab(paste0("Relative Frequency: ", interactionCells)) + rremove("legend") +
  stat_compare_means(comparisons = my_comparisons,label = "p.signif") + theme_bw() +
  theme(
    axis.ticks = element_line(size = 1),
    axis.line = element_line(size = 1),
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    legend.position = "none"
  ) + rotate_x_text(45)

saveRDS(toPlot, "D:/COVID_TISSUE_PROJECT/PlacentaPaper/REVIEWER_COMMENTS/Updated_Figures/Fig4/Code_for_plotting/fig4d_2.rds")

pdf("D:/COVID_TISSUE_PROJECT/PlacentaPaper/REVIEWER_COMMENTS/Updated_Figures/Fig4/interaction_cd8.pdf", width = 4, height = 6)

ggboxplot(toPlot, x = "Condition", y = "Likelihood",ncol = 6, palette = c("#D43F3A","#EEA236","#46B8DA","#F8DAAF","#AEE0EF","#CCCCCC"), fill = "Condition", color = "black", size = 1.05,
          bxp.errorbar = TRUE, bxp.errorbar.width = 0.3) +
  xlab("")   + ylab(paste0("Relative Frequency: ", interactionCells)) + ggdist::theme_ggdist() +
  geom_jitter(
    width = 0,
    color = "black",
    size = 2
  ) + theme(panel.grid.major = element_blank(),
            panel.grid.minor = element_blank(),
            legend.position = "none",
            axis.text.x = element_text(size = 4, colour = "black"),
            axis.text.y = element_text(size = 18, colour = "black"),
            axis.title.y = element_text(size = 20, colour = "black"),
            axis.title.x = element_text(size = 20, colour = "black"),
            axis.ticks.length=unit(.25, "cm"),
            axis.line.x = element_line(color = "black", linewidth = rel(0.5)),
            axis.line.y = element_line(color = "black", linewidth = rel(0.5)),
            axis.ticks = element_line(color="black"))  + rotate_x_text(60) +
  stat_compare_means(comparisons = my_comparisons, tip.length = 0) 

dev.off()







### UMAP visualisation

a <- downSampleVisualisation(clusters = covidPlacentaHighRes@meta.data$LowResolutionClustering, downsamplePercent = 0.01)

colnames(covidPlacentaHighRes@Expression$MultiIHC@counts)

covidPlacentaHighRes@current.identity$ActiveIdents <- "LowResolutionClustering"

orderMarkers = list(c(29),c(7),c(20),c(2),c(4),c(20,4),c(15),c(17),c(12),c(11,13),c(9,13),c(14))

orderMarkers = rev(orderMarkers)

levels(covidPlacentaHighRes@meta.data$LowResolutionClustering)

a <- downSampleVisualisation(clusters = covidPlacentaHighRes@meta.data$LowResolutionClustering, downsamplePercent = 0.01, countsMatrix = covidPlacentaHighRes@Expression$MultiIHC@counts, markers = orderMarkers)

a <- a[!(a$membership %in% c("Noise","Unclassified")),]

markers <- c(2,4,6,7,9,11,12,13,14,15,20,26,28)

reducedExp <- covidPlacentaHighRes@Expression$MultiIHC@counts[as.numeric(a$index),markers]

umap_out <- uwot::umap(reducedExp, ret_model = TRUE, verbose = TRUE,
                       n_threads = 14, n_neighbors = 15,
                       min_dist = 0.01, spread = 1, n_trees = 50,
                       n_epochs = NULL)


umap_plot <- data.frame(x = umap_out$embedding[,1],
                        y = umap_out$embedding[,2],
                        Clusters = a$membership)
umap_plot$Clusters <- as.factor(umap_plot$Clusters)
#umap_plot$Regions <- as.factor(umap_plot$Regions)


ggplot(umap_plot) + geom_point(aes(x=x, y=y, color=Clusters), alpha = 0.7) +
  guides(colour = guide_legend(override.aes = list(size=5))) +
  theme_bw() + labs(y= "UMAP 2", x = "UMAP 1") +
  theme_classic() + scale_color_manual(values = c("#ff566e",
                                                  "#a74621",
                                                  "#b66700",
                                                  "#deb200",
                                                  "#84aa00",
                                                  "#0ab807",
                                                  "#019258",
                                                  "#0072e4",
                                                  "#884c96",
                                                  "#ff86dd",
                                                  "#d3008d"))



a$membership[a$membership == "Fibroblasts"] <- "Stromal cells"


rtsne_out <- Rtsne::Rtsne(reducedExp, verbose = TRUE, perplexity = 30, num_threads = 14)

rtsne_plot <- data.frame(x = rtsne_out$Y[,1],
                         y = rtsne_out$Y[,2],
                         Clusters = a$membership)
rtsne_plot$Clusters <- as.factor(rtsne_plot$Clusters)

colorIn <- rtsne_plot$Clusters

centroids.data <- zeros(length(levels(colorIn)),2)
colnames(centroids.data) <- c("centroids_x","centroids_y")
levelsTo <- levels(colorIn)
for(i in 1:length(unique(colorIn))){
  x_coords <- rtsne_plot$x[colorIn == levelsTo[i]]
  y_coords <- rtsne_plot$y[colorIn == levelsTo[i]]
  centroids.data[i,1] <- mean(x_coords)
  centroids.data[i,2] <- mean(y_coords)
  rm(x_coords, y_coords)
}

qualpalette <- qualpalr::qualpal(10, "pretty")   # pretty, pretty_dark, rainbow, pastels
qualpalette <- qualpalr::qualpal(5, "pretty")

ggplot() + geom_point(data = rtsne_plot, aes(x=x, y=y, color=Clusters), alpha = 0.7) +
  guides(colour = guide_legend(override.aes = list(size=4))) +
  theme_bw() + labs(y= "tSNE 2", x = "tSNE 1") +
  theme_classic() + scale_color_manual(values = qualpalette$hex) +
  geom_text(data = as.data.frame(centroids.data), mapping = aes(x = centroids_x,
                                 y = centroids_y,
                                 label = 1:length(unique(colorIn))),
            color = "black", size = 6, fontface = 2)


# + scale_color_manual(values = brewer.pal(length(levels(rtsne_plot$Clusters)), "Paired"))

# scale_color_manual(values = brewer.pal(length(levels(umap_plot$Clusters)), "Set3"))

#+ theme(axis.line = element_line(colour = 'black', size = 2))

### Virus microenvironment ###

colnames(covidPlacentaHighRes@Expression$MultiIHC@counts)

funcMarkers <- c(1,3,5,6,8,10,19,21,22,25,27,30,31)

covidPlacentaHighRes@current.identity$ActiveIdents <- "InterResolutionClustering1"
covidPlacentaHighRes@current.identity$ActiveIdents <- "VMEclusters"
covidPlacentaHighRes@current.identity$ActiveIdents <- "VMEclusters"



unique(covidPlacentaHighRes@meta.data$InterResolutionClustering1)

covidPlacentaHighRes <- calculateUniqueNeighbourhood(MISSILeObject = covidPlacentaHighRes, functionalMarkers = funcMarkers, cellOfInterest = "SARS-CoV-2+ Epithelium")

covidPlacentaHighRes <- calculateUniqueNeighbourhood(MISSILeObject = covidPlacentaHighRes, functionalMarkers = funcMarkers, cellOfInterest = "Non-infected Epithelium")

covidPlacentaHighRes <- calculateUniqueNeighbourhood(MISSILeObject = covidPlacentaHighRes, functionalMarkers = funcMarkers, cellOfInterest = "Control Epithelium")

covidPlacentaHighRes <- calculateUniqueNeighbourhood(MISSILeObject = covidPlacentaHighRes, functionalMarkers = funcMarkers, cellOfInterest = "CHI Control Epithelium")


phenotypes <- grep('B cells', unique(covidPlacentaHighRes@meta.data$InterResolutionClustering1), value=TRUE)



plotContinuousUniqueNeighbourhood(MISSILeObject = covidPlacentaHighRes, phenotype = phenotypes,
                                  neighbourhoodName = "SARS-CoV-2+ Epithelium", colours = gg_color_hue(length(phenotypes)))


plotContinuousUniqueNeighbourhood(MISSILeObject = covidPlacentaHighRes, phenotype = phenotypes,
                                  neighbourhoodName = "Non-infected Epithelium", colours = gg_color_hue(length(phenotypes)))

# this one
  plotContinuousUniqueNeighbourhood(MISSILeObject = covidPlacentaHighRes, phenotype = "PDL1+ M2 Macrophages",
                                  neighbourhoodName = "SARS-CoV-2+ Epithelium", cellComparisons = c("SARS-CoV-2+ Epithelium","Non-infected Epithelium"),
                                  colours = c("#00e1ff","#ff8c00"))

plotContinuousUniqueNeighbourhood(MISSILeObject = covidPlacentaHighRes, phenotype = phenotypes[1],
                                  neighbourhoodName = "CHI Control Epithelium", cellComparisons = c("CHI Control Epithelium","Non-infected Epithelium"),
                                  colours = c("#99EDC3","#A45EE5"))




otherDisease <- c("COVID_Mother","VUE","CHI","Normal")

phenotypes <- grep('Epithelium', unique(covidPlacentaHighRes@meta.data$InterResolutionClustering1), value=TRUE)[c(1:5)]

# Rename ACE2+ epithelium to ACE2+ Epithelium
covidPlacentaHighRes@meta.data$InterResolutionClustering1[covidPlacentaHighRes@meta.data$InterResolutionClustering1 == "ACE2+ epithelium"] <- "ACE2+ Epithelium"
covidPlacentaHighRes@meta.data$HighResolutionClustering[covidPlacentaHighRes@meta.data$HighResolutionClustering == "ACE2+ epithelium"] <- "ACE2+ Epithelium"

# Create VME specific clustering slot in the metadata
covidPlacentaHighRes@meta.data$VMEclusters <- covidPlacentaHighRes@meta.data$InterResolutionClustering1

# How many SARS-CoV-2+ epithelium cells are there?
numVirusCells <- length(covidPlacentaHighRes@meta.data$VMEclusters[covidPlacentaHighRes@meta.data$VMEclusters == "SARS-CoV-2+ Epithelium"])

# Randomly sample control epithelium from the other diseases - and rename to "Control Epithelium"
numControlCell <- length(covidPlacentaHighRes@meta.data$VMEclusters[covidPlacentaHighRes@meta.data$VMEclusters %in% phenotypes & covidPlacentaHighRes@meta.data$Disease %in% "CHI"])
covidPlacentaHighRes@meta.data$VMEclusters[covidPlacentaHighRes@meta.data$VMEclusters %in% phenotypes & covidPlacentaHighRes@meta.data$Disease == "CHI"][sample(1:numControlCell, numVirusCells, replace=T)] <- "CHI Control Epithelium"

# Randomly sample control epithelium from the other diseases - and rename to "Non-infected Epithelium"
numNICell <- length(covidPlacentaHighRes@meta.data$VMEclusters[covidPlacentaHighRes@meta.data$VMEclusters %in% phenotypes & covidPlacentaHighRes@meta.data$Disease %in% c("COVID_Placenta")])
covidPlacentaHighRes@meta.data$VMEclusters[covidPlacentaHighRes@meta.data$VMEclusters %in% phenotypes & covidPlacentaHighRes@meta.data$Disease %in% c("COVID_Placenta")][sample(1:numNICell, numVirusCells, replace=T)] <- "Non-infected Epithelium"



covidPlacentaHighRes@neighbourhoods$uniqueNeighbourhood$`CHI Control Epithelium` <- covidPlacentaHighRes@neighbourhoods$uniqueNeighbourhood$`CHI Control Epithelium`$`CHI Control Epithelium`




saveRDS(covidPlacentaHighRes,"D:/COVID_TISSUE_PROJECT/PlacentaPaper/RDS_Objects/highResClusteringMISSILeObject_afterAllvirusME_05_10_2021.rds")
#covidPlacentaHighRes@neighbourhoods$uniqueNeighbourhood <- vector(mode = "list", length = 0)


#### Plot functional marker expression of virus microenvironment

colnames(covidPlacentaHighRes@Expression$MultiIHC@counts)[funcMarkers][9]

plotUniqueNeighbourhoodExpression(MISSILeObject = covidPlacentaHighRes, marker = colnames(covidPlacentaHighRes@Expression$MultiIHC@counts)[funcMarkers][6],neighbourhoodName = "SARS-CoV-2+ Epithelium",
                                  cellComparisons = c("SARS-CoV-2+ Epithelium","Non-infected Epithelium"), colours = c("#00e1ff","#ff8c00"))

#### Spatial likelihood ratio

covidPlacentaHighRes@current.identity$ActiveIdents <- "InterResolutionClustering1"
a <- localSpatialRelationships(MISSILeObject = covidPlacentaHighRes, cellType = cellType, microenvironmentSize = 200, maxDistance = 100)

# Plotting spatial likelihood ratio


#### Neighbourhoods

covidPlacentaHighRes@current.identity$ActiveIdents <- "ClusteringForNeighbourhoods"
covidPlacentaHighRes@current.identity$ActiveIdents <- "InterResolutionClustering1"


covidPlacentaHighRes@meta.data$LowResolutionClustering <- as.character(covidPlacentaHighRes@meta.data$LowResolutionClustering)
covidPlacentaHighRes@meta.data$LowResolutionClustering[is.na(covidPlacentaHighRes@meta.data$LowResolutionClustering)] = "Epithelium"
covidPlacentaHighRes@meta.data$LowResolutionClustering <- as.factor(covidPlacentaHighRes@meta.data$LowResolutionClustering)


covidPlacentaHighRes@meta.data$ClusteringForNeighbourhoods <- as.character(covidPlacentaHighRes@meta.data$LowResolutionClustering)
covidPlacentaHighRes@meta.data$ClusteringForNeighbourhoods[covidPlacentaHighRes@meta.data$InterResolutionClustering1 == "SARS-CoV-2+ Epithelium"] = "SARS-CoV-2+ Epithelium"
covidPlacentaHighRes@meta.data$ClusteringForNeighbourhoods <- as.factor(covidPlacentaHighRes@meta.data$ClusteringForNeighbourhoods)

covidPlacentaHighRes@meta.data$InterResolutionClustering1 <- as.factor(covidPlacentaHighRes@meta.data$InterResolutionClustering1)


covidPlacentaHighRes <- cellNeighbourhoods(MISSILeObject = covidPlacentaHighRes, numOfCells = 10, kMeans = 10, verbose = TRUE, ignoreClusters = c("Unclassified","Noise"),
                                           neighbourhoodName = "GeneralNeighbourhood1")

neighbourhoodEnrichment(MISSILeObject = covidPlacentaHighRes, ignoreClusters = c("Unclassified","Noise"))

a <- neighbourhoodFrequencies(MISSILeObject = covidPlacentaHighRes, ignoreClusters = c("Unclassified","Noise"), neighbourhoodNames = NULL)

ggboxplot(
  a, x = "Condition", y = "Frequency", fill = "Neighbourhood", color = "Neighbourhood",
  scales = "free", add = "jitter", facet.by = "Neighbourhood", ncol = 8, shape = 17, size = 1.1, alpha = 0.2
) + rotate_x_text(45) +
  stat_compare_means(comparisons = my_comparisons,label = "p.format") + xlab("")

covidPlacentaHighRes@current.identity$ActiveIdents <- "InterResolutionClustering1"
covidPlacentaHighRes@current.identity$ActiveIdents <- "LowResolutionClustering"

a <- neighbourhoodCellAbundance(MISSILeObject = covidPlacentaHighRes, cellType = c("CD4+ T cells"), ignoreClusters = c("Unclassified","Noise"))

#a <- as.data.frame(apply(a, MARGIN = 2, scale))

a$Disease <- metadataInput$Disease

library(tidyr)
longA <- a %>% gather(Neighbourhood, Enrichment, -c(Disease))

my_comparisons <- list( c("COVID_Placenta", "COVID_Mother"), c("COVID_Placenta", "VUE"), c("COVID_Placenta", "CHI"), c("COVID_Placenta", "Normal") )

ggboxplot(
  longA, x = "Disease", y = "Enrichment", fill = "Neighbourhood", color = "Neighbourhood",
  scales = "free", add = "jitter", facet.by = "Neighbourhood", ncol = 8, shape = 17, size = 1.1, alpha = 0.2
) + rotate_x_text(45) +
  stat_compare_means(comparisons = my_comparisons,label = "p.format") + xlab("")


b <- neighbourhoodMarkerExpression(MISSILeObject = covidPlacentaHighRes, marker = "PD1")
b <- as.data.frame(apply(b, MARGIN = 2, scale))

b$Disease <- metadataInput$Disease

longB <- b %>% gather(Marker, Enrichment, -c(Disease))

ggboxplot(
  longB, x = "Disease", y = "Enrichment", fill = "Marker", color = "Marker", palette = "nejm",
  scales = "free", add = "jitter", facet.by = "Marker", ncol = 8, shape = 17, size = 1.1, alpha = 0.2
) + rotate_x_text(45) +
  stat_compare_means(comparisons = my_comparisons,label = "p.format") + xlab("")

######

covidPlacentaHighRes <- readRDS("D:/COVID_TISSUE_PROJECT/PlacentaPaper/RDS_Objects/highResClusteringMISSILeObject_withNewVirus.rds")

covidPlacentaHighRes@current.identity$ActiveIdents <- "PhenographClustering"

plotClusterIntensities(MISSILeObject = covidPlacentaHighRes, markers = lineageMarkers, clusteringName = 'PhenographClustering', heatType = "hello", includeMetaData = c("size"),
                       enrichmentLimits = c(-2,2), expSet = "counts")


colnames(covidPlacentaHighRes@Expression$MultiIHC@counts)

lineageMarkers <- c(2,4,6,7,9,11,12,13,14,15,20,26,28)

colnames(covidPlacentaHighRes@Expression$MultiIHC@counts)[lineageMarkers]

###### Tidy up clusters - include some unclassified clusters after bubble plot ######

# Cluster 0 - M2 Macrophage
covidPlacentaHighRes@meta.data$LowResolutionClustering[covidPlacentaHighRes@meta.data$PhenographClustering == 0] <- "M2 Macrophage"
covidPlacentaHighRes@meta.data$InterResolutionClustering1[covidPlacentaHighRes@meta.data$PhenographClustering == 0] <- "M2 Macrophages"
covidPlacentaHighRes@meta.data$HighResolutionClustering[covidPlacentaHighRes@meta.data$PhenographClustering == 0] <- "M2 Macrophage"
covidPlacentaHighRes@meta.data$VMEclusters[covidPlacentaHighRes@meta.data$PhenographClustering == 0] <- "M2 Macrophage"
covidPlacentaHighRes@meta.data$ClusteringForNeighbourhoods[covidPlacentaHighRes@meta.data$PhenographClustering == 0] <- "M2 Macrophage"

# Cluster 2 - Smooth muscle
covidPlacentaHighRes@meta.data$LowResolutionClustering[covidPlacentaHighRes@meta.data$PhenographClustering == 2] <- "Smooth muscle"
covidPlacentaHighRes@meta.data$InterResolutionClustering1[covidPlacentaHighRes@meta.data$PhenographClustering == 2] <- "Smooth muscle"
covidPlacentaHighRes@meta.data$HighResolutionClustering[covidPlacentaHighRes@meta.data$PhenographClustering == 2] <- "Smooth muscle"
covidPlacentaHighRes@meta.data$VMEclusters[covidPlacentaHighRes@meta.data$PhenographClustering == 2] <- "Smooth muscle"
covidPlacentaHighRes@meta.data$ClusteringForNeighbourhoods[covidPlacentaHighRes@meta.data$PhenographClustering == 2] <- "Smooth muscle"

# Cluster 10 - Epithelium
covidPlacentaHighRes@meta.data$LowResolutionClustering[covidPlacentaHighRes@meta.data$PhenographClustering == 10] <- "Epithelium"
covidPlacentaHighRes@meta.data$InterResolutionClustering1[covidPlacentaHighRes@meta.data$PhenographClustering == 10] <- "Epithelium"
covidPlacentaHighRes@meta.data$HighResolutionClustering[covidPlacentaHighRes@meta.data$PhenographClustering == 10] <- "Epithelium"
covidPlacentaHighRes@meta.data$VMEclusters[covidPlacentaHighRes@meta.data$PhenographClustering == 10] <- "Epithelium"
covidPlacentaHighRes@meta.data$ClusteringForNeighbourhoods[covidPlacentaHighRes@meta.data$PhenographClustering == 10] <- "Epithelium"

# Cluster 19 - Epithelium
covidPlacentaHighRes@meta.data$LowResolutionClustering[covidPlacentaHighRes@meta.data$PhenographClustering == 19] <- "Epithelium"
covidPlacentaHighRes@meta.data$InterResolutionClustering1[covidPlacentaHighRes@meta.data$PhenographClustering == 19] <- "Epithelium"
covidPlacentaHighRes@meta.data$HighResolutionClustering[covidPlacentaHighRes@meta.data$PhenographClustering == 19] <- "Epithelium"
covidPlacentaHighRes@meta.data$VMEclusters[covidPlacentaHighRes@meta.data$PhenographClustering == 19] <- "Epithelium"
covidPlacentaHighRes@meta.data$ClusteringForNeighbourhoods[covidPlacentaHighRes@meta.data$PhenographClustering == 19] <- "Epithelium"

# Cluster 24 - Epithelium
covidPlacentaHighRes@meta.data$LowResolutionClustering[covidPlacentaHighRes@meta.data$PhenographClustering == 24] <- "Epithelium"
covidPlacentaHighRes@meta.data$InterResolutionClustering1[covidPlacentaHighRes@meta.data$PhenographClustering == 24] <- "Epithelium"
covidPlacentaHighRes@meta.data$HighResolutionClustering[covidPlacentaHighRes@meta.data$PhenographClustering == 24] <- "Epithelium"
covidPlacentaHighRes@meta.data$VMEclusters[covidPlacentaHighRes@meta.data$PhenographClustering == 24] <- "Epithelium"
covidPlacentaHighRes@meta.data$ClusteringForNeighbourhoods[covidPlacentaHighRes@meta.data$PhenographClustering == 24] <- "Epithelium"

# Cluster 37 - NKs?
covidPlacentaHighRes@meta.data$LowResolutionClustering[covidPlacentaHighRes@meta.data$PhenographClustering == 37] <- "NKs"
covidPlacentaHighRes@meta.data$InterResolutionClustering1[covidPlacentaHighRes@meta.data$PhenographClustering == 37] <- "NKs"
covidPlacentaHighRes@meta.data$HighResolutionClustering[covidPlacentaHighRes@meta.data$PhenographClustering == 37] <- "NKs"
covidPlacentaHighRes@meta.data$VMEclusters[covidPlacentaHighRes@meta.data$PhenographClustering == 37] <- "NKs"
covidPlacentaHighRes@meta.data$ClusteringForNeighbourhoods[covidPlacentaHighRes@meta.data$PhenographClustering == 37] <- "NKs"




# Export data for cluster validation





forClusterVal_inter <- cbind(covidPlacentaHighRes@meta.data$InterResolutionClustering1, covidPlacentaHighRes@meta.data$RegionIdents, covidPlacentaHighRes@Expression$MultiIHC@spatial.data)
colnames(forClusterVal_inter) <- c("InterCluster","Region","x","y")
write.csv(forClusterVal_inter, "D:/COVID_TISSUE_PROJECT/PlacentaPaper/ClusterVal/interResolution.csv")



forClusterVal_low <- cbind(covidPlacentaHighRes@meta.data$LowResolutionClustering, covidPlacentaHighRes@meta.data$RegionIdents, covidPlacentaHighRes@Expression$MultiIHC@spatial.data)
colnames(forClusterVal_low) <- c("LowCluster","Region","x","y")
write.csv(forClusterVal_low, "D:/COVID_TISSUE_PROJECT/PlacentaPaper/ClusterVal/lowResolution.csv")

### Check virus clustering

markers <- c(2,4,6,7,9,11,12,13,14,15,16,20,23,26)

covidPlacentaHighRes<- clusterMISSILe(MISSILeObject = covidPlacentaHighRes, markers = markers, includeMetaData = c("size"),
                                      numNeighbours = 11, expSet = "counts", clusteringName = "VirusClusterIdentification")


markers <- colnames(covidPlacentaHighRes@Expression$MultiIHC@counts)[c(2,4,6,7,9,11,12,13,14,15,16,20,23,26)]

BubblePlot(MISSILeObject = covidPlacentaHighRes, identities = "VirusClusterIdentification", markers = markers,
           threshold.percent = 0.5,
           colour.scale = c(0,1.5), transposePlot = TRUE)



forClusterVal_virus <- cbind(covidPlacentaHighRes@meta.data$VirusClusterIdentification, covidPlacentaHighRes@meta.data$RegionIdents, covidPlacentaHighRes@Expression$MultiIHC@spatial.data)
colnames(forClusterVal_virus) <- c("VirusCluster","Region","x","y")
write.csv(forClusterVal_virus, "D:/COVID_TISSUE_PROJECT/PlacentaPaper/ClusterVal/virusClusterResolution.csv")

#

virusCellRecluster <- covidPlacentaHighRes@Expression$MultiIHC@counts[covidPlacentaHighRes@meta.data$VirusClusterIdentification == 8, c(2,4,6,7,9,11,12,13,14,15,16,20,23,26)]

Rphenograph_out <- FastPG::fastCluster(data = as.matrix(virusCellRecluster), k = 400, num_threads = 14, verbose = TRUE)

Rphenograph_out$communities


# Dot plot

ignoreClusters <- NULL

identities <- Rphenograph_out$communities
uniqueClusters <- levels(as.factor(Rphenograph_out$communities))
plottingDF <- as.data.frame(phonTools::zeros(length(markers) * length(uniqueClusters),4))
colnames(plottingDF) <- c("Cluster","Protein","Percent.Exp","Mean.Exp")
plottingDF$Protein <- rep(markers, times = length(uniqueClusters))
for(i in 1:length(uniqueClusters)){
  if(i == 1){
    clusterColumn <- rep(uniqueClusters[i], times = length(markers))
  } else {
    clusterColumnTemp <- rep(uniqueClusters[i], times = length(markers))
    clusterColumn <- c(clusterColumn, clusterColumnTemp)
    rm(clusterColumnTemp)
  }
}
plottingDF$Cluster <- clusterColumn
rm(clusterColumn)
if(!is.null(ignoreClusters)){
  uniqueClusters <- uniqueClusters[! uniqueClusters %in% ignoreClusters]
}
dataToScale <- virusCellRecluster
#dataToScale <- covidPlacentaHighRes@Expression$MultiIHC@counts[!(covidPlacentaHighRes@meta.data[,identities] %in% ignoreClusters),markers]
scaledData <- as.data.frame(apply(dataToScale, MARGIN = 2, scale))

threshold <- 0.5

for(i in 1:length(uniqueClusters)){
  tempClusterData <- scaledData[Rphenograph_out$communities == uniqueClusters[i],]
  meanExpression <- phonTools::zeros(length(markers))
  aboveThreshold <- phonTools::zeros(length(markers))

  for(j in 1:length(markers)){

    meanExpression[j] <- mean(tempClusterData[,j])

    aboveThresholdTemp <- tempClusterData[,j] > threshold

    aboveThreshold[j] <- (length(aboveThresholdTemp[aboveThresholdTemp == TRUE]) / length(aboveThresholdTemp))*100

  }
  plottingDF$Percent.Exp[plottingDF$Cluster == uniqueClusters[i]] <- aboveThreshold
  plottingDF$Mean.Exp[plottingDF$Cluster == uniqueClusters[i]] <- meanExpression
}

plottingDF$Cluster <- as.numeric(plottingDF$Cluster)

plottingDF$Mean.Exp[plottingDF$Mean.Exp > 3] <- 3

plottingDF$Mean.Exp[plottingDF$Mean.Exp < 0] <- 0

# double for loop
# first loop through markers
# second loop through clusters

plot <- ggplot(data = plottingDF, mapping = aes_string(x = 'Protein', y = 'Cluster')) +
  geom_point(mapping = aes_string(size = 'Percent.Exp', color = 'Mean.Exp')) +
  guides(size = guide_legend(title = 'Percent Expressed')) + theme_classic() +
  xlab("") + ylab("") +
  theme(axis.text.x = element_text(angle = 45, hjust = 1, size = 12),
        axis.text.y = element_text(size = 12)) + scale_colour_gradient2(
          low = "#132B43",
          high = "#FF0000",
          space = "Lab",
          na.value = "grey50",
          guide = "colourbar",
          aesthetics = "colour"
        )

plot +
  scale_y_continuous(breaks = unique(plottingDF$Cluster))




table(covidPlacentaHighRes@meta.data$LowResolutionClustering[covidPlacentaHighRes@meta.data$VirusClusterIdentification == 8])
table(covidPlacentaHighRes@meta.data$InterResolutionClustering1[covidPlacentaHighRes@meta.data$VirusClusterIdentification == 8])

covidPlacentaHighRes@meta.data$InterResolutionClustering2 <- covidPlacentaHighRes@meta.data$InterResolutionClustering1


# Set all SARS-CoV-2+ back to regular
# Set SARS and cell from inter1 equal to inter2

covidPlacentaHighRes@meta.data$InterResolutionClustering2[covidPlacentaHighRes@meta.data$InterResolutionClustering2 == "SARS-CoV-2+ Macrophages"] <- "M1 Macrophages"
covidPlacentaHighRes@meta.data$InterResolutionClustering2[covidPlacentaHighRes@meta.data$InterResolutionClustering2 == "SARS-CoV-2+ Neutrophils"] <- "Neutrophils"
covidPlacentaHighRes@meta.data$InterResolutionClustering2[covidPlacentaHighRes@meta.data$InterResolutionClustering2 == "SARS-CoV-2+ Epithelium"] <- "Epithelium"

covidPlacentaHighRes@meta.data$InterResolutionClustering2[covidPlacentaHighRes@meta.data$VirusClusterIdentification == 8 & covidPlacentaHighRes@meta.data$InterResolutionClustering1 == "SARS-CoV-2+ Macrophages"] <- "SARS-CoV-2+ Macrophages"
covidPlacentaHighRes@meta.data$InterResolutionClustering2[covidPlacentaHighRes@meta.data$VirusClusterIdentification == 8 & covidPlacentaHighRes@meta.data$InterResolutionClustering1 == "SARS-CoV-2+ Epithelium"] <- "SARS-CoV-2+ Epithelium"
covidPlacentaHighRes@meta.data$InterResolutionClustering2[covidPlacentaHighRes@meta.data$VirusClusterIdentification == 8 & covidPlacentaHighRes@meta.data$InterResolutionClustering1 == "Unclassified"] <- "SARS-CoV-2+ Epithelium"

saveRDS(covidPlacentaHighRes,"D:/COVID_TISSUE_PROJECT/PlacentaPaper/RDS_Objects/highResClusteringMISSILeObject_withNewVirus.rds")

# Set up clusters for virus microenvironment (control epithelium in same samples)

covidPlacentaHighRes@meta.data$VMEclustersRedo <- covidPlacentaHighRes@meta.data$InterResolutionClustering2

numEpi <- length(covidPlacentaHighRes@meta.data$VMEclustersRedo[covidPlacentaHighRes@meta.data$VMEclustersRedo == "Epithelium" & covidPlacentaHighRes@meta.data$Disease == "COVID_Placenta"])
numVirusEpi <- length(covidPlacentaHighRes@meta.data$VMEclustersRedo[covidPlacentaHighRes@meta.data$VMEclustersRedo == "SARS-CoV-2+ Epithelium"])

randomNums <- sample(1:numEpi,numVirusEpi)


covidPlacentaHighRes@meta.data$VMEclustersRedo <- as.character(covidPlacentaHighRes@meta.data$VMEclustersRedo)
covidPlacentaHighRes@meta.data$VMEclustersRedo[covidPlacentaHighRes@meta.data$VMEclustersRedo == "Epithelium" & covidPlacentaHighRes@meta.data$Disease == "COVID_Placenta"][randomNums] <- "Non-infected epithelium"
covidPlacentaHighRes@meta.data$VMEclustersRedo <- as.factor(covidPlacentaHighRes@meta.data$VMEclustersRedo)


covidPlacentaHighRes@current.identity$ActiveIdents <- "VMEclustersRedo"

# VME redo

emptyList <- vector(mode = "list", length = 0)
covidPlacentaHighRes@neighbourhoods$uniqueNeighbourhood <- emptyList

funcMarkers <- c(1,3,5,6,8,10,19,21,22,25,27,30,31)

covidPlacentaHighRes <- calculateUniqueNeighbourhood(MISSILeObject = covidPlacentaHighRes, functionalMarkers = funcMarkers, cellOfInterest = "SARS-CoV-2+ Epithelium")

covidPlacentaHighRes <- calculateUniqueNeighbourhood(MISSILeObject = covidPlacentaHighRes, functionalMarkers = funcMarkers, cellOfInterest = "Non-infected epithelium")


saveRDS(covidPlacentaHighRes, "D:/COVID_TISSUE_PROJECT/PlacentaPaper/RDS_Objects/highResClusteringMISSILeObject_withNewVirus_andVME_Redo.rds")


bb <- unique(covidPlacentaHighRes@meta.data$VMEclustersRedo)
phenotypeToPlot <- grep('Epithelium', bb, value=TRUE)
phenotypeToPlot

phenotypeToPlot <- c("PDL1+ M2 Macrophages")

plotContinuousUniqueNeighbourhood(MISSILeObject = covidPlacentaHighRes, phenotype = phenotypeToPlot,
                                  neighbourhoodName = "SARS-CoV-2+ Epithelium", cellComparisons = c("SARS-CoV-2+ Epithelium","Non-infected epithelium"),
                                  colours = c("#00e1ff","#ff8c00"))



covidPlacentaHighRes <- readRDS("D:/COVID_TISSUE_PROJECT/PlacentaPaper/RDS_Objects/highResClusteringMISSILeObject_withNewVirus_andVME_Redo.rds")


# Neighbourhoods - these are from NeighbourhoodTesting.R

covidPlacentaHighRes@current.identity$ActiveIdents <- "LowResolutionClustering"

covidPlacentaHighRes <- cellNeighbourhoods(MISSILeObject = covidPlacentaHighRes,
                          numOfCells = 7, kMeans = 6,
                          verbose = TRUE,
                          ignoreClusters = NULL,
                          neighbourhoodName = "GeneralNeighbourhood1")

spatialNeighbourhood(MISSILeObject = covidPlacentaHighRes, region = 18, neighbourhoodName = "GeneralNeighbourhood1", pt.size = 0.9)

neighbourhoodEnrichment(MISSILeObject = covidPlacentaHighRes, neighbourhoodName = "GeneralNeighbourhood1", ignoreClusters = c("Noise"))

a <- neighbourhoodFrequencies(MISSILeObject = covidPlacentaHighRes, neighbourhoodName = "GeneralNeighbourhood1", ignoreClusters = c("Noise"), comparison = "Disease3")

a$Condition <- factor(a$Condition, levels = c("COVID_Placenta","COVID_Mother_CHI","COVID_Mother_CV","CHI","VUE","Normal"))

bxp <- ggboxplot(
  a, x = "Condition", y = "Frequency", fill = "Condition", color = "black", palette = "nejm",
  scales = "free", facet.by = "Neighbourhood", ncol = 8, shape = 17, size = 1.1, alpha = 1
)

my_comparisons <- list(c("COVID_Placenta", "COVID_Mother_CHI"), c("COVID_Placenta", "COVID_Mother_CV"), c("COVID_Mother_CHI", "COVID_Mother_CV"), c("COVID_Placenta", "CHI"), c("COVID_Placenta", "VUE"), c("COVID_Placenta", "Normal") )

bxp + xlab("") + ylab("Frequency [%]") + theme_bw() +
  theme(strip.text = element_text(colour = "black", face = "bold", size = 9.5)) +
  theme(strip.background = element_rect(fill=NA, colour = "black", linetype="solid"))  + theme(axis.text.x = element_text(color = "black", size = 9),
                                                                                               axis.text.y = element_text(color = "black", size = 11),
                                                                                               axis.title.x = element_text(color = "black", size = 12),
                                                                                               axis.title.y = element_text(color = "black", size = 13),
                                                                                               legend.position = "none",
                                                                                               panel.grid.major = element_blank(),
                                                                                               panel.grid.minor = element_blank()) + rotate_x_text(45) +
  stat_compare_means(comparisons = my_comparisons, label = "p.signif", hide.ns = TRUE)


saveRDS(covidPlacentaHighRes, "D:/COVID_TISSUE_PROJECT/PlacentaPaper/RDS_Objects/highResClusteringMISSILeObject_withNewVirus_andVME_andNeighbourhood_Redo.rds")

plotVoronoi(covidPlacentaHighRes@delaunay$deltri@cleaned.deldir[[17]], colScheme = "rainbow", polygonCentre = c(2016,1512) , polygonRadius = 1512,
            showPoints = FALSE, number = FALSE, backGround = "black", clusterMatrix = aaa,
            lineWidth = 3, edgeTiles = TRUE, fontRatio = 0.6, uniquePhenotypes = levels(covidPlacentaHighRes@meta.data$LowResolutionClustering))

#################### Extracting data for Noah ####################

covidPlacentaHighRes <- readRDS("D:/COVID_TISSUE_PROJECT/PlacentaPaper/RDS_Objects/highResClusteringMISSILeObject_withNewVirus_andVME_andNeighbourhood_Redo.rds")

region <- "J1"

phenotypeTable <- cbind(covidPlacentaHighRes@Expression$MultiIHC@spatial.data[covidPlacentaHighRes@meta.data$RegionIdents == region,],covidPlacentaHighRes@meta.data$LowResolutionClustering[covidPlacentaHighRes@meta.data$RegionIdents == region])
colnames(phenotypeTable) <- c("x","y","phenotype")

#phenotypeTable <- phenotypeTable[phenotypeTable$phenotype == "Epithelium",]

phenotypeTable$phenotype[phenotypeTable$phenotype == "Unclassified"] <- "Epithelium"

low_x <- 1000
high_x <- 3000
low_y <- 1000
high_y <- 3000

polygon <- data.frame(c(low_x, high_x, high_x, low_x),c(low_y, low_y, high_y, high_y))
inPolygonIDs <- secr::pointsInPolygon(phenotypeTable, polygon, logical = TRUE)

phenotypeTable <- phenotypeTable[inPolygonIDs,]
phenotypeTable$y <- abs(phenotypeTable$y - max(phenotypeTable$y))

#phenotypeTable$y <- phenotypeTable$y - low_y
phenotypeTable$x <- phenotypeTable$x - low_x

ggscatter(phenotypeTable, x = "x", y = "y", color = "phenotype")


write.csv(phenotypeTable, "D:/Noah_Greenwald/TestData/J1/phenotypeTable.csv", row.names = FALSE)
