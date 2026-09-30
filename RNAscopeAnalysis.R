library(stringr)

#sample5759 <- read.delim("F:/RNAscope/COVID-19/QuPathSegmentationOutput/5759.txt")

filelist = list.files(path = "F:/RNAscope/COVID-19/QuPathSegmentationOutput/", pattern = ".*.txt")

filelist <- stringr::str_c("F:/RNAscope/COVID-19/QuPathSegmentationOutput/",filelist)

rnaScopeList <- lapply(filelist, function(x)read.delim(x))

for(i in 1:length(rnaScopeList)){
  colnames(rnaScopeList[[i]]) <- colnames(rnaScopeList[[1]])
}

#rnaScopeList <- do.call("rbind", rnaScopeList)

rnaScopeDF <- do.call(rbind.data.frame, rnaScopeList)

rm(rnaScopeList)

# Reduce DF columns to ones needed

rnaScopeDFreduced <- rnaScopeDF[,c(1:8,14,62,72,76,80,84,88,92,96)]

# Plot all markers across slides

ggboxplot(rnaScopeDFreduced, x = "Image", y = "Cell..Opal.780.mean")

# Set positivity status on each channel
rnaScopeDFreduced$Macrophage <- rep("TRUE", nrow(rnaScopeDFreduced))
rnaScopeDFreduced$Macrophage[rnaScopeDFreduced$Cell..Opal.780.mean < 2] <- FALSE

rnaScopeDFreduced$IL6 <- rep("FALSE", nrow(rnaScopeDFreduced))
rnaScopeDFreduced$IL6[rnaScopeDFreduced$Cell..Opal.690.mean > 3] <- TRUE

rnaScopeDFreduced$IL8 <- rep("FALSE", nrow(rnaScopeDFreduced))
rnaScopeDFreduced$IL8[rnaScopeDFreduced$Cell..Opal.620.mean > 3] <- TRUE

rnaScopeDFreduced$CXCL10 <- rep("FALSE", nrow(rnaScopeDFreduced))
rnaScopeDFreduced$CXCL10[rnaScopeDFreduced$Cell..Opal.520.mean > 2.5] <- TRUE

# Create macrophage dataframe for plotting

rnaScopeDFreduced$SlideID <- sapply(strsplit(rnaScopeDFreduced$Image,"_"),"[[",1)

rnaScopeDFmacrophages <- rnaScopeDFreduced

  #rnaScopeDFreduced[rnaScopeDFreduced$Macrophage == "TRUE",]

# Organise samples, slides and disease states in the dataframe

# Slide ID
#rnaScopeDFmacrophages$SlideID <- sapply(strsplit(rnaScopeDFmacrophages$Image,"_"),"[[",1)

# Sample ID
rnaScopeDFmacrophages$SampleID <- rnaScopeDFmacrophages$SlideID

rnaScopeDFmacrophages$SampleID[rnaScopeDFmacrophages$SampleID == 3 & rnaScopeDFmacrophages$Parent == "TopRight"] <- "189"
rnaScopeDFmacrophages$SampleID[rnaScopeDFmacrophages$SampleID == 3 & rnaScopeDFmacrophages$Parent == "BottomRight"] <- "129"
rnaScopeDFmacrophages$SampleID[rnaScopeDFmacrophages$SampleID == 3 & rnaScopeDFmacrophages$Parent == "BottomLeft"] <- "6868"

rnaScopeDFmacrophages$SampleID[rnaScopeDFmacrophages$SampleID == 4 & rnaScopeDFmacrophages$Parent == "TopRight"] <- "5889"
rnaScopeDFmacrophages$SampleID[rnaScopeDFmacrophages$SampleID == 4 & rnaScopeDFmacrophages$Parent == "BottomRight"] <- "799"
rnaScopeDFmacrophages$SampleID[rnaScopeDFmacrophages$SampleID == 4 & rnaScopeDFmacrophages$Parent == "BottomLeft"] <- "835"
rnaScopeDFmacrophages$SampleID[rnaScopeDFmacrophages$SampleID == 4 & rnaScopeDFmacrophages$Parent == "TopLeft"] <- "6827"

rnaScopeDFmacrophages$SampleID[rnaScopeDFmacrophages$SampleID == 5 & rnaScopeDFmacrophages$Parent == "TopRight"] <- "1238"
rnaScopeDFmacrophages$SampleID[rnaScopeDFmacrophages$SampleID == 5 & rnaScopeDFmacrophages$Parent == "BottomRight"] <- "7363"
rnaScopeDFmacrophages$SampleID[rnaScopeDFmacrophages$SampleID == 5 & rnaScopeDFmacrophages$Parent == "BottomLeft"] <- "1005"
rnaScopeDFmacrophages$SampleID[rnaScopeDFmacrophages$SampleID == 5 & rnaScopeDFmacrophages$Parent == "TopLeft"] <- "1380"

# Disease state
rnaScopeDFmacrophages$DiseaseState <- NA

diseaseStates <- c("M+P+","M+P- CHI","Control CHI","M+P- CV","Control CHI","M+P+",
                   "Control CV","Control","Control CV","Control","M+P- CV","M+P+",
                   "M+P- CHI","M+P+","M+P- CV","Control CHI","M+P- CV","M+P- CHI",
                   "M+P+")

samples <- unique(rnaScopeDFmacrophages$SampleID)

for(i in 1:length(diseaseStates)){

  rnaScopeDFmacrophages$DiseaseState[rnaScopeDFmacrophages$SampleID == samples[i]] <- diseaseStates[i]

}

# Plot

library(phonTools)

# Target+ Macrophage abundance as a percent of total cells - Fig. 4g

samples <- unique(rnaScopeDFmacrophages$SampleID)
toPlot <- as.data.frame(zeros(length(samples),5))
colnames(toPlot) <- c("Sample","Disease","IL6+ Macrophage Abundance","CXCL10+ Macrophage Abundance","IL8+ Macrophage Abundance")
toPlot$Sample <- samples
toPlot$Disease <- diseaseStates

target <- c("IL6","CXCL10","IL8")

for(i in 1:(length(samples))){
  for(j in 1:length(target)){
    tempSample <- rnaScopeDFmacrophages[rnaScopeDFmacrophages[,"SampleID"] == samples[i],]
    toPlot[i,j+2] <- (nrow(tempSample[tempSample[,target[j]] == "TRUE" & tempSample[,"Macrophage"] == "TRUE",]) / nrow(tempSample))*100
  }
}

saveRDS(toPlot, "D:/COVID_TISSUE_PROJECT/PlacentaPaper/REVIEWER_COMMENTS/Updated_Figures/Fig4/Code_for_plotting/fig_4g_data.rds")

# Target+ Macrophage abundance as a percent of total macrophages

samples <- unique(rnaScopeDFmacrophages$SampleID)
toPlot <- as.data.frame(zeros(length(samples),5))
colnames(toPlot) <- c("Sample","Disease","IL6+ Macrophage Abundance","CXCL10+ Macrophage Abundance","IL8+ Macrophage Abundance")
toPlot$Sample <- samples
toPlot$Disease <- diseaseStates

target <- c("IL6","CXCL10","IL8")

for(i in 1:(length(samples))){
  for(j in 1:length(target)){
    tempSample <- rnaScopeDFmacrophages[rnaScopeDFmacrophages[,"SampleID"] == samples[i],]
    toPlot[i,j+2] <- (nrow(tempSample[tempSample[,target[j]] == "TRUE" & tempSample[,"Macrophage"] == "TRUE",]) / nrow(tempSample[tempSample[,"Macrophage"] == "TRUE",]))*100
  }
}

saveRDS(toPlot, "D:/COVID_TISSUE_PROJECT/PlacentaPaper/REVIEWER_COMMENTS/Updated_Figures/Fig4/Code_for_plotting/ex_4_a_data.rds")

## Just target positive cells as a percent of total cells

samples <- unique(rnaScopeDFmacrophages$SampleID)
toPlot <- as.data.frame(zeros(length(samples),5))
colnames(toPlot) <- c("Sample","Disease","IL6+ Macrophage Abundance","CXCL10+ Macrophage Abundance","IL8+ Macrophage Abundance")
toPlot$Sample <- samples
toPlot$Disease <- diseaseStates

target <- c("IL6","CXCL10","IL8")

for(i in 1:(length(samples))){
  for(j in 1:length(target)){
    tempSample <- rnaScopeDFmacrophages[rnaScopeDFmacrophages[,"SampleID"] == samples[i],]
    toPlot[i,j+2] <- (nrow(tempSample[tempSample[,target[j]] == "TRUE",]) / nrow(tempSample))*100
  }
}

saveRDS(toPlot, "D:/COVID_TISSUE_PROJECT/PlacentaPaper/REVIEWER_COMMENTS/Updated_Figures/Fig4/Code_for_plotting/target_cells_of_all_cells.rds")

toPlot <- readRDS("D:/COVID_TISSUE_PROJECT/PlacentaPaper/REVIEWER_COMMENTS/Updated_Figures/Fig4/Code_for_plotting/target_cells_of_all_cells.rds")

##### Plot

toPlot$Disease <- factor(toPlot$Disease, levels = c("M+P+", "M+P- CHI", "M+P- CV", "Control CHI", "Control CV", "Control"))

#my_comparisons <- list( c("M+P+", "M+P- CHI"), c("M+P+", "M+P- CV"), c("M+P+", "Control CHI"), c("M+P+", "Control CV"), c("M+P+", "Control"), c("M+P- CHI", "M+P- CV") )
#my_comparisons2 <- list( c("M+P+", "M+P- CHI"), c("M+P+", "M+P- CV"), c("M+P- CHI", "M+P- CV") )

toPlotLong <- tidyr::gather(toPlot, marker, abundance,
                            `IL6+ Macrophage Abundance`:`IL8+ Macrophage Abundance`,
                            factor_key = TRUE)

toPlotLongCut <- toPlotLong[toPlotLong$Disease %in% c("M+P+","M+P- CHI","M+P- CV"),]

toPlotLongCut$marker <- factor(toPlotLongCut$marker, levels = c("CXCL10+ Macrophage Abundance", "IL8+ Macrophage Abundance", "IL6+ Macrophage Abundance"))

toPlotLongCut$Disease <- factor(toPlotLongCut$Disease, levels = c("M+P+","M+P- CHI","M+P- CV"))

ggboxplot(toPlotLongCut, x = "marker", y = "abundance", fill = "Disease", size = 1.05, palette = "npg") + ylab("Abundance [%]") + xlab("") + theme_bw() +
  theme(panel.grid.major = element_blank(),
        panel.grid.minor = element_blank()) #+      # Add global p-value
  #stat_compare_means(comparisons = my_comparisons2)

# CXCL10 comparisons

wilcox.test(x = toPlotLongCut$abundance[toPlotLongCut$Disease == "M+P+" & toPlotLongCut$marker == "CXCL10+ Macrophage Abundance"],
            y = toPlotLongCut$abundance[toPlotLongCut$Disease == "M+P- CHI" & toPlotLongCut$marker == "CXCL10+ Macrophage Abundance"])
# p = 0.25

wilcox.test(x = toPlotLongCut$abundance[toPlotLongCut$Disease == "M+P+" & toPlotLongCut$marker == "CXCL10+ Macrophage Abundance"],
            y = toPlotLongCut$abundance[toPlotLongCut$Disease == "M+P- CV" & toPlotLongCut$marker == "CXCL10+ Macrophage Abundance"])
# p = 0.1111

wilcox.test(x = toPlotLongCut$abundance[toPlotLongCut$Disease == "M+P- CHI" & toPlotLongCut$marker == "CXCL10+ Macrophage Abundance"],
            y = toPlotLongCut$abundance[toPlotLongCut$Disease == "M+P- CV" & toPlotLongCut$marker == "CXCL10+ Macrophage Abundance"])
# p = 0.8571

# IL8 comparisons

wilcox.test(x = toPlotLongCut$abundance[toPlotLongCut$Disease == "M+P+" & toPlotLongCut$marker == "IL8+ Macrophage Abundance"],
            y = toPlotLongCut$abundance[toPlotLongCut$Disease == "M+P- CHI" & toPlotLongCut$marker == "IL8+ Macrophage Abundance"])
# p = 0.5714

wilcox.test(x = toPlotLongCut$abundance[toPlotLongCut$Disease == "M+P+" & toPlotLongCut$marker == "IL8+ Macrophage Abundance"],
            y = toPlotLongCut$abundance[toPlotLongCut$Disease == "M+P- CV" & toPlotLongCut$marker == "IL8+ Macrophage Abundance"])
# p = 0.01587

wilcox.test(x = toPlotLongCut$abundance[toPlotLongCut$Disease == "M+P- CHI" & toPlotLongCut$marker == "IL8+ Macrophage Abundance"],
            y = toPlotLongCut$abundance[toPlotLongCut$Disease == "M+P- CV" & toPlotLongCut$marker == "IL8+ Macrophage Abundance"])
# p = 0.05714

# IL6 comparisons

wilcox.test(x = toPlotLongCut$abundance[toPlotLongCut$Disease == "M+P+" & toPlotLongCut$marker == "IL6+ Macrophage Abundance"],
            y = toPlotLongCut$abundance[toPlotLongCut$Disease == "M+P- CHI" & toPlotLongCut$marker == "IL6+ Macrophage Abundance"])
# p = 0.03571

wilcox.test(x = toPlotLongCut$abundance[toPlotLongCut$Disease == "M+P+" & toPlotLongCut$marker == "IL6+ Macrophage Abundance"],
            y = toPlotLongCut$abundance[toPlotLongCut$Disease == "M+P- CV" & toPlotLongCut$marker == "IL6+ Macrophage Abundance"])
# p = 0.01587

wilcox.test(x = toPlotLongCut$abundance[toPlotLongCut$Disease == "M+P- CHI" & toPlotLongCut$marker == "IL6+ Macrophage Abundance"],
            y = toPlotLongCut$abundance[toPlotLongCut$Disease == "M+P- CV" & toPlotLongCut$marker == "IL6+ Macrophage Abundance"])
# p = 0.2286


##








pdf("D:/COVID_TISSUE_PROJECT/PlacentaPaper/REVIEWER_COMMENTS/Updated_Figures/Ex4/rnascope.pdf", width = 4.5, height = 6)

ggboxplot(toPlotLongCut, x = "marker", y = "abundance",ncol = 6, palette = c("#D43F3A","#EEA236","#46B8DA"), fill = "Disease", color = "black", size = 1.05, scales = "free",
          bxp.errorbar = TRUE, bxp.errorbar.width = 0.3) +
  xlab("")   + ylab("Abundance [% of total cells]") + ggdist::theme_ggdist() +
  geom_jitter(
    aes(x = marker, y = abundance, fill = Disease),
    position = position_jitterdodge(jitter.width = 0, dodge.width = 0.8),
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
            axis.ticks = element_line(color="black"))  + rotate_x_text(60) 

dev.off()









