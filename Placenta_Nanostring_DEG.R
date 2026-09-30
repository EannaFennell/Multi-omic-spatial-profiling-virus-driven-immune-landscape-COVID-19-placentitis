# To load all the data for this, see 'D:\COVID_TISSUE_PROJECT\PlacentaPaper\GeoMX\Placenta_Nanostring.R'

# Align this data with the annotation file

library(glmmSeq)

regions <- unique(annotationFile$Compartment)
whichRegion <- 4

readsTrophoblast <- reads[,annotationFile[,"Compartment"] == regions[whichRegion]]
annoTrophoblast <- annotationFile[annotationFile[,"Compartment"] == regions[whichRegion],]

# Virus positive trophoblast region vs. virus negative trophoblast region

diseases <- c("COVID+ Placenta") #,"COVID+ Mother")

readsTrophoblastModel <- readsTrophoblast[,annoTrophoblast$DiseaseState %in% diseases]
annoTrophoblastModel <- annoTrophoblast[annoTrophoblast$DiseaseState %in% diseases,]

annoTrophoblastModel$VirusPositive2 <- annoTrophoblastModel$`Nucleocapsid IHC`
annoTrophoblastModel$VirusPositive2[annoTrophoblastModel$`Nucleocapsid IHC` %in% c("Strong")] <- TRUE
annoTrophoblastModel$VirusPositive2[annoTrophoblastModel$`Nucleocapsid IHC` %in% c("Negative","Weak","NA")] <- FALSE

disp <- apply(readsTrophoblastModel, 1, function(x){
  (var(x, na.rm=TRUE)-mean(x, na.rm=TRUE))/(mean(x, na.rm=TRUE)**2)
})

head(disp)

colnames(annoTrophoblastModel)[c(19:20)] <- c("Gestational_Age","Maternal_Age")

statResults <- glmmSeq::glmmSeq(~ VirusPositive2 + (1|Number),
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
  foldChange[i,1] = log2(mean(as.numeric(readsTrophoblastModel[rownames(stats)[i],annoTrophoblastModel$VirusPositive2 == TRUE]))/mean(as.numeric(readsTrophoblastModel[rownames(stats)[i],annoTrophoblastModel$VirusPositive2 == FALSE])) )
}

stats <- cbind(stats,foldChange)

stats$adjustedP <- p.adjust(stats$P_VirusPositive2, method = "bonferroni")

write.table(stats, file = "D:/COVID_TISSUE_PROJECT/PlacentaPaper/GeoMX/DEG/LegitModel/glmmSeq_COVID_placenta_virus_positive_vs_virus_negative_trophoblast_LEGIT_adding_no_fixed_effects_age.csv", sep = ",", quote = FALSE)

stats <- read.table("D:/COVID_TISSUE_PROJECT/PlacentaPaper/GeoMX/DEG/LegitModel/Keep/glmmSeq_COVID_placenta_virus_positive_vs_virus_negative_trophoblast_LEGIT_adding_no_fixed_effects_age.csv", sep = ",")

colnames(stats) <- stats[1,]
rownames(stats) <- stats[,1]
stats <- stats[c(2:nrow(stats)), c(2:ncol(stats))]

stats$foldChange <- as.numeric(stats$foldChange)
stats$P_VirusPositive2 <- as.numeric(stats$P_VirusPositive2)
#stats$adjustedP <- as.numeric(stats$adjustedP)

keyvals.colour <- as.character(zeros(nrow(stats)) * NA)

keyvals.colour[stats$foldChange > 0.5] <- "red"
keyvals.colour[stats$foldChange < -0.5] <- "deepskyblue"
keyvals.colour[stats$P_VirusPositive2 > 0.05] <- "grey80"
keyvals.colour[is.na(keyvals.colour)] <- 'grey80'

names(keyvals.colour) <- rownames(stats)

EnhancedVolcano(stats,
                lab = rownames(stats),
                x = 'foldChange',
                y = 'P_VirusPositive2',
                cutoffLineType = 'blank',
                selectLab = c('VSIG8',"CCF1","PFN1","RAB1A","FUS","CCR1","CAMK2N2"),
                cutoffLineWidth = 0.8,
                colCustom = keyvals.colour,
                pointSize = 2.0,
                pCutoff = 0.05,
                FCcutoff = 0.75,
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
                border = "full",
                gridlines.major = FALSE,
                gridlines.minor = FALSE) + theme(legend.position = "none") + theme(
                  axis.title.x = element_text(size = 12),
                  axis.text.x = element_text(size = 12),
                  axis.title.y = element_text(size = 12),
                  axis.text.y = element_text(size = 12)) + xlab("Log2FC") + ylim(0, 25)






# Trophoblast region of Mother+Placenta+ vs. Mother+Placenta-

diseases <- c("COVID+ Placenta","COVID+ Mother")

readsTrophoblastModel <- readsTrophoblast[,annoTrophoblast$DiseaseState %in% diseases]
annoTrophoblastModel <- annoTrophoblast[annoTrophoblast$DiseaseState %in% diseases,]

annoTrophoblastModel$VirusPositive2 <- annoTrophoblastModel$`Nucleocapsid IHC`
annoTrophoblastModel$VirusPositive2[annoTrophoblastModel$`Nucleocapsid IHC` %in% c("Strong")] <- TRUE
annoTrophoblastModel$VirusPositive2[annoTrophoblastModel$`Nucleocapsid IHC` %in% c("Negative","Weak")] <- FALSE

disp <- apply(readsTrophoblastModel, 1, function(x){
  (var(x, na.rm=TRUE)-mean(x, na.rm=TRUE))/(mean(x, na.rm=TRUE)**2)
})

head(disp)

colnames(annoTrophoblastModel)[c(19:20)] <- c("Gestational_Age","Maternal_Age")

statResults <- glmmSeq::glmmSeq(~ DiseaseState + Gestational_Age + (1|Number),
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
  foldChange[i,1] = log2(mean(as.numeric(readsTrophoblastModel[rownames(stats)[i],annoTrophoblastModel$DiseaseState == diseases[1]]))/mean(as.numeric(readsTrophoblastModel[rownames(stats)[i],annoTrophoblastModel$DiseaseState == diseases[2]])) )
}

stats <- cbind(stats,foldChange)

stats$adjustedP <- p.adjust(stats$P_DiseaseState, method = "bonferroni")

write.table(stats, file = "D:/COVID_TISSUE_PROJECT/PlacentaPaper/GeoMX/DEG/LegitModel/glmmSeq_COVID_placenta_vs_COVID_mother_trophoblast_LEGIT_adding_no_maternal_age.csv", sep = ",", quote = FALSE)

keyvals.colour <- as.character(zeros(nrow(stats)) * NA)

keyvals.colour[stats$foldChange > 0] <- "red"
keyvals.colour[stats$foldChange < 0] <- "deepskyblue"
keyvals.colour[stats$P_DiseaseState > 0.05] <- "grey80"
keyvals.colour[is.na(keyvals.colour)] <- 'grey80'
keyvals.colour[rownames(stats) %in% c("HLA-B")] <- "green"

names(keyvals.colour) <- rownames(stats)


#stats$P_DiseaseState[stats$P_DiseaseState < 0.00001] <- 0.00001

EnhancedVolcano(stats,
                lab = rownames(stats),
                x = 'foldChange',
                y = 'P_DiseaseState',
                cutoffLineType = 'blank',
                selectLab = c('B2M','LYZ',"HLA-B","MSR1","PPT1","SEC22B","GPX3"),
                cutoffLineWidth = 0.8,
                colCustom = keyvals.colour,
                pointSize = 1.5,
                pCutoff = 0.05,
                FCcutoff = 0.5,
                labFace = 'italic',
                #boxedLabels = TRUE,
                #drawConnectors = TRUE,
                widthConnectors = 1.0,
                colConnectors = 'black',
                title = "",
                subtitle = "",
                caption = "",
                labSize = 4.0,
                colAlpha = 0.9,
                gridlines.major = FALSE,
                border = "full",
                gridlines.minor = FALSE) + theme(legend.position = "none") + theme(
                  axis.title.x = element_text(size = 11),
                  axis.text.x = element_text(size = 10),
                  axis.title.y = element_text(size = 11),
                  axis.text.y = element_text(size = 10)) + xlab("Log2FC")  + ylab(expression(-log[10](P~value))) # + ylim(0, max(abs(stats$foldChange)) + 0.5)

# Macrophage - COVID Placenta vs COVID Mother

diseases <- c("COVID+ Placenta","COVID+ Mother")

readsTrophoblastModel <- readsTrophoblast[,annoTrophoblast$DiseaseState %in% diseases]
annoTrophoblastModel <- annoTrophoblast[annoTrophoblast$DiseaseState %in% diseases,]

annoTrophoblastModel$VirusPositive2 <- annoTrophoblastModel$`Nucleocapsid IHC`
annoTrophoblastModel$VirusPositive2[annoTrophoblastModel$`Nucleocapsid IHC` %in% c("Strong")] <- TRUE
annoTrophoblastModel$VirusPositive2[annoTrophoblastModel$`Nucleocapsid IHC` %in% c("Negative","Weak")] <- FALSE

disp <- apply(readsTrophoblastModel, 1, function(x){
  (var(x, na.rm=TRUE)-mean(x, na.rm=TRUE))/(mean(x, na.rm=TRUE)**2)
})

head(disp)

colnames(annoTrophoblastModel)[c(19:20)] <- c("Gestational_Age","Maternal_Age")

statResults <- glmmSeq::glmmSeq(~ DiseaseState + Gestational_Age + (1|Number),
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
  foldChange[i,1] = log2(mean(as.numeric(readsTrophoblastModel[rownames(stats)[i],annoTrophoblastModel$DiseaseState == diseases[1]]))/mean(as.numeric(readsTrophoblastModel[rownames(stats)[i],annoTrophoblastModel$DiseaseState == diseases[2]])) )
}

stats <- cbind(stats,foldChange)

stats$adjustedP <- p.adjust(stats$P_DiseaseState, method = "bonferroni")

write.table(stats, file = "D:/COVID_TISSUE_PROJECT/PlacentaPaper/GeoMX/DEG/LegitModel/glmmSeq_COVID_placenta_vs_COVID_mother_macrophage_LEGIT_adding_no_maternal_age.csv", sep = ",", quote = FALSE)

stats <-read.table("D:/COVID_TISSUE_PROJECT/PlacentaPaper/GeoMX/DEG/LegitModel/Keep/glmmSeq_COVID_placenta_vs_COVID_mother_macrophage_LEGIT_adding_no_maternal_age.csv", sep = ",")

colnames(stats) <- stats[1,]
rownames(stats) <- stats[,1]
stats <- stats[c(2:nrow(stats)), c(2:ncol(stats))]

stats$foldChange <- as.numeric(stats$foldChange)
stats$P_DiseaseState <- as.numeric(stats$P_DiseaseState)
stats$adjustedP <- as.numeric(stats$adjustedP)

keyvals.colour <- as.character(zeros(nrow(stats)) * NA)

keyvals.colour[stats$foldChange > 0] <- "red"
keyvals.colour[stats$foldChange < 0] <- "deepskyblue"
keyvals.colour[stats$P_DiseaseState > 0.05] <- "grey80"
keyvals.colour[is.na(keyvals.colour)] <- 'grey80'
#keyvals.colour[rownames(stats) %in% c("HLA-B")] <- "green"

names(keyvals.colour) <- rownames(stats)


#stats$P_DiseaseState[stats$P_DiseaseState < 0.00001] <- 0.00001

EnhancedVolcano(stats,
                lab = rownames(stats),
                x = 'foldChange',
                y = 'adjustedP',
                cutoffLineType = 'blank',
                selectLab = c('LAPTM5','B2M',"IFI30","FCGR2C","OSCAR","CTSZ"),
                cutoffLineWidth = 0.8,
                colCustom = keyvals.colour,
                pointSize = 2,
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
                colAlpha = 0.7,
                gridlines.major = FALSE,
                border = "full",
                gridlines.minor = FALSE) + theme(legend.position = "none") + theme(
                  axis.title.x = element_text(size = 11),
                  axis.text.x = element_text(size = 10),
                  axis.title.y = element_text(size = 11),
                  axis.text.y = element_text(size = 10)) + xlab("Log2FC")  + ylab(expression(-log[10](Q~value))) # + ylim(0, max(abs(stats$foldChange)) + 0.5)

# Macrophage - COVID Placenta vs CHI

diseases <- c("COVID+ Placenta","CHI")

readsTrophoblastModel <- readsTrophoblast[,annoTrophoblast$DiseaseState %in% diseases]
annoTrophoblastModel <- annoTrophoblast[annoTrophoblast$DiseaseState %in% diseases,]

annoTrophoblastModel$VirusPositive2 <- annoTrophoblastModel$`Nucleocapsid IHC`
annoTrophoblastModel$VirusPositive2[annoTrophoblastModel$`Nucleocapsid IHC` %in% c("Strong")] <- TRUE
annoTrophoblastModel$VirusPositive2[annoTrophoblastModel$`Nucleocapsid IHC` %in% c("Negative","Weak")] <- FALSE

disp <- apply(readsTrophoblastModel, 1, function(x){
  (var(x, na.rm=TRUE)-mean(x, na.rm=TRUE))/(mean(x, na.rm=TRUE)**2)
})

head(disp)

colnames(annoTrophoblastModel)[c(19:20)] <- c("Gestational_Age","Maternal_Age")

statResults <- glmmSeq::glmmSeq(~ DiseaseState + Gestational_Age + (1|Number),
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
  foldChange[i,1] = log2(mean(as.numeric(readsTrophoblastModel[rownames(stats)[i],annoTrophoblastModel$DiseaseState == diseases[1]]))/mean(as.numeric(readsTrophoblastModel[rownames(stats)[i],annoTrophoblastModel$DiseaseState == diseases[2]])) )
}

stats <- cbind(stats,foldChange)

stats$adjustedP <- p.adjust(stats$P_DiseaseState, method = "bonferroni")

write.table(stats, file = "D:/COVID_TISSUE_PROJECT/PlacentaPaper/GeoMX/DEG/LegitModel/glmmSeq_COVID_placenta_vs_CHI_macrophage_LEGIT_adding_no_maternal_age.csv", sep = ",", quote = FALSE)

keyvals.colour <- as.character(zeros(nrow(stats)) * NA)

keyvals.colour[stats$foldChange > 0] <- "red"
keyvals.colour[stats$foldChange < 0] <- "deepskyblue"
keyvals.colour[stats$P_DiseaseState > 0.05] <- "grey80"
keyvals.colour[is.na(keyvals.colour)] <- 'grey80'
keyvals.colour[rownames(stats) %in% c("HLA-B")] <- "green"

names(keyvals.colour) <- rownames(stats)


#stats$P_DiseaseState[stats$P_DiseaseState < 0.00001] <- 0.00001

EnhancedVolcano(stats,
                lab = rownames(stats),
                x = 'foldChange',
                y = 'P_DiseaseState',
                cutoffLineType = 'blank',
                #selectLab = c('B2M','LYZ',"HLA-B","MSR1","PPT1","SEC22B","GPX3"),
                cutoffLineWidth = 0.8,
                colCustom = keyvals.colour,
                pointSize = 1.5,
                pCutoff = 0.05,
                FCcutoff = 0.5,
                labFace = 'italic',
                #boxedLabels = TRUE,
                #drawConnectors = TRUE,
                widthConnectors = 1.0,
                colConnectors = 'black',
                title = "",
                subtitle = "",
                caption = "",
                labSize = 4.0,
                colAlpha = 0.9,
                gridlines.major = FALSE,
                border = "full",
                gridlines.minor = FALSE) + theme(legend.position = "none") + theme(
                  axis.title.x = element_text(size = 11),
                  axis.text.x = element_text(size = 10),
                  axis.title.y = element_text(size = 11),
                  axis.text.y = element_text(size = 10)) + xlab("Log2FC")  + ylab(expression(-log[10](P~value))) # + ylim(0, max(abs(stats$foldChange)) + 0.5)

# Decidua - COVID Placenta vs COVID Mother

diseases <- c("COVID+ Placenta","COVID+ Mother")

readsTrophoblastModel <- readsTrophoblast[,annoTrophoblast$DiseaseState %in% diseases]
annoTrophoblastModel <- annoTrophoblast[annoTrophoblast$DiseaseState %in% diseases,]

annoTrophoblastModel$VirusPositive2 <- annoTrophoblastModel$`Nucleocapsid IHC`
annoTrophoblastModel$VirusPositive2[annoTrophoblastModel$`Nucleocapsid IHC` %in% c("Strong")] <- TRUE
annoTrophoblastModel$VirusPositive2[annoTrophoblastModel$`Nucleocapsid IHC` %in% c("Negative","Weak")] <- FALSE

disp <- apply(readsTrophoblastModel, 1, function(x){
  (var(x, na.rm=TRUE)-mean(x, na.rm=TRUE))/(mean(x, na.rm=TRUE)**2)
})

head(disp)

colnames(annoTrophoblastModel)[c(19:20)] <- c("Gestational_Age","Maternal_Age")

statResults <- glmmSeq::glmmSeq(~ DiseaseState + Gestational_Age + (1|Number),
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
  foldChange[i,1] = log2(mean(as.numeric(readsTrophoblastModel[rownames(stats)[i],annoTrophoblastModel$DiseaseState == diseases[1]]))/mean(as.numeric(readsTrophoblastModel[rownames(stats)[i],annoTrophoblastModel$DiseaseState == diseases[2]])) )
}

stats <- cbind(stats,foldChange)

stats$adjustedP <- p.adjust(stats$P_DiseaseState, method = "bonferroni")

write.table(stats, file = "D:/COVID_TISSUE_PROJECT/PlacentaPaper/GeoMX/DEG/LegitModel/glmmSeq_COVID_placenta_vs_COVID_mother_decidua_LEGIT_adding_no_maternal_age.csv", sep = ",", quote = FALSE)

keyvals.colour <- as.character(zeros(nrow(stats)) * NA)

keyvals.colour[stats$foldChange > 0] <- "red"
keyvals.colour[stats$foldChange < 0] <- "deepskyblue"
keyvals.colour[stats$P_DiseaseState > 0.05] <- "grey80"
keyvals.colour[is.na(keyvals.colour)] <- 'grey80'
keyvals.colour[rownames(stats) %in% c("HLA-B")] <- "green"

names(keyvals.colour) <- rownames(stats)


#stats$P_DiseaseState[stats$P_DiseaseState < 0.00001] <- 0.00001

EnhancedVolcano(stats,
                lab = rownames(stats),
                x = 'foldChange',
                y = 'P_DiseaseState',
                cutoffLineType = 'blank',
                #selectLab = c('B2M','LYZ',"HLA-B","MSR1","PPT1","SEC22B","GPX3"),
                cutoffLineWidth = 0.8,
                colCustom = keyvals.colour,
                pointSize = 1.5,
                pCutoff = 0.05,
                FCcutoff = 0.5,
                labFace = 'italic',
                #boxedLabels = TRUE,
                #drawConnectors = TRUE,
                widthConnectors = 1.0,
                colConnectors = 'black',
                title = "",
                subtitle = "",
                caption = "",
                labSize = 4.0,
                colAlpha = 0.9,
                gridlines.major = FALSE,
                border = "full",
                gridlines.minor = FALSE) + theme(legend.position = "none") + theme(
                  axis.title.x = element_text(size = 11),
                  axis.text.x = element_text(size = 10),
                  axis.title.y = element_text(size = 11),
                  axis.text.y = element_text(size = 10)) + xlab("Log2FC")  + ylab(expression(-log[10](P~value))) # + ylim(0, max(abs(stats$foldChange)) + 0.5)

# Decidua - COVID Placenta vs CHI

diseases <- c("COVID+ Placenta","CHI")

readsTrophoblastModel <- readsTrophoblast[,annoTrophoblast$DiseaseState %in% diseases]
annoTrophoblastModel <- annoTrophoblast[annoTrophoblast$DiseaseState %in% diseases,]

annoTrophoblastModel$VirusPositive2 <- annoTrophoblastModel$`Nucleocapsid IHC`
annoTrophoblastModel$VirusPositive2[annoTrophoblastModel$`Nucleocapsid IHC` %in% c("Strong")] <- TRUE
annoTrophoblastModel$VirusPositive2[annoTrophoblastModel$`Nucleocapsid IHC` %in% c("Negative","Weak")] <- FALSE

disp <- apply(readsTrophoblastModel, 1, function(x){
  (var(x, na.rm=TRUE)-mean(x, na.rm=TRUE))/(mean(x, na.rm=TRUE)**2)
})

head(disp)

colnames(annoTrophoblastModel)[c(19:20)] <- c("Gestational_Age","Maternal_Age")

statResults <- glmmSeq::glmmSeq(~ DiseaseState + Gestational_Age + (1|Number),
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
  foldChange[i,1] = log2(mean(as.numeric(readsTrophoblastModel[rownames(stats)[i],annoTrophoblastModel$DiseaseState == diseases[1]]))/mean(as.numeric(readsTrophoblastModel[rownames(stats)[i],annoTrophoblastModel$DiseaseState == diseases[2]])) )
}

stats <- cbind(stats,foldChange)

stats$adjustedP <- p.adjust(stats$P_DiseaseState, method = "bonferroni")

write.table(stats, file = "D:/COVID_TISSUE_PROJECT/PlacentaPaper/GeoMX/DEG/LegitModel/glmmSeq_COVID_placenta_vs_CHI_decidua_LEGIT_adding_no_maternal_age.csv", sep = ",", quote = FALSE)

keyvals.colour <- as.character(zeros(nrow(stats)) * NA)

keyvals.colour[stats$foldChange > 0] <- "red"
keyvals.colour[stats$foldChange < 0] <- "deepskyblue"
keyvals.colour[stats$P_DiseaseState > 0.05] <- "grey80"
keyvals.colour[is.na(keyvals.colour)] <- 'grey80'
keyvals.colour[rownames(stats) %in% c("HLA-B")] <- "green"

names(keyvals.colour) <- rownames(stats)


#stats$P_DiseaseState[stats$P_DiseaseState < 0.00001] <- 0.00001

EnhancedVolcano(stats,
                lab = rownames(stats),
                x = 'foldChange',
                y = 'P_DiseaseState',
                cutoffLineType = 'blank',
                #selectLab = c('B2M','LYZ',"HLA-B","MSR1","PPT1","SEC22B","GPX3"),
                cutoffLineWidth = 0.8,
                colCustom = keyvals.colour,
                pointSize = 1.5,
                pCutoff = 0.05,
                FCcutoff = 0.5,
                labFace = 'italic',
                #boxedLabels = TRUE,
                #drawConnectors = TRUE,
                widthConnectors = 1.0,
                colConnectors = 'black',
                title = "",
                subtitle = "",
                caption = "",
                labSize = 4.0,
                colAlpha = 0.9,
                gridlines.major = FALSE,
                border = "full",
                gridlines.minor = FALSE) + theme(legend.position = "none") + theme(
                  axis.title.x = element_text(size = 11),
                  axis.text.x = element_text(size = 10),
                  axis.title.y = element_text(size = 11),
                  axis.text.y = element_text(size = 10)) + xlab("Log2FC")  + ylab(expression(-log[10](P~value))) # + ylim(0, max(abs(stats$foldChange)) + 0.5)

# Villous stroma - COVID Placenta vs COVID Mother

diseases <- c("COVID+ Placenta","COVID+ Mother")

readsTrophoblastModel <- readsTrophoblast[,annoTrophoblast$DiseaseState %in% diseases]
annoTrophoblastModel <- annoTrophoblast[annoTrophoblast$DiseaseState %in% diseases,]

annoTrophoblastModel$VirusPositive2 <- annoTrophoblastModel$`Nucleocapsid IHC`
annoTrophoblastModel$VirusPositive2[annoTrophoblastModel$`Nucleocapsid IHC` %in% c("Strong")] <- TRUE
annoTrophoblastModel$VirusPositive2[annoTrophoblastModel$`Nucleocapsid IHC` %in% c("Negative","Weak")] <- FALSE

disp <- apply(readsTrophoblastModel, 1, function(x){
  (var(x, na.rm=TRUE)-mean(x, na.rm=TRUE))/(mean(x, na.rm=TRUE)**2)
})

head(disp)

colnames(annoTrophoblastModel)[c(19:20)] <- c("Gestational_Age","Maternal_Age")

statResults <- glmmSeq::glmmSeq(~ DiseaseState + Gestational_Age + (1|Number),
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
  foldChange[i,1] = log2(mean(as.numeric(readsTrophoblastModel[rownames(stats)[i],annoTrophoblastModel$DiseaseState == diseases[1]]))/mean(as.numeric(readsTrophoblastModel[rownames(stats)[i],annoTrophoblastModel$DiseaseState == diseases[2]])) )
}

stats <- cbind(stats,foldChange)

stats$adjustedP <- p.adjust(stats$P_DiseaseState, method = "bonferroni")

write.table(stats, file = "D:/COVID_TISSUE_PROJECT/PlacentaPaper/GeoMX/DEG/LegitModel/glmmSeq_COVID_placenta_vs_COVID_Mother_villous_stroma_LEGIT_adding_no_maternal_age.csv", sep = ",", quote = FALSE)

keyvals.colour <- as.character(zeros(nrow(stats)) * NA)

keyvals.colour[stats$foldChange > 0] <- "red"
keyvals.colour[stats$foldChange < 0] <- "deepskyblue"
keyvals.colour[stats$P_DiseaseState > 0.05] <- "grey80"
keyvals.colour[is.na(keyvals.colour)] <- 'grey80'
keyvals.colour[rownames(stats) %in% c("HLA-B")] <- "green"

names(keyvals.colour) <- rownames(stats)


#stats$P_DiseaseState[stats$P_DiseaseState < 0.00001] <- 0.00001

EnhancedVolcano(stats,
                lab = rownames(stats),
                x = 'foldChange',
                y = 'P_DiseaseState',
                cutoffLineType = 'blank',
                #selectLab = c('B2M','LYZ',"HLA-B","MSR1","PPT1","SEC22B","GPX3"),
                cutoffLineWidth = 0.8,
                colCustom = keyvals.colour,
                pointSize = 1.5,
                pCutoff = 0.05,
                FCcutoff = 0.5,
                labFace = 'italic',
                #boxedLabels = TRUE,
                #drawConnectors = TRUE,
                widthConnectors = 1.0,
                colConnectors = 'black',
                title = "",
                subtitle = "",
                caption = "",
                labSize = 4.0,
                colAlpha = 0.9,
                gridlines.major = FALSE,
                border = "full",
                gridlines.minor = FALSE) + theme(legend.position = "none") + theme(
                  axis.title.x = element_text(size = 11),
                  axis.text.x = element_text(size = 10),
                  axis.title.y = element_text(size = 11),
                  axis.text.y = element_text(size = 10)) + xlab("Log2FC")  + ylab(expression(-log[10](P~value))) # + ylim(0, max(abs(stats$foldChange)) + 0.5)

# Villous stroma - COVID Placenta vs CHI

diseases <- c("COVID+ Placenta","CHI")

readsTrophoblastModel <- readsTrophoblast[,annoTrophoblast$DiseaseState %in% diseases]
annoTrophoblastModel <- annoTrophoblast[annoTrophoblast$DiseaseState %in% diseases,]

annoTrophoblastModel$VirusPositive2 <- annoTrophoblastModel$`Nucleocapsid IHC`
annoTrophoblastModel$VirusPositive2[annoTrophoblastModel$`Nucleocapsid IHC` %in% c("Strong")] <- TRUE
annoTrophoblastModel$VirusPositive2[annoTrophoblastModel$`Nucleocapsid IHC` %in% c("Negative","Weak")] <- FALSE

disp <- apply(readsTrophoblastModel, 1, function(x){
  (var(x, na.rm=TRUE)-mean(x, na.rm=TRUE))/(mean(x, na.rm=TRUE)**2)
})

head(disp)

colnames(annoTrophoblastModel)[c(19:20)] <- c("Gestational_Age","Maternal_Age")

statResults <- glmmSeq::glmmSeq(~ DiseaseState + Gestational_Age + (1|Number),
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
  foldChange[i,1] = log2(mean(as.numeric(readsTrophoblastModel[rownames(stats)[i],annoTrophoblastModel$DiseaseState == diseases[1]]))/mean(as.numeric(readsTrophoblastModel[rownames(stats)[i],annoTrophoblastModel$DiseaseState == diseases[2]])) )
}

stats <- cbind(stats,foldChange)

stats$adjustedP <- p.adjust(stats$P_DiseaseState, method = "bonferroni")

write.table(stats, file = "D:/COVID_TISSUE_PROJECT/PlacentaPaper/GeoMX/DEG/LegitModel/glmmSeq_COVID_placenta_vs_CHI_villous_stroma_LEGIT_adding_no_maternal_age.csv", sep = ",", quote = FALSE)

keyvals.colour <- as.character(zeros(nrow(stats)) * NA)

keyvals.colour[stats$foldChange > 0] <- "red"
keyvals.colour[stats$foldChange < 0] <- "deepskyblue"
keyvals.colour[stats$P_DiseaseState > 0.05] <- "grey80"
keyvals.colour[is.na(keyvals.colour)] <- 'grey80'
keyvals.colour[rownames(stats) %in% c("HLA-B")] <- "green"

names(keyvals.colour) <- rownames(stats)


#stats$P_DiseaseState[stats$P_DiseaseState < 0.00001] <- 0.00001

EnhancedVolcano(stats,
                lab = rownames(stats),
                x = 'foldChange',
                y = 'P_DiseaseState',
                cutoffLineType = 'blank',
                #selectLab = c('B2M','LYZ',"HLA-B","MSR1","PPT1","SEC22B","GPX3"),
                cutoffLineWidth = 0.8,
                colCustom = keyvals.colour,
                pointSize = 1.5,
                pCutoff = 0.05,
                FCcutoff = 0.5,
                labFace = 'italic',
                #boxedLabels = TRUE,
                #drawConnectors = TRUE,
                widthConnectors = 1.0,
                colConnectors = 'black',
                title = "",
                subtitle = "",
                caption = "",
                labSize = 4.0,
                colAlpha = 0.9,
                gridlines.major = FALSE,
                border = "full",
                gridlines.minor = FALSE) + theme(legend.position = "none") + theme(
                  axis.title.x = element_text(size = 11),
                  axis.text.x = element_text(size = 10),
                  axis.title.y = element_text(size = 11),
                  axis.text.y = element_text(size = 10)) + xlab("Log2FC")  + ylab(expression(-log[10](P~value))) # + ylim(0, max(abs(stats$foldChange)) + 0.5)

# Plotting


stats <- read.table(file = "D:/COVID_TISSUE_PROJECT/PlacentaPaper/GeoMX/DEG/LegitModel/Keep/glmmSeq_COVID_placenta_vs_COVID_Mother_villous_stroma_LEGIT_adding_no_maternal_age.csv", sep = ",")

colnames(stats) <- stats[1,]
stats <- stats[c(2:nrow(stats)),]

genes <- stats[,1]
stats <- stats[,c(2:ncol(stats))]



genes <- rownames(stats)
stats <- as.data.frame(apply(stats,MARGIN = 2, as.numeric))
rownames(stats) <- genes

keyvals.colour <- as.character(zeros(nrow(stats)) * NA)

keyvals.colour[stats$foldChange > 0] <- "red"
keyvals.colour[stats$foldChange < 0] <- "deepskyblue"
keyvals.colour[stats$P_DiseaseState > 0.05] <- "grey80"
keyvals.colour[is.na(keyvals.colour)] <- 'grey80'
keyvals.colour[rownames(stats) %in% c("HLA-B")] <- "green"

names(keyvals.colour) <- rownames(stats)


#stats$P_DiseaseState[stats$P_DiseaseState < 0.00001] <- 0.00001

EnhancedVolcano(stats,
                lab = rownames(stats),
                x = 'foldChange',
                y = 'P_DiseaseState',
                cutoffLineType = 'blank',
                #selectLab = c('B2M','LYZ',"HLA-B","MSR1","PPT1","SEC22B","GPX3"),
                cutoffLineWidth = 0.8,
                colCustom = keyvals.colour,
                pointSize = 1.5,
                pCutoff = 0.05,
                FCcutoff = 0.5,
                labFace = 'italic',
                #boxedLabels = TRUE,
                #drawConnectors = TRUE,
                widthConnectors = 1.0,
                colConnectors = 'black',
                title = "",
                subtitle = "",
                caption = "",
                labSize = 4.0,
                colAlpha = 0.9,
                gridlines.major = FALSE,
                border = "full",
                gridlines.minor = FALSE) + theme(legend.position = "none") + theme(
                  axis.title.x = element_text(size = 11),
                  axis.text.x = element_text(size = 10),
                  axis.title.y = element_text(size = 11),
                  axis.text.y = element_text(size = 10)) + xlab("Log2FC") # + ylab(expression(-log[10](P~value))) + ylim(0, 2)


# Trophoblast - COVID Placenta vs Normal

diseases <- c("COVID+ Placenta","Normal")

readsTrophoblastModel <- readsTrophoblast[,annoTrophoblast$DiseaseState %in% diseases]
annoTrophoblastModel <- annoTrophoblast[annoTrophoblast$DiseaseState %in% diseases,]

annoTrophoblastModel$VirusPositive2 <- annoTrophoblastModel$`Nucleocapsid IHC`
annoTrophoblastModel$VirusPositive2[annoTrophoblastModel$`Nucleocapsid IHC` %in% c("Strong")] <- TRUE
annoTrophoblastModel$VirusPositive2[annoTrophoblastModel$`Nucleocapsid IHC` %in% c("Negative","Weak")] <- FALSE

disp <- apply(readsTrophoblastModel, 1, function(x){
  (var(x, na.rm=TRUE)-mean(x, na.rm=TRUE))/(mean(x, na.rm=TRUE)**2)
})

head(disp)

colnames(annoTrophoblastModel)[c(19:20)] <- c("Gestational_Age","Maternal_Age")

statResults <- glmmSeq::glmmSeq(~ DiseaseState + Gestational_Age + (1|Number),
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
  foldChange[i,1] = log2(mean(as.numeric(readsTrophoblastModel[rownames(stats)[i],annoTrophoblastModel$DiseaseState == diseases[1]]))/mean(as.numeric(readsTrophoblastModel[rownames(stats)[i],annoTrophoblastModel$DiseaseState == diseases[2]])) )
}

stats <- cbind(stats,foldChange)

stats$adjustedP <- p.adjust(stats$P_DiseaseState, method = "bonferroni")

write.table(stats, file = "D:/COVID_TISSUE_PROJECT/PlacentaPaper/GeoMX/DEG/LegitModel/glmmSeq_COVID_placenta_vs_Normal_trophoblast_LEGIT_adding_no_maternal_age.csv", sep = ",", quote = FALSE)

keyvals.colour <- as.character(zeros(nrow(stats)) * NA)

keyvals.colour[stats$foldChange > 0] <- "red"
keyvals.colour[stats$foldChange < 0] <- "deepskyblue"
keyvals.colour[stats$P_DiseaseState > 0.05] <- "grey80"
keyvals.colour[is.na(keyvals.colour)] <- 'grey80'
keyvals.colour[rownames(stats) %in% c("HLA-B")] <- "green"

names(keyvals.colour) <- rownames(stats)


#stats$P_DiseaseState[stats$P_DiseaseState < 0.00001] <- 0.00001

EnhancedVolcano(stats,
                lab = rownames(stats),
                x = 'foldChange',
                y = 'P_DiseaseState',
                cutoffLineType = 'blank',
                #selectLab = c('B2M','LYZ',"HLA-B","MSR1","PPT1","SEC22B","GPX3"),
                cutoffLineWidth = 0.8,
                colCustom = keyvals.colour,
                pointSize = 1.5,
                pCutoff = 0.05,
                FCcutoff = 0.5,
                labFace = 'italic',
                #boxedLabels = TRUE,
                #drawConnectors = TRUE,
                widthConnectors = 1.0,
                colConnectors = 'black',
                title = "",
                subtitle = "",
                caption = "",
                labSize = 4.0,
                colAlpha = 0.9,
                gridlines.major = FALSE,
                border = "full",
                gridlines.minor = FALSE) + theme(legend.position = "none") + theme(
                  axis.title.x = element_text(size = 11),
                  axis.text.x = element_text(size = 10),
                  axis.title.y = element_text(size = 11),
                  axis.text.y = element_text(size = 10)) + xlab("Log2FC")  + ylab(expression(-log[10](P~value))) # + ylim(0, max(abs(stats$foldChange)) + 0.5)

# Decidua - COVID Placenta vs Normal

diseases <- c("COVID+ Placenta","Normal")

readsTrophoblastModel <- readsTrophoblast[,annoTrophoblast$DiseaseState %in% diseases]
annoTrophoblastModel <- annoTrophoblast[annoTrophoblast$DiseaseState %in% diseases,]

annoTrophoblastModel$VirusPositive2 <- annoTrophoblastModel$`Nucleocapsid IHC`
annoTrophoblastModel$VirusPositive2[annoTrophoblastModel$`Nucleocapsid IHC` %in% c("Strong")] <- TRUE
annoTrophoblastModel$VirusPositive2[annoTrophoblastModel$`Nucleocapsid IHC` %in% c("Negative","Weak")] <- FALSE

disp <- apply(readsTrophoblastModel, 1, function(x){
  (var(x, na.rm=TRUE)-mean(x, na.rm=TRUE))/(mean(x, na.rm=TRUE)**2)
})

head(disp)

colnames(annoTrophoblastModel)[c(19:20)] <- c("Gestational_Age","Maternal_Age")

statResults <- glmmSeq::glmmSeq(~ DiseaseState + Gestational_Age + (1|Number),
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
  foldChange[i,1] = log2(mean(as.numeric(readsTrophoblastModel[rownames(stats)[i],annoTrophoblastModel$DiseaseState == diseases[1]]))/mean(as.numeric(readsTrophoblastModel[rownames(stats)[i],annoTrophoblastModel$DiseaseState == diseases[2]])) )
}

stats <- cbind(stats,foldChange)

stats$adjustedP <- p.adjust(stats$P_DiseaseState, method = "bonferroni")

write.table(stats, file = "D:/COVID_TISSUE_PROJECT/PlacentaPaper/GeoMX/DEG/LegitModel/glmmSeq_COVID_placenta_vs_Normal_decidua_LEGIT_adding_no_maternal_age.csv", sep = ",", quote = FALSE)

keyvals.colour <- as.character(zeros(nrow(stats)) * NA)

keyvals.colour[stats$foldChange > 0] <- "red"
keyvals.colour[stats$foldChange < 0] <- "deepskyblue"
keyvals.colour[stats$P_DiseaseState > 0.05] <- "grey80"
keyvals.colour[is.na(keyvals.colour)] <- 'grey80'
keyvals.colour[rownames(stats) %in% c("HLA-B")] <- "green"

names(keyvals.colour) <- rownames(stats)


#stats$P_DiseaseState[stats$P_DiseaseState < 0.00001] <- 0.00001

EnhancedVolcano(stats,
                lab = rownames(stats),
                x = 'foldChange',
                y = 'P_DiseaseState',
                cutoffLineType = 'blank',
                #selectLab = c('B2M','LYZ',"HLA-B","MSR1","PPT1","SEC22B","GPX3"),
                cutoffLineWidth = 0.8,
                colCustom = keyvals.colour,
                pointSize = 1.5,
                pCutoff = 0.05,
                FCcutoff = 0.5,
                labFace = 'italic',
                #boxedLabels = TRUE,
                #drawConnectors = TRUE,
                widthConnectors = 1.0,
                colConnectors = 'black',
                title = "",
                subtitle = "",
                caption = "",
                labSize = 4.0,
                colAlpha = 0.9,
                gridlines.major = FALSE,
                border = "full",
                gridlines.minor = FALSE) + theme(legend.position = "none") + theme(
                  axis.title.x = element_text(size = 11),
                  axis.text.x = element_text(size = 10),
                  axis.title.y = element_text(size = 11),
                  axis.text.y = element_text(size = 10)) + xlab("Log2FC")  + ylab(expression(-log[10](P~value))) # + ylim(0, max(abs(stats$foldChange)) + 0.5)

# Villous stroma - COVID Placenta vs Normal

diseases <- c("COVID+ Placenta","Normal")

readsTrophoblastModel <- readsTrophoblast[,annoTrophoblast$DiseaseState %in% diseases]
annoTrophoblastModel <- annoTrophoblast[annoTrophoblast$DiseaseState %in% diseases,]

annoTrophoblastModel$VirusPositive2 <- annoTrophoblastModel$`Nucleocapsid IHC`
annoTrophoblastModel$VirusPositive2[annoTrophoblastModel$`Nucleocapsid IHC` %in% c("Strong")] <- TRUE
annoTrophoblastModel$VirusPositive2[annoTrophoblastModel$`Nucleocapsid IHC` %in% c("Negative","Weak")] <- FALSE

disp <- apply(readsTrophoblastModel, 1, function(x){
  (var(x, na.rm=TRUE)-mean(x, na.rm=TRUE))/(mean(x, na.rm=TRUE)**2)
})

head(disp)

colnames(annoTrophoblastModel)[c(19:20)] <- c("Gestational_Age","Maternal_Age")

statResults <- glmmSeq::glmmSeq(~ DiseaseState + Gestational_Age + (1|Number),
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
  foldChange[i,1] = log2(mean(as.numeric(readsTrophoblastModel[rownames(stats)[i],annoTrophoblastModel$DiseaseState == diseases[1]]))/mean(as.numeric(readsTrophoblastModel[rownames(stats)[i],annoTrophoblastModel$DiseaseState == diseases[2]])) )
}

stats <- cbind(stats,foldChange)

stats$adjustedP <- p.adjust(stats$P_DiseaseState, method = "bonferroni")

write.table(stats, file = "D:/COVID_TISSUE_PROJECT/PlacentaPaper/GeoMX/DEG/LegitModel/glmmSeq_COVID_placenta_vs_Normal_villous_stroma_LEGIT_adding_no_maternal_age.csv", sep = ",", quote = FALSE)

keyvals.colour <- as.character(zeros(nrow(stats)) * NA)

keyvals.colour[stats$foldChange > 0] <- "red"
keyvals.colour[stats$foldChange < 0] <- "deepskyblue"
keyvals.colour[stats$P_DiseaseState > 0.05] <- "grey80"
keyvals.colour[is.na(keyvals.colour)] <- 'grey80'
keyvals.colour[rownames(stats) %in% c("HLA-B")] <- "green"

names(keyvals.colour) <- rownames(stats)


#stats$P_DiseaseState[stats$P_DiseaseState < 0.00001] <- 0.00001

EnhancedVolcano(stats,
                lab = rownames(stats),
                x = 'foldChange',
                y = 'P_DiseaseState',
                cutoffLineType = 'blank',
                #selectLab = c('B2M','LYZ',"HLA-B","MSR1","PPT1","SEC22B","GPX3"),
                cutoffLineWidth = 0.8,
                colCustom = keyvals.colour,
                pointSize = 1.5,
                pCutoff = 0.05,
                FCcutoff = 0.5,
                labFace = 'italic',
                #boxedLabels = TRUE,
                #drawConnectors = TRUE,
                widthConnectors = 1.0,
                colConnectors = 'black',
                title = "",
                subtitle = "",
                caption = "",
                labSize = 4.0,
                colAlpha = 0.9,
                gridlines.major = FALSE,
                border = "full",
                gridlines.minor = FALSE) + theme(legend.position = "none") + theme(
                  axis.title.x = element_text(size = 11),
                  axis.text.x = element_text(size = 10),
                  axis.title.y = element_text(size = 11),
                  axis.text.y = element_text(size = 10)) + xlab("Log2FC")  + ylab(expression(-log[10](P~value))) # + ylim(0, max(abs(stats$foldChange)) + 0.5)

# Decidua - COVID_Placenta + COVID_Mother vs Normal


annoTrophoblast$DiseaseState2 <- annoTrophoblast$DiseaseState
annoTrophoblast$DiseaseState2[annoTrophoblast$DiseaseState2 %in% c("COVID+ Placenta","COVID+ Mother")] <- "COVID"

diseases <- c("COVID","Normal")

readsTrophoblastModel <- readsTrophoblast[,annoTrophoblast$DiseaseState2 %in% diseases]
annoTrophoblastModel <- annoTrophoblast[annoTrophoblast$DiseaseState2 %in% diseases,]

annoTrophoblastModel$VirusPositive2 <- annoTrophoblastModel$`Nucleocapsid IHC`
annoTrophoblastModel$VirusPositive2[annoTrophoblastModel$`Nucleocapsid IHC` %in% c("Strong")] <- TRUE
annoTrophoblastModel$VirusPositive2[annoTrophoblastModel$`Nucleocapsid IHC` %in% c("Negative","Weak")] <- FALSE

disp <- apply(readsTrophoblastModel, 1, function(x){
  (var(x, na.rm=TRUE)-mean(x, na.rm=TRUE))/(mean(x, na.rm=TRUE)**2)
})

head(disp)

colnames(annoTrophoblastModel)[c(19:20)] <- c("Gestational_Age","Maternal_Age")

statResults <- glmmSeq::glmmSeq(~ DiseaseState + Gestational_Age + (1|Number),
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
  foldChange[i,1] = log2(mean(as.numeric(readsTrophoblastModel[rownames(stats)[i],annoTrophoblastModel$DiseaseState2 == diseases[1]]))/mean(as.numeric(readsTrophoblastModel[rownames(stats)[i],annoTrophoblastModel$DiseaseState2 == diseases[2]])) )
}

stats <- cbind(stats,foldChange)

stats$adjustedP <- p.adjust(stats$P_DiseaseState, method = "bonferroni")

write.table(stats, file = "D:/COVID_TISSUE_PROJECT/PlacentaPaper/GeoMX/DEG/LegitModel/glmmSeq_COVID_TOTAL_vs_Normal_decidua_LEGIT_adding_no_maternal_age_proper.csv", sep = ",", quote = FALSE)

keyvals.colour <- as.character(zeros(nrow(stats)) * NA)

keyvals.colour[stats$foldChange > 0] <- "red"
keyvals.colour[stats$foldChange < 0] <- "deepskyblue"
keyvals.colour[stats$P_DiseaseState > 0.05] <- "grey80"
keyvals.colour[is.na(keyvals.colour)] <- 'grey80'
keyvals.colour[rownames(stats) %in% c("HLA-B")] <- "green"

names(keyvals.colour) <- rownames(stats)


#stats$P_DiseaseState[stats$P_DiseaseState < 0.00001] <- 0.00001

EnhancedVolcano(stats,
                lab = rownames(stats),
                x = 'foldChange',
                y = 'P_DiseaseState',
                cutoffLineType = 'blank',
                #selectLab = c('B2M','LYZ',"HLA-B","MSR1","PPT1","SEC22B","GPX3"),
                cutoffLineWidth = 0.8,
                colCustom = keyvals.colour,
                pointSize = 1.5,
                pCutoff = 0.05,
                FCcutoff = 0.5,
                labFace = 'italic',
                #boxedLabels = TRUE,
                #drawConnectors = TRUE,
                widthConnectors = 1.0,
                colConnectors = 'black',
                title = "",
                subtitle = "",
                caption = "",
                labSize = 4.0,
                colAlpha = 0.9,
                gridlines.major = FALSE,
                border = "full",
                gridlines.minor = FALSE) + theme(legend.position = "none") + theme(
                  axis.title.x = element_text(size = 11),
                  axis.text.x = element_text(size = 10),
                  axis.title.y = element_text(size = 11),
                  axis.text.y = element_text(size = 10)) + xlab("Log2FC")  + ylab(expression(-log[10](P~value))) # + ylim(0, max(abs(stats$foldChange)) + 0.5)

###
# Virus positive trophoblast region vs. virus negative villous stroma region

diseases <- c("COVID+ Placenta") #,"COVID+ Mother")

readsTrophoblastModel <- readsTrophoblast[,annoTrophoblast$DiseaseState %in% diseases]
annoTrophoblastModel <- annoTrophoblast[annoTrophoblast$DiseaseState %in% diseases,]

annoTrophoblastModel$VirusPositive2 <- annoTrophoblastModel$`Nucleocapsid IHC`
annoTrophoblastModel$VirusPositive2[annoTrophoblastModel$`Nucleocapsid IHC` %in% c("Strong")] <- TRUE
annoTrophoblastModel$VirusPositive2[annoTrophoblastModel$`Nucleocapsid IHC` %in% c("Negative","Weak","NA")] <- FALSE

disp <- apply(readsTrophoblastModel, 1, function(x){
  (var(x, na.rm=TRUE)-mean(x, na.rm=TRUE))/(mean(x, na.rm=TRUE)**2)
})

head(disp)

colnames(annoTrophoblastModel)[c(19:20)] <- c("Gestational_Age","Maternal_Age")

statResults <- glmmSeq::glmmSeq(~ VirusPositive2 + (1|Number),
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
  foldChange[i,1] = log2(mean(as.numeric(readsTrophoblastModel[rownames(stats)[i],annoTrophoblastModel$VirusPositive2 == TRUE]))/mean(as.numeric(readsTrophoblastModel[rownames(stats)[i],annoTrophoblastModel$VirusPositive2 == FALSE])) )
}

stats <- cbind(stats,foldChange)

stats$adjustedP <- p.adjust(stats$P_VirusPositive2, method = "bonferroni")

write.table(stats, file = "D:/COVID_TISSUE_PROJECT/PlacentaPaper/GeoMX/DEG/LegitModel/glmmSeq_COVID_placenta_virus_positive_vs_virus_negative_VS_LEGIT_adding_no_fixed_effects_age.csv", sep = ",", quote = FALSE)

stats <- read.table("D:/COVID_TISSUE_PROJECT/PlacentaPaper/GeoMX/DEG/LegitModel/Keep/glmmSeq_COVID_placenta_virus_positive_vs_virus_negative_VS_LEGIT_adding_no_fixed_effects_age.csv", sep = ",")

colnames(stats) <- stats[1,]
rownames(stats) <- stats[,1]
stats <- stats[c(2:nrow(stats)), c(2:ncol(stats))]

stats$foldChange <- as.numeric(stats$foldChange)
stats$P_VirusPositive2 <- as.numeric(stats$P_VirusPositive2)
#stats$adjustedP <- as.numeric(stats$adjustedP)

keyvals.colour <- as.character(zeros(nrow(stats)) * NA)

keyvals.colour[stats$foldChange > 0.5] <- "red"
keyvals.colour[stats$foldChange < -0.5] <- "deepskyblue"
keyvals.colour[stats$P_VirusPositive2 > 0.05] <- "grey80"
keyvals.colour[is.na(keyvals.colour)] <- 'grey80'

names(keyvals.colour) <- rownames(stats)

EnhancedVolcano(stats,
                lab = rownames(stats),
                x = 'foldChange',
                y = 'P_VirusPositive2',
                cutoffLineType = 'blank',
                selectLab = c('VSIG8',"CCF1","PFN1","RAB1A","FUS","CCR1","CAMK2N2"),
                cutoffLineWidth = 0.8,
                colCustom = keyvals.colour,
                pointSize = 2.0,
                pCutoff = 0.05,
                FCcutoff = 0.75,
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
                border = "full",
                gridlines.major = FALSE,
                gridlines.minor = FALSE) + theme(legend.position = "none") + theme(
                  axis.title.x = element_text(size = 12),
                  axis.text.x = element_text(size = 12),
                  axis.title.y = element_text(size = 12),
                  axis.text.y = element_text(size = 12)) + xlab("Log2FC") + ylim(0, 25)

#### New split in diseases

diseases <- c("COVID+ Placenta","COVID+ Mother CHI","CHI")

readsTrophoblastModel <- readsTrophoblast[,annoTrophoblast$DiseaseState2 %in% diseases]
annoTrophoblastModel <- annoTrophoblast[annoTrophoblast$DiseaseState2 %in% diseases,]

annoTrophoblastModel$DiseaseState2[annoTrophoblastModel$DiseaseState2 %in% c("COVID+ Placenta","COVID+ Mother CHI")] <- "COVID CHI"
diseases <- c("COVID CHI","CHI")

annoTrophoblastModel$VirusPositive2 <- annoTrophoblastModel$`Nucleocapsid IHC`
annoTrophoblastModel$VirusPositive2[annoTrophoblastModel$`Nucleocapsid IHC` %in% c("Strong")] <- TRUE
annoTrophoblastModel$VirusPositive2[annoTrophoblastModel$`Nucleocapsid IHC` %in% c("Negative","Weak")] <- FALSE

disp <- apply(readsTrophoblastModel, 1, function(x){
  (var(x, na.rm=TRUE)-mean(x, na.rm=TRUE))/(mean(x, na.rm=TRUE)**2)
})

head(disp)

colnames(annoTrophoblastModel)[c(20:21)] <- c("Gestational_Age","Maternal_Age")

statResults <- glmmSeq::glmmSeq(~ DiseaseState2 + Gestational_Age + (1|Number),
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
  foldChange[i,1] = log2(mean(as.numeric(readsTrophoblastModel[rownames(stats)[i],annoTrophoblastModel$DiseaseState2 == diseases[1]]))/mean(as.numeric(readsTrophoblastModel[rownames(stats)[i],annoTrophoblastModel$DiseaseState2 == diseases[2]])) )
}

stats <- cbind(stats,foldChange)

stats$adjustedP <- p.adjust(stats$P_DiseaseState, method = "bonferroni")

#stats <- read.table("D:/COVID_TISSUE_PROJECT/PlacentaPaper/GeoMX/DEG/LegitModel/glmmSeq_COVID_placenta_vs_COVID_mother_CV_trophoblast_LEGIT_adding_no_maternal_age_new_split_v2.csv", sep = ",")
#colnames(stats) <- stats[1,]
#rownames(stats) <- stats[,1]
#stats <- stats[c(2:nrow(stats)), c(2:ncol(stats))]
#stats$foldChange <- as.numeric(stats$foldChange)
#stats$P_VirusPositive2 <- as.numeric(stats$P_DiseaseState)


write.table(stats, file = "D:/COVID_TISSUE_PROJECT/PlacentaPaper/GeoMX/DEG/LegitModel/glmmSeq_COVID_all_CHI_vs_Control_CHI_macrophage_LEGIT_adding_no_maternal_age_new_split_v2.csv", sep = ",", quote = FALSE)


stats <- read.table(file = "D:/COVID_TISSUE_PROJECT/PlacentaPaper/GeoMX/DEG/LegitModel/glmmSeq_COVID_placenta_vs_COVID_mother_CHI_macrophage_LEGIT_adding_no_maternal_age_new_split_v2.csv", sep = ",",
                    header = TRUE, row.names = 1)


stats2 <- stats
stats <- stats[!is.na(stats$P_DiseaseState2),]
notSignif <- stats$P_DiseaseState2 > 0.05
perc.70 <- round(sum(notSignif) * 0.85)
button.5 <- which(notSignif == TRUE)
sampled.70 <- sample(button.5, perc.70)
stats <- stats[-sampled.70, ]


keyvals.colour <- as.character(zeros(nrow(stats)) * NA)

keyvals.colour[stats$foldChange > 0] <- "red"
keyvals.colour[stats$foldChange < 0] <- "deepskyblue"
keyvals.colour[stats$P_DiseaseState > 0.05] <- "grey80"
keyvals.colour[is.na(keyvals.colour)] <- 'grey80'
keyvals.colour[rownames(stats) %in% c("HLA-B")] <- "green"

names(keyvals.colour) <- rownames(stats)


#stats$P_DiseaseState[stats$P_DiseaseState < 0.00001] <- 0.00001

EnhancedVolcano(stats,
                lab = rownames(stats),
                x = 'foldChange',
                y = 'P_DiseaseState2',
                cutoffLineType = 'blank',
                selectLab = c("FN1","B2M","FBN2","IFI30","SEC22B","SND1",
                              "CLIC3","FBLN1","IFNE","HSPA9","HBA2"),
                cutoffLineWidth = 0.8,
                colCustom = keyvals.colour,
                pointSize = 1.5,
                pCutoff = 0.05,
                FCcutoff = 0.5,
                labFace = 'italic',
                #boxedLabels = TRUE,
                #drawConnectors = TRUE,
                widthConnectors = 1.0,
                colConnectors = 'black',
                title = "",
                subtitle = "",
                caption = "",
                labSize = 4.0,
                colAlpha = 0.9,
                gridlines.major = FALSE,
                border = "full",
                gridlines.minor = FALSE) + theme_bw() + theme(
                  legend.position = "none",
                  panel.grid.major = element_blank(),
                  panel.grid.minor = element_blank(),
                  axis.title.x = element_text(size = 11),
                  axis.text.x = element_text(size = 10),
                  axis.title.y = element_text(size = 11),
                  axis.text.y = element_text(size = 10)) + xlab("Log2FC")  + ylab(expression(-log[10](P))) # + ylim(0, max(abs(stats$foldChange)) + 0.5)

### Use top and bottom genes from M+P+ CHI vs M+P- CHI DEG to input into heatmap

stats <- read.table(file = "D:/COVID_TISSUE_PROJECT/PlacentaPaper/GeoMX/DEG/LegitModel/glmmSeq_COVID_placenta_vs_COVID_mother_CHI_macrophage_LEGIT_adding_no_maternal_age_new_split_v2.csv", sep = ",",
                    header = TRUE, row.names = 1)

numberOf <- 20

upReg <- stats[stats$foldChange > 0,]
upReg$P_DiseaseState <- as.numeric(upReg$P_DiseaseState)
upReg <- upReg[order(upReg$P_DiseaseState),]  # decreasing = TRUE

top30 <- rownames(upReg)[c(1:numberOf)]

downReg <- stats[stats$foldChange < 0,]
downReg$P_DiseaseState <- as.numeric(downReg$P_DiseaseState)
downReg <- downReg[order(downReg$P_DiseaseState),]  # decreasing = TRUE

bottom30 <- rownames(downReg)[c(1:numberOf)]

topBottom <- c(top30,bottom30)


# heatmap

genes <- topBottom

ROI <- c("Macrophage")

readsHM <- reads[,annotationFile[,"Compartment"] %in% ROI]
annoHM <- annotationFile[annotationFile[,"Compartment"] %in% ROI,]

logTransformedData <- log2(readsHM+1)
scaledData <- as.data.frame(t(apply(na.omit(logTransformedData[genes,]), MARGIN = 1, scale)))

#scaledData <- na.omit(logTransformedData[genes,])

upperLimit <- 3
lowerLimit <- -3
scaledData[scaledData > upperLimit] <- upperLimit
scaledData[scaledData < lowerLimit] <- lowerLimit

annoHM$DiseaseState2 <- factor(annoHM$DiseaseState2, levels = c("COVID+ Placenta","COVID+ Mother CHI", "CHI", "COVID+ Mother CV"))

ha = HeatmapAnnotation(
  Patient = as.factor(annoHM$Number)
)
split = annoHM$DiseaseState2

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

#### DEG Macrophage Region clusters

diseases <- c("COVID+ Placenta","COVID+ Mother CHI")

readsTrophoblastModel <- readsTrophoblast[,annoTrophoblast$DiseaseState2 %in% diseases]
annoTrophoblastModel <- annoTrophoblast[annoTrophoblast$DiseaseState2 %in% diseases,]

annoTrophoblastModel$TrophoblastClustersCOVIDCHI[annoTrophoblastModel$TrophoblastClustersCOVIDCHI %in% c("0","1")] <- "Rest"

#readsTrophoblastModel <- readsTrophoblastModel[,annoTrophoblastModel$MacrophageClustersCOVIDCHI %in% c("0","1")]
#annoTrophoblastModel <- annoTrophoblastModel[annoTrophoblastModel$MacrophageClustersCOVIDCHI %in% c("0","1"),]

annoTrophoblastModel$VirusPositive2 <- annoTrophoblastModel$`Nucleocapsid IHC`
annoTrophoblastModel$VirusPositive2[annoTrophoblastModel$`Nucleocapsid IHC` %in% c("Strong")] <- TRUE
annoTrophoblastModel$VirusPositive2[annoTrophoblastModel$`Nucleocapsid IHC` %in% c("Negative","Weak")] <- FALSE

disp <- apply(readsTrophoblastModel, 1, function(x){
  (var(x, na.rm=TRUE)-mean(x, na.rm=TRUE))/(mean(x, na.rm=TRUE)**2)
})

head(disp)

colnames(annoTrophoblastModel)[c(20:21)] <- c("Gestational_Age","Maternal_Age")

statResults <- glmmSeq::glmmSeq(~ TrophoblastClustersCOVIDCHI + (1|Number),
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
  foldChange[i,1] = log2(mean(as.numeric(readsTrophoblastModel[rownames(stats)[i],annoTrophoblastModel$TrophoblastClustersCOVIDCHI == "2"]))/mean(as.numeric(readsTrophoblastModel[rownames(stats)[i],annoTrophoblastModel$TrophoblastClustersCOVIDCHI == "Rest"])) )
}

stats <- cbind(stats,foldChange)

stats$adjustedP <- p.adjust(stats$P_TrophoblastClustersCOVIDCHI, method = "bonferroni")

#stats <- read.table("D:/COVID_TISSUE_PROJECT/PlacentaPaper/GeoMX/DEG/LegitModel/glmmSeq_MacrophageClusters_cluster0_vs_cluster1_LEGIT_adding_no_maternal_age_new_split.csv", sep = ",")
#colnames(stats) <- stats[1,]
#rownames(stats) <- stats[,1]
#stats <- stats[c(2:nrow(stats)), c(2:ncol(stats))]
#stats$foldChange <- as.numeric(stats$foldChange)
#stats$P_VirusPositive2 <- as.numeric(stats$P_DiseaseState)


write.table(stats, file = "D:/COVID_TISSUE_PROJECT/PlacentaPaper/GeoMX/DEG/LegitModel/glmmSeq_trophoblastClusters_COVID_CHI_cluster2_vs_rest_LEGIT_clustering_noGestationalAge.csv", sep = ",", quote = FALSE)

keyvals.colour <- as.character(zeros(nrow(stats)) * NA)

keyvals.colour[stats$foldChange > 0] <- "red"
keyvals.colour[stats$foldChange < 0] <- "deepskyblue"
keyvals.colour[stats$P_MacrophageClustersCOVIDCHI > 0.05] <- "grey80"
keyvals.colour[is.na(keyvals.colour)] <- 'grey80'
keyvals.colour[rownames(stats) %in% c("HLA-B")] <- "green"

names(keyvals.colour) <- rownames(stats)


#stats$P_DiseaseState[stats$P_DiseaseState < 0.00001] <- 0.00001

EnhancedVolcano(stats,
                lab = rownames(stats),
                x = 'foldChange',
                y = 'P_TrophoblastClustersCOVIDCHI',
                cutoffLineType = 'blank',
                cutoffLineWidth = 0.8,
                colCustom = keyvals.colour,
                pointSize = 1.5,
                pCutoff = 0.05,
                FCcutoff = 0.5,
                labFace = 'italic',
                #boxedLabels = TRUE,
                #drawConnectors = TRUE,
                widthConnectors = 1.0,
                colConnectors = 'black',
                title = "",
                subtitle = "",
                caption = "",
                labSize = 4.0,
                colAlpha = 0.9,
                gridlines.major = FALSE,
                border = "full",
                gridlines.minor = FALSE) + theme(legend.position = "none") + theme(
                  axis.title.x = element_text(size = 11),
                  axis.text.x = element_text(size = 10),
                  axis.title.y = element_text(size = 11),
                  axis.text.y = element_text(size = 10)) + xlab("Log2FC")  + ylab(expression(-log[10](P~value))) # + ylim(0, max(abs(stats$foldChange)) + 0.5)

#### DEG Macrophage Virus Hi vs. Virus Lo

diseases <- c("COVID+ Placenta")

readsTrophoblastModel <- readsTrophoblast[,annoTrophoblast$DiseaseState2 %in% diseases]
annoTrophoblastModel <- annoTrophoblast[annoTrophoblast$DiseaseState2 %in% diseases,]

#annoTrophoblastModel$TrophoblastClustersCOVIDCHI[annoTrophoblastModel$TrophoblastClustersCOVIDCHI %in% c("0","1")] <- "Rest"

#readsTrophoblastModel <- readsTrophoblastModel[,annoTrophoblastModel$MacrophageClustersCOVIDCHI %in% c("0","1")]
#annoTrophoblastModel <- annoTrophoblastModel[annoTrophoblastModel$MacrophageClustersCOVIDCHI %in% c("0","1"),]

annoTrophoblastModel$VirusPositive2 <- annoTrophoblastModel$`Nucleocapsid IHC`
annoTrophoblastModel$VirusPositive2[annoTrophoblastModel$`Nucleocapsid IHC` %in% c("Strong")] <- TRUE
annoTrophoblastModel$VirusPositive2[annoTrophoblastModel$`Nucleocapsid IHC` %in% c("Negative","Weak")] <- FALSE

disp <- apply(readsTrophoblastModel, 1, function(x){
  (var(x, na.rm=TRUE)-mean(x, na.rm=TRUE))/(mean(x, na.rm=TRUE)**2)
})

head(disp)

colnames(annoTrophoblastModel)[c(20:21)] <- c("Gestational_Age","Maternal_Age")

statResults <- glmmSeq::glmmSeq(~ VirusPositive2 + (1|Number),
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
  foldChange[i,1] = log2(mean(as.numeric(readsTrophoblastModel[rownames(stats)[i],annoTrophoblastModel$VirusPositive2 == TRUE]))/mean(as.numeric(readsTrophoblastModel[rownames(stats)[i],annoTrophoblastModel$VirusPositive2 == FALSE])) )
}

stats <- cbind(stats,foldChange)

stats$adjustedP <- p.adjust(stats$P_VirusPositive2, method = "bonferroni")

write.table(stats, file = "D:/COVID_TISSUE_PROJECT/PlacentaPaper/GeoMX/DEG/LegitModel/glmmSeq_macrophage_virus_high_vs_virus_low.csv", sep = ",", quote = FALSE)

# reduce number of dots that are not significant

stats2 <- stats
stats <- stats[!is.na(stats$P_VirusPositive2),]
notSignif <- stats$P_VirusPositive2 > 0.05
perc.70 <- round(sum(notSignif) * 0.8)
button.5 <- which(notSignif == TRUE)
sampled.70 <- sample(button.5, perc.70)
stats <- stats[-sampled.70, ]


keyvals.colour <- as.character(zeros(nrow(stats)) * NA)

keyvals.colour[stats$foldChange > 0] <- "red"
keyvals.colour[stats$foldChange < 0] <- "deepskyblue"
keyvals.colour[stats$P_MacrophageClustersCOVIDCHI > 0.05] <- "grey80"
keyvals.colour[is.na(keyvals.colour)] <- 'grey80'
keyvals.colour[rownames(stats) %in% c("HLA-B")] <- "green"

names(keyvals.colour) <- rownames(stats)

#stats$P_DiseaseState[stats$P_DiseaseState < 0.00001] <- 0.00001

EnhancedVolcano(stats,
                lab = rownames(stats),
                x = 'foldChange',
                y = 'P_VirusPositive2',
                cutoffLineType = 'blank',
                cutoffLineWidth = 0.8,
                colCustom = keyvals.colour,
                selectLab = c("C1QC","VAMP8","NMI","INHBB","IL17D","HSPA5","IGF2","SAT1","TGM2","HLA-A"),
                pointSize = 1.5,
                pCutoff = 0.05,
                FCcutoff = 0.5,
                labFace = 'italic',
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
                border = "full",
                gridlines.minor = FALSE) + theme(legend.position = "none") + theme(
                  axis.title.x = element_text(size = 11),
                  axis.text.x = element_text(size = 10),
                  axis.title.y = element_text(size = 11),
                  axis.text.y = element_text(size = 10)) + xlab("Log2FC")  + ylab(expression(-log[10](P~value))) # + ylim(0, max(abs(stats$foldChange)) + 0.5)

## Plotting imported tables

stats <- read.table("D:/COVID_TISSUE_PROJECT/PlacentaPaper/GeoMX/DEG/LegitModel/Macrophage_Clustering_COVID_CHI/glmmSeq_MacrophageClusters_COVID_CHI_cluster2_vs_rest_LEGIT_second_clustering.csv", sep = ",")
colnames(stats) <- stats[1,]
rownames(stats) <- stats[,1]
stats <- stats[c(2:nrow(stats)), c(2:ncol(stats))]
stats$foldChange <- as.numeric(stats$foldChange)
stats$P_MacrophageClustersCOVIDCHI <- as.numeric(stats$P_MacrophageClustersCOVIDCHI)

# Reduce number of non-significant dots for easier visualisation

stats2 <- stats
stats <- stats[!is.na(stats$P_MacrophageClustersCOVIDCHI),]
notSignif <- stats$P_MacrophageClustersCOVIDCHI > 0.05
perc.70 <- round(sum(notSignif) * 0.8)
button.5 <- which(notSignif == TRUE)
sampled.70 <- sample(button.5, perc.70)
stats <- stats[-sampled.70, ]

keyvals.colour <- as.character(zeros(nrow(stats)) * NA)

keyvals.colour[stats$foldChange > 0] <- "red"
keyvals.colour[stats$foldChange < 0] <- "deepskyblue"
keyvals.colour[stats$P_MacrophageClustersCOVIDCHI > 0.05] <- "grey80"
keyvals.colour[is.na(keyvals.colour)] <- 'grey80'
#keyvals.colour[rownames(stats) %in% c("HLA-B")] <- "green"

names(keyvals.colour) <- rownames(stats)

EnhancedVolcano(stats,
                lab = rownames(stats),
                x = 'foldChange',
                y = 'P_MacrophageClustersCOVIDCHI',
                cutoffLineType = 'blank',
                cutoffLineWidth = 0.8,
                colCustom = keyvals.colour,
                #selectLab = c("C1QC","VAMP8","NMI","INHBB","IL17D","HSPA5","IGF2","SAT1","TGM2","HLA-A"),
                pointSize = 2,
                pCutoff = 0.05,
                FCcutoff = 0.5,
                labFace = 'italic',
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
                border = "full",
                gridlines.minor = FALSE) + theme(legend.position = "none") + theme(
                  axis.title.x = element_text(size = 11),
                  axis.text.x = element_text(size = 10),
                  axis.title.y = element_text(size = 11),
                  axis.text.y = element_text(size = 10)) + xlab("Log2FC")  + ylab(expression(-log[10](P~value))) # + ylim(0, max(abs(stats$foldChange)) + 0.5)




#####
