library(readxl)
library(reshape)

clinicalInfo <- as.data.frame(readxl::read_xlsx("D:/COVID_TISSUE_PROJECT/PlacentaPaper/clinicalDataForTimelinePlot.xlsx",
                                                sheet = "Sheet2"))

DF.SYMBOL <- as.data.frame(readxl::read_xlsx("D:/COVID_TISSUE_PROJECT/PlacentaPaper/clinicalDataForTimelinePlot.xlsx",
                                             sheet = "Symbols"))

rownames(clinicalInfo) <- clinicalInfo[,1]

clinicalInfo <- clinicalInfo[,c(2:ncol(clinicalInfo))]

clinicalInfoTransformed <- melt(cbind(clinicalInfo, Phase = rownames(clinicalInfo)), id.vars = c('Phase'))

clinicalInfoTransformed$Phase <- factor(clinicalInfoTransformed$Phase,levels=c("Symptoms","ICU","Hospital"))

clinicalInfoTransformed$variable <- factor(clinicalInfoTransformed$variable, levels=c(paste0("Patient ",c(2:10)),paste0("Patient ",c(1,11:15))))

ggplot() +
  geom_col(data = clinicalInfoTransformed, aes(x = variable, y = ifelse(Phase %in% "Symptoms", -value, value),
                                        fill = Phase),colour="black", width = 0.8) +
  geom_point(data = DF.SYMBOL, aes(x = Patient, y = Time, shape = Event), size = 2)  +
  coord_flip() + ylab("Days") + xlab("") +
  scale_y_continuous(limits=c(-13, 15), breaks = c(-13, -10,-5,0,5,10,15)) + theme_classic() +
  theme(axis.text.x = element_text(color = "black", size = 10, face = "plain"),
        axis.text.y = element_text(color = "black", size = 10, face = "plain"),
        axis.title.x = element_text(color = "black", size = 11, face = "plain"),
        axis.title.y = element_text(color = "black", size = 11, face = "plain")) +  scale_fill_npg() +
  scale_colour_discrete(guide=guide_legend(override.aes=list(size=5))) +
  theme(legend.title=element_blank(),
        legend.key.size = unit(0.4, 'cm')) +
  geom_line(data = data.frame(x = c(4.6,15.5), y = c(-12,-12)),
            aes(x = x, y = y), colour = "red", size = 1.2) +
  geom_line(data = data.frame(x = c(4.4,0.5), y = c(-12,-12)),
            aes(x = x, y = y), colour = "blue", size = 1.2) +
  geom_text() +
  annotate(
    "text", label = "SARS-CoV-2+",
    x = 10, y = -12.8, size = 5, colour = "red", angle = 90
  ) +
  geom_text() +
  annotate(
    "text", label = "SARS-CoV-2-",
    x = 2.5, y = -12.8, size = 5, colour = "blue", angle = 90
  )


# geom_point(data = DF.SYMBOL, aes(x = TIME, fill = EVENT, shape = EVENT), size = )


#############

p1 <- ggplot(clinicalInfo,aes(x=Start, y=Patient, color=Phase)) +
  geom_segment(aes(x=Start,xend=End,yend=Patient),size=10) +
  scale_colour_discrete(guide=guide_legend(override.aes=list(size=7)))

p1 + xlab("Days") + ylab("") + theme_classic() +
  theme(axis.text.x = element_text(color = "black", size = 10, face = "plain"),
        axis.text.y = element_text(color = "black", size = 10, face = "plain"),
        axis.title.x = element_text(color = "black", size = 11, face = "plain"),
        axis.title.y = element_text(color = "black", size = 14, face = "plain"))


ggplot(clinicalInfo,aes(x=End, y=Patient, fill=Phase)) +
  geom_col()


dat <- read.table(text = "    ONE TWO THREE
                  1   23  234 324
                  2   34  534 12
                  3   56  324 124
                  4   34  234 124
                  5   123 534 654",sep = "",header = TRUE)
datm <- melt(cbind(dat, ind = rownames(dat)), id.vars = c('ind'))
ggplot(datm, aes(x = variable, y = ifelse(ind %in% 1:2, -value, value), fill = ind)) +
  geom_col(colour="black", width = 0.5) +
  coord_flip() + ylab("Days") + xlab("") + theme_classic() +
  theme(axis.text.x = element_text(color = "black", size = 10, face = "plain"),
        axis.text.y = element_text(color = "black", size = 10, face = "plain"),
        axis.title.x = element_text(color = "black", size = 11, face = "plain"),
        axis.title.y = element_text(color = "black", size = 11, face = "plain"))


ID<-rep(c(1,2),each=2)
EVENT <- rep(c("TBR","PBR"))
TIME <- c(90, 220,120,200)
DF.SYMBOL<-data.frame(ID,EVENT,TIME)

