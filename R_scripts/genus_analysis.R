############ Genus ##########################################


setwd("C:/Users/HASI9S/OneDrive/Documents/Alignments/KrakenAlignments/Kraken2")

#Kraken2 alignments

AllK2Files<-list.files()
GenusFileList3<-grep("_genus_abundance.txt", AllK2Files)
K2GenusFiles<-AllK2Files[GenusFileList3]
K2FileList<-gsub("_genus_abundance.txt", "", K2GenusFiles)


NICUFiles<-as.character(AllSampleKey$Sample)
K2NICUFiles<-K2FileList[K2FileList %in% NICUFiles]
ALLNICUFiles<-K2NICUFiles



MissingNICU<-subset(NICUFiles, ! NICUFiles %in% c(K2NICUFiles))



# Get all the species counts
# get the files
for(f in 1:length(K2NICUFiles)){
  fnr = paste(K2NICUFiles[f], "_genus_abundance.txt", sep = "")
  filename<-K2NICUFiles[f]
  
  assign(filename, read.csv(paste(fnr), header = TRUE, stringsAsFactors = FALSE))
}

# ok getting Bracken reassigned.

NewGenusTable<-read.csv(paste(K2NICUFiles[1], "_genus_abundance.txt", sep = ""), sep = "\t", header = TRUE)
names(NewGenusTable)[1]<-"Genus"
NewGenusTable$Genus<-gsub(" ", ".", NewGenusTable$Genus, fixed = TRUE)
NewGenusTable$Genus<-gsub("..", ".", NewGenusTable$Genus, fixed = TRUE)
NewGenusTable$Genus<-gsub("_", ".", NewGenusTable$Genus, fixed = TRUE)
NewGenusTable$Genus<-gsub("-", ".", NewGenusTable$Genus, fixed = TRUE)
NewGenusTable$Genus<-gsub("/", ".", NewGenusTable$Genus, fixed = TRUE)
NewGenusTable$Genus<-gsub("X,", "", NewGenusTable$Genus, fixed = TRUE)
NewGenusTable$Genus<-gsub("[", "", NewGenusTable$Genus, fixed = TRUE)
NewGenusTable$Genus<-gsub("]", "", NewGenusTable$Genus, fixed = TRUE)
NewGenusTable<-as.data.frame(NewGenusTable[,c(1,6)])
names(NewGenusTable)<-c("Genus",paste(paste(K2NICUFiles[1])))
NewGenusTable<-subset(NewGenusTable, ! duplicated(NewGenusTable$Genus))

for ( i in 2:length(K2NICUFiles)){
  x <- read.csv(paste(K2NICUFiles[i], "_genus_abundance.txt", sep = ""), sep = "\t", header = TRUE)
  x<-as.data.frame(x[,c(1,6)])
  names( x ) <- c("Genus", paste(K2NICUFiles[i]))
  x$Genus<-gsub(" ", ".", x$Genus, fixed = TRUE)
  x$Genus<-gsub("..", ".", x$Genus, fixed = TRUE)
  x$Genus<-gsub("_", ".", x$Genus, fixed = TRUE)
  x$Genus<-gsub("-", ".", x$Genus, fixed = TRUE)
  x$Genus<-gsub("/", ".", x$Genus, fixed = TRUE)
  x$Genus<-gsub("X,", "", x$Genus, fixed = TRUE)
  x$Genus<-gsub("[", "", x$Genus, fixed = TRUE)
  x$Genus<-gsub("]", "", x$Genus, fixed = TRUE)
  x<-subset(x, ! duplicated(x$Genus))
  NewGenusTable<-merge(x, NewGenusTable, by = "Genus", all = TRUE)
}

row.names(NewGenusTable)<-NewGenusTable$Genus
NewGenusTable$Genus<-NULL
NewGenusTable[is.na(NewGenusTable)]<-0

rm(list = K2NICUFiles)


AllGenusRaw<-t(NewGenusTable)
Genusdf<-as.data.frame(AllGenusRaw)
Genusdf$Sample<-row.names(AllGenusRaw)
GenusTable<-as.data.frame(t(NewGenusTable))


# Now filter species to be present in at least 10% of samples, account for > 0.005% of counts,
# And get rid of samples that have < 1 million reads

NICUGenus<-as.data.frame(AllGenusRaw)
SaliniCol<-grep("Salinibacter", names(NICUGenus))
NICUGenus<-NICUGenus[,-SaliniCol]
# Human crap data
Homo<-grep("Homo", names(NICUGenus))
NICUGenus<-NICUGenus[,-Homo]
EpichloeCol<-grep("Epichloe", names(NICUGenus))
NICUGenus<-NICUGenus[,-EpichloeCol]
DermatoCol<-grep("Dermatophagoides", names(NICUGenus))
NICUGenus<-NICUGenus[,-DermatoCol]
ToxoCol<-grep("Toxoplasma", names(NICUGenus))
NICUGenus<-NICUGenus[,-ToxoCol]
MyriCol<-grep("Myriosclerotinia", names(NICUGenus))
NICUGenus<-NICUGenus[,-MyriCol]
MelCol<-grep("Melampsora", names(NICUGenus))
NICUGenus<-NICUGenus[,-MelCol]
MycoCol<-grep("Mycobacteroides", names(NICUGenus))
NICUGenus<-NICUGenus[,-MycoCol]
Zasmidium<-grep("Zasmidium", names(NICUGenus))
NICUGenus<-NICUGenus[,-Zasmidium]
Yamadazyma<-grep("Yamadazyma", names(NICUGenus))
NICUGenus<-NICUGenus[,-Yamadazyma]



TenPercentCutoff<-floor(nrow(NICUGenus)/20)

NonZeroCounts<-list()
for (i in 1:ncol(NICUGenus)){
  NonZeroCounts[i]<-length(which(NICUGenus[,i] > 0))
  
}

TenPercentNotZero<-which(NonZeroCounts >= TenPercentCutoff)
NICUGenusReduced<-NICUGenus[,TenPercentNotZero]
LowSamples<-which(rowSums(NICUGenusReduced) <= 500000)
NICUGenusReduced<-NICUGenusReduced[-LowSamples,]

NICUGenusNR<-as.data.frame(t(noise.removal(t(NICUGenusReduced), 0.01)))
minCount<-min(rowSums(NICUGenusNR))
NICUGenusNR<-data.frame(rrarefy(NICUGenusNR, minCount))
NICUGenusNR$Sample<-row.names(NICUGenusNR)


NICUGenusNR$Sample<-row.names(NICUGenusNR)
# NICUGenusReduced has all samples and all species except the three with low read counts

NICUGenusReduced$Sample<-row.names(NICUGenusReduced) 
NICUGenusReduced<-subset(NICUGenusReduced, ! duplicated(NICUGenusReduced$Sample))

NonRarefiedGenus<-merge(AllSampleKey, NICUGenusReduced, by = "Sample", all.y = TRUE)
NonRarefiedGenus<-subset(NonRarefiedGenus, ! duplicated(NonRarefiedGenus$Sample))

# Get diversity of all samples:
NonRarefiedGenus<-subset(NonRarefiedGenus, ! is.na(NonRarefiedGenus$SampleType))
row.names(NonRarefiedGenus)<-NonRarefiedGenus$Sample
# NonRarefiedGenus<-subset(NonRarefiedGenus, NonRarefiedGenus$Location == "Cincinnati" & NonRarefiedGenus$SampleType == "Stool")
SampleData<-NonRarefiedGenus[, 1:17]
NonRarefiedGenus<-NonRarefiedGenus[, 18:ncol(NonRarefiedGenus)]

DiversityGenus<-NonRarefiedGenus

# Calculate sample diversity using 'diversity' function in Vegan
H <- data.frame(diversity(DiversityGenus))
simpson <- data.frame(diversity(DiversityGenus, "simpson"))
shannon<-data.frame(diversity(DiversityGenus, "shannon"))
invsimp <- data.frame(diversity(DiversityGenus, "inv"))
alpha <- data.frame(fisher.alpha(DiversityGenus))

#BergerParker 
library(diverse) # this is going to over ride vegan package
bergerparker<-data.frame(diversity(as.matrix(DiversityGenus), type =  "berger-parker"))
library(vegan)
## Genus richness (S) and Pielou's evenness (J):
S <- data.frame(specnumber(DiversityGenus))
J <- data.frame(H/log(S))
Diversity<-cbind(simpson, shannon, invsimp, alpha, S, J, bergerparker)
Diversity$Sample<-row.names(Diversity)
names(Diversity)<-c("Simpson", "Shannon",  "InvSimpson", "Alpha", "GenusNo", "Evenness", "Sample", "berger.parker.D", "berger.parker.I")
Diversity<-Diversity[,c(7,1,2,3,4, 5, 6)]
pairs(cbind(simpson, shannon, alpha, J, S), pch="+", col="blue")


Diversity<-merge(SampleData, Diversity, by = "Sample", all.y  = TRUE)
Diversity<-subset(Diversity, ! duplicated(Diversity$Sample))

setwd("C:/Users/dbhas/OneDrive/Documents/Code/Metagenomics/Yanping/NICU_Microbiome/Hangzhou/NoHumanDNA20220929")
setwd("C:/Users/HASI9S/OneDrive/Documents/Code/Metagenomics/Yanping/NICU_Microbiome/Hangzhou/NoHumanDNA20220929")
setwd("/home/david/Documents/Code/Metagenomics/Yanping/NICU_Microbiome/Hangzhou/NoHumanDNA20220929")


## Back to rarefied species

NICUGenusNR<-merge(AllSampleKey, NICUGenusNR, by = "Sample", all.y = TRUE)
NICUGenusNR<-subset(NICUGenusNR, ! duplicated(NICUGenusNR$Sample))


# What's the mean gestational age in each cohort

GATable<-NICUGenusNR[,c(1,2,4,5,6,7,8,9,10)]
GATable<-subset(GATable, GATable$SampleType == "Stool")
GATable<-subset(GATable, GATable$SampleCollectionWeek == "Week.1")
GATable$Weeks<- as.integer(sapply(strsplit(GATable$GestationalAge, " "), `[`, 1))
GATable$Days<- sapply(strsplit(GATable$GestationalAge, " "), `[`, 2)
GATable$Days[is.na(GATable$Days)]<-"0/7"
GATable$Fraction<- as.numeric(sapply(strsplit(GATable$Days, "/"), `[`, 1))/7
GATable$GA<-GATable$Weeks + GATable$Fraction

GAAbx<-ggplot(GATable,
              aes(x=PostNatalAntibiotics, y=as.numeric(GA), fill=PostNatalAntibiotics)) + geom_boxplot(lwd=1,aes(color=factor(PostNatalAntibiotics),fill = NA), outlier.size = 3)  +
  stat_summary(fun=mean, geom="point", shape=5, size=8) + scale_fill_manual(values = col) + scale_colour_manual(values = col ) +
  geom_point(size=4,aes(color = factor(PostNatalAntibiotics))) + xlab(NULL) + 
  ylab("Gestational Age \n") + paramsBox() + scale_y_log10() + xlab(NULL)  + scale_color_tableau() +
  ylim(24,37)
ggsave(GAAbx, file = "GestataionalAgeAbx.pdf", width = 6, height = 8, limitsize = FALSE)

pairwise.wilcox.test(GATable$GA, GATable$PostNatalAntibiotics)

GATableMeanAge<-as.data.frame(GATable %>%
                                group_by(GestationCohort, PostNatalAntibiotics) %>%
                                summarize(mean(GA)))

GATableSubjectCount<-as.data.frame(GATable %>%
                                     group_by(GestationCohort, PostNatalAntibiotics) %>%
                                     summarize(n()))

# What about just the 28-32 weeks infants

GATable2<-NICUGenusNR[,c(1,2,4,5,6,7,8,9,10)]
GATable2<-subset(GATable2, GATable2$SampleType == "Stool" & GATable2$GestationCohort == "28-32 Weeks")
GATable2<-subset(GATable2, GATable2$SampleCollectionWeek == "Week.1")
GATable2$Weeks<- as.integer(sapply(strsplit(GATable2$GestationalAge, " "), `[`, 1))
GATable2$Days<- sapply(strsplit(GATable2$GestationalAge, " "), `[`, 2)
GATable2$Days[is.na(GATable2$Days)]<-"0/7"
GATable2$Fraction<- as.numeric(sapply(strsplit(GATable2$Days, "/"), `[`, 1))/7
GATable2$GA<-GATable2$Weeks + GATable2$Fraction

GATable3<-subset(GATable2, GATable2$SampleType == "Stool" & GATable2$GA < 31)


GAAbx3<-ggplot(GATable3,
               aes(x=PostNatalAntibiotics, y=as.numeric(GA), fill=PostNatalAntibiotics)) + geom_boxplot(lwd=1,aes(color=factor(PostNatalAntibiotics),fill = NA), outlier.size = 3)  +
  stat_summary(fun=mean, geom="point", shape=5, size=8) + scale_fill_manual(values = col) + scale_colour_manual(values = col ) +
  geom_point(size=4,aes(color = factor(PostNatalAntibiotics))) + xlab(NULL) + 
  ylab("Gestational Age \n") + paramsBox() + scale_y_log10() + xlab(NULL)  + scale_color_tableau() +
  ylim(24,34)
ggsave(GAAbx3, file = "GALessThan31.pdf", width = 6, height = 8, limitsize = FALSE)

pairwise.wilcox.test(GATable3$GA, GATable3$PostNatalAntibiotics)

save.image(file = "NICUData20250209")

##### Generalized Linear Mixed Effects Models ####


GLMGenus<-NICUGenusNR
GLMGenus$PatientID<-paste(GLMGenus$PatientID, GLMGenus$Location, sep = "_")

Location=factor(GLMGenus$Location)
Abx = factor(GLMGenus$PostNatalAntibiotics); table(Abx)
MatAbx = factor(GLMGenus$MaternalAntibiotics)
Gestation = GLMGenus$GestationTime
Subject = GLMGenus[, "PatientID"]; table(Subject)
Milk = factor(GLMGenus$AnyMilk)
SampleType = factor(GLMGenus$SampleType)
Wk1Abx=factor(GLMGenus$Week1Abx)
Wk3Abx=factor(GLMGenus$Week3Abx)
Week<-as.integer(ifelse(GLMGenus$SampleCollectionWeek == "Week.1", 7, ifelse(GLMGenus$SampleCollectionWeek == "Week.3", 21, NA)))

GLMGenusLocation<-glmer(GLMGenus[,18] ~ Location+Week+SampleType+Wk1Abx+Wk3Abx+Abx+Gestation+MatAbx+Milk+(1|Subject), data = NICUGenusNR[,c(18:ncol(NICUGenusNR))], family = poisson)
GLMGenusLocationZINB <- update(GLMGenusLocation,family=nbinom2)
pvalues<-as.data.frame(summary(GLMGenusLocationZINB)$coeff[-1,4])
names(pvalues)<-names(GLMGenus)[18]
estimates<-as.data.frame(summary(GLMGenusLocationZINB)$coeff[-1,1])
names(estimates)<-names(GLMGenus)[18]


for(i in 19:ncol(GLMGenus)){
  tryCatch({
    GLMGenusLocation<-glmer(GLMGenus[,i] ~ Location+Week+SampleType+Wk1Abx+Wk3Abx+Abx+Gestation+MatAbx+Milk+(1|Subject), data = NICUGenusNR, family = poisson)
    GLMGenusLocationZINB <- update(GLMGenusLocation,family=nbinom2)
  }, error=function(e){cat("ERROR :",conditionMessage(e), "\n")})
  pvaluesNew<-as.data.frame(summary(GLMGenusLocationZINB)$coeff[-1,4])
  names(pvaluesNew)<-names(GLMGenus)[i]
  estimatesNew<-as.data.frame(summary(GLMGenusLocationZINB)$coeff[-1,1])
  names(estimatesNew)<-names(GLMGenus)[i]
  pvalues<-cbind(pvalues, pvaluesNew)
  estimates<-cbind(estimates, estimatesNew)
  
}

tpvalues<-as.data.frame(t(pvalues))
tpvalues$Location.FDR<-p.adjust(tpvalues$LocationHangzhou, method = "fdr")
tpvalues$Week.FDR<-p.adjust(tpvalues$Week, method = "fdr")
tpvalues$SampleTypeGroin.FDR<-p.adjust(tpvalues$SampleTypeGroin, method = "fdr")
tpvalues$SampleTypeStool.FDR<-p.adjust(tpvalues$SampleTypeStool, method = "fdr")
tpvalues$Wk1AbxHigh.FDR<-p.adjust(tpvalues$Wk1AbxHigh.Infant.Abx, method = "fdr")
tpvalues$Wk1AbxLow.FDR<-p.adjust(tpvalues$Wk1AbxLow.Infant.Abx, method = "fdr")
tpvalues$Wk1AbxNo.FDR<-p.adjust(tpvalues$Wk1AbxNo.Infant.Abx, method = "fdr")
tpvalues$Wk3AbxHigh.FDR<-p.adjust(tpvalues$Wk3AbxHigh.Infant.Abx, method = "fdr")
tpvalues$Wk3AbxLow.FDR<-p.adjust(tpvalues$Wk3AbxLow.Infant.Abx, method = "fdr")
tpvalues$Wk3AbxNo.FDR<-p.adjust(tpvalues$Wk3AbxNo.Infant.Abx, method = "fdr")
tpvalues$AbxInfant.FDR<-p.adjust(tpvalues$AbxInfant.Abx, method = "fdr")
tpvalues$AbxMat.FDR<-p.adjust(tpvalues$MatAbxMat.Abx, method = "fdr")
tpvalues$Gestation.FDR<-p.adjust(tpvalues$Gestation, method = "fdr")
tpvalues$MilkNo.FDR<-p.adjust(tpvalues$MilkNo.Milk, method = "fdr")
tpvalues$MilkMother.FDR<-p.adjust(tpvalues$MilkMother, method = "fdr")

HangzhouCincinnatiGLMMGenus<-tpvalues
write.csv(tpvalues, file = "HangzhouCincinnatiGLMMGenus.csv")

HangzhouCincinnatiGLMMGenusEstimates<-t(estimates)
write.csv(HangzhouCincinnatiGLMMGenusEstimates, file = "HangzhouCincinnatiGLMMGenusEstimates.csv")


# here we can try to plot the estimates, p-values, etc

library(lme4)
library(lmerTest)
library(broom.mixed)
library(ggplot2)

# Fit a GLMM model (example)
glmm_model <- lmer(response ~ predictor1 + predictor2 + (1 | group), data = dataset, family = "binomial")

# Extract model summary
model_summary <- summary(glmm_model)

# Extract estimates and confidence intervals
tidy_model <- tidy(glmm_model, effects = "fixed", conf.int = TRUE)

# Plotting
ggplot(tidy_model, aes(x = term, y = estimate)) +
  geom_point() +
  geom_errorbar(aes(ymin = conf.low, ymax = conf.high), width = 0.2) +
  geom_text(aes(label = format.pval(p.value, digits = 3)), vjust = -1.5)

# or use forest plots, like here: https://www.mihiretukebede.com/posts/2020-09-30-2020-09-30-plotting-model-coefficients-in-a-forest-plot/



save.image(file = "NICUData20250209")

Week1AbxGenusTable<-as.data.frame(tpvalues[,20])
row.names(Week1AbxGenusTable)<-row.names(tpvalues)
names(Week1AbxGenusTable)<-"AbxInfant.FDR"
Week1AbxGenusTable<-subset(Week1AbxGenusTable, Week1AbxGenusTable$AbxInfant.FDR < 0.1)

Week3AbxGenusTable<-as.data.frame(tpvalues[,23])
row.names(Week3AbxGenusTable)<-row.names(tpvalues)
names(Week3AbxGenusTable)<-"AbxInfant.FDR"
Week3AbxGenusTable<-subset(Week3AbxGenusTable, Week3AbxGenusTable$AbxInfant.FDR < 0.1)


MatAbxGenusTable<-as.data.frame(tpvalues[,27])
row.names(MatAbxGenusTable)<-row.names(tpvalues)
names(MatAbxGenusTable)<-"MatAbx.FDR"
SigGenusMatAbxTable<-subset(MatAbxGenusTable, MatAbxGenusTable$MatAbx.FDR < 0.1)



GestationGenusTable<-as.data.frame(tpvalues[,28])
row.names(GestationGenusTable)<-row.names(tpvalues)
names(GestationGenusTable)<-"Gestation.FDR"
GestationGenusTable<-subset(GestationGenusTable, GestationGenusTable$Gestation.FDR < 0.1)

MilkGenusTable<-as.data.frame(tpvalues[,30])
row.names(MilkGenusTable)<-row.names(tpvalues)
names(MilkGenusTable)<-"Milk.FDR"
MilkGenusTable<-subset(MilkGenusTable, MilkGenusTable$Milk.FDR < 0.1)

NoMilkGenusTable<-as.data.frame(tpvalues[,29])
row.names(NoMilkGenusTable)<-row.names(tpvalues)
names(NoMilkGenusTable)<-"NoMilk.FDR"
NoMilkGenusTable<-subset(NoMilkGenusTable, NoMilkGenusTable$NoMilk.FDR < 0.1)

# DeliveryGenusTable<-as.data.frame(tpvalues[,10])
# row.names(DeliveryGenusTable)<-row.names(tpvalues)
# names(DeliveryGenusTable)<-"Delivery.FDR"
# DeliveryGenusTable<-subset(DeliveryGenusTable, DeliveryGenusTable$Delivery.FDR < 0.1)

LocationGenusTable<-as.data.frame(tpvalues[,16])
row.names(LocationGenusTable)<-row.names(tpvalues)
names(LocationGenusTable)<-"Location.FDR"
LocationGenusTable<-subset(LocationGenusTable, LocationGenusTable$Location.FDR < 0.1)

PostnatalGenusTable<-as.data.frame(tpvalues[,17])
row.names(PostnatalGenusTable)<-row.names(tpvalues)
names(PostnatalGenusTable)<-"Postnatal.FDR"
PostnatalGenusTable<-subset(PostnatalGenusTable, PostnatalGenusTable$Postnatal.FDR < 0.1)

SampleTypeGroinGenusTable<-as.data.frame(tpvalues[,18])
row.names(SampleTypeGroinGenusTable)<-row.names(tpvalues)
names(SampleTypeGroinGenusTable)<-"SampleTypeGroin.FDR"
SampleTypeGroinGenusTable<-subset(SampleTypeGroinGenusTable, SampleTypeGroinGenusTable$SampleTypeGroin.FDR < 0.1)

save.image(file = "NICUData20250209")


DummyGenusTable<-NICUGenusNR

# col=tableau.colors
# DummyGenusTable$Clostridium[DummyGenusTable$Clostridium==0]<-1
# ClostridiumGestation<-ggplot(subset(DummyGenusTable, DummyGenusTable$SampleType %in% c("Stool")),
#                              aes(x=GestationCohort, y=as.numeric(Clostridium), fill=GestationCohort)) + geom_boxplot(lwd=1,aes(color=factor(GestationCohort),fill = NA), outlier.size = 3)  +
#   stat_summary(fun=mean, geom="point", shape=5, size=8) + scale_fill_manual(values = col) + scale_colour_manual(values = col ) +
#   geom_point(size=4,aes(color = factor(GestationCohort))) + xlab(NULL) + facet_grid(rows = . ~ PostNatalAntibiotics) +
#   ylab("Clostridium \n") + paramsAngled() + scale_y_log10() + xlab(NULL)  + scale_color_tableau() +
#   theme(axis.title.y = element_text(face = "italic"))
# ggsave(ClostridiumGestation, file = "ClostridiumGestation.pdf", width = 12, height = 8, limitsize = FALSE)

### Location graphs
Clostridium<-ggplot(DummyGenusTable, 
                    aes(x=Location, y=as.numeric(Clostridium), fill=Location)) + geom_boxplot(lwd=1,aes(color=factor(Location),fill = NA), outlier.size = 3)  +
  stat_summary(fun=mean, geom="point", shape=5, size=8) + scale_fill_manual(values = col) + scale_colour_manual(values = col ) +
  geom_point(size=4,aes(color = factor(Location))) + xlab(NULL) + facet_grid(rows =  SampleCollectionWeek ~ SampleType) +
  ylab("Clostridium \n") + paramsAngled() + scale_y_log10() + xlab(NULL)  + scale_color_tableau() +
  theme(axis.title.y = element_text(face = "italic"))
ggsave(Clostridium, file = "Clostridium.pdf", width = 14, height = 12, limitsize = FALSE)


Staphylococcus<-ggplot(DummyGenusTable, 
                       aes(x=Location, y=as.numeric(Staphylococcus), fill=Location)) + geom_boxplot(lwd=1,aes(color=factor(Location),fill = NA), outlier.size = 3)  +
  stat_summary(fun=mean, geom="point", shape=5, size=8) + scale_fill_manual(values = col) + scale_colour_manual(values = col ) +
  geom_point(size=4,aes(color = factor(Location))) + xlab(NULL) + facet_grid(rows =  SampleCollectionWeek ~ SampleType) +
  ylab("Staphylococcus \n") + paramsAngled() + scale_y_log10() + xlab(NULL)  + scale_color_tableau() +
  theme(axis.title.y = element_text(face = "italic"))
ggsave(Staphylococcus, file = "Staphylococcus_new.pdf", width = 12, height = 8, limitsize = FALSE)

Bacteroides<-ggplot(DummyGenusTable, 
                    aes(x=Location, y=as.numeric(Bacteroides), fill=Location)) + geom_boxplot(lwd=1,aes(color=factor(Location),fill = NA), outlier.size = 3)  +
  stat_summary(fun=mean, geom="point", shape=5, size=8) + scale_fill_manual(values = col) + scale_colour_manual(values = col ) +
  geom_point(size=4,aes(color = factor(Location))) + xlab(NULL) + facet_grid(rows =  SampleCollectionWeek ~ SampleType) +
  ylab("Bacteroides \n") + paramsAngled() + scale_y_log10() + xlab(NULL)  + scale_color_tableau() +
  theme(axis.title.y = element_text(face = "italic"))
ggsave(Bacteroides, file = "Bacteroides.pdf", width = 12, height = 8, limitsize = FALSE)


Klebsiella<-ggplot(DummyGenusTable, 
                   aes(x=Location, y=as.numeric(Klebsiella), fill=Location)) + geom_boxplot(lwd=1,aes(color=factor(Location),fill = NA), outlier.size = 3)  +
  stat_summary(fun=mean, geom="point", shape=5, size=8) + scale_fill_manual(values = col) + scale_colour_manual(values = col ) +
  geom_point(size=4,aes(color = factor(Location))) + xlab(NULL) + 
  ylab("Klebsiella \n") + paramsAngled() + scale_y_log10() + xlab(NULL)  + scale_color_tableau() +
  theme(axis.title.y = element_text(face = "italic")) + facet_grid(rows =  SampleCollectionWeek ~ SampleType) 
ggsave(Klebsiella, file = "Klebsiella_new.pdf", width = 12, height = 8, limitsize = FALSE)

Streptococcus<-ggplot(DummyGenusTable, 
                      aes(x=Location, y=as.numeric(Streptococcus), fill=Location)) + geom_boxplot(lwd=1,aes(color=factor(Location),fill = NA), outlier.size = 3)  +
  stat_summary(fun=mean, geom="point", shape=5, size=8) + scale_fill_manual(values = col) + scale_colour_manual(values = col ) +
  geom_point(size=4,aes(color = factor(Location))) + xlab(NULL) + 
  ylab("Streptococcus \n") + paramsAngled() + scale_y_log10() + xlab(NULL)  + scale_color_tableau() +
  theme(axis.title.y = element_text(face = "italic")) + facet_grid(rows =  SampleCollectionWeek ~ SampleType) 
ggsave(Streptococcus, file = "Streptococcus_new.pdf", width = 12, height = 8, limitsize = FALSE)

Veillonella<-ggplot(DummyGenusTable, 
                    aes(x=Location, y=as.numeric(Veillonella), fill=Location)) + geom_boxplot(lwd=1,aes(color=factor(Location),fill = NA), outlier.size = 3)  +
  stat_summary(fun=mean, geom="point", shape=5, size=8) + scale_fill_manual(values = col) + scale_colour_manual(values = col ) +
  geom_point(size=4,aes(color = factor(Location))) + xlab(NULL) + 
  ylab("Veillonella \n") + paramsAngled() + scale_y_log10() + xlab(NULL)  + scale_color_tableau() +
  theme(axis.title.y = element_text(face = "italic")) + facet_grid(rows =  SampleCollectionWeek ~ SampleType) 
ggsave(Veillonella, file = "Veillonella.pdf", width = 12, height = 8, limitsize = FALSE)

Enterobacter<-ggplot(DummyGenusTable, 
                     aes(x=Location, y=as.numeric(Enterobacter), fill=Location)) + geom_boxplot(lwd=1,aes(color=factor(Location),fill = NA), outlier.size = 3)  +
  stat_summary(fun=mean, geom="point", shape=5, size=8) + scale_fill_manual(values = col) + scale_colour_manual(values = col ) +
  geom_point(size=4,aes(color = factor(Location))) + xlab(NULL) + 
  ylab("Enterobacter \n") + paramsAngled() + scale_y_log10() + xlab(NULL)  + scale_color_tableau() +
  theme(axis.title.y = element_text(face = "italic")) + facet_grid(rows =  SampleCollectionWeek ~ SampleType) 
ggsave(Enterobacter, file = "Enterobacter_new.pdf", width = 12, height = 8, limitsize = FALSE)

Acinetobacter<-ggplot(DummyGenusTable, 
                      aes(x=Location, y=as.numeric(Acinetobacter), fill=Location)) + geom_boxplot(lwd=1,aes(color=factor(Location),fill = NA), outlier.size = 3)  +
  stat_summary(fun=mean, geom="point", shape=5, size=8) + scale_fill_manual(values = col) + scale_colour_manual(values = col ) +
  geom_point(size=4,aes(color = factor(Location))) + xlab(NULL) + 
  ylab("Acinetobacter \n") + paramsAngled() + scale_y_log10() + xlab(NULL)  + scale_color_tableau() +
  theme(axis.title.y = element_text(face = "italic")) + facet_grid(rows =  SampleCollectionWeek ~ SampleType) 
ggsave(Acinetobacter, file = "Acinetobacter_new.pdf", width = 12, height = 8, limitsize = FALSE)

Escherichia<-ggplot(DummyGenusTable, 
                    aes(x=Location, y=as.numeric(Escherichia), fill=Location)) + geom_boxplot(lwd=1,aes(color=factor(Location),fill = NA), outlier.size = 3)  +
  stat_summary(fun=mean, geom="point", shape=5, size=8) + scale_fill_manual(values = col) + scale_colour_manual(values = col ) +
  geom_point(size=4,aes(color = factor(Location))) + xlab(NULL) + 
  ylab("Escherichia \n") + paramsAngled() + scale_y_log10() + xlab(NULL)  + scale_color_tableau() +
  theme(axis.title.y = element_text(face = "italic")) + facet_grid(rows =  SampleCollectionWeek ~ SampleType) 
ggsave(Escherichia, file = "Escherichia_new.pdf", width = 12, height = 8, limitsize = FALSE)

Escherichia<-ggplot(DummyGenusTable, 
                    aes(x=Location, y=as.numeric(Escherichia), fill=Location)) + geom_boxplot(lwd=1,aes(color=factor(Location),fill = NA), outlier.size = 3)  +
  stat_summary(fun=mean, geom="point", shape=5, size=8) + scale_fill_manual(values = col) + scale_colour_manual(values = col ) +
  geom_point(size=4,aes(color = factor(Location))) + xlab(NULL) + 
  ylab("Escherichia \n") + paramsAngled() + scale_y_log10() + xlab(NULL)  + scale_color_tableau() +
  theme(axis.title.y = element_text(face = "italic")) + facet_grid(rows =  SampleCollectionWeek) 
ggsave(Escherichia, file = "Escherichia_new_location.pdf", width = 12, height = 8, limitsize = FALSE)

Acinetobacter<-ggplot(DummyGenusTable, 
                      aes(x=Location, y=as.numeric(Acinetobacter), fill=Location)) + geom_boxplot(lwd=1,aes(color=factor(Location),fill = NA), outlier.size = 3)  +
  stat_summary(fun=mean, geom="point", shape=5, size=8) + scale_fill_manual(values = col) + scale_colour_manual(values = col ) +
  geom_point(size=4,aes(color = factor(Location))) + xlab(NULL) + 
  ylab("Acinetobacter \n") + paramsAngled() + scale_y_log10() + xlab(NULL)  + scale_color_tableau() +
  theme(axis.title.y = element_text(face = "italic")) + facet_grid(rows =  SampleCollectionWeek) 
ggsave(Acinetobacter, file = "Acinetobacter_new_location.pdf", width = 12, height = 8, limitsize = FALSE)

Klebsiella<-ggplot(DummyGenusTable, 
                   aes(x=Location, y=as.numeric(Klebsiella), fill=Location)) + geom_boxplot(lwd=1,aes(color=factor(Location),fill = NA), outlier.size = 3)  +
  stat_summary(fun=mean, geom="point", shape=5, size=8) + scale_fill_manual(values = col) + scale_colour_manual(values = col ) +
  geom_point(size=4,aes(color = factor(Location))) + xlab(NULL) + 
  ylab("Klebsiella \n") + paramsAngled() + scale_y_log10() + xlab(NULL)  + scale_color_tableau() +
  theme(axis.title.y = element_text(face = "italic")) + facet_grid(facets =  SampleCollectionWeek) 
ggsave(Klebsiella, file = "Klebsiella_new_location.pdf", width = 12, height = 8, limitsize = FALSE)


save.image(file = "NICUData20250209")


DummyStoolTable<-subset(DummyGenusTable, DummyGenusTable$SampleType == "Stool")

fit1 <- lm(Clostridium ~ GestationTime, data = DummyStoolTable)
plot1<-ggplotRegression(fit1) + scale_y_log10()
pdf("ClostridiumGestation.pdf")
print(plot1)
dev.off()

fit2 <- lm(Klebsiella ~ GestationTime, data = DummyStoolTable)
plot2<-ggplotRegression(fit2) + scale_y_log10()
pdf("KlebsiellaGestation.pdf")
print(plot2)
dev.off()

fit3 <- lm(Clostridium ~ GestationTime, data = DummyStoolTable)
plot3<-ggplotRegression(fit3) + scale_y_log10()
pdf("ClostridiumGestation.pdf")
print(plot3)
dev.off()

# Ok this is actually both week 1 and week 3

CorrelationTable<-as.data.frame(cor(DummyStoolTable[,18:ncol(DummyStoolTable)], y = DummyStoolTable[,8],
                                    method = c("spearman")), use = "pairwise.complete.obs")
names(CorrelationTable)<-"GestationCorrelation"
write.csv(CorrelationTable, file = "GenusGestationCorrelation.csv")


fit = lm(DummyStoolTable[,8] ~ DummyStoolTable[,18])
RegressionTable<-as.data.frame(summary(fit)$coefficients[,4])[2,]
names(RegressionTable)<-names(DummyStoolTable[18])

for (i in 19:ncol(DummyStoolTable)){
  fit = lm(DummyStoolTable[,8] ~ DummyStoolTable[,i])
  x<-as.data.frame(summary(fit)$coefficients[,4])[2,]
  names(x)<-names(DummyStoolTable[i])
  RegressionTable<-rbind(RegressionTable, x)
}

row.names(RegressionTable)<-colnames(DummyStoolTable[18:ncol(DummyStoolTable)])
RegressionTable<-as.data.frame(RegressionTable)
colnames(RegressionTable)<-"Regression_pValue"
write.csv(RegressionTable, file = "Week1GestationRegressionTable.csv")

SigCorrStoolWeek1Gestation<-intersect(row.names(CorrelationTable)[abs(CorrelationTable$GestationCorrelation) > 0.2],
                                      row.names(RegressionTable)[RegressionTable$Regression_pValue < 0.05])

# do it again with week 3

StoolWeek3Data<-subset(DummyStoolTable, DummyStoolTable$SampleCollectionWeek == "Week.3")

CorrelationTableWeek3<-as.data.frame(cor(StoolWeek3Data[,18:ncol(StoolWeek3Data)], y = StoolWeek3Data[,8],
                                         method = c("spearman")), use = "pairwise.complete.obs")
names(CorrelationTableWeek3)<-c("GestationCorrelation")
SigCorrelationTableWeek3<-subset(CorrelationTableWeek3, abs(CorrelationTableWeek3$GestationCorrelation)>0.2)


fit = lm(StoolWeek3Data[,8] ~ StoolWeek3Data[,18])
RegressionTable<-as.data.frame(summary(fit)$coefficients[,4])[2,]
names(RegressionTable)<-names(StoolWeek3Data[18])

for (i in 19:ncol(StoolWeek3Data)){
  fit = lm(StoolWeek3Data[,8] ~ StoolWeek3Data[,i])
  x<-as.data.frame(summary(fit)$coefficients[,4])[2,]
  names(x)<-names(StoolWeek3Data[i])
  RegressionTable<-rbind(RegressionTable, x)
}

row.names(RegressionTable)<-colnames(StoolWeek3Data[18:ncol(StoolWeek3Data)])
RegressionTable<-as.data.frame(RegressionTable)
colnames(RegressionTable)<-"Regression_pValue"
write.csv(RegressionTable, file = "Week3GenusGestationRegressionTable.csv")

SigCorrStoolWeek3Gestation<-intersect(row.names(CorrelationTable)[abs(CorrelationTableWeek3$GestationCorrelation) > 0.2],
                                      row.names(RegressionTable)[RegressionTable$Regression_pValue < 0.05])


fitStaph <- lm(Staphylococcus ~ GestationTime, data = StoolWeek3Data)
plotStaph<-ggplotRegression(fitStaph) + scale_y_log10()
pdf("StaphWeek3Gestation.pdf")
print(plotStaph)
dev.off()

fitStrep <- lm(Streptococcus ~ GestationTime, data = StoolWeek3Data)
plotStrep<-ggplotRegression(fitStrep) + scale_y_log10()
pdf("StrepWeek3Gestation.pdf")
print(plotStrep)
dev.off()

fitBacteroides <- lm(Bacteroides ~ GestationTime, data = StoolWeek3Data)
plotBacteroides<-ggplotRegression(fitBacteroides) + scale_y_log10()
pdf("BacteroidesWeek3Gestation.pdf")
print(plotBacteroides)
dev.off()

fitBifido <- lm(Bifidobacterium ~ GestationTime, data = StoolWeek3Data)
plotBifido<-ggplotRegression(fitBifido) #+ scale_y_log10()
pdf("BifidoWeek3Gestation.pdf")
print(plotBifido)
dev.off()


# Alright we're going to look at the two largest gestation cohorts and calcuate effect size for the genera between those two gestations at week 1 and week 3

StoolGenusTable<-subset(NICUGenusNR, NICUGenusNR$SampleType == "Stool")

# StoolWeek1Data$GestationCohort<-gsub("28-32 Weeks", "Cohort.1", StoolWeek1Data$GestationCohort)
# StoolWeek1Data$GestationCohort<-gsub("33-36 Weeks", "Cohort.2", StoolWeek1Data$GestationCohort)


metadata <- StoolGenusTable[, 1:17]
cts <- as.matrix(StoolGenusTable[, -(1:17)])
rownames(cts) <- metadata$Sample

cts_l2 <- glog2(cts)
grps <- StoolGenusTable$Location


res <- compute_ef_Location(d = cts_l2, g = StoolGenusTable$Location, min_shrink = 0.3)

# res <- subset(res, res$Species %in% SigCorrStoolWeek1Gestation)

resLocation<-res
resLocation$lfdr



write.table(resGestationWeek1, file = "SignificantStoolGenusLocatoin.csv", sep = ",", row.names = FALSE, col.names = TRUE,  quote = FALSE)


Gestation_Effect<- ggplot(filter(resGestationWeek1, abs(resGestationWeek1$ef_shrunk) > 0.3),
                          aes(x = reorder(Species, ef_shrunk), y = ef_shrunk, fill = HigherIn)) +
  geom_bar(colour="black", stat="identity") +
  #guides(fill=FALSE) + scale_fill_hue(l=40) +
  guides(fill=guide_legend(override.aes=list(colour=NA))) +
  coord_flip() +
  labs( x = NULL) +
  ylab(expression(atop("Effect Size"))) +
  guides(fill=guide_legend(title="More Abundant In:", )) +
  #geom_hline(yintercept = c(-0.2, 0.2), color = "firebrick", lty = 2) +
  theme(axis.text.y = element_text(size= 12, color="black")) +
  theme(legend.title = element_text(size=12)) +
  theme(legend.text = element_text(size = 11))

Gestation_Effect<-Gestation_Effect + scale_fill_manual(values = tableau.colors[c(1,2)])


ggsave(filename = "GestationWeek1_Effect.pdf", plot = Gestation_Effect, width = 10,
       height = 12, limitsize = FALSE)

# Week 3 gestation effect size

StoolWeek3GestationCohortTable<-subset(StoolWeek3Data, StoolWeek3Data$GestationCohort %in% c("28-32 Weeks", "33-36 Weeks"))

StoolWeek3Data$GestationCohort<-gsub("28-32 Weeks", "Cohort.1", StoolWeek3Data$GestationCohort)
StoolWeek3Data$GestationCohort<-gsub("33-36 Weeks", "Cohort.2", StoolWeek3Data$GestationCohort)


metadata <- StoolWeek3Data[, 1:17]
cts <- as.matrix(StoolWeek3Data[, -(1:17)])
rownames(cts) <- metadata$Sample

cts_l2 <- glog2(cts)
grps <- StoolWeek3Data$GestationCohort


res <- compute_ef_Gestation(d = cts_l2, g = StoolWeek3Data$GestationCohort, min_shrink = 0.3)

res <- subset(res, res$Species %in% SigCorrStoolWeek3Gestation)

resGestationWeek3<-res

resGestationWeek3$HigherIn<-ifelse(res$ef_shrunk > 0, "Cohort.1", "Cohort.2")
resGestationWeek3$HigherIn<-factor(resGestationWeek3$HigherIn, levels = c("Cohort.1", "Cohort.2"))


write.table(resGestationWeek3, file = "SignificantWeek3Gestation.csv", sep = ",", row.names = FALSE, col.names = TRUE,  quote = FALSE)


Gestation_Effect<- ggplot(filter(resGestationWeek3, abs(resGestationWeek3$ef_shrunk) > 0.3),
                          aes(x = reorder(Species, ef_shrunk), y = ef_shrunk, fill = HigherIn)) +
  geom_bar(colour="black", stat="identity") +
  #guides(fill=FALSE) + scale_fill_hue(l=40) +
  guides(fill=guide_legend(override.aes=list(colour=NA))) +
  coord_flip() +
  labs( x = NULL) +
  ylab(expression(atop("Effect Size"))) +
  guides(fill=guide_legend(title="More Abundant In:", )) +
  #geom_hline(yintercept = c(-0.2, 0.2), color = "firebrick", lty = 2) +
  theme(axis.text.y = element_text(size= 12, color="black")) +
  theme(legend.title = element_text(size=12)) +
  theme(legend.text = element_text(size = 11))

Gestation_Effect<-Gestation_Effect + scale_fill_manual(values = tableau.colors[c(1,2)])

ggsave(filename = "GestationWeek3_Effect.pdf", plot = Gestation_Effect, width = 10,
       height = 12, limitsize = FALSE)


DummyStoolTable$Clostridioides[DummyStoolTable$Clostridioides==0]<-1
ClostridioidesGestation<-ggplot(subset(DummyStoolTable, DummyStoolTable$SampleType %in% c("Stool")),
                                aes(x=GestationCohort, y=as.numeric(Clostridioides), fill=GestationCohort)) + geom_boxplot(lwd=1,aes(color=factor(GestationCohort),fill = NA), outlier.size = 3)  +
  stat_summary(fun=mean, geom="point", shape=5, size=8) + scale_fill_manual(values = col) + scale_colour_manual(values = col ) +
  geom_point(size=4,aes(color = factor(GestationCohort))) + xlab(NULL) + #facet_grid(rows = . ~ PostNatalAntibiotics) +
  ylab("Clostridioides \n") + paramsAngled() + scale_y_log10() + xlab(NULL)  + scale_color_tableau() +
  theme(axis.title.y = element_text(face = "italic"))
ggsave(ClostridioidesGestation, file = "ClostridioidesGestationWeek1.pdf", width = 6, height = 8, limitsize = FALSE)

DummyStoolTable2$Clostridioides[DummyStoolTable2$Clostridioides==0]<-1
ClostridioidesGestation<-ggplot(subset(DummyStoolTable2, DummyStoolTable2$SampleType %in% c("Stool")),
                                aes(x=GestationCohort, y=as.numeric(Clostridioides), fill=GestationCohort)) + geom_boxplot(lwd=1,aes(color=factor(GestationCohort),fill = NA), outlier.size = 3)  +
  stat_summary(fun=mean, geom="point", shape=5, size=8) + scale_fill_manual(values = col) + scale_colour_manual(values = col ) +
  geom_point(size=4,aes(color = factor(GestationCohort))) + xlab(NULL) + #facet_grid(rows = . ~ PostNatalAntibiotics) +
  ylab("Clostridioides \n") + paramsAngled() + scale_y_log10() + xlab(NULL)  + scale_color_tableau() +
  theme(axis.title.y = element_text(face = "italic"))
ggsave(ClostridioidesGestation, file = "ClostridioidesGestationWeek3.pdf", width = 6, height = 8, limitsize = FALSE)


DummyStoolTable$Klebsiella[DummyStoolTable$Klebsiella==0]<-1
KlebsiellaGestation<-ggplot(subset(DummyStoolTable, DummyStoolTable$SampleType %in% c("Stool")),
                            aes(x=GestationCohort, y=as.numeric(Klebsiella), fill=GestationCohort)) + geom_boxplot(lwd=1,aes(color=factor(GestationCohort),fill = NA), outlier.size = 3)  +
  stat_summary(fun=mean, geom="point", shape=5, size=8) + scale_fill_manual(values = col) + scale_colour_manual(values = col ) +
  geom_point(size=4,aes(color = factor(GestationCohort))) + xlab(NULL) + #facet_grid(rows = . ~ PostNatalAntibiotics) +
  ylab("Klebsiella \n") + paramsAngled() + scale_y_log10() + xlab(NULL)  + scale_color_tableau() +
  theme(axis.title.y = element_text(face = "italic"))
ggsave(KlebsiellaGestation, file = "KlebsiellaGestation.pdf", width = 6, height = 8, limitsize = FALSE)

DummyGenusNR<-NICUGenusNR

DummyGenusNR$Bacteroides[DummyGenusNR$Bacteroides == 0]<-1
BacteroidesMatAbx<-ggplot(DummyGenusNR,
                          aes(x=MaternalAntibiotics, y=as.numeric(Bacteroides), fill=MaternalAntibiotics)) + geom_boxplot(lwd=1,aes(color=factor(MaternalAntibiotics),fill = NA), outlier.size = 3)  +
  stat_summary(fun=mean, geom="point", shape=5, size=8) + scale_fill_manual(values = col) + scale_colour_manual(values = col ) +
  geom_point(size=4,aes(color = factor(MaternalAntibiotics))) + xlab(NULL) + facet_grid( rows = . ~ SampleCollectionWeek) +
  ylab("Bacteroides \n") + paramsBox() + scale_y_log10() + xlab(NULL)  + scale_color_tableau() +
  theme(axis.title.y = element_text(face = "italic"))
ggsave(BacteroidesMatAbx, file = "BacteroidesMatAbx.pdf", width = 8, height = 8, limitsize = FALSE)

BacteroidesMilk<-ggplot(subset(DummyGenusNR, DummyGenusNR$AnyMilk %in% c("No.Milk", "Mother")),
                        aes(x=AnyMilk, y=as.numeric(Bacteroides), fill=AnyMilk)) + geom_boxplot(lwd=1,aes(color=factor(AnyMilk),fill = NA), outlier.size = 3)  +
  stat_summary(fun=mean, geom="point", shape=5, size=8) + scale_fill_manual(values = col) + scale_colour_manual(values = col ) +
  geom_point(size=4,aes(color = factor(AnyMilk))) + xlab(NULL) + facet_grid( rows = . ~ SampleCollectionWeek) +
  ylab("Bacteroides \n") + paramsBox() + scale_y_log10() + xlab(NULL)  + scale_color_tableau() +
  theme(axis.title.y = element_text(face = "italic"))
ggsave(BacteroidesMilk, file = "BacteroidesMilk.pdf", width = 8, height = 8, limitsize = FALSE)

pairwise.wilcox.test(DummyStoolTable$Bacteroides, DummyStoolTable$AnyMilk)
pairwise.wilcox.test(DummyStoolTable$Bacteroides, DummyStoolTable$MaternalAntibiotics)
pairwise.wilcox.test(DummyStoolTable$Bacteroides, DummyStoolTable$SampleCollectionWeek)


DummyGenusNR$Lactobacillus[DummyGenusNR$Lactobacillus == 0]<-1
LactobacillusMatAbx<-ggplot(DummyGenusNR,
                            aes(x=MaternalAntibiotics, y=as.numeric(Lactobacillus), fill=MaternalAntibiotics)) + geom_boxplot(lwd=1,aes(color=factor(MaternalAntibiotics),fill = NA), outlier.size = 3)  +
  stat_summary(fun=mean, geom="point", shape=5, size=8) + scale_fill_manual(values = col) + scale_colour_manual(values = col ) +
  geom_point(size=4,aes(color = factor(MaternalAntibiotics))) + xlab(NULL) + facet_grid( rows = . ~ SampleCollectionWeek) +
  ylab("Lactobacillus \n") + paramsBox() + scale_y_log10() + xlab(NULL)  + scale_color_tableau() +
  theme(axis.title.y = element_text(face = "italic"))
ggsave(LactobacillusMatAbx, file = "LactobacillusMatAbx.pdf", width = 8, height = 8, limitsize = FALSE)

## This is week 3

StoolWeek3Data<-subset(NICUGenusNR, NICUGenusNR$SampleType == "Stool" & NICUGenusNR$SampleCollectionWeek == "Week.3")

Abx = factor(StoolWeek3Data$PostNatalAntibiotics); table(Abx)
MatAbx = factor(StoolWeek3Data$MaternalAntibiotics)
Gestation = StoolWeek3Data$GestationTime
Subject = StoolWeek3Data[, "PatientID"]; table(Subject)
Milk = factor(StoolWeek3Data$AnyMilk)
Delivery = factor(StoolWeek3Data$Delivery)

GLMGenus<-glmer(StoolWeek3Data[,19] ~ Abx+Gestation+MatAbx+Milk+Delivery+(1|Subject), data = StoolWeek3Data, family = poisson)
GLMGenusZINB <- update(GLMGenus,family=nbinom2)
pvalues<-as.data.frame(summary(GLMGenusZINB)$coeff[-1,4])
names(pvalues)<-names(StoolWeek3Data)[19]
estimates<-as.data.frame(summary(GLMGenusZINB)$coeff[-1,1])
names(estimates)<-names(StoolWeek3Data)[19]

for(i in 20:ncol(StoolWeek3Data)){
  tryCatch({
    GLMGenus<-glmer(StoolWeek3Data[,i] ~ Abx+Gestation+MatAbx+Milk+Delivery+(1|Subject), data = StoolWeek3Data, family = poisson)
    GLMGenusZINB <- update(GLMGenus,family=nbinom2)
  }, error=function(e){cat("ERROR :",conditionMessage(e), "\n")})
  pvaluesNew<-as.data.frame(summary(GLMGenusZINB)$coeff[-1,4])
  names(pvaluesNew)<-names(StoolWeek3Data)[i]
  estimatesNew<-as.data.frame(summary(GLMGenusZINB)$coeff[-1,1])
  names(estimatesNew)<-names(StoolWeek3Data)[i]
  pvalues<-cbind(pvalues, pvaluesNew)
  estimates<-cbind(estimates, estimatesNew)
  
}

tpvalues<-as.data.frame(t(pvalues))
tpvalues$AbxInfant.FDR<-p.adjust(tpvalues$AbxInfant.Abx, method = "fdr")
tpvalues$AbxMat.FDR<-p.adjust(tpvalues$MatAbxMat.Abx, method = "fdr")
tpvalues$Gestation.FDR<-p.adjust(tpvalues$Gestation, method = "fdr")
tpvalues$Milk.FDR<-p.adjust(tpvalues$MilkMilk, method = "fdr")
tpvalues$Delivery.FDR<-p.adjust(tpvalues$DeliveryVaginal, method = "fdr")



StoolWeek3AbxTable<-as.data.frame(tpvalues[,6])
row.names(StoolWeek3AbxTable)<-row.names(tpvalues)
names(StoolWeek3AbxTable)<-"AbxInfant.FDR"
SigStoolWeek3AbxTable<-subset(StoolWeek3AbxTable, StoolWeek3AbxTable$AbxInfant.FDR < 0.1)

StoolWeek3MatAbxTable<-as.data.frame(tpvalues[,7])
row.names(StoolWeek3MatAbxTable)<-row.names(tpvalues)
names(StoolWeek3MatAbxTable)<-"MatAbx.FDR"
SigStoolWeek3MatAbxTable<-subset(StoolWeek3MatAbxTable, StoolWeek3MatAbxTable$MatAbx.FDR < 0.1)

StoolWeek3GestationTable<-as.data.frame(tpvalues[,8])
row.names(StoolWeek3GestationTable)<-row.names(tpvalues)
names(StoolWeek3GestationTable)<-"Gestation.FDR"
SigStoolWeek3GestationTable<-subset(StoolWeek3GestationTable, StoolWeek3GestationTable$Gestation.FDR < 0.1)

StoolWeek3MilkTable<-as.data.frame(tpvalues[,9])
row.names(StoolWeek3MilkTable)<-row.names(tpvalues)
names(StoolWeek3MilkTable)<-"Milk.FDR"
StoolWeek3MilkTable<-subset(StoolWeek3MilkTable, StoolWeek3MilkTable$Milk.FDR < 0.1)

StoolWeek3DeliveryTable<-as.data.frame(tpvalues[,10])
row.names(StoolWeek3DeliveryTable)<-row.names(tpvalues)
names(StoolWeek3DeliveryTable)<-"Delivery.FDR"
StoolWeek3DeliveryTable<-subset(StoolWeek3DeliveryTable, StoolWeek3DeliveryTable$Delivery.FDR < 0.1)

DummyStoolTable2<-StoolWeek3Data
# DummyStoolTable2<-subset(DummyStoolTable2, DummyStoolTable2$PostNatalAntibiotics == "No.Infant.Abx")

# DummyStoolTable2$Clostridium[DummyStoolTable2$Clostridium==0]<-1
fit3 <- lm(Clostridium ~ GestationTime, data = DummyStoolTable2)
plot3<-ggplotRegression(fit3) + scale_y_log10()+ xlim(27,35)
pdf("ClostridiumGestationWeek3.pdf")
print(plot3)
dev.off()

# DummyStoolTable2$Rothia[DummyStoolTable2$Rothia==0]<-1
fit4 <- lm(Rothia ~ GestationTime, data = DummyStoolTable2)
plot4<-ggplotRegression(fit4) + scale_y_log10() + xlim(27.8,35)
pdf("RothiaGestationWeek3.pdf")
print(plot4)
dev.off()

# Ok we can do a Venn diagram of species impacted by Postnatal age, maternal Abx, Postnatal Abx, and Gestation


IntersectingAgeAbxGenus<-intersect(row.names(SigStoolWeek1Week3WeekTable), row.names(SigStoolWeek3AbxTable))
area1<-length(row.names(SigStoolWeek1Week3WeekTable))
area2<-length(row.names(SigStoolWeek3AbxTable))
area3<-length(row.names(SigStoolWeek3GestationTable))
area4<-length(row.names(SigStoolWeek3MatAbxTable))


L12=  intersect(row.names(SigStoolWeek1Week3WeekTable), row.names(SigStoolWeek3AbxTable))
L13=  intersect(row.names(SigStoolWeek1Week3WeekTable), row.names(SigStoolWeek3GestationTable))
L14=  intersect(row.names(SigStoolWeek1Week3WeekTable), row.names(SigStoolWeek3MatAbxTable))
L23=  intersect(row.names(SigStoolWeek3AbxTable), row.names(SigStoolWeek3GestationTable))
L24= intersect(row.names(SigStoolWeek3AbxTable), row.names(SigStoolWeek3MatAbxTable))
L34=  intersect(row.names(SigStoolWeek3GestationTable), row.names(SigStoolWeek3MatAbxTable))

n12=  length(intersect(row.names(SigStoolWeek1Week3WeekTable), row.names(SigStoolWeek3AbxTable)))
n13=  length(intersect(row.names(SigStoolWeek1Week3WeekTable), row.names(SigStoolWeek3GestationTable)))
n14=  length(intersect(row.names(SigStoolWeek1Week3WeekTable), row.names(SigStoolWeek3MatAbxTable)))
n23=  length(intersect(row.names(SigStoolWeek3AbxTable), row.names(SigStoolWeek3GestationTable)))
n24= length(intersect(row.names(SigStoolWeek3AbxTable), row.names(SigStoolWeek3MatAbxTable)))
n34=  length(intersect(row.names(SigStoolWeek3GestationTable), row.names(SigStoolWeek3MatAbxTable)))

L123=  intersect(L12, row.names(SigStoolWeek3GestationTable))
L124=  intersect(L24, row.names(SigStoolWeek1Week3WeekTable))
L134=  intersect(L34, row.names(SigStoolWeek1Week3WeekTable))
L234=  intersect(L23, row.names(SigStoolWeek3MatAbxTable))
L1234= intersect(row.names(SigStoolWeek1Week3WeekTable), L234)

n134= length(L134)
n123=  length(L123)
n124=  length(L124)
n234=  length(L234)
n1234= length(L1234)

write.csv(L12, file = "PostnatalAgePostnatalAbxIntersectingGenusAtWeek3.csv")

VennDiagram<-draw.pairwise.venn(area1 = area1, area2 = area2, area3 = area3, area4=area4,
                                cross.area = length(IntersectingAgeAbxGenus), category = c("Postnatal Age",  "Postnatal Antibiotics "), alpha = rep(0.2, 2),
                                label.col = "black", 
                                fontface = "plain", fontfamily = "serif", 
                                cat.just =
                                  rep(list(c(0.5, 0.5)), 2), rotation.degree = 0,
                                rotation.centre = c(0.5, 0.5), ind = TRUE, cex.prop =
                                  NULL, print.mode = "raw", sigdigs = 3, direct.area =
                                  FALSE, area.vector = 0,   fill = c("blue", "red"),
                                lty = "dashed",
                                cex = 2,
                                cat.cex = 2,
                                cat.col = c("blue", "red"))

VennDiagramAll<-draw.quad.venn(area1, area2, area3, area4, n12, n13, n14, n23, n24,
                               n34, n123, n124, n134, n234, n1234, category =c("Postnatal Age",  "Postnatal Antibiotics", "Gestational Age", "Maternal Antibiotics"),
                               alpha = rep(0.2, 4),
                               label.col = rep("black", 15), 
                               fontface = "plain", fontfamily = rep("serif",
                                                                    15), cat.pos = c(-15, 15, 0, 0), cat.dist = c(0.22,
                                                                                                                  0.22, 0.11, 0.11),  cat.fontface = rep("plain", 4),
                               cat.just =
                                 rep(list(c(0.5, 0.5)), 4), rotation.degree = 0,
                               rotation.centre = c(0.5, 0.5), ind = TRUE, cex.prop =
                                 NULL, print.mode = "raw", sigdigs = 3, direct.area =
                                 FALSE, area.vector = 0,   fill = c("blue", "red", "green", "orange"),
                               lty = "dashed",
                               cex = 2,
                               cat.cex = 2,
                               cat.col = c("blue", "red", "green", "orange"))

# ggsave(filename = "AbxAgeVennDiagram.pdf", plot = VennDiagram, width = 6,
#        height = 6, limitsize = FALSE)
# 
# ggsave(filename = "VennDiagramAll.pdf", plot = VennDiagramAll, width = 6,
#        height = 6, limitsize = FALSE)



# Do PCA for US samples: SampleType and Week

## Next do axilla, groin, and stool. (US only) to see what changes with antibiotics or from week 1 to week 3



AxillaGroinStoolGenus<-subset(NICUGenusNR, NICUGenusNR$SampleType %in% c("Axilla", "Groin", "Stool") )

# This is going to be only babies less than 31 weeks 
LessThan31WksSamples<-NICUGenusNR$Sample[NICUGenusNR$PatientID %in% GATable3$PatientID]


AxillaGroinStoolGenusNoEarly<-subset(AxillaGroinStoolGenus, AxillaGroinStoolGenus$Sample %in% LessThan31WksSamples)

metadata <- AxillaGroinStoolGenus[, 1:17]
cts <- as.matrix(AxillaGroinStoolGenus[, -(1:17)])
rownames(cts) <- metadata$Sample
cts_l2 <- glog2(cts)
grps <- AxillaGroinStoolGenus$TypeWeek

col=tableau.colors
newPCA<-PCA(cts_l2, scale.unit = FALSE, ncp = 5, graph = FALSE)
pcaPlot<-fviz_pca_ind(newPCA,
                      geom.ind = "point", 
                      pointsize = 1.5,
                      point.alph = 0.1,
                      title = "Unsupervised Principal Coordinate Analysis",
                      subtitle = "Samples Colored by Sample Type and Collection Week : All Subjects",
                      col.ind = grps, 
                      addEllipses = TRUE, 
                      ellipse.alpha = 0.2,
                      ellipse.type = "confidence",
                      ellipse.level = 0.95,
                      legend.title = "Sample", 
                      legend.size = 11,
                      mean.point = TRUE,
                      palette = col,
                      axes.linetype = "blank"
)
SampleTypeWeekPCA<- pcaPlot + xlim(-23,36) + ylim(-24,33)

pdf("SampleTypeGenusWeekPCAAllSubjects.pdf") 
print(SampleTypeWeekPCA)
dev.off() 

NewAxillaGroinStoolGenus<-AxillaGroinStoolGenus
NewAxillaGroinStoolGenus<-subset(NewAxillaGroinStoolGenus, NewAxillaGroinStoolGenus$SampleType %in% c("Stool", "Groin"))
NewAxillaGroinStoolGenus$TypeWeekMatAbx<-paste(NewAxillaGroinStoolGenus$MaternalAntibiotics,NewAxillaGroinStoolGenus$TypeWeek , sep = "-")
TypeWeekMatAbxCol<-grep("TypeWeekMatAbx", names(NewAxillaGroinStoolGenus))
metadata <- NewAxillaGroinStoolGenus[, 1:17, TypeWeekMatAbxCol]
cts <- as.matrix(NewAxillaGroinStoolGenus[, -c(1:17,TypeWeekMatAbxCol)])
rownames(cts) <- metadata$Sample
cts_l2 <- glog2(cts)
grps <- NewAxillaGroinStoolGenus$TypeWeekMatAbx

col=rep(tableau.colors,2)
newPCA<-PCA(cts_l2, scale.unit = FALSE, ncp = 5, graph = FALSE)
pcaPlot<-fviz_pca_ind(newPCA,
                      geom.ind = "point", 
                      pointsize = 1.5,
                      point.alph = 0.1,
                      title = "Unsupervised Principal Coordinate Analysis",
                      subtitle = "Samples Colored by Sample Type and Collection Week : All Subjects",
                      col.ind = grps, 
                      addEllipses = TRUE, 
                      ellipse.alpha = 0.2,
                      ellipse.type = "confidence",
                      ellipse.level = 0.95,
                      legend.title = "Sample", 
                      legend.size = 11,
                      mean.point = TRUE,
                      palette = col,
                      axes.linetype = "blank"
)
SampleTypeWeekPCAMatAbx<- pcaPlot 

pdf("MatAbxSampleType.pdf") 
print(SampleTypeWeekPCAMatAbx)
dev.off() 

#MRPP of maternal antibiotics (this includes all sample types and times)



#this is < 31 weeks infants
metadata <- AxillaGroinStoolGenusNoEarly[, 1:17]
cts <- as.matrix(AxillaGroinStoolGenusNoEarly[, -(1:17)])
rownames(cts) <- metadata$Sample
cts_l2 <- glog2(cts)
grps <- AxillaGroinStoolGenusNoEarly$TypeWeek

col=tableau.colors
newPCA<-PCA(cts_l2, scale.unit = FALSE, ncp = 5, graph = FALSE)
pcaPlot<-fviz_pca_ind(newPCA,
                      geom.ind = "point",
                      pointsize = 1.5,
                      point.alph = 0.1,
                      title = "Unsupervised Principal Coordinate Analysis",
                      subtitle = "Samples Colored by Sample Type and Collection Week : All Subjects",
                      col.ind = grps,
                      addEllipses = TRUE,
                      ellipse.alpha = 0.2,
                      ellipse.type = "confidence",
                      ellipse.level = 0.95,
                      legend.title = "Sample",
                      legend.size = 11,
                      mean.point = TRUE,
                      palette = col,
                      axes.linetype = "blank"
)
SampleTypeWeekPCA2<- pcaPlot + xlim(-23,36) + ylim(-24,33)

pdf("SampleTypeGenusWeekPCAAllSubjectsLessThan31.pdf")
print(SampleTypeWeekPCA2)
dev.off()


AbxSampleIndAllWks<-AxillaGroinStoolGenus$Sample[AxillaGroinStoolGenus$PostNatalAntibiotics == "Infant.Abx"]
NoAbxSampleIndAllWks<-AxillaGroinStoolGenus$Sample[AxillaGroinStoolGenus$PostNatalAntibiotics == "No.Infant.Abx"]

AbxSampleIndUnder31Wks<-AxillaGroinStoolGenusNoEarly$Sample[AxillaGroinStoolGenusNoEarly$PostNatalAntibiotics == "Infant.Abx"]
NoAbxSampleIndUnder31Wks<-AxillaGroinStoolGenusNoEarly$Sample[AxillaGroinStoolGenusNoEarly$PostNatalAntibiotics == "No.Infant.Abx"]


col=tableau.colors
# newPCA<-PCA(cts_l2, scale.unit = FALSE, ncp = 5, graph = FALSE)
pcaPlot<-fviz_pca_ind(newPCA,
                      geom.ind = "point", 
                      pointsize = 1.5,
                      point.alph = 0.1,
                      title = "Unsupervised Principal Coordinate Analysis",
                      subtitle = "Samples Colored by Sample Type and Collection Week : No Antibiotics",
                      col.ind = grps, 
                      addEllipses = TRUE, 
                      ellipse.alpha = 0.2,
                      ellipse.type = "confidence",
                      ellipse.level = 0.95,
                      legend.title = "Sample", 
                      legend.size = 11,
                      mean.point = TRUE,
                      palette = col,
                      axes.linetype = "blank",
                      select.ind = list(name = NoAbxSampleIndAllWks)
)
SampleTypeWeekPCANOAbx<-  pcaPlot + xlim(-23,36) + ylim(-24,33)

pdf("SampleTypeGenusWeekPCANoAbx.pdf") 
print(SampleTypeWeekPCANOAbx)
dev.off() 

pcaPlot<-fviz_pca_ind(newPCA,
                      geom.ind = "point", 
                      pointsize = 1.5,
                      point.alph = 0.1,
                      title = "Unsupervised Principal Coordinate Analysis",
                      subtitle = "Samples Colored by Sample Type and Collection Week : Antibiotics Exposed",
                      col.ind = grps, 
                      addEllipses = TRUE, 
                      ellipse.alpha = 0.2,
                      ellipse.type = "confidence",
                      ellipse.level = 0.95,
                      legend.title = "Sample", 
                      legend.size = 11,
                      mean.point = TRUE,
                      palette = col,
                      axes.linetype = "blank",
                      select.ind = list(name = AbxSampleIndAllWks)
)
SampleTypeWeekPCAAbx<- pcaPlot + xlim(-23,36) + ylim(-24,33)

pdf("SampleTypeGenusWeekPCAAbx.pdf") 
print(SampleTypeWeekPCAAbx)
dev.off() 

# Do PCA for Gestation age in non-antibiotic treated group # this includes all babies and going to make two week cohorts

AxillaGroinStoolGenusGA<-AxillaGroinStoolGenus
AxillaGroinStoolGenusGA$GA2Wk<-ifelse(AxillaGroinStoolGenusGA$GestationTime >= 24 & AxillaGroinStoolGenusGA$GestationTime <= 25.99, "24-26", 
                                      ifelse(AxillaGroinStoolGenusGA$GestationTime >= 26 & AxillaGroinStoolGenusGA$GestationTime <=27.99, "26-28",
                                             ifelse(AxillaGroinStoolGenusGA$GestationTime >= 28 & AxillaGroinStoolGenusGA$GestationTime <= 29.99, "28-30", 
                                                    ifelse(AxillaGroinStoolGenusGA$GestationTime >= 30 & AxillaGroinStoolGenusGA$GestationTime <=31.99, "30-32",
                                                           ifelse(AxillaGroinStoolGenusGA$GestationTime >= 30 & AxillaGroinStoolGenusGA$GestationTime <= 31.99, "30-32", 
                                                                  ifelse(AxillaGroinStoolGenusGA$GestationTime >= 32 & AxillaGroinStoolGenusGA$GestationTime <= 33.99, "32-34",
                                                                         ifelse(AxillaGroinStoolGenusGA$GestationTime >= 34 & AxillaGroinStoolGenusGA$GestationTime <= 35.99, "34-36", 
                                                                                ifelse(AxillaGroinStoolGenusGA$GestationTime >= 36 & AxillaGroinStoolGenusGA$GestationTime <= 37.99, "36-38", "NA"))))))))                                                  ))))
AxillaGroinStoolGenusGA<-subset(AxillaGroinStoolGenusGA, ! is.na(AxillaGroinStoolGenusGA$GA2Wk))

AxillaGroinStoolGenusGA$AgeTypeWeek<-paste(AxillaGroinStoolGenusGA$SampleType, AxillaGroinStoolGenusGA$SampleCollectionWeek,AxillaGroinStoolGenusGA$GA2Wk,  sep = "-")


AxillaGroinStoolGenusGA<-AxillaGroinStoolGenusGA[,c(1,78,79,2:ncol(AxillaGroinStoolGenusGA)-2)]
AxillaGroinStoolGenusGANoAbx<-subset(AxillaGroinStoolGenusGA, AxillaGroinStoolGenusGA$PostNatalAntibiotics == "No.Infant.Abx")

AxillaGroinStoolGenusGANoAbx<-subset(AxillaGroinStoolGenusGANoAbx, AxillaGroinStoolGenusGANoAbx$SampleType == "Stool")

metadata <- AxillaGroinStoolGenusGANoAbx[, 1:20]
cts <- as.matrix(AxillaGroinStoolGenusGANoAbx[, -(1:20)])
rownames(cts) <- metadata$Sample
cts_l2 <- glog2(cts)
grps <- AxillaGroinStoolGenusGANoAbx$AgeTypeWeek

Fig1IndvList<-which(AxillaGroinStoolGenusGANoAbx$AgeTypeWeek %in% c("Stool-Week.1-28-30", "Stool-Week.1-30-32", "Stool-Week.3-28-30"))
Fig1Indv<-AxillaGroinStoolGenusGANoAbx$Sample[Fig1IndvList]

Week1IndNew<-subset(AxillaGroinStoolGenusGANoAbx$Sample, AxillaGroinStoolGenusGANoAbx$SampleCollectionWeek == "Week.1")
Week3IndNew<-subset(AxillaGroinStoolGenusGANoAbx$Sample, AxillaGroinStoolGenusGANoAbx$SampleCollectionWeek == "Week.3")

col=rep(tableau.colors,2)
newPCA<-PCA(cts_l2, scale.unit = FALSE, ncp = 5, graph = FALSE)
pcaPlot<-fviz_pca_ind(newPCA,
                      geom.ind = "point", 
                      pointsize = 1.5,
                      point.alph = 0.1,
                      title = "Unsupervised Principal Coordinate Analysis",
                      subtitle = "Samples Colored by Sample Type, Gestational Age and Collection Week : No Antibiotics",
                      col.ind = grps, 
                      addEllipses = TRUE, 
                      ellipse.alpha = 0.2,
                      ellipse.type = "confidence",
                      ellipse.level = 0.95,
                      legend.title = "Sample", 
                      legend.size = 11,
                      mean.point = TRUE,
                      palette = col,
                      axes.linetype = "blank"
)

SampleTypeWeek1GestationPCANOAbx<-  pcaPlot + xlim(-34,20) + ylim(-20,24)

pdf("SampleTypeWeek13GestationPCANOAbx.pdf") 
print(SampleTypeWeek1GestationPCANOAbx)
dev.off() 

## Alright split this into Week 1 and Week 3 for Fig 5A and 5B

col=rep(tableau.colors,2)
col=col[1:4]
# newPCA<-PCA(cts_l2, scale.unit = FALSE, ncp = 5, graph = FALSE)
pcaPlot<-fviz_pca_ind(newPCA,
                      geom.ind = "point", 
                      pointsize = 1.5,
                      point.alph = 0.1,
                      title = "Unsupervised Principal Coordinate Analysis",
                      subtitle = "Samples Colored by Sample Type, Gestational Age and Collection Week : No Antibiotics",
                      col.ind = grps, 
                      addEllipses = TRUE, 
                      ellipse.alpha = 0.2,
                      ellipse.type = "confidence",
                      ellipse.level = 0.95,
                      legend.title = "Sample", 
                      legend.size = 11,
                      mean.point = TRUE,
                      palette = col,
                      axes.linetype = "blank",
                      select.ind = list(names=Week1IndNew)
)

Week1SampleTypeWeek1GestationPCANOAbx<-  pcaPlot + xlim(-34,20) + ylim(-20,24)

pdf("Week1SampleTypeWeek1GestationPCANOAbx.pdf") 
print(Week1SampleTypeWeek1GestationPCANOAbx)
dev.off() 

col=rep(tableau.colors,2)
col=col[5:8]
# newPCA<-PCA(cts_l2, scale.unit = FALSE, ncp = 5, graph = FALSE)
pcaPlot<-fviz_pca_ind(newPCA,
                      geom.ind = "point", 
                      pointsize = 1.5,
                      point.alph = 0.1,
                      title = "Unsupervised Principal Coordinate Analysis",
                      subtitle = "Samples Colored by Sample Type, Gestational Age and Collection Week : No Antibiotics",
                      col.ind = grps, 
                      addEllipses = TRUE, 
                      ellipse.alpha = 0.2,
                      ellipse.type = "confidence",
                      ellipse.level = 0.95,
                      legend.title = "Sample", 
                      legend.size = 11,
                      mean.point = TRUE,
                      palette = col,
                      axes.linetype = "blank",
                      select.ind = list(names=Week3IndNew)
)

Week3SampleTypeWeek3GestationPCANOAbx<-  pcaPlot + xlim(-34,20) + ylim(-20,24)

pdf("Week3SampleTypeWeek3GestationPCANOAbx.pdf") 
print(Week3SampleTypeWeek3GestationPCANOAbx)
dev.off() 

# Now split up the gestational ages

New28to30Ind<-subset(AxillaGroinStoolGenusGANoAbx$Sample, AxillaGroinStoolGenusGANoAbx$GestationTime  >= 28 & AxillaGroinStoolGenusGANoAbx$GestationTime <= 29.99)
New30to32Ind<-subset(AxillaGroinStoolGenusGANoAbx$Sample, AxillaGroinStoolGenusGANoAbx$GestationTime  >= 30 & AxillaGroinStoolGenusGANoAbx$GestationTime <= 31.99)
New32to34Ind<-subset(AxillaGroinStoolGenusGANoAbx$Sample, AxillaGroinStoolGenusGANoAbx$GestationTime  >= 32 & AxillaGroinStoolGenusGANoAbx$GestationTime <= 33.99)
New34to36Ind<-subset(AxillaGroinStoolGenusGANoAbx$Sample, AxillaGroinStoolGenusGANoAbx$GestationTime  >= 34 & AxillaGroinStoolGenusGANoAbx$GestationTime <= 35.99)

col=rep(tableau.colors,2)
col=col[c(1,5)]
# newPCA<-PCA(cts_l2, scale.unit = FALSE, ncp = 5, graph = FALSE)
pcaPlot<-fviz_pca_ind(newPCA,
                      geom.ind = "point", 
                      pointsize = 1.5,
                      point.alph = 0.1,
                      title = "Unsupervised Principal Coordinate Analysis",
                      subtitle = "Samples Colored by Sample Type, Gestational Age and Collection Week : No Antibiotics",
                      col.ind = grps, 
                      addEllipses = TRUE, 
                      ellipse.alpha = 0.2,
                      ellipse.type = "confidence",
                      ellipse.level = 0.95,
                      legend.title = "Sample", 
                      legend.size = 11,
                      mean.point = TRUE,
                      palette = col,
                      axes.linetype = "blank",
                      select.ind = list(names=)
)

SampleTypeWeek28to30PCANOAbx<-  pcaPlot + xlim(-34,20) + ylim(-20,24)

pdf("SampleTypeWeek28to30PCANOAbx.pdf") 
print(SampleTypeWeek28to30PCANOAbx)
dev.off() 

col=rep(tableau.colors,2)
col=col[c(2,6)]
# newPCA<-PCA(cts_l2, scale.unit = FALSE, ncp = 5, graph = FALSE)
pcaPlot<-fviz_pca_ind(newPCA,
                      geom.ind = "point", 
                      pointsize = 1.5,
                      point.alph = 0.1,
                      title = "Unsupervised Principal Coordinate Analysis",
                      subtitle = "Samples Colored by Sample Type, Gestational Age and Collection Week : No Antibiotics",
                      col.ind = grps, 
                      addEllipses = TRUE, 
                      ellipse.alpha = 0.2,
                      ellipse.type = "confidence",
                      ellipse.level = 0.95,
                      legend.title = "Sample", 
                      legend.size = 11,
                      mean.point = TRUE,
                      palette = col,
                      axes.linetype = "blank",
                      select.ind = list(names=New30to32Ind)
)

SampleTypeWeek30to32PCANOAbx<-  pcaPlot + xlim(-34,20) + ylim(-20,24)

pdf("SampleTypeWeek30to32PCANOAbx.pdf") 
print(SampleTypeWeek30to32PCANOAbx)
dev.off() 

col=rep(tableau.colors,2)
col=col[c(3,7)]
# newPCA<-PCA(cts_l2, scale.unit = FALSE, ncp = 5, graph = FALSE)
pcaPlot<-fviz_pca_ind(newPCA,
                      geom.ind = "point", 
                      pointsize = 1.5,
                      point.alph = 0.1,
                      title = "Unsupervised Principal Coordinate Analysis",
                      subtitle = "Samples Colored by Sample Type, Gestational Age and Collection Week : No Antibiotics",
                      col.ind = grps, 
                      addEllipses = TRUE, 
                      ellipse.alpha = 0.2,
                      ellipse.type = "confidence",
                      ellipse.level = 0.95,
                      legend.title = "Sample", 
                      legend.size = 11,
                      mean.point = TRUE,
                      palette = col,
                      axes.linetype = "blank",
                      select.ind = list(names=New32to34Ind)
)

SampleTypeWeek32to34PCANOAbx<-  pcaPlot + xlim(-34,20) + ylim(-20,24)

pdf("SampleTypeWeek32to34PCANOAbx.pdf") 
print(SampleTypeWeek32to34PCANOAbx)
dev.off() 

col=rep(tableau.colors,2)
col=col[c(4,8)]
# newPCA<-PCA(cts_l2, scale.unit = FALSE, ncp = 5, graph = FALSE)
pcaPlot<-fviz_pca_ind(newPCA,
                      geom.ind = "point", 
                      pointsize = 1.5,
                      point.alph = 0.1,
                      title = "Unsupervised Principal Coordinate Analysis",
                      subtitle = "Samples Colored by Sample Type, Gestational Age and Collection Week : No Antibiotics",
                      col.ind = grps, 
                      addEllipses = TRUE, 
                      ellipse.alpha = 0.2,
                      ellipse.type = "confidence",
                      ellipse.level = 0.95,
                      legend.title = "Sample", 
                      legend.size = 11,
                      mean.point = TRUE,
                      palette = col,
                      axes.linetype = "blank",
                      select.ind = list(names=New34to36Ind)
)

SampleTypeWeek34to36PCANOAbx<-  pcaPlot + xlim(-34,20) + ylim(-20,24)

pdf("SampleTypeWeek34to36PCANOAbx.pdf") 
print(SampleTypeWeek34to36PCANOAbx)
dev.off() 


pcaBiplot<-fviz_pca_biplot(
  newPCA,
  axes = c(1, 2),
  geom = c("point", "text"),
  geom.ind = "point",
  geom.var = c("arrow", "text"),
  col.ind = grps,
  addEllipses = TRUE, 
  ellipse.alpha = 0.2,
  ellipse.type = "confidence",
  ellipse.level = 0.95,
  legend.title = "Sample", 
  pointsize = 0.5,
  point.alph = 0.1,
  legend.size = 11,
  label = "all",
  invisible = "none",
  repel = TRUE,
  habillage = "none",
  palette = col,
  title = "Unsupervised Principal Coordinate Analysis",
  subtitle = "Samples Colored by Sample Type, Gestational Age and Collection Week : No Antibiotics",
  axes.linetype = "blank",
  select.var = list(name = NULL, cos2 = 40)
)

SampleTypeWeek1GestationBiplotNOAbx<-  pcaBiplot #+ xlim(-30,70) + ylim(-20,70)

pdf("SampleTypeWeek13GestationBiplotNOAbx.pdf") 
print(SampleTypeWeek1GestationBiplotNOAbx)
dev.off() 

pcaBiplot<-fviz_pca_var(
  newPCA,
  axes = c(1, 2),
  geom = c("point", "text"),
  geom.ind = "point",
  geom.var = c("arrow", "text"),
  col.ind = grps,
  addEllipses = TRUE, 
  ellipse.alpha = 0.2,
  ellipse.type = "confidence",
  ellipse.level = 0.95,
  legend.title = "Sample", 
  legend.size = 11,
  label = "all",
  invisible = "none",
  repel = TRUE,
  habillage = "none",
  palette = col,
  title = "Unsupervised Principal Coordinate Analysis",
  subtitle = "Samples Colored by Sample Type, Gestational Age and Collection Week : No Antibiotics",
  axes.linetype = "blank"
)




col=tableau.colors[c(1,2,5)]
pcaPlot<-fviz_pca_ind(newPCA,
                      geom.ind = "point", 
                      pointsize = 1.5,
                      point.alph = 0.1,
                      title = "Unsupervised Principal Coordinate Analysis",
                      subtitle = "Samples Colored by Sample Type, Gestational Age and Collection Week : No Antibiotics",
                      col.ind = grps, 
                      addEllipses = TRUE, 
                      ellipse.alpha = 0.2,
                      ellipse.type = "confidence",
                      ellipse.level = 0.95,
                      legend.title = "Sample", 
                      legend.size = 11,
                      mean.point = TRUE,
                      palette = col,
                      axes.linetype = "blank",
                      select.ind = list(name = Fig1Indv)
)

SampleTypeWeekFig1GestationPCANOAbx<- pcaPlot + xlim(-34,20) + ylim(-20,24)

pdf("SampleTypeWeekFig1GestationPCANOAbx.pdf") 
print(SampleTypeWeekFig1GestationPCANOAbx)
dev.off() 

Fig2IndvList<-which(AxillaGroinStoolGenusGANoAbx$AgeTypeWeek %in% c("Stool-Week.1-30-32", "Stool-Week.1-32-34", "Stool-Week.3-30-32"))
Fig2Indv<-AxillaGroinStoolGenusGANoAbx$Sample[Fig2IndvList]

col=tableau.colors[c(2,4,6)]
pcaPlot<-fviz_pca_ind(newPCA,
                      geom.ind = "point", 
                      pointsize = 1.5,
                      point.alph = 0.1,
                      title = "Unsupervised Principal Coordinate Analysis",
                      subtitle = "Samples Colored by Sample Type, Gestational Age and Collection Week : No Antibiotics",
                      col.ind = grps, 
                      addEllipses = TRUE, 
                      ellipse.alpha = 0.2,
                      ellipse.type = "confidence",
                      ellipse.level = 0.95,
                      legend.title = "Sample", 
                      legend.size = 11,
                      mean.point = TRUE,
                      palette = col,
                      axes.linetype = "blank",
                      select.ind = list(name = Fig2Indv)
)

SampleTypeWeekFig2GestationPCANOAbx<- pcaPlot + xlim(-34,20) + ylim(-20,24)

pdf("SampleTypeWeekFig2GestationPCANOAbx.pdf") 
print(SampleTypeWeekFig2GestationPCANOAbx)
dev.off() 

Fig3IndvList<-which(AxillaGroinStoolGenusGANoAbx$AgeTypeWeek %in% c("Stool-Week.1-32-34", "Stool-Week.1-34-36", "Stool-Week.3-32-34"))
Fig3Indv<-AxillaGroinStoolGenusGANoAbx$Sample[Fig3IndvList]

col=tableau.colors[c(4,8,7)]
pcaPlot<-fviz_pca_ind(newPCA,
                      geom.ind = "point", 
                      pointsize = 1.5,
                      point.alph = 0.1,
                      title = "Unsupervised Principal Coordinate Analysis",
                      subtitle = "Samples Colored by Sample Type, Gestational Age and Collection Week : No Antibiotics",
                      col.ind = grps, 
                      addEllipses = TRUE, 
                      ellipse.alpha = 0.2,
                      ellipse.type = "confidence",
                      ellipse.level = 0.95,
                      legend.title = "Sample", 
                      legend.size = 11,
                      mean.point = TRUE,
                      palette = col,
                      axes.linetype = "blank",
                      select.ind = list(name = Fig3Indv)
)

SampleTypeWeekFig3GestationPCANOAbx<- pcaPlot + xlim(-34,20) + ylim(-20,24)

pdf("SampleTypeWeekFig3GestationPCANOAbx.pdf") 
print(SampleTypeWeekFig3GestationPCANOAbx)
dev.off() 

Fig4IndvList<-which(AxillaGroinStoolGenusGANoAbx$AgeTypeWeek %in% c("Stool-Week.1-34-36", "Stool-Week.1-36-38", "Stool-Week.3-34-36"))
Fig4Indv<-AxillaGroinStoolGenusGANoAbx$Sample[Fig4IndvList]

col=tableau.colors[c(8,3)]
pcaPlot<-fviz_pca_ind(newPCA,
                      geom.ind = "point", 
                      pointsize = 1.5,
                      point.alph = 0.1,
                      title = "Unsupervised Principal Coordinate Analysis",
                      subtitle = "Samples Colored by Sample Type, Gestational Age and Collection Week : No Antibiotics",
                      col.ind = grps, 
                      addEllipses = TRUE, 
                      ellipse.alpha = 0.2,
                      ellipse.type = "confidence",
                      ellipse.level = 0.95,
                      legend.title = "Sample", 
                      legend.size = 11,
                      mean.point = TRUE,
                      palette = col,
                      axes.linetype = "blank",
                      select.ind = list(name = Fig4Indv)
)

SampleTypeWeekFig4GestationPCANOAbx<- pcaPlot + xlim(-34,20) + ylim(-20,24)

pdf("SampleTypeWeekFig4GestationPCANOAbx.pdf") 
print(SampleTypeWeekFig4GestationPCANOAbx)
dev.off() 


# ok get MRPP for stool samples at various gestational age:


New28to30IndforMRPP<-subset(AxillaGroinStoolGenusGANoAbx, AxillaGroinStoolGenusGANoAbx$GestationTime >= 28 & AxillaGroinStoolGenusGANoAbx$GestationTime <= 29.99 & AxillaGroinStoolGenusGANoAbx$SampleType == "Stool")
New30to32IndforMRPP<-subset(AxillaGroinStoolGenusGANoAbx, AxillaGroinStoolGenusGANoAbx$GestationTime >= 30 & AxillaGroinStoolGenusGANoAbx$GestationTime <= 31.99 & AxillaGroinStoolGenusGANoAbx$SampleType == "Stool")
New32to34IndforMRPP<-subset(AxillaGroinStoolGenusGANoAbx, AxillaGroinStoolGenusGANoAbx$GestationTime >= 32 & AxillaGroinStoolGenusGANoAbx$GestationTime <= 33.99 & AxillaGroinStoolGenusGANoAbx$SampleType == "Stool")
New34to36IndforMRPP<-subset(AxillaGroinStoolGenusGANoAbx, AxillaGroinStoolGenusGANoAbx$GestationTime >= 34 & AxillaGroinStoolGenusGANoAbx$GestationTime <= 35.99 & AxillaGroinStoolGenusGANoAbx$SampleType == "Stool")

MRPP28to30.mrpp<-mrpp(New28to30IndforMRPP[,21:ncol(New28to30IndforMRPP)], New28to30IndforMRPP$SampleCollectionWeek, distance = "bray")
MRPP30to32.mrpp<-mrpp(New30to32IndforMRPP[,21:ncol(New30to32IndforMRPP)], New30to32IndforMRPP$SampleCollectionWeek, distance = "bray")
MRPP32to34.mrpp<-mrpp(New32to34IndforMRPP[,21:ncol(New32to34IndforMRPP)], New32to34IndforMRPP$SampleCollectionWeek, distance = "bray")
MRPP34to36.mrpp<-mrpp(New34to36IndforMRPP[,21:ncol(New34to36IndforMRPP)], New34to36IndforMRPP$SampleCollectionWeek, distance = "bray")

## Alright now the p values for gestational age at week 1 and week 3
Week1NoAbxSamples28vs30<-subset(AxillaGroinStoolGenusGANoAbx, AxillaGroinStoolGenusGANoAbx$GA2Wk %in% c("28-30", "30-32") & AxillaGroinStoolGenusGANoAbx$SampleType == "Stool" & AxillaGroinStoolGenusGANoAbx$SampleCollectionWeek == "Week.1")
Week1NoAbxSamples28vs30.mrpp<-mrpp(Week1NoAbxSamples28vs30[,21:ncol(Week1NoAbxSamples28vs30)], Week1NoAbxSamples28vs30$GA2Wk, distance = "bray")
Week1NoAbxSamples30vs32<-subset(AxillaGroinStoolGenusGANoAbx, AxillaGroinStoolGenusGANoAbx$GA2Wk %in% c("32-34", "30-32") & AxillaGroinStoolGenusGANoAbx$SampleType == "Stool" & AxillaGroinStoolGenusGANoAbx$SampleCollectionWeek == "Week.1")
Week1NoAbxSamples30vs32.mrpp<-mrpp(Week1NoAbxSamples30vs32[,21:ncol(Week1NoAbxSamples30vs32)], Week1NoAbxSamples30vs32$GA2Wk, distance = "bray")
Week1NoAbxSamples32vs34<-subset(AxillaGroinStoolGenusGANoAbx, AxillaGroinStoolGenusGANoAbx$GA2Wk %in% c("32-34", "34-36") & AxillaGroinStoolGenusGANoAbx$SampleType == "Stool" & AxillaGroinStoolGenusGANoAbx$SampleCollectionWeek == "Week.1")
Week1NoAbxSamples32vs34.mrpp<-mrpp(Week1NoAbxSamples32vs34[,21:ncol(Week1NoAbxSamples32vs34)], Week1NoAbxSamples32vs34$GA2Wk, distance = "bray")


Week3NoAbxSamples28vs30<-subset(AxillaGroinStoolGenusGANoAbx, AxillaGroinStoolGenusGANoAbx$GA2Wk %in% c("28-30", "30-32") & AxillaGroinStoolGenusGANoAbx$SampleType == "Stool" & AxillaGroinStoolGenusGANoAbx$SampleCollectionWeek == "Week.3")
Week3NoAbxSamples28vs30.mrpp<-mrpp(Week3NoAbxSamples28vs30[,21:ncol(Week3NoAbxSamples28vs30)], Week3NoAbxSamples28vs30$GA2Wk, distance = "bray")
Week3NoAbxSamples30vs32<-subset(AxillaGroinStoolGenusGANoAbx, AxillaGroinStoolGenusGANoAbx$GA2Wk %in% c("32-34", "30-32") & AxillaGroinStoolGenusGANoAbx$SampleType == "Stool" & AxillaGroinStoolGenusGANoAbx$SampleCollectionWeek == "Week.3")
Week3NoAbxSamples30vs32.mrpp<-mrpp(Week3NoAbxSamples30vs32[,21:ncol(Week3NoAbxSamples30vs32)], Week3NoAbxSamples30vs32$GA2Wk, distance = "bray")
Week3NoAbxSamples32vs34<-subset(AxillaGroinStoolGenusGANoAbx, AxillaGroinStoolGenusGANoAbx$GA2Wk %in% c("32-34", "34-36") & AxillaGroinStoolGenusGANoAbx$SampleType == "Stool" & AxillaGroinStoolGenusGANoAbx$SampleCollectionWeek == "Week.3")
Week3NoAbxSamples32vs34.mrpp<-mrpp(Week3NoAbxSamples32vs34[,21:ncol(Week3NoAbxSamples32vs34)], Week3NoAbxSamples32vs34$GA2Wk, distance = "bray")

### ok now I'm going to try to calculate the Bray curtis distance between these groups:
BCMatrix<-vegdist(AxillaGroinStoolGenusGA[,21:ncol(AxillaGroinStoolGenusGA)], method = "bray")

BCMatrix<-as.matrix(BCMatrix)
BCMatrix<-as.data.frame(BCMatrix)
row.names(BCMatrix)<-AxillaGroinStoolGenusGA$Sample
BCDistanceTable<-as.data.frame(BCMatrix)
names(BCDistanceTable)<-AxillaGroinStoolGenusGA$Sample
BCDistanceTable$Sample<-row.names(BCDistanceTable)
BCDistanceTable<-melt(BCDistanceTable, id.vars = "Sample")
names(BCDistanceTable)<-c("Sample", "SampleID", "BCDistance")
BCDistanceTable<-merge(BCDistanceTable, AxillaGroinStoolGenusGA[,c(1,3,12)], by = "Sample", all.x = TRUE)
names(BCDistanceTable)<-c("Sample", "SampleID", "BCDistance", "SampleGroup", "SampleAbx")
BCDistanceTable<-merge(BCDistanceTable, AxillaGroinStoolGenusGA[,c(6,3,12)], by = "SampleID", all.x = TRUE)
names(BCDistanceTable)<-c("Sample", "SampleID", "BCDistance", "SampleGroup", "SampleAbx", "ComparitorGroup", "ComparitorAbx")
BCDistanceTable$ComparisonGroups<-paste(BCDistanceTable$SampleGroup, BCDistanceTable$SampleAbx, BCDistanceTable$ComparitorGroup, BCDistanceTable$ComparitorAbx, sep = "-")

MeanDistanceTable<-BCDistanceTable %>%
  group_by(ComparisonGroups) %>%
  summarize(mean(BCDistance))

names(MeanDistanceTable)<-c("ComparisonGroups", "MeanDistance")
MeanDistanceTable<-MeanDistanceTable[order(MeanDistanceTable$MeanDistance, decreasing = FALSE),]
write.csv(MeanDistanceTable, file = "MeanDistanceTable.csv")

stoolRows<-grep("Stool", BCDistanceTable$SampleGroup)
BCStoolDistanceTable<-BCDistanceTable[stoolRows,]
stoolRows2<-grep("Stool", BCStoolDistanceTable$ComparitorGroup)
BCStoolDistanceTable<-BCStoolDistanceTable[stoolRows2,]

write.csv(BCStoolDistanceTable, file = "BCDistanceTableStoolSamples.csv")

MeanStoolDistanceTable<-BCStoolDistanceTable %>%
  group_by(ComparisonGroups) %>%
  summarize(mean(BCDistance))


MeanStoolDistanceTable$ComparisonGroups<-gsub("Stool-", "", MeanStoolDistanceTable$ComparisonGroups)
MeanStoolDistanceTable$SampleWeek<-sapply(strsplit(MeanStoolDistanceTable$ComparisonGroups, "-"),`[`, 1)
MeanStoolDistanceTable$ComparitorWeek<-sapply(strsplit(MeanStoolDistanceTable$ComparisonGroups, "-"),`[`, 5)
MeanStoolDistanceTable$SampleAbx<-sapply(strsplit(MeanStoolDistanceTable$ComparisonGroups, "-"),`[`, 4)
MeanStoolDistanceTable$ComparitorAbx<-sapply(strsplit(MeanStoolDistanceTable$ComparisonGroups, "-"),`[`, 8)

write.csv(MeanStoolDistanceTable, file = "MeanStoolDistanceTable.csv")


names(MeanStoolDistanceTable)[2]<-"BCDistance"
col<-rep(tableau.colors,2)
NoAbxBCDistance<-ggplot(subset(BCDistanceTable, BCDistanceTable$SampleGroup %in% c("Stool-Week.1-28-30", "Stool-Week.1-30-32", "Stool-Week.1-32-34", "Stool-Week.1-34-36") & BCDistanceTable$ComparitorGroup %in% c("Stool-Week.1-28-30", "Stool-Week.1-30-32", "Stool-Week.1-32-34", "Stool-Week.1-34-36")
                               &   BCDistanceTable$SampleAbx == "No.Infant.Abx"  & BCDistanceTable$ComparitorAbx == "No.Infant.Abx" ),
                        aes(x=ComparisonGroups, y=as.numeric(BCDistance), fill=ComparisonGroups)) + geom_boxplot(lwd=1,aes(color=factor(ComparisonGroups),fill = NA), outlier.size = 3)  +
  stat_summary(fun=mean, geom="point", shape=5, size=8) + scale_fill_manual(values = col) + scale_colour_manual(values = col ) +
  geom_point(size=4,aes(color = factor(ComparisonGroups))) + xlab(NULL)  +
  ylab("Bray Curtis Distance \n") + paramsAngled()  + xlab(NULL)  + scale_color_tableau()
ggsave(NoAbxBCDistance, file = "NoAbxBCDistance.pdf", width = 12, height = 8, limitsize = FALSE)


NoAbxDistanceTable<-subset(BCDistanceTable, BCDistanceTable$SampleAbx == "No.Infant.Abx"  & BCDistanceTable$ComparitorAbx == "No.Infant.Abx")
NoAbxStoolDistanceTable<-subset(NoAbxDistanceTable,
                                NoAbxDistanceTable$SampleGroup %in% c("Stool-Week.1-28-30", "Stool-Week.1-30-32", "Stool-Week.1-32-34", "Stool-Week.1-34-36", "Stool-Week.3-28-30", "Stool-Week.3-30-32", "Stool-Week.3-32-34", "Stool-Week.3-34-36") & 
                                  NoAbxDistanceTable$ComparitorGroup %in% c("Stool-Week.1-28-30", "Stool-Week.1-30-32", "Stool-Week.1-32-34", "Stool-Week.1-34-36", "Stool-Week.3-28-30", "Stool-Week.3-30-32", "Stool-Week.3-32-34", "Stool-Week.3-34-36")) 


col<-"blue"
AllAbxBCDistance<-ggplot(NoAbxStoolDistanceTable,
                         aes(x=ComparisonGroups, y=as.numeric(BCDistance), fill=ComparisonGroups)) + geom_boxplot(lwd=1,aes(color=factor(ComparisonGroups),fill = NA), outlier.size = 3)  +
  stat_summary(fun=mean, geom="point", shape=5, size=8) + scale_fill_manual(values = col) + scale_colour_manual(values = col ) +
  geom_point(size=4,aes(color = factor(ComparisonGroups))) + xlab(NULL)  +
  ylab("Bray Curtis Distance \n") + paramsAngled()  + xlab(NULL)  #+ scale_color_tableau()
ggsave(NoAbxBCDistance, file = "NoAbxBCDistance.pdf", width = 8, height = 8, limitsize = FALSE)


# Get only the no infant antibiotics:
AbxRows<-grep("-Infant.Abx", MeanStoolDistanceTable$ComparisonGroups)
NoAbxMeanStoolDistanceTable<-MeanStoolDistanceTable[-AbxRows,]


names(MeanDistanceTable)<-c("ComparisonGroups", "MeanDistance")
MeanDistanceTable<-MeanDistanceTable[order(MeanDistanceTable$MeanDistance, decreasing = FALSE),]

### alright now don't separate by gestational age and just look week 1 to week 3 and antibiotics or no antibiotics:

BCDistanceTableNoGA<-as.data.frame(BCMatrix)
names(BCDistanceTableNoGA)<-AxillaGroinStoolGenusGA$Sample
BCDistanceTableNoGA$Sample<-row.names(BCDistanceTableNoGA)
BCDistanceTableNoGA<-melt(BCDistanceTableNoGA, id.vars = "Sample")
names(BCDistanceTableNoGA)<-c("Sample", "SampleID", "BCDistance")
BCDistanceTableNoGA<-merge(BCDistanceTableNoGA, AxillaGroinStoolGenusGA[,c(1,7,11,12)], by = "Sample", all.x = TRUE)
names(BCDistanceTableNoGA)<-c("Sample", "SampleID", "BCDistance", "SampleType", "SampleWeek", "SampleAbx")
BCDistanceTableNoGA$TypeWeekAbx<-paste(BCDistanceTableNoGA$SampleType, BCDistanceTableNoGA$SampleWeek, BCDistanceTableNoGA$SampleAbx, sep = "-")
StoolBCDistanceTableNoGA<-subset(BCDistanceTableNoGA, BCDistanceTableNoGA$SampleType == "Stool")
StoolBCDistanceTableNoGA$SampleType<-NULL
StoolBCDistanceTableNoGA$SampleWeek<-NULL
StoolBCDistanceTableNoGA$SampleAbx<-NULL

StoolBCDistanceTableNoGA<-merge(StoolBCDistanceTableNoGA, AxillaGroinStoolGenusGA[,c(6,7,11,12)], by = "SampleID", all.x = TRUE)
StoolBCDistanceTableNoGA$ComparitorTypeWeekAbx<-paste(StoolBCDistanceTableNoGA$SampleType, StoolBCDistanceTableNoGA$SampleCollectionWeek, StoolBCDistanceTableNoGA$PostNatalAntibiotics, sep = "-")
StoolBCDistanceTableNoGA$SampleType<-NULL
StoolBCDistanceTableNoGA$SampleCollectionWeek<-NULL
StoolBCDistanceTableNoGA$PostNatalAntibiotics<-NULL
StoolComparitors<-grep("Stool", StoolBCDistanceTableNoGA$ComparitorTypeWeekAbx)
StoolBCDistanceTableNoGA<-StoolBCDistanceTableNoGA[StoolComparitors,]



names(StoolBCDistanceTableNoGA)<-c("SampleID", "Sample", "BCDistance", "TypeWeekAbx", "ComparitorTypeWeekAbx")
StoolBCDistanceTableNoGA$ComparisonGroups<-paste(StoolBCDistanceTableNoGA$TypeWeekAbx, StoolBCDistanceTableNoGA$ComparitorTypeWeekAbx, sep = "-")

StoolMeanDistanceTableNoGA<-StoolBCDistanceTableNoGA %>%
  group_by(ComparisonGroups) %>%
  summarize(mean(BCDistance))

names(StoolMeanDistanceTableNoGA)<-c("ComparisonGroups", "MeanDistance")
StoolMeanDistanceTableNoGA<-StoolMeanDistanceTableNoGA[order(StoolMeanDistanceTableNoGA$MeanDistance, decreasing = FALSE),]

StoolMeanDistanceTableNoGA$ComparisonGroups<-gsub("Stool-", "", StoolMeanDistanceTableNoGA$ComparisonGroups)
StoolMeanDistanceTableNoGA$SampleWeek<-sapply(strsplit(StoolMeanDistanceTableNoGA$ComparisonGroups, "-"),`[`, 1)
StoolMeanDistanceTableNoGA$ComparitorWeek<-sapply(strsplit(StoolMeanDistanceTableNoGA$ComparisonGroups, "-"),`[`, 3)
StoolMeanDistanceTableNoGA$SampleAbx<-sapply(strsplit(StoolMeanDistanceTableNoGA$ComparisonGroups, "-"),`[`, 2)
StoolMeanDistanceTableNoGA$ComparitorAbx<-sapply(strsplit(StoolMeanDistanceTableNoGA$ComparisonGroups, "-"),`[`, 4)

write.csv(StoolMeanDistanceTableNoGA, file = "StoolMeanDistanceTableNoGA.csv")

## whoops in this case we don't want to subset to stool. The figure has all sample types

BCDistanceTableNoGA$SampleType<-NULL
BCDistanceTableNoGA$SampleWeek<-NULL
BCDistanceTableNoGA$SampleAbx<-NULL

BCDistanceTableNoGA<-merge(BCDistanceTableNoGA, AxillaGroinStoolGenusGA[,c(6,7,11,12)], by = "SampleID", all.x = TRUE)
BCDistanceTableNoGA$ComparitorTypeWeekAbx<-paste(BCDistanceTableNoGA$SampleType, BCDistanceTableNoGA$SampleCollectionWeek, BCDistanceTableNoGA$PostNatalAntibiotics, sep = "-")
BCDistanceTableNoGA$SampleType<-NULL
BCDistanceTableNoGA$SampleCollectionWeek<-NULL
BCDistanceTableNoGA$PostNatalAntibiotics<-NULL


names(BCDistanceTableNoGA)<-c("SampleID", "Sample", "BCDistance", "TypeWeekAbx", "ComparitorTypeWeekAbx")
BCDistanceTableNoGA$ComparisonGroups<-paste(BCDistanceTableNoGA$TypeWeekAbx, BCDistanceTableNoGA$ComparitorTypeWeekAbx, sep = "-")

MeanDistanceTableNoGA<-BCDistanceTableNoGA %>%
  group_by(ComparisonGroups) %>%
  summarize(mean(BCDistance))

names(MeanDistanceTableNoGA)<-c("ComparisonGroups", "MeanDistance")
MeanDistanceTableNoGA<-MeanDistanceTableNoGA[order(MeanDistanceTableNoGA$MeanDistance, decreasing = FALSE),]

MeanDistanceTableNoGA$SampleType<-sapply(strsplit(MeanDistanceTableNoGA$ComparisonGroups, "-"),`[`, 1)
MeanDistanceTableNoGA$ComparitorType<-sapply(strsplit(MeanDistanceTableNoGA$ComparisonGroups, "-"),`[`, 4)
MeanDistanceTableNoGA$SampleWeek<-sapply(strsplit(MeanDistanceTableNoGA$ComparisonGroups, "-"),`[`, 2)
MeanDistanceTableNoGA$ComparitorWeek<-sapply(strsplit(MeanDistanceTableNoGA$ComparisonGroups, "-"),`[`, 5)
MeanDistanceTableNoGA$SampleAbx<-sapply(strsplit(MeanDistanceTableNoGA$ComparisonGroups, "-"),`[`, 3)
MeanDistanceTableNoGA$ComparitorAbx<-sapply(strsplit(MeanDistanceTableNoGA$ComparisonGroups, "-"),`[`, 6)

## mark whether sample type and antibiotic exposure are the same
MeanDistanceTableNoGA$SameType<-ifelse(MeanDistanceTableNoGA$SampleType == MeanDistanceTableNoGA$ComparitorType, "Same.Type", "Different.Type")
MeanDistanceTableNoGA$SameAbx<-ifelse(MeanDistanceTableNoGA$SampleAbx == MeanDistanceTableNoGA$ComparitorAbx, "Same.Abx", "Different.Abx")


write.csv(MeanDistanceTableNoGA, file = "AllSamplesMeanDistanceTableNoGA.csv")

save.image(file = "NICUData20250209")

# We should get all values to do box plots and statistics. And we need to compare one subject to him/herself
BCDistanceTableNoGA$SampleType<-sapply(strsplit(BCDistanceTableNoGA$ComparisonGroups, "-"),`[`, 1)
BCDistanceTableNoGA$ComparitorType<-sapply(strsplit(BCDistanceTableNoGA$ComparisonGroups, "-"),`[`, 4)
BCDistanceTableNoGA$SampleWeek<-sapply(strsplit(BCDistanceTableNoGA$ComparisonGroups, "-"),`[`, 2)
BCDistanceTableNoGA$ComparitorWeek<-sapply(strsplit(BCDistanceTableNoGA$ComparisonGroups, "-"),`[`, 5)
BCDistanceTableNoGA$SampleAbx<-sapply(strsplit(BCDistanceTableNoGA$ComparisonGroups, "-"),`[`, 3)
BCDistanceTableNoGA$ComparitorAbx<-sapply(strsplit(BCDistanceTableNoGA$ComparisonGroups, "-"),`[`, 6)
BCDistanceTableNoGA$SampleSubject<-sapply(strsplit(as.character(BCDistanceTableNoGA$SampleID), "_"),`[`, 1)
BCDistanceTableNoGA$ComparitorSubject<-sapply(strsplit(as.character(BCDistanceTableNoGA$Sample), "_"),`[`, 1)

## mark whether sample type and antibiotic exposure are the same
BCDistanceTableNoGA$SameType<-ifelse(BCDistanceTableNoGA$SampleType == BCDistanceTableNoGA$ComparitorType, "Same.Type", "Different.Type")
BCDistanceTableNoGA$SameAbx<-ifelse(BCDistanceTableNoGA$SampleAbx == BCDistanceTableNoGA$ComparitorAbx, "Same.Abx", "Different.Abx")
BCDistanceTableNoGA$SameSubject<-ifelse(BCDistanceTableNoGA$SampleSubject == BCDistanceTableNoGA$ComparitorSubject, "Same.Subject", "Different.Subject")

BCDistanceTableNoGA<-subset(BCDistanceTableNoGA, BCDistanceTableNoGA$SameSubject=="Same.Subject")
BCDistanceTableNoGA<-subset(BCDistanceTableNoGA, BCDistanceTableNoGA$SameAbx=="Same.Abx")
BCDistanceTableNoGA<-subset(BCDistanceTableNoGA, BCDistanceTableNoGA$SameType=="Same.Type")

write.csv(BCDistanceTableNoGA, file = "AllSamplesBCDistanceTableNoGA.csv")


# so it turns out that all the 23 to 27 weekers got antibiotics. We'll need to exclude them from the prior analysis
## Next do axilla, groin, and stool. (US only) to see what changes with antibiotics or from week 1 to week 3

Week1AxillaGroinStoolGenus<-subset(NICUGenusNR, NICUGenusNR$SampleType %in% c("Axilla", "Groin", "Stool") &
                                     NICUGenusNR$SampleCollectionWeek == "Week.1")
AbxSampleInd<-Week1AxillaGroinStoolGenus$Sample[Week1AxillaGroinStoolGenus$PostNatalAntibiotics == "Infant.Abx"]
NoAbxSampleInd<-Week1AxillaGroinStoolGenus$Sample[Week1AxillaGroinStoolGenus$PostNatalAntibiotics == "No.Infant.Abx"]


AxillaGroinStoolGenus<-subset(NICUGenusNR, NICUGenusNR$SampleType %in% c("Axilla", "Groin", "Stool"))

Week1Ind<-AxillaGroinStoolGenus$Sample[AxillaGroinStoolGenus$SampleCollectionWeek == "Week.1"]
Week3Ind<-AxillaGroinStoolGenus$Sample[AxillaGroinStoolGenus$SampleCollectionWeek == "Week.3"]


metadata <- AxillaGroinStoolGenus[, 1:17]
cts <- as.matrix(AxillaGroinStoolGenus[, -(1:17)])
rownames(cts) <- metadata$Sample
cts_l2 <- glog2(cts)
grps <- AxillaGroinStoolGenus$TypeAbx

col=tableau.colors
newPCA<-PCA(cts_l2, scale.unit = FALSE, ncp = 5, graph = FALSE)
pcaPlot<-fviz_pca_ind(newPCA,
                      geom.ind = "point", 
                      pointsize = 1.5,
                      point.alph = 0.1,
                      title = "Unsupervised Principal Coordinate Analysis",
                      subtitle = "Samples Colored by Sample Type and Abx Exposure : Weeks 1 and 3",
                      col.ind = grps, 
                      addEllipses = TRUE, 
                      ellipse.alpha = 0.2,
                      ellipse.type = "confidence",
                      ellipse.level = 0.95,
                      legend.title = "Sample", 
                      legend.size = 11,
                      mean.point = TRUE,
                      palette = col,
                      axes.linetype = "blank"
)
SampleTypeAbxPCA<- pcaPlot + xlim(-25,34) + ylim(-25,33)

pdf("SampleTypeAbxBothWeeksGenusPCA.pdf") 
print(SampleTypeAbxWk1PCA)
dev.off() 



# This should be same figure as above but week 1

col<-tableau.colors
pcaPlot<-fviz_pca_ind(newPCA,
                      geom.ind = "point", 
                      pointsize = 1.5,
                      point.alph = 0.1,
                      title = "Unsupervised Principal Coordinate Analysis",
                      subtitle = "Week 1 Samples : All Subjects",
                      col.ind = grps, 
                      addEllipses = TRUE, 
                      ellipse.alpha = 0.2,
                      ellipse.type = "confidence",
                      ellipse.level = 0.95,
                      legend.title = "Week 1 Sample", 
                      legend.size = 11,
                      mean.point = TRUE,
                      palette = col,
                      axes.linetype = "blank",
                      select.ind = list(name = Week1Ind)
)
SampleTypeWk1PCA<- pcaPlot + xlim(-25,34) + ylim(-25,33)

pdf("SampleTypeWk1GenusPCA.pdf") 
print(SampleTypeWk1PCA)
dev.off() 

# col<-tableau.colors[c(1,3,5, 6,7,8)]
pcaPlot<-fviz_pca_ind(newPCA,
                      geom.ind = "point", 
                      pointsize = 1.5,
                      point.alph = 0.1,
                      title = "Unsupervised Principal Coordinate Analysis",
                      subtitle = "Week 3 Samples : All Subjects",
                      col.ind = grps, 
                      addEllipses = TRUE, 
                      ellipse.alpha = 0.2,
                      ellipse.type = "confidence",
                      ellipse.level = 0.95,
                      legend.title = "Week 3 Sample", 
                      legend.size = 11,
                      mean.point = TRUE,
                      palette = col,
                      axes.linetype = "blank",
                      select.ind = list(name = Week3Ind)
)
SampleTypeWk3PCA<- pcaPlot + xlim(-25,34) + ylim(-25,33)

pdf("SampleTypeWk3GenusPCA.pdf") 
print(SampleTypeWk3PCA)
dev.off() 

## Ok do a bunch of MRPP for different times and sites to see if abx have significant effect

Week1StoolGenus<-subset(NICUGenusNR, NICUGenusNR$SampleType %in% c("Stool") &
                          NICUGenusNR$SampleCollectionWeek == "Week.1")
#p = 0.001
Wk1StoolGenus.mrpp<-mrpp(Week1StoolGenus[,18:ncol(Week1StoolGenus)], Week1StoolGenus$PostNatalAntibiotics, distance = "bray")
Wk1StoolGenus.mrpp

# what about GLMM to account for GA, maternal antibiotics etc.



Week1AxillaGenus<-subset(NICUGenusNR, NICUGenusNR$SampleType %in% c("Axilla") &
                           NICUGenusNR$SampleCollectionWeek == "Week.1")
#p = 0.005
Wk1AxillaGenus.mrpp<-mrpp(Week1AxillaGenus[,18:ncol(Week1AxillaGenus)], Week1AxillaGenus$PostNatalAntibiotics, distance = "bray")
Wk1AxillaGenus.mrpp

Week1GroinGenus<-subset(NICUGenusNR, NICUGenusNR$SampleType %in% c("Groin") &
                          NICUGenusNR$SampleCollectionWeek == "Week.1")
#p = 0.152
Wk1GroinGenus.mrpp<-mrpp(Week1GroinGenus[,18:ncol(Week1GroinGenus)], Week1GroinGenus$PostNatalAntibiotics, distance = "bray")
Wk1GroinGenus.mrpp

Week3StoolGenus<-subset(NICUGenusNR, NICUGenusNR$SampleType %in% c("Stool") &
                          NICUGenusNR$SampleCollectionWeek == "Week.3")
#p = 0.048
Wk3StoolGenus.mrpp<-mrpp(Week3StoolGenus[,18:ncol(Week3StoolGenus)], Week3StoolGenus$PostNatalAntibiotics, distance = "bray")
Wk3StoolGenus.mrpp

Week3AxillaGenus<-subset(NICUGenusNR, NICUGenusNR$SampleType %in% c("Axilla") &
                           NICUGenusNR$SampleCollectionWeek == "Week.3")
#p = 0.486
Wk3AxillaGenus.mrpp<-mrpp(Week3AxillaGenus[,18:ncol(Week3AxillaGenus)], Week3AxillaGenus$PostNatalAntibiotics, distance = "bray")
Wk3AxillaGenus.mrpp

Week3GroinGenus<-subset(NICUGenusNR, NICUGenusNR$SampleType %in% c("Groin") &
                          NICUGenusNR$SampleCollectionWeek == "Week.3")
#p = 0.051
Wk3GroinGenus.mrpp<-mrpp(Week3GroinGenus[,18:ncol(Week3GroinGenus)], Week3GroinGenus$PostNatalAntibiotics, distance = "bray")
Wk3GroinGenus.mrpp

Week1AxillaGroinStoolGenus<-subset(NICUGenusNR, NICUGenusNR$SampleType %in% c("Axilla", "Groin", "Stool") &
                                     NICUGenusNR$SampleCollectionWeek == "Week.1")

metadata <- Week1AxillaGroinStoolGenus[, 1:17]
cts <- as.matrix(Week1AxillaGroinStoolGenus[, -(1:17)])
rownames(cts) <- metadata$Sample
cts_l2 <- glog2(cts)
grps <- Week1AxillaGroinStoolGenus$TypeAbx

col=tableau.colors
newPCA<-PCA(cts_l2, scale.unit = FALSE, ncp = 5, graph = FALSE)
pcaPlot<-fviz_pca_ind(newPCA,
                      geom.ind = "point", 
                      pointsize = 1.5,
                      point.alph = 0.1,
                      title = "Unsupervised Principal Coordinate Analysis",
                      subtitle = "Week 1 Samples Grouped by Sample Type and Abx Exposure",
                      col.ind = grps, 
                      addEllipses = TRUE, 
                      ellipse.alpha = 0.2,
                      ellipse.type = "confidence",
                      ellipse.level = 0.95,
                      legend.title = "Week 1 Sample", 
                      legend.size = 11,
                      mean.point = TRUE,
                      palette = col,
                      axes.linetype = "blank"
)
SampleTypeAbxWk1PCA<- pcaPlot + ylim(-25, 36) + xlim(-26,34)

pdf("SampleTypeAbxWk1PCA.pdf") 
print(SampleTypeAbxWk1PCA)
dev.off() 


Week3AxillaGroinStoolGenus<-subset(NICUGenusNR, NICUGenusNR$SampleType %in% c("Axilla", "Groin", "Stool") &
                                     NICUGenusNR$SampleCollectionWeek == "Week.3")

metadata <- Week3AxillaGroinStoolGenus[, 1:17]
cts <- as.matrix(Week3AxillaGroinStoolGenus[, -(1:17)])
rownames(cts) <- metadata$Sample
cts_l2 <- glog2(cts)
grps <- Week3AxillaGroinStoolGenus$TypeAbx

col=tableau.colors
newPCA<-PCA(cts_l2, scale.unit = FALSE, ncp = 5, graph = FALSE)
pcaPlot<-fviz_pca_ind(newPCA,
                      geom.ind = "point", 
                      pointsize = 1.5,
                      point.alph = 0.1,
                      title = "Unsupervised Principal Coordinate Analysis",
                      subtitle = "Week 3 Samples Grouped by Sample Type and Abx Exposure",
                      col.ind = grps, 
                      addEllipses = TRUE, 
                      ellipse.alpha = 0.2,
                      ellipse.type = "confidence",
                      ellipse.level = 0.95,
                      legend.title = "Week 3 Sample", 
                      legend.size = 11,
                      mean.point = TRUE,
                      palette = col,
                      axes.linetype = "blank"
)
SampleTypeAbxWk3PCA<-  pcaPlot + ylim(-25, 36) + xlim(-26,34)

pdf("SampleTypeAbxWk3GenusPCA.pdf") 
print(SampleTypeAbxWk3PCA)
dev.off() 

# OK look at week 1 stools to see what's different:

Week1Stools<-subset(NICUGenusNR, NICUGenusNR$SampleType == "Stool" & NICUGenusNR$SampleCollectionWeek =="Week.1")
Week1Stools<-Week1Stools[,c(1,10,18:ncol(Week1Stools))]


# p = 0.017 if No.Abx versus Abx
Week1Stools.mrpp<-mrpp(Week1Stools[,3:ncol(Week1Stools)], Week1Stools$PostNatalAntibiotics, distance = "bray")
Week1Stools.mrpp


# get all pairwise.wilcox for every species and aggregate into a table:

# ok need to combine high and low antibiotics for this part


Wilcox<-pairwise.wilcox.test(Week1Stools[,3], Week1Stools[,2])
WilcoxTable<-as.data.frame(Wilcox$p.value)
names(WilcoxTable)<-names(Week1Stools)[3]
row.names(WilcoxTable)<-"Wilcox"


for (i in 4:length(colnames(Week1Stools))){
  x<-pairwise.wilcox.test(Week1Stools[,i], Week1Stools[,2])
  y<-as.data.frame(x$p.value)
  names(y)<-names(Week1Stools)[i]
  WilcoxTable<-cbind(WilcoxTable, y)
}

WilcoxTable<-as.data.frame(t(WilcoxTable))
names(WilcoxTable)<-c("Unadjusted_p")
WilcoxTable$FDR<-p.adjust(WilcoxTable$Unadjusted_p, method = "fdr")
WilcoxTable$Genus<-row.names(WilcoxTable)
WilcoxTable<-WilcoxTable[order(WilcoxTable$Unadjusted_p),]



# Now let's figure out relative contribution of each species, fold difference between groups, and direction of difference

GenusCounts<-as.data.frame(colSums(Week1Stools[3:length(colnames(Week1Stools))]))
names(GenusCounts)<-"GenusTotal"
GenusCountsBigTable<-as.data.frame(colSums(Week1Stools[3:length(colnames(Week1Stools))]))
names(GenusCountsBigTable)<-"GenusTotal"
TotalCountsAll<-sum(GenusCountsBigTable$GenusTotal)
GenusCounts$Genus<-row.names(GenusCounts)

for (i in 1:length(rownames(GenusCounts))){
  GenusCounts$Fraction[i]<-GenusCounts[i,1] / TotalCountsAll
}



GenusGroupMeans<-aggregate(.~Week1Stools$PostNatalAntibiotics, mean, data=Week1Stools[3:ncol(Week1Stools)])
row.names(GenusGroupMeans)<-GenusGroupMeans$`Week1Stools$PostNatalAntibiotics`
GenusGroupMeans<-as.data.frame(t(GenusGroupMeans[,2:ncol(GenusGroupMeans)]))
GenusGroupMeans$Genus<-row.names(GenusGroupMeans)
OverallMeans<-as.data.frame(colMeans(Week1Stools[3:ncol(Week1Stools)]))
OverallMeans$Genus<-row.names(OverallMeans)
names(OverallMeans)<-c("OverallMean", "Genus")
GenusGroupMeans<-merge(GenusGroupMeans, OverallMeans, by = "Genus", all.x = TRUE)
GenusGroupMeans$No.Infant.Abx<-as.numeric(GenusGroupMeans$No.Infant.Abx)
GenusGroupMeans$Infant.Abx<-as.numeric(GenusGroupMeans$Infant.Abx)

GenusGroupMeans$FoldChange<-foldchange(GenusGroupMeans$No.Infant.Abx, GenusGroupMeans$Infant.Abx)
GenusGroupMeans$logFC<-foldchange2logratio(GenusGroupMeans$FoldChange, base = 2)
GenusGroupMeans$logFC[is.na(GenusGroupMeans$logFC)]<-0
GenusGroupMeans$logFC[GenusGroupMeans$logFC == -Inf]<-0

# GenusGroupMeans$FoldChange<-(GenusGroupMeans$Infant.Abx -GenusGroupMeans$No.Infant.Abx) / GenusGroupMeans$OverallMean
# GenusGroupMeans$Log2<-log2(GenusGroupMeans$Infant.Abx/GenusGroupMeans$OverallMean)
GenusGroupTable<-merge(GenusGroupMeans, GenusCounts, by = "Genus", all.x = TRUE)
GenusGroupTable<-merge(GenusGroupTable, WilcoxTable, by = "Genus", all.x = TRUE)
GenusGroupTable<-GenusGroupTable[order(GenusGroupTable$Unadjusted_p),]
GenusGroupTable<-GenusGroupTable[,c(1,8,4, 2,3,5,9,10)]
names(GenusGroupTable)<-c("Genus", "Abundance", "OverallMean", "No.Infant.Abx.Mean", "Infant.Abx.Mean", "Fold.Change", "p_Unadjusted", "FDR")
GenusGroupTable$Abundance<-GenusGroupTable$Abundance * 100
GenusGroupTable$Fold.Change[is.na(GenusGroupTable$Fold.Change)]<-0
GenusGroupTable<-GenusGroupTable[order(GenusGroupTable$Fold.Change, decreasing = TRUE),]

SigGenusGroupTable<-subset(GenusGroupTable, GenusGroupTable$p_Unadjusted < 0.05)

metadata <- Week1Stools[, 1:2]
cts <- as.matrix(Week1Stools[, -(1:2)])
rownames(cts) <- metadata$Sample

cts_l2 <- glog2(cts)
grps <- Week1Stools$PostNatalAntibiotics


res <- compute_ef_Abx(d = cts_l2, g = Week1Stools$PostNatalAntibiotics, min_shrink = 0.3)
names(res)[1]<-"Genus"

res <- arrange(left_join(res, SigGenusGroupTable), desc(abs(ef_shrunk))) %>%
  mutate(Genus = as_factor(Genus))

res<-subset(res, ! is.na(res$FDR)) 
resWk1Abx<-res

resWk1Abx$HigherIn<-ifelse(res$ef_shrunk > 0, "Infant.Abx", "No.Infant.Abx")
resWk1Abx$HigherIn<-factor(resWk1Abx$HigherIn, levels = c("Infant.Abx", "No.Infant.Abx"))


write.table(resWk1Abx, file = "SignificantWk1AbxGenusTableNew.csv", sep = ",", row.names = FALSE, col.names = TRUE,  quote = FALSE)
write.table(GenusGroupTable, file = "Wk1GenusTable.csv", sep = ",", row.names = FALSE, col.names = TRUE,  quote = FALSE)


Week1Abx_Effect<- ggplot(filter(resWk1Abx, abs(resWk1Abx$ef_shrunk) > 0.3),
                         aes(x = reorder(Genus, ef_shrunk), y = ef_shrunk, fill = HigherIn)) +
  geom_bar(colour="black", stat="identity") +
  #guides(fill=FALSE) + scale_fill_hue(l=40) +
  guides(fill=guide_legend(override.aes=list(colour=NA))) + 
  coord_flip() +
  labs( x = NULL) +
  ylab(expression(atop("Effect Size"))) +
  guides(fill=guide_legend(title="More Abundant In:", )) +
  #geom_hline(yintercept = c(-0.2, 0.2), color = "firebrick", lty = 2) +
  theme(axis.text.y = element_text(size= 12, color="black")) +
  theme(legend.title = element_text(size=12)) +
  theme(legend.text = element_text(size = 11)) 

Week1Abx_Effect<-Week1Abx_Effect + scale_fill_manual(values = tableau.colors[c(1,2)])

ggsave(filename = "Week1AbxGenus_Effect.pdf", plot = Week1Abx_Effect, width = 10, 
       height = 12, limitsize = FALSE)

# OK look at week 3 stools to see what's different:

Week3Stools<-subset(NICUGenusNR, NICUGenusNR$SampleType == "Stool" & NICUGenusNR$SampleCollectionWeek =="Week.3")
Week3Stools<-Week3Stools[,c(1,9,18:ncol(Week3Stools))]


# p = 0.048 if No.Abx versus Abx
Week3Stools.mrpp<-mrpp(Week3Stools[,3:ncol(Week3Stools)], Week3Stools$PostNatalAntibiotics, distance = "bray")
Week3Stools.mrpp


# get all pairwise.wilcox for every species and aggregate into a table:

# ok need to combine high and low antibiotics for this part


Wilcox<-pairwise.wilcox.test(Week3Stools[,3], Week3Stools[,2])
WilcoxTable<-as.data.frame(Wilcox$p.value)
names(WilcoxTable)<-names(Week3Stools)[3]
row.names(WilcoxTable)<-"Wilcox"


for (i in 4:length(colnames(Week3Stools))){
  x<-pairwise.wilcox.test(Week3Stools[,i], Week3Stools[,2])
  y<-as.data.frame(x$p.value)
  names(y)<-names(Week3Stools)[i]
  WilcoxTable<-cbind(WilcoxTable, y)
}

WilcoxTable<-as.data.frame(t(WilcoxTable))
names(WilcoxTable)<-c("Unadjusted_p")
WilcoxTable$FDR<-p.adjust(WilcoxTable$Unadjusted_p, method = "fdr")
WilcoxTable$Genus<-row.names(WilcoxTable)
WilcoxTable<-WilcoxTable[order(WilcoxTable$Unadjusted_p),]



# Now let's figure out relative contribution of each species, fold difference between groups, and direction of difference

GenusCounts<-as.data.frame(colSums(Week3Stools[3:length(colnames(Week3Stools))]))
names(GenusCounts)<-"GenusTotal"
GenusCountsBigTable<-as.data.frame(colSums(Week3Stools[3:length(colnames(Week3Stools))]))
names(GenusCountsBigTable)<-"GenusTotal"
TotalCountsAll<-sum(GenusCountsBigTable$GenusTotal)
GenusCounts$Genus<-row.names(GenusCounts)

for (i in 1:length(rownames(GenusCounts))){
  GenusCounts$Fraction[i]<-GenusCounts[i,1] / TotalCountsAll
}



GenusGroupMeans<-aggregate(.~Week3Stools$PostNatalAntibiotics, mean, data=Week3Stools[3:ncol(Week3Stools)])
row.names(GenusGroupMeans)<-GenusGroupMeans$`Week3Stools$PostNatalAntibiotics`
GenusGroupMeans<-as.data.frame(t(GenusGroupMeans[,2:ncol(GenusGroupMeans)]))
GenusGroupMeans$Genus<-row.names(GenusGroupMeans)
OverallMeans<-as.data.frame(colMeans(Week3Stools[3:ncol(Week3Stools)]))
OverallMeans$Genus<-row.names(OverallMeans)
names(OverallMeans)<-c("OverallMean", "Genus")
GenusGroupMeans<-merge(GenusGroupMeans, OverallMeans, by = "Genus", all.x = TRUE)
GenusGroupMeans$No.Infant.Abx<-as.numeric(GenusGroupMeans$No.Infant.Abx)
GenusGroupMeans$Infant.Abx<-as.numeric(GenusGroupMeans$Infant.Abx)
GenusGroupMeans$FoldChange<-(GenusGroupMeans$Infant.Abx -GenusGroupMeans$No.Infant.Abx) / GenusGroupMeans$OverallMean
GenusGroupMeans$Log2<-log2(GenusGroupMeans$Infant.Abx/GenusGroupMeans$OverallMean)
GenusGroupTable<-merge(GenusGroupMeans, GenusCounts, by = "Genus", all.x = TRUE)
GenusGroupTable<-merge(GenusGroupTable, WilcoxTable, by = "Genus", all.x = TRUE)
GenusGroupTable<-GenusGroupTable[order(GenusGroupTable$Unadjusted_p),]
GenusGroupTable<-GenusGroupTable[,c(1,8,4, 2,3,5,9,10)]
names(GenusGroupTable)<-c("Genus", "Abundance", "OverallMean", "No.Infant.Abx.Mean", "Infant.Abx.Mean", "Fold.Change", "p_Unadjusted", "FDR")
GenusGroupTable$Abundance<-GenusGroupTable$Abundance * 100
GenusGroupTable$Fold.Change[is.na(GenusGroupTable$Fold.Change)]<-0
GenusGroupTable<-GenusGroupTable[order(GenusGroupTable$Fold.Change, decreasing = TRUE),]

SigGenusGroupTable<-subset(GenusGroupTable, GenusGroupTable$p_Unadjusted < 0.05)

# Ok try out new effect size calculation from Klaus
# there's probably a way to shorten this part but we'll use his approach and first split out class and counts:

metadata <- Week3Stools[, 1:2]
cts <- as.matrix(Week3Stools[, -(1:2)])
rownames(cts) <- metadata$Sample

cts_l2 <- glog2(cts)
grps <- Week3Stools$PostNatalAntibiotics


res <- compute_ef_Abx(d = cts_l2, g = Week3Stools$PostNatalAntibiotics, min_shrink = 0.3)
names(res)[1]<-"Genus"

res <- arrange(left_join(res, SigGenusGroupTable), desc(abs(ef_shrunk))) %>%
  mutate(Genus = as_factor(Genus))

res<-subset(res, ! is.na(res$FDR)) 
resWk3Abx<-res

resWk3Abx$HigherIn<-ifelse(res$ef_shrunk > 0, "Infant.Abx", "No.Infant.Abx")
resWk3Abx$HigherIn<-factor(resWk3Abx$HigherIn, levels = c("Infant.Abx", "No.Infant.Abx"))


write.table(resWk3Abx, file = "SignificantWk3AbxGenusTableNew.csv", sep = ",", row.names = FALSE, col.names = TRUE,  quote = FALSE)
write.table(GenusGroupTable, file = "Wk3GenusTable.csv", sep = ",", row.names = FALSE, col.names = TRUE,  quote = FALSE)


Week3Abx_Effect<- ggplot(filter(resWk3Abx, abs(resWk3Abx$ef_shrunk) > 0.3),
                         aes(x = reorder(Genus, ef_shrunk), y = ef_shrunk, fill = HigherIn)) +
  geom_bar(colour="black", stat="identity") +
  #guides(fill=FALSE) + scale_fill_hue(l=40) +
  guides(fill=guide_legend(override.aes=list(colour=NA))) + 
  coord_flip() +
  labs( x = NULL) +
  ylab(expression(atop("Effect Size"))) +
  guides(fill=guide_legend(title="More Abundant In:", )) +
  #geom_hline(yintercept = c(-0.2, 0.2), color = "firebrick", lty = 2) +
  theme(axis.text.y = element_text(size= 12, color="black")) +
  theme(legend.title = element_text(size=12)) +
  theme(legend.text = element_text(size = 11)) 

Week3Abx_Effect<-Week3Abx_Effect + scale_fill_manual(values = tableau.colors[c(2)])

ggsave(filename = "Week3AbxGenus_Effect.pdf", plot = Week3Abx_Effect, width = 10, 
       height = 12, limitsize = FALSE)




# Make a dummy table so can set 0 counts to 1

DummyGenusNR<-NICUGenusNR


DummyGenusNR$Bifidobacterium[DummyGenusNR$Bifidobacterium==0]<-1
BifidoAbx<-ggplot(subset(DummyGenusNR, DummyGenusNR$SampleType %in% c("Axilla", "Groin", "Stool")),
                  aes(x=PostNatalAntibiotics, y=as.numeric(Bifidobacterium), fill=PostNatalAntibiotics)) + geom_boxplot(lwd=1,aes(color=factor(PostNatalAntibiotics),fill = NA), outlier.size = 3)  +
  stat_summary(fun=mean, geom="point", shape=5, size=8) + scale_fill_manual(values = col) + scale_colour_manual(values = col ) +
  geom_point(size=4,aes(color = factor(PostNatalAntibiotics))) + xlab(NULL) + facet_grid(rows = SampleType ~ SampleCollectionWeek) +
  ylab("Bifidobacterium \n") + paramsAngled() + scale_y_log10() + xlab(NULL)  + scale_color_tableau() +
  theme(axis.title.y = element_text(face = "italic"))
ggsave(BifidoAbx, file = "BifidoGenusAbx.pdf", width = 8, height = 8, limitsize = FALSE)

DummyGenusNR$Bifidobacterium[DummyGenusNR$Bifidobacterium==0]<-1
BifidoGestation<-ggplot(subset(DummyGenusNR, DummyGenusNR$SampleType %in% c("Axilla", "Groin", "Stool") & DummyGenusNR$PostNatalAntibiotics == "No.Infant.Abx"),
                        aes(x=GestationCohort, y=as.numeric(Bifidobacterium), fill=GestationCohort)) + geom_boxplot(lwd=1,aes(color=factor(GestationCohort),fill = NA), outlier.size = 3)  +
  stat_summary(fun=mean, geom="point", shape=5, size=8) + scale_fill_manual(values = col) + scale_colour_manual(values = col ) +
  geom_point(size=4,aes(color = factor(GestationCohort))) + xlab(NULL) + facet_grid(rows = SampleType ~ SampleCollectionWeek) +
  ylab("Bifidobacterium \n") + paramsAngled() + scale_y_log10() + xlab(NULL)  + scale_color_tableau() +
  theme(axis.title.y = element_text(face = "italic"))
ggsave(BifidoGestation, file = "BifidoGenusGestationNoAbx.pdf", width = 8, height = 8, limitsize = FALSE)


DummyGenusNR$Enterobacter[DummyGenusNR$Enterobacter==0]<-1
EnterobacterAbx<-ggplot(subset(DummyGenusNR, DummyGenusNR$SampleType %in% c("Axilla", "Groin", "Stool")),
                        aes(x=PostNatalAntibiotics, y=as.numeric(Enterobacter), fill=PostNatalAntibiotics)) + geom_boxplot(lwd=1,aes(color=factor(PostNatalAntibiotics),fill = NA), outlier.size = 3)  +
  stat_summary(fun=mean, geom="point", shape=5, size=8) + scale_fill_manual(values = col) + scale_colour_manual(values = col ) +
  geom_point(size=4,aes(color = factor(PostNatalAntibiotics))) + xlab(NULL) + facet_grid(rows = SampleType ~ SampleCollectionWeek) +
  ylab("Enterobacter \n") + paramsAngled() + scale_y_log10() + xlab(NULL)  + scale_color_tableau() +
  theme(axis.title.y = element_text(face = "italic"))
ggsave(EnterobacterAbx, file = "EnterobacterAbx.pdf", width = 8, height = 8, limitsize = FALSE)

DummyGenusNR$Escherichia[DummyGenusNR$Escherichia==0]<-1
EscherichiaAbx<-ggplot(subset(DummyGenusNR, DummyGenusNR$SampleType %in% c("Axilla", "Groin", "Stool")),
                       aes(x=PostNatalAntibiotics, y=as.numeric(Escherichia), fill=PostNatalAntibiotics)) + geom_boxplot(lwd=1,aes(color=factor(PostNatalAntibiotics),fill = NA), outlier.size = 3)  +
  stat_summary(fun=mean, geom="point", shape=5, size=8) + scale_fill_manual(values = col) + scale_colour_manual(values = col ) +
  geom_point(size=4,aes(color = factor(PostNatalAntibiotics))) + xlab(NULL) + facet_grid(rows = SampleType ~ SampleCollectionWeek) +
  ylab("Escherichia \n") + paramsAngled() + scale_y_log10() + xlab(NULL)  + scale_color_tableau() +
  theme(axis.title.y = element_text(face = "italic"))
ggsave(EscherichiaAbx, file = "EscherichiaAbx.pdf", width = 8, height = 8, limitsize = FALSE)


DummyGenusNR$Streptococcus[DummyGenusNR$Streptococcus==0]<-1
StreptococcusAbx<-ggplot(subset(DummyGenusNR, DummyGenusNR$SampleType %in% c("Axilla", "Groin", "Stool")),
                         aes(x=PostNatalAntibiotics, y=as.numeric(Streptococcus), fill=PostNatalAntibiotics)) + geom_boxplot(lwd=1,aes(color=factor(PostNatalAntibiotics),fill = NA), outlier.size = 3)  +
  stat_summary(fun=mean, geom="point", shape=5, size=8) + scale_fill_manual(values = col) + scale_colour_manual(values = col ) +
  geom_point(size=4,aes(color = factor(PostNatalAntibiotics))) + xlab(NULL) + facet_grid(rows = SampleType ~ SampleCollectionWeek) +
  ylab("Streptococcus \n") + paramsAngled() + scale_y_log10() + xlab(NULL)  + scale_color_tableau() +
  theme(axis.title.y = element_text(face = "italic"))
ggsave(StreptococcusAbx, file = "StreptococcusAbx.pdf", width = 8, height = 8, limitsize = FALSE)


StreptococcusStoolAbx<-ggplot(subset(DummyGenusNR, DummyGenusNR$SampleType %in% c("Stool")),
                              aes(x=PostNatalAntibiotics, y=as.numeric(Streptococcus), fill=PostNatalAntibiotics)) + geom_boxplot(lwd=1,aes(color=factor(PostNatalAntibiotics),fill = NA), outlier.size = 3)  +
  stat_summary(fun=mean, geom="point", shape=5, size=8) + scale_fill_manual(values = col) + scale_colour_manual(values = col ) +
  geom_point(size=4,aes(color = factor(PostNatalAntibiotics))) + xlab(NULL) + facet_grid(rows = SampleCollectionWeek ~ .) +
  ylab("Streptococcus \n") + paramsAngled() + scale_y_log10() + xlab(NULL)  + scale_color_tableau() +
  theme(axis.title.y = element_text(face = "italic"))
ggsave(StreptococcusStoolAbx, file = "StreptococcusStoolAbx.pdf", width = 5, height = 8, limitsize = FALSE)

DummyGenusNR$Clostridium[DummyGenusNR$Clostridium==0]<-1
ClostridiumStoolAbx<-ggplot(subset(DummyGenusNR, DummyGenusNR$SampleType %in% c("Stool")),
                            aes(x=PostNatalAntibiotics, y=as.numeric(Clostridium), fill=PostNatalAntibiotics)) + geom_boxplot(lwd=1,aes(color=factor(PostNatalAntibiotics),fill = NA), outlier.size = 3)  +
  stat_summary(fun=mean, geom="point", shape=5, size=8) + scale_fill_manual(values = col) + scale_colour_manual(values = col ) +
  geom_point(size=4,aes(color = factor(PostNatalAntibiotics))) + xlab(NULL) + facet_grid(rows = SampleCollectionWeek ~ .) +
  ylab("Clostridium \n") + paramsAngled() + scale_y_log10() + xlab(NULL)  + scale_color_tableau() +
  theme(axis.title.y = element_text(face = "italic"))
ggsave(ClostridiumStoolAbx, file = "ClostridiumStoolAbx.pdf", width = 5, height = 8, limitsize = FALSE)

DummyGenusNR$Candida[DummyGenusNR$Candida==0]<-1
CandidaStoolAbx<-ggplot(subset(DummyGenusNR, DummyGenusNR$SampleType %in% c("Stool")),
                        aes(x=PostNatalAntibiotics, y=as.numeric(Candida), fill=PostNatalAntibiotics)) + geom_boxplot(lwd=1,aes(color=factor(PostNatalAntibiotics),fill = NA), outlier.size = 3)  +
  stat_summary(fun=mean, geom="point", shape=5, size=8) + scale_fill_manual(values = col) + scale_colour_manual(values = col ) +
  geom_point(size=4,aes(color = factor(PostNatalAntibiotics))) + xlab(NULL) + facet_grid(rows = SampleCollectionWeek ~ .) +
  ylab("Candida \n") + paramsAngled() + scale_y_log10() + xlab(NULL)  + scale_color_tableau() +
  theme(axis.title.y = element_text(face = "italic"))
ggsave(CandidaStoolAbx, file = "CandidaStoolAbx.pdf", width = 5, height = 8, limitsize = FALSE)

DummyGenusNR$Staphylococcus[DummyGenusNR$Staphylococcus==0]<-1
StaphylococcusStoolAbx<-ggplot(subset(DummyGenusNR, DummyGenusNR$SampleType %in% c("Stool")),
                               aes(x=PostNatalAntibiotics, y=as.numeric(Staphylococcus), fill=PostNatalAntibiotics)) + geom_boxplot(lwd=1,aes(color=factor(PostNatalAntibiotics),fill = NA), outlier.size = 3)  +
  stat_summary(fun=mean, geom="point", shape=5, size=8) + scale_fill_manual(values = col) + scale_colour_manual(values = col ) +
  geom_point(size=4,aes(color = factor(PostNatalAntibiotics))) + xlab(NULL) + facet_grid(rows = SampleCollectionWeek ~ .) +
  ylab("Staphylococcus \n") + paramsAngled() + scale_y_log10() + xlab(NULL)  + scale_color_tableau() +
  theme(axis.title.y = element_text(face = "italic"))
ggsave(StaphylococcusStoolAbx, file = "StaphylococcusStoolAbx.pdf", width = 5, height = 8, limitsize = FALSE)

DummyGenusNR$Staphylococcus[DummyGenusNR$Staphylococcus==0]<-1
StaphylococcusAbx<-ggplot(subset(DummyGenusNR, DummyGenusNR$SampleType %in% c("Axilla", "Groin", "Stool")),
                          aes(x=PostNatalAntibiotics, y=as.numeric(Staphylococcus), fill=PostNatalAntibiotics)) + geom_boxplot(lwd=1,aes(color=factor(PostNatalAntibiotics),fill = NA), outlier.size = 3)  +
  stat_summary(fun=mean, geom="point", shape=5, size=8) + scale_fill_manual(values = col) + scale_colour_manual(values = col ) +
  geom_point(size=4,aes(color = factor(PostNatalAntibiotics))) + xlab(NULL) + facet_grid(rows = SampleType ~ SampleCollectionWeek) +
  ylab("Staphylococcus \n") + paramsAngled() + scale_y_log10() + xlab(NULL)  + scale_color_tableau() +
  theme(axis.title.y = element_text(face = "italic"))
ggsave(StaphylococcusAbx, file = "StaphylococcusAbx.pdf", width = 8, height = 8, limitsize = FALSE)
# same for Staph, 

DummyGenusNR$Candida[DummyGenusNR$Candida==0]<-1
CandidaAbx<-ggplot(subset(DummyGenusNR, DummyGenusNR$SampleType %in% c("Axilla", "Groin", "Stool")),
                   aes(x=PostNatalAntibiotics, y=as.numeric(Candida), fill=PostNatalAntibiotics)) + geom_boxplot(lwd=1,aes(color=factor(PostNatalAntibiotics),fill = NA), outlier.size = 3)  +
  stat_summary(fun=mean, geom="point", shape=5, size=8) + scale_fill_manual(values = col) + scale_colour_manual(values = col ) +
  geom_point(size=4,aes(color = factor(PostNatalAntibiotics))) + xlab(NULL) + facet_grid(rows = SampleType ~ SampleCollectionWeek) +
  ylab("Candida \n") + paramsAngled() + scale_y_log10() + xlab(NULL)  + scale_color_tableau() +
  theme(axis.title.y = element_text(face = "italic"))
ggsave(CandidaAbx, file = "CandidaAbx.pdf", width = 8, height = 8, limitsize = FALSE)


DummyGenusNR$Pseudomonas[DummyGenusNR$Pseudomonas==0]<-1
PseudomonasAbx<-ggplot(subset(DummyGenusNR, DummyGenusNR$SampleType %in% c("Axilla", "Groin", "Stool")),
                       aes(x=PostNatalAntibiotics, y=as.numeric(Pseudomonas), fill=PostNatalAntibiotics)) + geom_boxplot(lwd=1,aes(color=factor(PostNatalAntibiotics),fill = NA), outlier.size = 3)  +
  stat_summary(fun=mean, geom="point", shape=5, size=8) + scale_fill_manual(values = col) + scale_colour_manual(values = col ) +
  geom_point(size=4,aes(color = factor(PostNatalAntibiotics))) + xlab(NULL) + facet_grid(rows = SampleType ~ SampleCollectionWeek) +
  ylab("Pseudomonas \n") + paramsAngled() + scale_y_log10() + xlab(NULL)  + scale_color_tableau() +
  theme(axis.title.y = element_text(face = "italic"))
ggsave(PseudomonasAbx, file = "PseudomonasAbx.pdf", width = 8, height = 8, limitsize = FALSE)

DummyGenusNR$Klebsiella[DummyGenusNR$Klebsiella==0]<-1
KlebsiellaAbx<-ggplot(subset(DummyGenusNR, DummyGenusNR$SampleType %in% c("Axilla", "Groin", "Stool")),
                      aes(x=PostNatalAntibiotics, y=as.numeric(Klebsiella), fill=PostNatalAntibiotics)) + geom_boxplot(lwd=1,aes(color=factor(PostNatalAntibiotics),fill = NA), outlier.size = 3)  +
  stat_summary(fun=mean, geom="point", shape=5, size=8) + scale_fill_manual(values = col) + scale_colour_manual(values = col ) +
  geom_point(size=4,aes(color = factor(PostNatalAntibiotics))) + xlab(NULL) + facet_grid(rows = SampleType ~ SampleCollectionWeek) +
  ylab("Klebsiella \n") + paramsAngled() + scale_y_log10() + xlab(NULL)  + scale_color_tableau() +
  theme(axis.title.y = element_text(face = "italic"))
ggsave(KlebsiellaAbx, file = "KlebsiellaAbx.pdf", width = 8, height = 8, limitsize = FALSE)


DummyGenusNR$Enterococcus[DummyGenusNR$Enterococcus==0]<-1
EnterococcusAbx<-ggplot(subset(DummyGenusNR, DummyGenusNR$SampleType %in% c("Axilla", "Groin", "Stool")),
                        aes(x=PostNatalAntibiotics, y=as.numeric(Enterococcus), fill=PostNatalAntibiotics)) + geom_boxplot(lwd=1,aes(color=factor(PostNatalAntibiotics),fill = NA), outlier.size = 3)  +
  stat_summary(fun=mean, geom="point", shape=5, size=8) + scale_fill_manual(values = col) + scale_colour_manual(values = col ) +
  geom_point(size=4,aes(color = factor(PostNatalAntibiotics))) + xlab(NULL) + facet_grid(rows = SampleType ~ SampleCollectionWeek) +
  ylab("Enterococcus \n") + paramsAngled() + scale_y_log10() + xlab(NULL)  + scale_color_tableau() +
  theme(axis.title.y = element_text(face = "italic"))
ggsave(EnterococcusAbx, file = "EnterococcusAbx.pdf", width = 8, height = 8, limitsize = FALSE)
#Enterococcus decreased at Week 1 in stool with or without 23-27 wks infants

DummyGenusNR$Clostridioides[DummyGenusNR$Clostridioides==0]<-1
ClostridioidesAbx<-ggplot(subset(DummyGenusNR, DummyGenusNR$SampleType %in% c("Axilla", "Groin", "Stool")),
                          aes(x=PostNatalAntibiotics, y=as.numeric(Clostridioides), fill=PostNatalAntibiotics)) + geom_boxplot(lwd=1,aes(color=factor(PostNatalAntibiotics),fill = NA), outlier.size = 3)  +
  stat_summary(fun=mean, geom="point", shape=5, size=8) + scale_fill_manual(values = col) + scale_colour_manual(values = col ) +
  geom_point(size=4,aes(color = factor(PostNatalAntibiotics))) + xlab(NULL) + facet_grid(rows = SampleType ~ SampleCollectionWeek) +
  ylab("Clostridioides \n") + paramsAngled() + scale_y_log10() + xlab(NULL)  + scale_color_tableau() +
  theme(axis.title.y = element_text(face = "italic"))
ggsave(ClostridioidesAbx, file = "ClostridioidesAbx.pdf", width = 8, height = 8, limitsize = FALSE)
# Large impact of antibiotics on Clostridiodes in stool with our without 23-27 weeks 


pairwise.wilcox.test(DummyGenusNR$Streptococcus[DummyGenusNR$SampleType== "Stool" & DummyGenusNR$SampleCollectionWeek=="Week.1"], DummyGenusNR$PostNatalAntibiotics[DummyGenusNR$SampleType== "Stool" & DummyGenusNR$SampleCollectionWeek=="Week.1"])
pairwise.wilcox.test(DummyGenusNR$Streptococcus[DummyGenusNR$SampleType== "Stool" & DummyGenusNR$SampleCollectionWeek=="Week.3"], DummyGenusNR$PostNatalAntibiotics[DummyGenusNR$SampleType== "Stool" & DummyGenusNR$SampleCollectionWeek=="Week.3"])

pairwise.wilcox.test(DummyGenusNR$Staphylococcus[DummyGenusNR$SampleType== "Stool" & DummyGenusNR$SampleCollectionWeek=="Week.1"], DummyGenusNR$PostNatalAntibiotics[DummyGenusNR$SampleType== "Stool" & DummyGenusNR$SampleCollectionWeek=="Week.1"])
pairwise.wilcox.test(DummyGenusNR$Staphylococcus[DummyGenusNR$SampleType== "Stool" & DummyGenusNR$SampleCollectionWeek=="Week.3"], DummyGenusNR$PostNatalAntibiotics[DummyGenusNR$SampleType== "Stool" & DummyGenusNR$SampleCollectionWeek=="Week.3"])
pairwise.wilcox.test(DummyGenusNR$Candida[DummyGenusNR$SampleType== "Stool" & DummyGenusNR$SampleCollectionWeek=="Week.1"], DummyGenusNR$PostNatalAntibiotics[DummyGenusNR$SampleType== "Stool" & DummyGenusNR$SampleCollectionWeek=="Week.1"])
pairwise.wilcox.test(DummyGenusNR$Candida[DummyGenusNR$SampleType== "Stool" & DummyGenusNR$SampleCollectionWeek=="Week.3"], DummyGenusNR$PostNatalAntibiotics[DummyGenusNR$SampleType== "Stool" & DummyGenusNR$SampleCollectionWeek=="Week.3"])
pairwise.wilcox.test(DummyGenusNR$Clostridium[DummyGenusNR$SampleType== "Stool" & DummyGenusNR$SampleCollectionWeek=="Week.1"], DummyGenusNR$PostNatalAntibiotics[DummyGenusNR$SampleType== "Stool" & DummyGenusNR$SampleCollectionWeek=="Week.1"])
pairwise.wilcox.test(DummyGenusNR$Clostridium[DummyGenusNR$SampleType== "Stool" & DummyGenusNR$SampleCollectionWeek=="Week.3"], DummyGenusNR$PostNatalAntibiotics[DummyGenusNR$SampleType== "Stool" & DummyGenusNR$SampleCollectionWeek=="Week.3"])



Week1StoolGenus<-subset(NICUGenusNR, NICUGenusNR$SampleCollectionWeek == "Week.1"& NICUGenusNR$SampleType == "Stool")
Week1AxillaGenus<-subset(NICUGenusNR, NICUGenusNR$SampleCollectionWeek == "Week.1" & NICUGenusNR$SampleType == "Axilla")
Week1GroinGenus<-subset(NICUGenusNR, NICUGenusNR$SampleCollectionWeek == "Week.1" & NICUGenusNR$SampleType == "Groin")
Week3StoolGenus<-subset(NICUGenusNR, NICUGenusNR$SampleCollectionWeek == "Week.3" & NICUGenusNR$SampleType == "Stool")
Week3AxillaGenus<-subset(NICUGenusNR, NICUGenusNR$SampleCollectionWeek == "Week.3" & NICUGenusNR$SampleType == "Axilla")
Week3GroinGenus<-subset(NICUGenusNR, NICUGenusNR$SampleCollectionWeek == "Week.3" & NICUGenusNR$SampleType == "Groin")

pairwise.wilcox.test(Week1StoolGenus$Staphylococcus, Week1StoolGenus$Location)
pairwise.wilcox.test(Week3StoolGenus$Staphylococcus, Week3StoolGenus$Location)

pairwise.wilcox.test(Week1AxillaGenus$Staphylococcus, Week1AxillaGenus$Location)
pairwise.wilcox.test(Week3AxillaGenus$Staphylococcus, Week3AxillaGenus$Location)

pairwise.wilcox.test(Week1GroinGenus$Staphylococcus, Week1GroinGenus$Location)
pairwise.wilcox.test(Week3GroinGenus$Staphylococcus, Week3GroinGenus$Location)

# check them all:



for (i in 18:ncol(Week1StoolGenus)){
  
  tryCatch({
    
    a = pairwise.wilcox.test(Week1StoolGenus[,i], Week1StoolGenus$Location)[3]
    b = pairwise.wilcox.test(Week3StoolGenus[,i], Week3StoolGenus$Location)[3]
    
    c = pairwise.wilcox.test(Week1AxillaGenus[,i], Week1AxillaGenus$Location)[3]
    d = pairwise.wilcox.test(Week3AxillaGenus[,i], Week3AxillaGenus$Location)[3]
    
    e = pairwise.wilcox.test(Week1GroinGenus[,i], Week1GroinGenus$Location)[3]
    f = pairwise.wilcox.test(Week3GroinGenus[,i], Week3GroinGenus$Location)[3]
    
    head = (colnames(Week1StoolGenus[i]))
    
    
    GenusPvalueTableRow<-cbind(a, b, c, d, e, f)
    row.names(GenusPvalueTableRow) = head
    
    GenusPvalueTable<-rbind(GenusPvalueTable, GenusPvalueTableRow)        
    
    
  }, error=function(e){cat("ERROR :",conditionMessage(e), "\n")})}

GenusPvalueTable<-as.data.frame(GenusPvalueTable)

names(GenusPvalueTable)<-c("Week1.Stool", "Week3.Stool", "Week1.Axilla", "Week3.Axilla", "Week1.Groin", "Week3.Groin")

GenusPvalueTableMatrix<-as.matrix(GenusPvalueTable)

write.csv(GenusPvalueTableMatrix, file = "GenusPValueTable.csv")

# combine body sites

Week1Genus<-subset(NICUGenusNR, NICUGenusNR$SampleCollectionWeek == "Week.1")
Week3Genus<-subset(NICUGenusNR, NICUGenusNR$SampleCollectionWeek == "Week.3")


# check them all:



for (i in 18:ncol(Week1Genus)){
  
  tryCatch({
    
    a = pairwise.wilcox.test(Week1Genus[,i], Week1Genus$Location)[3]
    b = pairwise.wilcox.test(Week3Genus[,i], Week3Genus$Location)[3]
    
    head = (colnames(Week1Genus[i]))
    
    
    GenusPvalueTableRow<-cbind(a, b)
    row.names(GenusPvalueTableRow) = head
    
    GenusPvalueTableNew<-rbind(GenusPvalueTableNew, GenusPvalueTableRow)        
    
    
  }, error=function(e){cat("ERROR :",conditionMessage(e), "\n")})}

GenusPvalueTableNew<-as.data.frame(GenusPvalueTableNew)

names(GenusPvalueTableNew)<-c("Week1", "Week3")

GenusPvalueTableMatrixNew<-as.matrix(GenusPvalueTableNew)

write.csv(GenusPvalueTableMatrixNew, file = "GenusPValueTableWk1Wk3.csv")



# Do species too: Big table



Week1StoolSpecies<-subset(NICUSpeciesNR, NICUSpeciesNR$SampleCollectionWeek == "Week.1"& NICUSpeciesNR$SampleType == "Stool")
Week1AxillaSpecies<-subset(NICUSpeciesNR, NICUSpeciesNR$SampleCollectionWeek == "Week.1" & NICUSpeciesNR$SampleType == "Axilla")
Week1GroinSpecies<-subset(NICUSpeciesNR, NICUSpeciesNR$SampleCollectionWeek == "Week.1" & NICUSpeciesNR$SampleType == "Groin")
Week3StoolSpecies<-subset(NICUSpeciesNR, NICUSpeciesNR$SampleCollectionWeek == "Week.3" & NICUSpeciesNR$SampleType == "Stool")
Week3AxillaSpecies<-subset(NICUSpeciesNR, NICUSpeciesNR$SampleCollectionWeek == "Week.3" & NICUSpeciesNR$SampleType == "Axilla")
Week3GroinSpecies<-subset(NICUSpeciesNR, NICUSpeciesNR$SampleCollectionWeek == "Week.3" & NICUSpeciesNR$SampleType == "Groin")


Week1Species<-subset(NICUSpeciesNR, NICUSpeciesNR$SampleCollectionWeek == "Week.1")
Week3Species<-subset(NICUSpeciesNR, NICUSpeciesNR$SampleCollectionWeek == "Week.3")


# check them all:



for (i in 18:ncol(Week1StoolSpecies)){
  
  tryCatch({
    
    a = pairwise.wilcox.test(Week1StoolSpecies[,i], Week1StoolSpecies$Location)[3]
    b = pairwise.wilcox.test(Week3StoolSpecies[,i], Week3StoolSpecies$Location)[3]
    
    c = pairwise.wilcox.test(Week1AxillaSpecies[,i], Week1AxillaSpecies$Location)[3]
    d = pairwise.wilcox.test(Week3AxillaSpecies[,i], Week3AxillaSpecies$Location)[3]
    
    e = pairwise.wilcox.test(Week1GroinSpecies[,i], Week1GroinSpecies$Location)[3]
    f = pairwise.wilcox.test(Week3GroinSpecies[,i], Week3GroinSpecies$Location)[3]
    
    head = (colnames(Week1StoolSpecies[i]))
    
    
    SpeciesPvalueTableRow<-cbind(a, b, c, d, e, f)
    row.names(SpeciesPvalueTableRow) = head
    
    SpeciesPvalueTable<-rbind(SpeciesPvalueTable, SpeciesPvalueTableRow)        
    
    
  }, error=function(e){cat("ERROR :",conditionMessage(e), "\n")})}

SpeciesPvalueTable<-as.data.frame(SpeciesPvalueTable)

names(SpeciesPvalueTable)<-c("Week1.Stool", "Week3.Stool", "Week1.Axilla", "Week3.Axilla", "Week1.Groin", "Week3.Groin")

SpeciesPvalueTableMatrix<-as.matrix(SpeciesPvalueTable)

write.csv(SpeciesPvalueTableMatrix, file = "SpeciesPValueTable.csv")


pairwise.wilcox.test(Week1StoolSpecies$Streptococcus.agalactiae, Week1StoolSpecies$Location)
pairwise.wilcox.test(Week3StoolSpecies$Streptococcus.agalactiae, Week3StoolSpecies$Location)

pairwise.wilcox.test(Week1AxillaSpecies$Streptococcus.agalactiae, Week1AxillaSpecies$Location)
pairwise.wilcox.test(Week3AxillaSpecies$Streptococcus.agalactiae, Week3AxillaSpecies$Location)

pairwise.wilcox.test(Week1GroinSpecies$Streptococcus.agalactiae, Week1GroinSpecies$Location)
pairwise.wilcox.test(Week3GroinSpecies$Streptococcus.agalactiae, Week3GroinSpecies$Location)


# combine sample types:


for (i in 18:ncol(Week1Species)){
  
  tryCatch({
    
    a = pairwise.wilcox.test(Week1Species[,i], Week1Species$Location)[3]
    b = pairwise.wilcox.test(Week3Species[,i], Week3Species$Location)[3]
    
    head = (colnames(Week3Species[i]))
    
    
    SpeciesPvalueTableRow<-cbind(a, b)
    row.names(SpeciesPvalueTableRow) = head
    
    SpeciesPvalueTableNew<-rbind(SpeciesPvalueTableNew, SpeciesPvalueTableRow)        
    
    
  }, error=function(e){cat("ERROR :",conditionMessage(e), "\n")})}

SpeciesPvalueTableNew<-as.data.frame(SpeciesPvalueTableNew)

names(SpeciesPvalueTableNew)<-c("Week.1", "Week3")

SpeciesPvalueTableMatrixNew<-as.matrix(SpeciesPvalueTableNew)

write.csv(SpeciesPvalueTableMatrixNew, file = "SpeciesPValueTableWk1Wk3.csv")

save.image(file = "NICUData20250209")

##### which genus are most changed from week 1 to week 3 in no antibiotics grouppairwise.wilcox.test(DummyGenusNR$Streptococcus[DummyGenusNR$SampleType== "Stool" & DummyGenusNR$SampleCollectionWeek=="Week.1"], DummyGenusNR$PostNatalAntibiotics[DummyGenusNR$SampleType== "Stool" & DummyGenusNR$SampleCollectionWeek=="Week.1"])

NoAbxGenus<-subset(NICUGenusNR, NICUGenusNR$PostNatalAntibiotics == "No.Infant.Abx" & NICUGenusNR$SampleType == "Stool")
NoAbxGenus<-NoAbxGenus[, c(1,8,18:ncol(NoAbxGenus))]


Wilcox<-pairwise.wilcox.test(NoAbxGenus[,3], NoAbxGenus[,2])
WilcoxTable<-as.data.frame(Wilcox$p.value)
names(WilcoxTable)<-names(NoAbxGenus)[3]
row.names(WilcoxTable)<-"Wilcox"


for (i in 4:length(colnames(NoAbxGenus))){
  x<-pairwise.wilcox.test(NoAbxGenus[,i], NoAbxGenus[,2])
  y<-as.data.frame(x$p.value)
  names(y)<-names(NoAbxGenus)[i]
  WilcoxTable<-cbind(WilcoxTable, y)
}

WilcoxTable<-as.data.frame(t(WilcoxTable))
names(WilcoxTable)<-c("Unadjusted_p")
WilcoxTable$FDR<-p.adjust(WilcoxTable$Unadjusted_p, method = "fdr")
WilcoxTable$Genus<-row.names(WilcoxTable)
WilcoxTable<-WilcoxTable[order(WilcoxTable$Unadjusted_p),]



# Now let's figure out relative contribution of each species, fold difference between groups, and direction of difference

GenusCounts<-as.data.frame(colSums(NoAbxGenus[3:length(colnames(NoAbxGenus))]))
names(GenusCounts)<-"GenusTotal"
GenusCountsBigTable<-as.data.frame(colSums(NoAbxGenus[3:length(colnames(NoAbxGenus))]))
names(GenusCountsBigTable)<-"GenusTotal"
TotalCountsAll<-sum(GenusCountsBigTable$GenusTotal)
GenusCounts$Genus<-row.names(GenusCounts)

for (i in 1:length(rownames(GenusCounts))){
  GenusCounts$Fraction[i]<-GenusCounts[i,1] / TotalCountsAll
}

GenusGroupMeans<-aggregate(.~NoAbxGenus$SampleCollectionWeek, mean, data=NoAbxGenus[3:ncol(NoAbxGenus)])
row.names(GenusGroupMeans)<-GenusGroupMeans$`NoAbxGenus`
GenusGroupMeans<-as.data.frame(t(GenusGroupMeans[,2:ncol(GenusGroupMeans)]))
GenusGroupMeans$Genus<-row.names(GenusGroupMeans)
OverallMeans<-as.data.frame(colMeans(NoAbxGenus[3:ncol(NoAbxGenus)]))
OverallMeans$Genus<-row.names(OverallMeans)
names(OverallMeans)<-c("OverallMean", "Genus")
GenusGroupMeans<-merge(GenusGroupMeans, OverallMeans, by = "Genus", all.x = TRUE)
GenusGroupMeans$Week.1<-as.numeric(GenusGroupMeans$Week.1)
GenusGroupMeans$Week.3<-as.numeric(GenusGroupMeans$Week.3)
# GenusGroupMeans$FoldChange<-(GenusGroupMeans$Week.3 -GenusGroupMeans$Week.1) / GenusGroupMeans$OverallMean
# GenusGroupMeans$Log2<-log2(GenusGroupMeans$Week.3/GenusGroupMeans$OverallMean)
GenusGroupMeans$FC<-foldchange(GenusGroupMeans$Week.3, GenusGroupMeans$Week.1)
GenusGroupMeans$Log2<-foldchange2logratio(GenusGroupMeans$FC)

GenusGroupTable<-merge(GenusGroupMeans, GenusCounts, by = "Genus", all.x = TRUE)
StoolWeek1Week3WeekTable$Genus<-row.names(StoolWeek1Week3WeekTable)

# ok now we are fetching significant genus from ZINB-GLM
GenusGroupTable<-merge(GenusGroupTable, StoolWeek1Week3WeekTable, by = "Genus", all.x = TRUE)
GenusGroupTable<-GenusGroupTable[order(GenusGroupTable$Week.FDR),]
GenusGroupTable<-GenusGroupTable[,c(1,8,4, 2,3,5,6,9)]
names(GenusGroupTable)<-c("Genus", "Abundance", "OverallMean", "Week.1.Mean", "Week.3.Mean", "Fold.Change", "Log2FC", "FDR")
GenusGroupTable$Abundance<-GenusGroupTable$Abundance * 100
GenusGroupTable$Fold.Change[is.na(GenusGroupTable$Fold.Change)]<-0
GenusGroupTable<-GenusGroupTable[order(GenusGroupTable$Fold.Change, decreasing = TRUE),]

# Ok now we're using NBZI-GLM with FDR < 0.1 to find the genera that are significantly different from week 1 to 3
SigGenusGroupTable<-subset(GenusGroupTable, GenusGroupTable$FDR < 0.1)


metadata <- NoAbxGenus[, 1:2]
cts <- as.matrix(NoAbxGenus[, -(1:2)])
rownames(cts) <- metadata$Sample

cts_l2 <- glog2(cts)
grps <- NoAbxGenus$SampleCollectionWeek

res <- compute_ef_Week(d = cts_l2, g = NoAbxGenus$SampleCollectionWeek, min_shrink = 0.3)
names(res)[1]<-"Genus"

res <- arrange(left_join(res, SigGenusGroupTable), desc(abs(ef_shrunk))) %>%
  mutate(Genus = as_factor(Genus))

res<-subset(res, ! is.na(res$FDR)) 
resWks<-res

resWks$HigherIn<-ifelse(res$ef_shrunk > 0, "Week.1", "Week.3")
resWks$HigherIn<-factor(resWks$HigherIn, levels = c("Week.3", "Week.1"))

resWks$HigherFC<-ifelse(res$Fold.Change > 0, "Week.1", "Week.3")
resWks$HigherFC<-factor(resWks$HigherIn, levels = c("Week.3", "Week.1"))

resWks<-subset(resWks, ! resWks$Fold.Change== -Inf)

write.table(resWks, file = "SignificantWk1toWk3GenusTableNoAbx.csv", sep = ",", row.names = FALSE, col.names = TRUE,  quote = FALSE)
write.table(GenusGroupTable, file = "Wk1toWk3GenusTableNoAbx.csv", sep = ",", row.names = FALSE, col.names = TRUE,  quote = FALSE)

Week1toWk3NoAbx_Effect<- ggplot(filter(resWks, abs(resWks$ef_shrunk) > 0.3),
                                aes(x = reorder(Genus, ef_shrunk), y = ef_shrunk, fill = HigherIn)) +
  geom_bar(colour="black", stat="identity") +
  #guides(fill=FALSE) + scale_fill_hue(l=40) +
  guides(fill=guide_legend(override.aes=list(colour=NA))) + 
  coord_flip() +
  labs( x = NULL) +
  ylab(expression(atop("Effect Size"))) +
  guides(fill=guide_legend(title="More Abundant In:", )) +
  #geom_hline(yintercept = c(-0.2, 0.2), color = "firebrick", lty = 2) +
  theme(axis.text.y = element_text(size= 12, color="black")) +
  theme(legend.title = element_text(size=12)) +
  theme(legend.text = element_text(size = 11)) 

Week1toWk3NoAbx_Effect<-Week1toWk3NoAbx_Effect + scale_fill_manual(values = tableau.colors[c(1,3)])

ggsave(filename = "Week1toWk3GenusNoAbx_Effect.pdf", plot = Week1toWk3NoAbx_Effect, width = 10, 
       height = 12, limitsize = FALSE)

Week1toWk3NoAbx_FC<- ggplot(filter(resWks, abs(resWks$Fold.Change) >= 2 & abs(resWks$Fold.Change) <= 100),
                            aes(x = reorder(Genus, Fold.Change), y = Fold.Change, fill = HigherFC)) +
  geom_bar(colour="black", stat="identity") +
  #guides(fill=FALSE) + scale_fill_hue(l=40) +
  guides(fill=guide_legend(override.aes=list(colour=NA))) + 
  coord_flip() +
  labs( x = NULL) +
  ylab(expression(atop("Fold Change"))) +
  guides(fill=guide_legend(title="More Abundant In:", )) +
  #geom_hline(yintercept = c(-0.2, 0.2), color = "firebrick", lty = 2) +
  theme(axis.text.y = element_text(size= 12, color="black")) +
  theme(legend.title = element_text(size=12)) +
  theme(legend.text = element_text(size = 11)) 

Week1toWk3NoAbx_FC<-Week1toWk3NoAbx_FC + scale_fill_manual(values = tableau.colors[c(1,3)])
ggsave(filename = "Week1toWk3GenusNoAbx_FC.pdf", plot = Week1toWk3NoAbx_FC, width = 10, 
       height = 12, limitsize = FALSE)


FoldChangeGenusNoAbx<-subset(GenusGroupTable, abs(GenusGroupTable$Fold.Change) > 2 & GenusGroupTable$FDR < 0.15)
NOGood<-grep("-Inf", FoldChangeGenusAbx$Fold.Change)

matrixtable<-FoldChangeGenusNoAbx[,6:7]
row.names(matrixtable)<-FoldChangeGenusNoAbx$Genus
matrixtable<-as.matrix(matrixtable)

# matrixtable<-matrixtable[-NOGood,]
heatmap.2(t(matrixtable), Rowv = NULL, Colv = FALSE, col = diverging.colors(100), trace = "none", 
          tracecol= "black",hline = NULL, vline = NULL, margins=c(20,3),
          xlab = NULL, main = NULL)

ggplot(FoldChangeGenusNoAbx, aes(x=foldchange2logratio(Fold.Change), y = log2(p_Unadjusted))) + geom_point()




##### which genus are most changed from week 1 to week 3 in  antibiotics group

AbxGenus<-subset(NICUGenusNR, NICUGenusNR$PostNatalAntibiotics == "Infant.Abx" & NICUGenusNR$SampleType == "Stool")
AbxGenus<-AbxGenus[, c(1,8,18:ncol(AbxGenus))]


Wilcox<-pairwise.wilcox.test(AbxGenus[,3], AbxGenus[,2])
WilcoxTable<-as.data.frame(Wilcox$p.value)
names(WilcoxTable)<-names(AbxGenus)[3]
row.names(WilcoxTable)<-"Wilcox"


for (i in 4:length(colnames(AbxGenus))){
  x<-pairwise.wilcox.test(AbxGenus[,i], AbxGenus[,2])
  y<-as.data.frame(x$p.value)
  names(y)<-names(AbxGenus)[i]
  WilcoxTable<-cbind(WilcoxTable, y)
}

WilcoxTable<-as.data.frame(t(WilcoxTable))
names(WilcoxTable)<-c("Unadjusted_p")
WilcoxTable$FDR<-p.adjust(WilcoxTable$Unadjusted_p, method = "fdr")
WilcoxTable$Genus<-row.names(WilcoxTable)
WilcoxTable<-WilcoxTable[order(WilcoxTable$Unadjusted_p),]



# Now let's figure out relative contribution of each species, fold difference between groups, and direction of difference

GenusCounts<-as.data.frame(colSums(AbxGenus[3:length(colnames(AbxGenus))]))
names(GenusCounts)<-"GenusTotal"
GenusCountsBigTable<-as.data.frame(colSums(AbxGenus[3:length(colnames(AbxGenus))]))
names(GenusCountsBigTable)<-"GenusTotal"
TotalCountsAll<-sum(GenusCountsBigTable$GenusTotal)
GenusCounts$Genus<-row.names(GenusCounts)

for (i in 1:length(rownames(GenusCounts))){
  GenusCounts$Fraction[i]<-GenusCounts[i,1] / TotalCountsAll
}



GenusGroupMeans<-aggregate(.~AbxGenus$SampleCollectionWeek, mean, data=AbxGenus[3:ncol(AbxGenus)])
row.names(GenusGroupMeans)<-GenusGroupMeans$`AbxGenus`
GenusGroupMeans<-as.data.frame(t(GenusGroupMeans[,2:ncol(GenusGroupMeans)]))
GenusGroupMeans$Genus<-row.names(GenusGroupMeans)
OverallMeans<-as.data.frame(colMeans(AbxGenus[3:ncol(AbxGenus)]))
OverallMeans$Genus<-row.names(OverallMeans)
names(OverallMeans)<-c("OverallMean", "Genus")
GenusGroupMeans<-merge(GenusGroupMeans, OverallMeans, by = "Genus", all.x = TRUE)
GenusGroupMeans$Week.1<-as.numeric(GenusGroupMeans$Week.1)
GenusGroupMeans$Week.3<-as.numeric(GenusGroupMeans$Week.3)
# GenusGroupMeans$FoldChange<-(GenusGroupMeans$Week.3 -GenusGroupMeans$Week.1) / GenusGroupMeans$OverallMean
# GenusGroupMeans$Log2<-log2(GenusGroupMeans$Week.3/GenusGroupMeans$OverallMean)
GenusGroupMeans$FC<-foldchange(GenusGroupMeans$Week.3, GenusGroupMeans$Week.1)
GenusGroupMeans$Log2<-foldchange2logratio(GenusGroupMeans$FC)

GenusGroupTable<-merge(GenusGroupMeans, GenusCounts, by = "Genus", all.x = TRUE)
GenusGroupTable<-merge(GenusGroupTable, WilcoxTable, by = "Genus", all.x = TRUE)
GenusGroupTable<-GenusGroupTable[order(GenusGroupTable$Unadjusted_p),]
GenusGroupTable<-GenusGroupTable[,c(1,8,4, 2,3,5,9,10)]
names(GenusGroupTable)<-c("Genus", "Abundance", "OverallMean", "Week.1.Mean", "Week.3.Mean", "Fold.Change", "p_Unadjusted", "FDR")
GenusGroupTable$Abundance<-GenusGroupTable$Abundance * 100
GenusGroupTable$Fold.Change[is.na(GenusGroupTable$Fold.Change)]<-0
GenusGroupTable<-GenusGroupTable[order(GenusGroupTable$Fold.Change, decreasing = TRUE),]

SigGenusGroupTable<-subset(GenusGroupTable, GenusGroupTable$FDR < 0.1)

metadata <- AbxGenus[, 1:2]
cts <- as.matrix(AbxGenus[, -(1:2)])
rownames(cts) <- metadata$Sample

cts_l2 <- glog2(cts)
grps <- AbxGenus$SampleCollectionWeek


res <- compute_ef_Week(d = cts_l2, g = AbxGenus$SampleCollectionWeek, min_shrink = 0.3)
names(res)[1]<-"Genus"

res <- arrange(left_join(res, SigGenusGroupTable), desc(abs(ef_shrunk))) %>%
  mutate(Genus = as_factor(Genus))

res<-subset(res, ! is.na(res$FDR)) 
resWks<-res

resWks$HigherIn<-ifelse(res$ef_shrunk > 0, "Week.1", "Week.3")
resWks$HigherIn<-factor(resWks$HigherIn, levels = c("Week.3", "Week.1"))

resWks$HigherFC<-ifelse(res$Fold.Change > 0, "Week.1", "Week.3")
resWks$HigherFC<-factor(resWks$HigherIn, levels = c("Week.3", "Week.1"))

resWks<-subset(resWks, ! resWks$Fold.Change== -Inf)

write.table(resWks, file = "SignificantWk1toWk3GenusTableNewAbxExposed.csv", sep = ",", row.names = FALSE, col.names = TRUE,  quote = FALSE)
write.table(GenusGroupTable, file = "Wk1toWk3GenusTableAbxExposed.csv", sep = ",", row.names = FALSE, col.names = TRUE,  quote = FALSE)



Week1toWk3Abx_Effect<- ggplot(filter(resWks, abs(resWks$ef_shrunk) > 0.3),
                              aes(x = reorder(Genus, ef_shrunk), y = ef_shrunk, fill = HigherIn)) +
  geom_bar(colour="black", stat="identity") +
  #guides(fill=FALSE) + scale_fill_hue(l=40) +
  guides(fill=guide_legend(override.aes=list(colour=NA))) + 
  coord_flip() +
  labs( x = NULL) +
  ylab(expression(atop("Effect Size"))) +
  guides(fill=guide_legend(title="More Abundant In:", )) +
  #geom_hline(yintercept = c(-0.2, 0.2), color = "firebrick", lty = 2) +
  theme(axis.text.y = element_text(size= 12, color="black")) +
  theme(legend.title = element_text(size=12)) +
  theme(legend.text = element_text(size = 11)) 

Week1toWk3Abx_Effect<-Week1toWk3Abx_Effect + scale_fill_manual(values = tableau.colors[c(1,3)])

ggsave(filename = "Week1toWk3GenusABXOnly_Effect.pdf", plot = Week1toWk3Abx_Effect, width = 10, 
       height = 12, limitsize = FALSE)



Week1toWk3Abx_FC<- ggplot(filter(resWks, abs(resWks$Fold.Change) > 2),
                          aes(x = reorder(Genus, Fold.Change), y = Fold.Change, fill = HigherFC)) +
  geom_bar(colour="black", stat="identity") +
  #guides(fill=FALSE) + scale_fill_hue(l=40) +
  guides(fill=guide_legend(override.aes=list(colour=NA))) + 
  coord_flip() +
  labs( x = NULL) +
  ylab(expression(atop("Fold Change"))) +
  guides(fill=guide_legend(title="More Abundant In:", )) +
  #geom_hline(yintercept = c(-0.2, 0.2), color = "firebrick", lty = 2) +
  theme(axis.text.y = element_text(size= 12, color="black")) +
  theme(legend.title = element_text(size=12)) +
  theme(legend.text = element_text(size = 11)) 

Week1toWk3Abx_FC<-Week1toWk3Abx_FC + scale_fill_manual(values = tableau.colors[c(1,3)])

ggsave(filename = "Week1toWk3GenusABXOnly_FC.pdf", plot = Week1toWk3Abx_FC, width = 10, 
       height = 12, limitsize = FALSE)



## ok how about heatmaps like with the pathways. Select for FC > 2

FoldChangeGenusAbx<-subset(GenusGroupTable, abs(GenusGroupTable$Fold.Change) > 2 & GenusGroupTable$FDR < 0.1)
# NOGood<-grep("-Inf", FoldChangeGenusAbx$Fold.Change)
FoldChangeGenusAbx$Abx<-"Infant.Abx"

FoldChangeGenusNoAbx$Abx<-"No.Infant.Abx"

FoldChangeGenusAbxAndNoAbx<-rbind(FoldChangeGenusNoAbx, FoldChangeGenusAbx)
NOGood<-grep("-Inf", FoldChangeGenusAbxAndNoAbx$Fold.Change)
# FoldChangeGenusAbxAndNoAbx<-FoldChangeGenusAbxAndNoAbx[-NOGood,]

FoldChangeGenusAbxAndNoAbx$Fold.Change<-foldchange2logratio(FoldChangeGenusAbxAndNoAbx$Fold.Change)

matrixtable<-FoldChangeGenusAbxAndNoAbx[,6:7]
row.names(matrixtable)<-paste(FoldChangeGenusAbxAndNoAbx$Genus, FoldChangeGenusAbxAndNoAbx$Abx, sep = ".")
matrixtable<-as.matrix(matrixtable)
heatmap.2(t(matrixtable), Rowv = NULL, Colv = FALSE, col = diverging.colors(100), trace = "none", 
          tracecol= "black",hline = NULL, vline = NULL, margins=c(20,3),
          xlab = NULL, main = NULL)
