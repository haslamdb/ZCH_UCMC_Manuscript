# Here's the updated SampleKey format
NewSampleKey<-read.csv("AllNICUSampleKeyRevised20250206_for_HangzhouCincinnatiSamples.csv", header = TRUE, stringsAsFactors = FALSE)
NewSampleKey$PostNatalAntibiotics<-factor(NewSampleKey$PostNatalAntibiotics, levels = c("No.Infant.Abx", "Infant.Abx"))
NewSampleKey$PostNatalAntibioticsNew<-factor(NewSampleKey$PostNatalAntibioticsNew, levels = c("No.Infant.Abx", "Low.Infant.Abx", "High.Infant.Abx"))
NewSampleKey$Sample<-NewSampleKey$SampleID

NewSampleKeyAbx<-read.csv("SampleKeyAbx20231106.csv", header = FALSE)


AllSampleKey<-NewSampleKey


AllSampleKey$GestationCohort<-factor(AllSampleKey$GestationCohort, levels =  c("23-27 Weeks", "28-32 Weeks", "33-36 Weeks", "37-42 Weeks"))

AllSampleKey$SampleCollectionWeek<-factor(AllSampleKey$SampleCollectionWeek, levels = c("Week.1", "Week.3"))
AllSampleKey$SampleType<-factor(AllSampleKey$SampleType, levels = c("Nares", "Axilla", "Groin", "Stool"))
AllSampleKey$MaternalAntibiotics<-factor(AllSampleKey$MaternalAntibiotics, levels = c("No.Mat.Abx", "Mat.Abx"))
AllSampleKey$TypeWeek<-paste(AllSampleKey$SampleType, AllSampleKey$SampleCollectionWeek, sep = "-")
AllSampleKey$TypeAbx<-paste(AllSampleKey$SampleType, AllSampleKey$PostNatalAntibiotics, sep = "-")
AllSampleKey$TypeAbx<-factor(AllSampleKey$TypeAbx, levels = c("Axilla-No.Infant.Abx", "Axilla-Infant.Abx", 
                                                              "Groin-No.Infant.Abx", "Groin-Infant.Abx", "Stool-No.Infant.Abx", "Stool-Infant.Abx", "Nares-No.Infant.Abx", 
                                                              "Nares-Infant.Abx"))

AllSampleKey<-subset(AllSampleKey, ! duplicated(AllSampleKey$SampleID) & AllSampleKey$SampleType %in% c("Stool", "Groin", "Axilla"))




HumanReactiveSpecies<-read.csv("../HumanReactiveKraken2.csv")
HumanReactiveSpecies$SubSpecies<-gsub("_", ".", HumanReactiveSpecies$SubSpecies)
HumanReactiveSpecies$SubSpecies<-gsub("-", ".", HumanReactiveSpecies$SubSpecies)
HumanReactiveSpecies$SubSpecies<-gsub("/", ".", HumanReactiveSpecies$SubSpecies, fixed = TRUE)
HumanReactiveSpecies$SubSpecies<-gsub("X,", "", HumanReactiveSpecies$SubSpecies, fixed = TRUE)
HumanReactiveSpecies$SubSpecies<-gsub("[", "", HumanReactiveSpecies$SubSpecies, fixed = TRUE)
HumanReactiveSpecies$SubSpecies<-gsub("]", "", HumanReactiveSpecies$SubSpecies, fixed = TRUE)
HumanReactiveSpecies$SubSpecies<-gsub("..", ".", HumanReactiveSpecies$SubSpecies, fixed = TRUE)
HumanReactiveSpecies$SubSpecies<-gsub("..", ".", HumanReactiveSpecies$SubSpecies, fixed = TRUE)
HumanReactiveSpecies$SubSpecies<-gsub(" ", ".", HumanReactiveSpecies$SubSpecies, fixed = TRUE)

#Kraken2 alignments
setwd("C:/Users/dbhas/OneDrive/Documents/Alignments/KrakenAlignments/Kraken2")
setwd("C:/Users/HASI9S/OneDrive/Documents/Alignments/KrakenAlignments/Kraken2")
# setwd("P:/Documents/Alignments/KrakenAlignments/Kraken2")
# setwd("~/pCloudDrive/Documents/Alignments/KrakenAlignments/Kraken2")

#Kraken2 alignments

AllK2Files<-list.files()
SpeciesFileList3<-grep("_species_abundance.txt", AllK2Files)
K2SpeciesFiles<-AllK2Files[SpeciesFileList3]
K2FileList<-gsub("_species_abundance.txt", "", K2SpeciesFiles)


NICUFiles<-as.character(AllSampleKey$SampleID)
K2NICUFiles<-K2FileList[K2FileList %in% NICUFiles]
ALLNICUFiles<-K2NICUFiles



MissingNICU<-subset(NICUFiles, ! NICUFiles %in% c(K2NICUFiles))
# write.csv(MissingNICU, file = "MissingNICUFiles20210730.txt")


# Get all the species counts
# get the files
for(f in 1:length(K2NICUFiles)){
  fnr = paste(K2NICUFiles[f], "_species_abundance.txt", sep = "")
  filename<-K2NICUFiles[f]
  
  assign(filename, read.csv(paste(fnr, sep=""),
                            sep="\t", header = TRUE, stringsAsFactors = FALSE))
}

# Bracken counts

NewSpeciesTable<-read.csv(paste(K2NICUFiles[1], "_species_abundance.txt", sep = ""), sep = "\t", header = TRUE)
names(NewSpeciesTable)[1]<-"Species"
NewSpeciesTable$Species<-gsub(" ", ".", NewSpeciesTable$Species, fixed = TRUE)
NewSpeciesTable$Species<-gsub("..", ".", NewSpeciesTable$Species, fixed = TRUE)
NewSpeciesTable$Species<-gsub("_", ".", NewSpeciesTable$Species, fixed = TRUE)
NewSpeciesTable$Species<-gsub("-", ".", NewSpeciesTable$Species, fixed = TRUE)
NewSpeciesTable$Species<-gsub("/", ".", NewSpeciesTable$Species, fixed = TRUE)
NewSpeciesTable$Species<-gsub("X,", "", NewSpeciesTable$Species, fixed = TRUE)
NewSpeciesTable$Species<-gsub("[", "", NewSpeciesTable$Species, fixed = TRUE)
NewSpeciesTable$Species<-gsub("]", "", NewSpeciesTable$Species, fixed = TRUE)
NewSpeciesTable<-as.data.frame(NewSpeciesTable[,c(1,6)])
names(NewSpeciesTable)<-c("Species",paste(paste(K2NICUFiles[1])))
NewSpeciesTable<-subset(NewSpeciesTable, ! duplicated(NewSpeciesTable$Species))


for ( i in 2:length(K2NICUFiles)){
  x <- get( K2NICUFiles[i] )
  x<-as.data.frame(x[,c(1,6)])
  names( x ) <- c("Species", paste(K2NICUFiles[i]))
  x$Species<-gsub(" ", ".", x$Species, fixed = TRUE)
  x$Species<-gsub("..", ".", x$Species, fixed = TRUE)
  x$Species<-gsub("_", ".", x$Species, fixed = TRUE)
  x$Species<-gsub("-", ".", x$Species, fixed = TRUE)
  x$Species<-gsub("/", ".", x$Species, fixed = TRUE)
  x$Species<-gsub("X,", "", x$Species, fixed = TRUE)
  x$Species<-gsub("[", "", x$Species, fixed = TRUE)
  x$Species<-gsub("]", "", x$Species, fixed = TRUE)
  x<-subset(x, ! duplicated(x$Species))
  NewSpeciesTable<-merge(x, NewSpeciesTable, by = "Species", all = TRUE)
  NewSpeciesTable[is.na(NewSpeciesTable)]<-0
}

row.names(NewSpeciesTable)<-NewSpeciesTable$Species

## Ok, here we're going to add the BSI count data

### Alright now we're going to add the bloodstream infection data and then calculate Bray Curtis dissimilarity

setwd("C:/Users/dbhas/OneDrive/Documents/Code/Metagenomics/Yanping/NICU_Microbiome/Hangzhou/NoHumanDNA20220929")
setwd("C:/Users/HASI9S/Documents/Code/Metagenomics/Yanping/NICU_Microbiome/Hangzhou/NoHumanDNA20220929")
setwd("/home/david/Documents/Code/Metagenomics/Yanping/NICU_Microbiome/Hangzhou/NoHumanDNA20220929")


NewSpeciesTable$Species<-NULL
NewSpeciesTable[is.na(NewSpeciesTable)]<-0

rm(list = K2NICUFiles)

AllSpeciesRaw<-t(NewSpeciesTable)
Speciesdf<-as.data.frame(AllSpeciesRaw)
Speciesdf$SampleID<-row.names(AllSpeciesRaw)
SpeciesTable<-as.data.frame(t(NewSpeciesTable))


# setwd("C:/Users/dbhas/OneDrive/Documents/Code/Metagenomics/Yanping/NICU_Microbiome/Hangzhou/NoHumanDNA20220929")
# setwd("C:/Users/HASI9S/OneDrive/Documents/Code/Metagenomics/Yanping/NICU_Microbiome/Hangzhou/NoHumanDNA20220929")
# setwd("/home/david/Documents/Code/Metagenomics/Yanping/NICU_Microbiome/Hangzhou/NoHumanDNA20220929")


# Now filter species to be present in at least 10% of samples, account for > 0.005% of counts,
# And get rid of samples that have < 1 million reads

library(vegan)

NICUSpecies<- as.data.frame(AllSpeciesRaw)
SaliniCol<-grep("Salinibacter.ruber", names(NICUSpecies))
names(NICUSpecies)[SaliniCol]
NICUSpecies<-NICUSpecies[,-SaliniCol]
HumanCol<-grep("Homo.sapiens", names(NICUSpecies))
names(NICUSpecies)[HumanCol]
NICUSpecies<-NICUSpecies[,-HumanCol]
HumanReactiveCols<-which(colnames(NICUSpecies) %in% HumanReactiveSpecies$SubSpecies)
NICUSpecies<-NICUSpecies[-HumanReactiveCols]


TenPercentCutoff<-floor(nrow(NICUSpecies)/20)

NonZeroCounts<-list()
for (i in 1:ncol(NICUSpecies)){
  NonZeroCounts[i]<-length(which(NICUSpecies[,i] > 0))
  
}

TenPercentNotZero<-which(NonZeroCounts >= TenPercentCutoff)
NICUSpeciesReduced<-NICUSpecies[,TenPercentNotZero]
LowSamples<-which(rowSums(NICUSpeciesReduced) <= 100000)
NICUSpeciesReduced<-NICUSpeciesReduced[-LowSamples,]

#make sure the metadata does not get merged twice (only run it once)
#check with head(name(NICUSpeciesNR), 40)

NICUSpeciesNR<-as.data.frame(t(noise.removal(t(NICUSpeciesReduced), 0.001)))
minCount<-min(rowSums(NICUSpeciesNR))
NICUSpeciesNR<-data.frame(rrarefy(NICUSpeciesReduced, minCount))
NICUSpeciesNR$SampleID<-row.names(NICUSpeciesNR)

# NICUSpeciesReduced and NICUSpeciesNR have all samples except the three with low read counts and all species except Human Reactive

# setwd("C:/Users/HASI9S/OneDrive/Documents/Code/Metagenomics/Yanping/NICU_Microbiome/Hangzhou/NoHumanDNA20220929")

NICUSpeciesNR<-merge(AllSampleKey, NICUSpeciesNR, by = "SampleID", all.x = TRUE)

## There's a lot of human DNA in Hangzhou samples
HumanDNA<-ggplot(NICUSpeciesNR,
                 aes(x=Location, y=as.numeric(Homo.sapiens ), fill=Location)) + geom_boxplot(lwd=1,aes(color=factor(Location),fill = NA), outlier.size = 3)  +
  stat_summary(fun=mean, geom="point", shape=5, size=8) + scale_fill_manual(values = col) + scale_colour_manual(values = col ) +
  geom_point(size=4,aes(color = factor(Location))) + xlab(NULL) + facet_grid(rows = SampleType ~ SampleCollectionWeek) +
  ylab("Homo.sapiens \n") + paramsBoxWide() + scale_y_log10() + xlab(NULL)  +
  theme(axis.title.y = element_text(face = "italic"))
ggsave(HumanDNA, file = "HumanDNA.pdf", width = 12, height = 8, limitsize = FALSE)

LongSpeciesTable<-melt(NICUSpeciesNR, id.vars = c("SampleID", "PatientID", "Sample", "SampleType", "Location", "GestationalAge", "GestationCohort", "GestationTime", "SampleCollectionWeek", "PostNatalAntibiotics", "MaternalAntibiotics", "Week1Abx", "Week3Abx", "PostNatalAntibioticsNew", "TypeWeek", "TypeAbx", "AnyMilk"))
names(LongSpeciesTable)[c(18,19)]<-c("Species", "Count")
LongSpeciesTable<-subset(LongSpeciesTable, ! LongSpeciesTable$Count == 0)
write.csv(LongSpeciesTable, file = "LongSpeciesTableNoHuman.csv")

NICUSpeciesReduced$SampleID<-row.names(NICUSpeciesReduced) 
NICUSpeciesReduced<-subset(NICUSpeciesReduced, ! duplicated(NICUSpeciesReduced$SampleID))

NonRarefiedSpecies<-merge(AllSampleKey, NICUSpeciesReduced, by = "SampleID", all.y = TRUE)
NonRarefiedSpecies<-subset(NonRarefiedSpecies, ! duplicated(NonRarefiedSpecies$SampleID))


# Get diversity of all samples:
NonRarefiedSpecies<-subset(NonRarefiedSpecies, ! is.na(NonRarefiedSpecies$SampleType))
row.names(NonRarefiedSpecies)<-NonRarefiedSpecies$SampleID
# NonRarefiedSpecies<-subset(NonRarefiedSpecies, NonRarefiedSpecies$Location == "Cincinnati" & NonRarefiedSpecies$SampleType == "Stool")
SampleData<-NonRarefiedSpecies[, 1:17]
NonRarefiedSpecies<-NonRarefiedSpecies[, 18:ncol(NonRarefiedSpecies)]

DiversitySpecies<-NonRarefiedSpecies

#where initial issues are:
# Calculate sample diversity using 'diversity' function in Vegan
H <- data.frame(diversity(DiversitySpecies))
simpson <- data.frame(diversity(DiversitySpecies, "simpson"))
shannon<-data.frame(diversity(DiversitySpecies, "shannon"))
invsimp <- data.frame(diversity(DiversitySpecies, "inv"))
alpha <- data.frame(fisher.alpha(DiversitySpecies))
## Species richness (S) and Pielou's evenness (J):
S <- data.frame(specnumber(DiversitySpecies))
J <- data.frame(H/log(S))
Diversity<-cbind(simpson, shannon, invsimp, alpha, S, J)
Diversity$SampleID<-row.names(Diversity)
names(Diversity)<-c("Simpson", "Shannon",  "InvSimpson", "Alpha", "SpeciesNo", "Evenness", "SampleID")
Diversity<-Diversity[,c(7,1,2,3,4, 5, 6)]
pairs(cbind(simpson, shannon, alpha, J, S), pch="+", col="blue")


Diversity<-merge(SampleData, Diversity, by = "SampleID", all.y  = TRUE)
Diversity<-subset(Diversity, ! duplicated(Diversity$SampleID))


ShannonLocation<-ggplot(subset(Diversity, Diversity$Location %in% c("Cincinnati", "Hangzhou") & Diversity$SampleCollectionWeek == "Week.1" &
                                 Diversity$SampleType == "Stool"),
                        aes(x=Location, y=Shannon, fill=Location)) + geom_boxplot(lwd=1,aes(color=factor(Location),fill = NA), outlier.size = 3)  +
  stat_summary(fun=mean, geom="point", shape=5, size=8) + scale_fill_manual(values = col) + 
  geom_point(size=4,aes(color = factor(Location))) + ylab("Shannon Diversity Index : Stool\n") + xlab(NULL)  +  paramsBox() + labs( caption = "*p = 0.036")#+ facet_grid(. ~ SampleType)
ggsave(filename = "ShannonStoolLocation.pdf", plot = ShannonLocation, width = 5,
       height = 8, limitsize = FALSE)

StoolDiversity<-subset(Diversity, Diversity$Location %in% c("Cincinnati", "Hangzhou") & Diversity$SampleCollectionWeek == "Week.1" & 
                         Diversity$SampleType == "Stool")
pairwise.wilcox.test(StoolDiversity$Shannon, StoolDiversity$Location)

col<-tableau.colors
ShannonCollectionWeek<-ggplot(subset(Diversity, Diversity$SampleType %in% c("Axilla", "Groin", "Stool")),
                              aes(x=SampleCollectionWeek, y=Shannon, fill=TypeWeek)) + geom_boxplot(lwd=1,aes(color=factor(TypeWeek),fill = NA), outlier.size = 3)  +
  stat_summary(fun=mean, geom="point", shape=5, size=8) + scale_fill_manual(values = col) + scale_colour_tableau() +
  geom_point(size=4,aes(color = factor(TypeWeek))) + ylab("Shannon Diversity Index\n")  + xlab(NULL) +  paramsBox()  + facet_grid(Location ~ SampleType) # labs( caption = "*p = 0.034")#
ShannonCollectionWeek<-ShannonCollectionWeek + scale_fill_manual(values = col)
ggsave(filename = "ShannonCollectionWeek.pdf", plot = ShannonCollectionWeek, width = 14, 
       height = 12, limitsize = FALSE)

DiversityAxilla<-subset(Diversity, Diversity$SampleType == "Axilla")
pairwise.wilcox.test(DiversityAxilla$Shannon, DiversityAxilla$SampleCollectionWeek)
DiversityGroin<-subset(Diversity, Diversity$SampleType == "Groin")
pairwise.wilcox.test(DiversityGroin$Shannon, DiversityGroin$SampleCollectionWeek)
DiversityStool<-subset(Diversity, Diversity$SampleType == "Stool")
pairwise.wilcox.test(DiversityStool$Shannon, DiversityStool$SampleCollectionWeek)

kruskal.test(DiversityStool$Shannon, DiversityStool$SampleCollectionWeek)
Dunn2<-dunn.test(DiversityStool$Shannon, DiversityStool$SampleCollectionWeek, method = "bonferroni")



col=tableau.colors[c(4:6)]
DiversityAgeLocation<-ggplot(subset(Diversity, Diversity$SampleType %in% c("Stool") &
                                      ! is.na(Diversity$GestationCohort)),
                             aes(x=GestationCohort, y=Shannon, fill=GestationCohort)) + geom_boxplot(lwd=1,aes(color=factor(GestationCohort),fill = NA), outlier.size = 3)  +
  stat_summary(fun=mean, geom="point", shape=5, size=8) + scale_fill_manual(values = col) + scale_colour_tableau() +
  geom_point(size=4,aes(color = factor(GestationCohort))) + ylab("Shannon Diversity Index : Stool\n")  +  paramsAngled() +
  xlab(NULL) + facet_grid(rows = Location ~ SampleCollectionWeek) 
ggsave(DiversityAgeLocation, file = "ShannonGestationalAge.pdf", width = 7, height = 12, limitsize = FALSE)


DiversityLocationSample<-ggplot(subset(Diversity, 
                                       ! is.na(Diversity$GestationCohort)),
                                aes(x=Location, y=Shannon, fill=Location)) + geom_boxplot(lwd=1,aes(color=factor(Location),fill = NA), outlier.size = 3)  +
  stat_summary(fun=mean, geom="point", shape=5, size=8) + scale_fill_manual(values = col) + scale_colour_tableau() +
  geom_point(size=4,aes(color = factor(Location))) + ylab("Shannon Diversity Index\n")  +  paramsAngled() +
  xlab(NULL) + facet_grid(rows = SampleType ~ SampleCollectionWeek) 
ggsave(DiversityLocationSample, file = "ShannonLocationType.pdf", width = 7, height = 12, limitsize = FALSE)

col<-c(tableau.colors[1], tableau.colors[1],tableau.colors[2],tableau.colors[2],tableau.colors[3],tableau.colors[3])
DiversityAgeSampleType<-ggplot(subset(Diversity, Diversity$SampleType %in% c("Axilla", "Groin", "Stool") &
                                        ! is.na(Diversity$GestationCohort)),
                               aes(x=GestationCohort, y=Shannon, fill=TypeWeek)) + geom_boxplot(lwd=1,aes(color=factor(TypeWeek),fill = NA), outlier.size = 3)  +
  stat_summary(fun=mean, geom="point", shape=5, size=8) + scale_fill_manual(values = col) + scale_colour_tableau() +
  geom_point(size=4,aes(color = factor(TypeWeek))) + ylab("Shannon Diversity Index\n")  +  paramsAngled() +
  xlab(NULL) + facet_grid(. ~ SampleType  + SampleCollectionWeek + Location )
ggsave(DiversityAgeSampleType, file = "ShannonAgeSampleWeek.pdf", width = 20, height = 8, limitsize = FALSE)


pairwise.wilcox.test(DiversityAxilla$Shannon, DiversityAxilla$GestationCohort)
pairwise.wilcox.test(DiversityGroin$Shannon, DiversityGroin$GestationCohort)
pairwise.wilcox.test(DiversityStool$Shannon, DiversityStool$GestationCohort)

DiversityStoolAbx<-ggplot(subset(Diversity, Diversity$SampleType %in% c("Stool")),
                          aes(x=PostNatalAntibioticsNew, y=Shannon, fill=PostNatalAntibioticsNew)) + geom_boxplot(lwd=1,aes(color=factor(PostNatalAntibioticsNew),fill = NA), outlier.size = 3)  +
  stat_summary(fun=mean, geom="point", shape=5, size=8) + scale_fill_manual(values = col) + scale_colour_tableau() +
  geom_point(size=4,aes(color = factor(PostNatalAntibioticsNew))) + ylab("Shannon Diversity Index : Stool\n")  +  paramsBoxWide() +
  xlab(NULL) + facet_grid(. ~ SampleCollectionWeek + Location)
ggsave(DiversityStoolAbx, file = "ShannonAgeAntibiotics.pdf", width = 20, height = 8, limitsize = FALSE)


DiversityStoolAbxMatAbx<-ggplot(subset(Diversity, Diversity$SampleType %in% c("Stool")),
                                aes(x=MaternalAntibiotics, y=Shannon, fill=MaternalAntibiotics)) + geom_boxplot(lwd=1,aes(color=factor(MaternalAntibiotics),fill = NA), outlier.size = 3)  +
  stat_summary(fun=mean, geom="point", shape=5, size=8) + scale_fill_manual(values = col) + scale_colour_tableau() +
  geom_point(size=4,aes(color = factor(MaternalAntibiotics))) + ylab("Shannon Diversity Index : Stool\n")  +  paramsBoxWide() +
  xlab(NULL) + facet_grid(SampleCollectionWeek ~ PostNatalAntibioticsNew + Location)
ggsave(DiversityStoolAbxMatAbx, file = "ShannonMatInfantAntibiotics.pdf", width = 20, height = 8, limitsize = FALSE)

DiversityStoolAbxMatAbx<-ggplot(subset(Diversity, Diversity$SampleType %in% c("Stool")),
                                aes(x=MaternalAntibiotics, y=Shannon, fill=MaternalAntibiotics)) + geom_boxplot(lwd=1,aes(color=factor(MaternalAntibiotics),fill = NA), outlier.size = 3)  +
  stat_summary(fun=mean, geom="point", shape=5, size=8) + scale_fill_manual(values = col) + scale_colour_tableau() +
  geom_point(size=4,aes(color = factor(MaternalAntibiotics))) + ylab("Shannon Diversity Index : Stool\n")  +  paramsBoxWide() +
  xlab(NULL) + facet_grid(SampleCollectionWeek ~ PostNatalAntibiotics + Location)
ggsave(DiversityStoolAbxMatAbx, file = "ShannonMatInfantAntibiotics.pdf", width = 20, height = 8, limitsize = FALSE)

DiversityStoolMatAbx<-ggplot(subset(Diversity, Diversity$SampleType %in% c("Stool")& Diversity$PostNatalAntibiotics == "No.Infant.Abx"),
                             aes(x=MaternalAntibiotics, y=Shannon, fill=MaternalAntibiotics)) + geom_boxplot(lwd=1,aes(color=factor(MaternalAntibiotics),fill = NA), outlier.size = 3)  +
  stat_summary(fun=mean, geom="point", shape=5, size=8) + scale_fill_manual(values = col) + scale_colour_tableau() +
  geom_point(size=4,aes(color = factor(MaternalAntibiotics))) + ylab("Shannon Diversity Index : Stool\n")  +  paramsBoxWide() +
  xlab(NULL) + facet_grid(. ~ SampleCollectionWeek + GestationCohort + Location)
ggsave(DiversityStoolMatAbx, file = "ShannonAgeMatAntibiotics.pdf", width = 20, height = 8, limitsize = FALSE)

DiversityStoolMatAbx<-ggplot(subset(Diversity, Diversity$SampleType %in% c("Stool") & Diversity$PostNatalAntibiotics == "No.Infant.Abx"),
                             aes(x=MaternalAntibiotics, y=Shannon, fill=MaternalAntibiotics)) + geom_boxplot(lwd=1,aes(color=factor(MaternalAntibiotics),fill = NA), outlier.size = 3)  +
  stat_summary(fun=mean, geom="point", shape=5, size=8) + scale_fill_manual(values = col) + scale_colour_tableau() +
  geom_point(size=4,aes(color = factor(MaternalAntibiotics))) + ylab("Shannon Diversity Index : Stool\n")  +  paramsBoxWide() +
  xlab(NULL) + facet_grid(PostNatalAntibioticsNew ~ SampleCollectionWeek + Location)
ggsave(DiversityStoolMatAbx, file = "ShannonStoolMatAntibioticsNoInfantAbx.pdf", width = 18, height = 8, limitsize = FALSE)

DiversityAxillaMatAbx<-ggplot(subset(Diversity, Diversity$SampleType %in% c("Axilla") & Diversity$PostNatalAntibiotics == "No.Infant.Abx"),
                              aes(x=MaternalAntibiotics, y=Shannon, fill=MaternalAntibiotics)) + geom_boxplot(lwd=1,aes(color=factor(MaternalAntibiotics),fill = NA), outlier.size = 3)  +
  stat_summary(fun=mean, geom="point", shape=5, size=8) + scale_fill_manual(values = col) + scale_colour_tableau() +
  geom_point(size=4,aes(color = factor(MaternalAntibiotics))) + ylab("Shannon Diversity Index : Axilla\n")  +  paramsBoxWide() +
  xlab(NULL) + facet_grid( PostNatalAntibiotics ~ SampleCollectionWeek + Location)
ggsave(DiversityAxillaMatAbx, file = "ShannonAxillaMatAntibiotics.pdf", width = 18, height = 8, limitsize = FALSE)

DiversityGroinaMatAbx<-ggplot(subset(Diversity, Diversity$SampleType %in% c("Groin") & Diversity$PostNatalAntibiotics == "No.Infant.Abx"),
                              aes(x=MaternalAntibiotics, y=Shannon, fill=MaternalAntibiotics)) + geom_boxplot(lwd=1,aes(color=factor(MaternalAntibiotics),fill = NA), outlier.size = 3)  +
  stat_summary(fun=mean, geom="point", shape=5, size=8) + scale_fill_manual(values = col) + scale_colour_tableau() +
  geom_point(size=4,aes(color = factor(MaternalAntibiotics))) + ylab("Shannon Diversity Index : Groin\n")  +  paramsBoxWide() +
  xlab(NULL) + facet_grid( PostNatalAntibiotics ~ SampleCollectionWeek + Location)
ggsave(DiversityGroinaMatAbx, file = "ShannonGroinMatAntibiotics.pdf", width = 18, height = 8, limitsize = FALSE)

col<-tableau.colors
DiversityAllTypesAbx<-ggplot(subset(Diversity, Diversity$SampleType %in% c("Axilla", "Groin", "Stool")),
                             aes(x=TypeAbx, y=Shannon, fill=TypeAbx)) + geom_boxplot(lwd=1,aes(color=factor(TypeAbx),fill = NA), outlier.size = 3)  +
  stat_summary(fun=mean, geom="point", shape=5, size=8) + scale_fill_manual(values = col) + scale_colour_tableau() +
  geom_point(size=4,aes(color = factor(TypeAbx))) + ylab("Shannon Diversity Index\n")  +  paramsAngled() +
  xlab(NULL) + facet_grid(SampleCollectionWeek ~ Location) + paramsAngled()
ggsave(DiversityAllTypesAbx, file = "DiversityAllTypesAbx.pdf", width = 20, height = 12, limitsize = FALSE)


DiversityWeek1Axilla<-subset(DiversityAxilla, DiversityAxilla$SampleCollectionWeek == "Week.1")
DiversityWeek3Axilla<-subset(DiversityAxilla, DiversityAxilla$SampleCollectionWeek == "Week.3")
DiversityWeek1Groin<-subset(DiversityGroin, DiversityGroin$SampleCollectionWeek == "Week.1")
DiversityWeek3Groin<-subset(DiversityGroin, DiversityGroin$SampleCollectionWeek == "Week.3")
DiversityWeek1Stool<-subset(DiversityStool, DiversityStool$SampleCollectionWeek == "Week.1")
DiversityWeek3Stool<-subset(DiversityStool, DiversityStool$SampleCollectionWeek == "Week.3")

pairwise.wilcox.test(DiversityWeek1Axilla$Shannon, DiversityWeek1Axilla$Location)
pairwise.wilcox.test(DiversityWeek3Axilla$Shannon, DiversityWeek3Axilla$Location)
pairwise.wilcox.test(DiversityWeek1Groin$Shannon, DiversityWeek1Groin$Location)
pairwise.wilcox.test(DiversityWeek3Groin$Shannon, DiversityWeek3Groin$Location)
pairwise.wilcox.test(DiversityWeek1Stool$Shannon, DiversityWeek1Stool$Location)
pairwise.wilcox.test(DiversityWeek3Stool$Shannon, DiversityWeek3Stool$Location)


DiversityCincinnatiAxilla<-subset(Diversity, Diversity$SampleType == "Axilla" & Diversity$Location == "Cincinnati")
DiversityHangzhouAxilla<-subset(Diversity, Diversity$SampleType == "Axilla" & Diversity$Location == "Hangzhou")
DiversityCincinnatiGroin<-subset(Diversity, Diversity$SampleType == "Groin" & Diversity$Location == "Cincinnati")
DiversityHangzhouGroin<-subset(Diversity, Diversity$SampleType == "Groin" & Diversity$Location == "Hangzhou")
DiversityCincinnatiStool<-subset(Diversity, Diversity$SampleType == "Stool" & Diversity$Location == "Cincinnati")
DiversityHangzhouStool<-subset(Diversity, Diversity$SampleType == "Stool" & Diversity$Location == "Hangzhou")


pairwise.wilcox.test(DiversityCincinnatiAxilla$Shannon, DiversityCincinnatiAxilla$SampleCollectionWeek)
pairwise.wilcox.test(DiversityHangzhouAxilla$Shannon, DiversityHangzhouAxilla$SampleCollectionWeek)
pairwise.wilcox.test(DiversityCincinnatiGroin$Shannon, DiversityCincinnatiGroin$SampleCollectionWeek)
pairwise.wilcox.test(DiversityHangzhouGroin$Shannon, DiversityHangzhouGroin$SampleCollectionWeek)
pairwise.wilcox.test(DiversityCincinnatiStool$Shannon, DiversityCincinnatiStool$SampleCollectionWeek)
pairwise.wilcox.test(DiversityHangzhouStool$Shannon, DiversityHangzhouStool$SampleCollectionWeek)


save.image(file = "NICUData20250209")

## Back to rarefied species

# NICUSpeciesNR<-merge(AllSampleKey, NICUSpeciesNR, by = "SampleID", all.y = TRUE) # we merged it above
NICUSpeciesNR<-subset(NICUSpeciesNR, ! duplicated(NICUSpeciesNR$SampleID))


# Analyze all Species samples and compare location:

LocationWeek1Species<-subset(NICUSpeciesNR, NICUSpeciesNR$Location %in% c("Cincinnati", "Hangzhou") & NICUSpeciesNR$SampleType == "Stool" &
                               NICUSpeciesNR$SampleCollectionWeek == "Week.1")

LocationWeek1Species<-LocationWeek1Species[,c(1,5,18:ncol(LocationWeek1Species))]
LocationWeek1Species<-na.omit(LocationWeek1Species)
#
LocationWeek1Species.mrpp<-mrpp(LocationWeek1Species[,3:ncol(LocationWeek1Species)], LocationWeek1Species$Location, distance = "bray")
LocationWeek1Species.mrpp

# MRPP for maternal antibiotics week 1
MaternalWeek1StoolSpecies<-subset(NICUSpeciesNR, NICUSpeciesNR$Location %in% c("Cincinnati", "Hangzhou") & NICUSpeciesNR$SampleType == "Stool" &
                                    NICUSpeciesNR$SampleCollectionWeek == "Week.1")
MaternalWeek1StoolSpecies<-MaternalWeek1StoolSpecies[,c(1,11,18:ncol(MaternalWeek1StoolSpecies))]
MaternalWeek1StoolSpecies<-na.omit(MaternalWeek1StoolSpecies)
MaternalAbxWeek1Species.mrpp<-mrpp(MaternalWeek1StoolSpecies[,3:ncol(MaternalWeek1StoolSpecies)], MaternalWeek1StoolSpecies$MaternalAntibiotics, distance = "bray")
MaternalAbxWeek1Species.mrpp

# MRPP for maternal antibiotics week 3
MaternalWeek3StoolSpecies<-subset(NICUSpeciesNR, NICUSpeciesNR$Location %in% c("Cincinnati", "Hangzhou") & NICUSpeciesNR$SampleType == "Stool" &
                                    NICUSpeciesNR$SampleCollectionWeek == "Week.3")
MaternalWeek3StoolSpecies<-MaternalWeek3StoolSpecies[,c(1,11,18:ncol(MaternalWeek1StoolSpecies))]
MaternalWeek3StoolSpecies<-na.omit(MaternalWeek3StoolSpecies)
MaternalAbxWeek3Species.mrpp<-mrpp(MaternalWeek3StoolSpecies[,3:ncol(MaternalWeek3StoolSpecies)], MaternalWeek3StoolSpecies$MaternalAntibiotics, distance = "bray")
MaternalAbxWeek3Species.mrpp



# 
# 
# 
# # get all pairwise.wilcox for every species and aggregate into a table:
# 
# 
Wilcox<-pairwise.wilcox.test(LocationWeek1Species[,3], LocationWeek1Species[,2])
WilcoxTable<-as.data.frame(Wilcox$p.value)
names(WilcoxTable)<-names(LocationWeek1Species)[3]
row.names(WilcoxTable)<-"Wilcox"

library(dunn.test)
Dunn<-dunn.test(LocationWeek1Species[,5], LocationWeek1Species[,2], method = "bonferroni")


for (i in 4:length(colnames(LocationWeek1Species))){
  x<-pairwise.wilcox.test(LocationWeek1Species[,i], LocationWeek1Species[,2])
  y<-as.data.frame(x$p.value)
  names(y)<-names(LocationWeek1Species)[i]
  WilcoxTable<-cbind(WilcoxTable, y)
}

WilcoxTable<-as.data.frame(t(WilcoxTable))
names(WilcoxTable)<-c("Unadjusted_p")
WilcoxTable$FDR<-p.adjust(WilcoxTable$Unadjusted_p, method = "fdr")
WilcoxTable$Species<-row.names(WilcoxTable)
WilcoxTable<-WilcoxTable[order(WilcoxTable$Unadjusted_p),]



SpeciesCounts<-as.data.frame(colSums(LocationWeek1Species[3:length(colnames(LocationWeek1Species))]))
names(SpeciesCounts)<-"SpeciesTotal"
SpeciesCountsBigTable<-as.data.frame(colSums(LocationWeek1Species[3:length(colnames(LocationWeek1Species))]))
names(SpeciesCountsBigTable)<-"SpeciesTotal"
TotalCountsAll<-sum(SpeciesCountsBigTable$SpeciesTotal)
SpeciesCounts$Species<-row.names(SpeciesCounts)

for (i in 1:length(rownames(SpeciesCounts))){
  SpeciesCounts$Fraction[i]<-SpeciesCounts[i,1] / TotalCountsAll
}


SpeciesGroupMeans<-aggregate(.~LocationWeek1Species$Location, mean, data=LocationWeek1Species[3:ncol(LocationWeek1Species)])
row.names(SpeciesGroupMeans)<-SpeciesGroupMeans$`LocationWeek1Species$Location`
SpeciesGroupMeans<-as.data.frame(t(SpeciesGroupMeans[,2:ncol(SpeciesGroupMeans)]))
SpeciesGroupMeans$Species<-row.names(SpeciesGroupMeans)
OverallMeans<-as.data.frame(colMeans(LocationWeek1Species[3:ncol(LocationWeek1Species)]))
OverallMeans$Species<-row.names(OverallMeans)
names(OverallMeans)<-c("OverallMean", "Species")
SpeciesGroupMeans<-merge(SpeciesGroupMeans, OverallMeans, by = "Species", all.x = TRUE)
SpeciesGroupMeans$Cincinnati<-as.numeric(SpeciesGroupMeans$Cincinnati)
SpeciesGroupMeans$Hangzhou<-as.numeric(SpeciesGroupMeans$Hangzhou)
# SpeciesGroupMeans$FoldChange<-(SpeciesGroupMeans$Hangzhou -SpeciesGroupMeans$Cincinnati) / SpeciesGroupMeans$OverallMean
# SpeciesGroupMeans$Log2<-log2(SpeciesGroupMeans$Hangzhou/SpeciesGroupMeans$OverallMean)
SpeciesGroupMeans$FC<-foldchange(SpeciesGroupMeans$Hangzhou, SpeciesGroupMeans$Cincinnati)
SpeciesGroupMeans$logratio<-foldchange2logratio(SpeciesGroupMeans$FC)

SpeciesGroupTable<-merge(SpeciesGroupMeans, SpeciesCounts, by = "Species", all.x = TRUE)
SpeciesGroupTable<-merge(SpeciesGroupTable, WilcoxTable, by = "Species", all.x = TRUE)
SpeciesGroupTable<-SpeciesGroupTable[order(SpeciesGroupTable$Unadjusted_p),]
SpeciesGroupTable<-SpeciesGroupTable[,c(1,8,4, 2,3,6,9,10)]
names(SpeciesGroupTable)<-c("Species", "Abundance", "OverallMean", "Cincinnati.Mean", "Hangzhou.Mean", "LogFC", "p_Unadjusted", "FDR")
SpeciesGroupTable$Abundance<-SpeciesGroupTable$Abundance * 100
SpeciesGroupTable$LogFC[is.na(SpeciesGroupTable$LogFC)]<-0



SpeciesGroupTable<-SpeciesGroupTable[order(SpeciesGroupTable$LogFC, decreasing = TRUE),]

SigSpeciesGroupTable<-subset(SpeciesGroupTable, SpeciesGroupTable$FDR < 0.1)

SigSpeciesGroupTable$LogFC[SigSpeciesGroupTable$LogFC == Inf]<-0
SigSpeciesGroupTable$LogFC[SigSpeciesGroupTable$LogFC == -Inf]<-0


SigSpeciesGroupTable$HigherIn<-ifelse(SigSpeciesGroupTable$LogFC > 0, "Hangzhou", "Cincinnati")
SigSpeciesGroupTable$HigherIn<-factor(SigSpeciesGroupTable$HigherIn, levels = c("Hangzhou", "Cincinnati"))

# SigSpeciesGroupTable$LogFC[SigSpeciesGroupTable$LogFC == Inf]<-10
# SigSpeciesGroupTable$LogFC[SigSpeciesGroupTable$LogFC == -Inf]<--10


FilteredSigSpeciesGroupTable<-subset(SigSpeciesGroupTable, SigSpeciesGroupTable$LogFC >= 8 | SigSpeciesGroupTable$LogFC <= -3)
LocationSpecies_LogFC<- ggplot(FilteredSigSpeciesGroupTable, 
                               aes(x = reorder(Species, LogFC), y = LogFC, fill = HigherIn)) +
  geom_bar(colour="black", stat="identity") +
  #guides(fill=FALSE) + scale_fill_hue(l=40) +
  guides(fill=guide_legend(override.aes=list(colour=NA))) +
  coord_flip() +
  labs( x = NULL) +
  ylab(expression(atop("Log Fold Change"))) +
  guides(fill=guide_legend(title="More Abundant In:")) +
  #geom_hline(yintercept = c(-0.2, 0.2), color = "firebrick", lty = 2) +
  theme(axis.text.y = element_text(size= 12, color="black")) +
  theme(legend.title = element_text(size=12)) +
  theme(legend.text = element_text(size = 11)) #+ ylim(-1.8, 1.5)
#labs( caption = "*FDR < 0.05 and SDA Effect Size > 0.5 \n All Samples")

LocationSpecies_LogFC<-LocationSpecies_LogFC + scale_fill_manual(values=col)

ggsave(filename = "LocationSpecies_LogFC_Week3.pdf", plot = LocationSpecies_LogFC, width = 10,
       height = 14, limitsize = FALSE)



SpeciesVolcFC<-  EnhancedVolcano(as.data.frame(SigSpeciesGroupTable),
                                 lab = SigSpeciesGroupTable$Species,
                                 x = 'LogFC',
                                 y = 'FDR',
                                 pCutoff = 10e-3,
                                 FCcutoff = 2,
                                 pointSize = 2,
                                 labSize = 3.5, 
                                 title = NULL,
                                 subtitle = NULL)

SpeciesVolcFC_Week3<-SpeciesVolcFC + ylim(0,10)

ggsave(filename = "SpeciesVolcanoLogFCPlotWeek3.pdf", plot = SpeciesVolcFC_Week3, limitsize = FALSE)

save.image(file = "NICUData20250209")


library(sda)

metadata <- LocationWeek1Species[, 1:2]
cts <- as.matrix(LocationWeek1Species[, -(1:2)])
rownames(cts) <- metadata$Sample

cts_l2 <- glog2(cts)
grps <- LocationWeek1Species$Location


res <- compute_ef_Location(d = cts_l2, g = LocationWeek1Species$Location, min_shrink = 0.3)

res <- arrange(left_join(res, SigSpeciesGroupTable), desc(abs(ef_shrunk))) %>%
  mutate(Species = as_factor(Species))

res<-subset(res, ! is.na(res$FDR))
resLocation<-res

resLocation$HigherIn<-ifelse(res$ef_shrunk > 0, "Hangzhou", "Cincinnati")
resLocation$HigherIn<-factor(resLocation$HigherIn, levels = c("Hangzhou", "Cincinnati"))


write.table(resLocation, file = "SignificantLocationTableNew.csv", sep = ",", row.names = FALSE, col.names = TRUE,  quote = FALSE)
write.table(SpeciesGroupTable, file = "LocationSpeciesTable.csv", sep = ",", row.names = FALSE, col.names = TRUE,  quote = FALSE)


Location_Effect<- ggplot(filter(resLocation, abs(resLocation$ef_shrunk) > 0.5),
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

Location_Effect<-Location_Effect + scale_fill_manual(values = tableau.colors[c(1,2)])


ggsave(filename = "Location_Effect_Week3.pdf", plot = Location_Effect, width = 10,
       height = 12, limitsize = FALSE)


LocationWeek1SpeciesTransformed<-as.data.frame(t(LocationWeek1Species))
LocationWeek1LocationRow<-LocationWeek1SpeciesTransformed[2,]
LocationWeek1SpeciesTransformedNoLocation<-as.data.frame(LocationWeek1SpeciesTransformed[-2,])
colnames(LocationWeek1SpeciesTransformedNoLocation)<-row.names(LocationWeek1SpeciesTransformedNoLocation)
CultureDataTransformed<-as.data.frame(t(CultureData))
MergedSpeciesData<-merge(LocationWeek1SpeciesTransformed, CultureDataTransformed, by = "")

col=tableau.colors[c(2,1)]
newPCA<-PCA(cts_l2, scale.unit = FALSE, ncp = 5, graph = FALSE)
pcaPlot<-fviz_pca_ind(newPCA,
                      geom.ind = "text",
                      pointsize = 1.5,
                      point.alph = 0.1,
                      title = "Unsupervised Principal Coordinate Analysis",
                      subtitle = "Week 3 Stool Samples Grouped by Location",
                      col.ind = grps,
                      addEllipses = TRUE,
                      ellipse.alpha = 0.2,
                      ellipse.type = "confidence",
                      ellipse.level = 0.95,
                      legend.title = "Location",
                      legend.size = 11,
                      mean.point = TRUE,
                      palette = col,
                      axes.linetype = "blank"
)
LocationPCA_week1<- pcaPlot + scale_fill_manual(values = tableau.colors[c(2,1)])

pdf("LocationPCAWeek1.pdf")
print(LocationPCA_week1)
dev.off()
# 
SigSpecies<-subset(resLocation, abs(resLocation$ef_shrunk) > 0.5)
SigSpeciesList<-as.list(as.character(SigSpecies$Species))

LocationBiplot<-fviz_pca_biplot(newPCA,label="var", habillage=factor(grps),
                                addEllipses=TRUE, ellipse.level=0.05 , select.var = list(name = c(SigSpeciesList)),
                                pointsize = 1, labelsize = 2, col.var = "gray24", palette  = col,
                                repel = TRUE, alpha.var = 0.1, title = "Location : PCA Biplot")

ggsave(plot = LocationBiplot, file = "LocationBiplot.pdf",  limitsize = FALSE)

## Location week 3 ###



LocationWeek3Species<-subset(NICUSpeciesNR, NICUSpeciesNR$Location %in% c("Cincinnati", "Hangzhou") & NICUSpeciesNR$SampleType == "Stool" &
                               NICUSpeciesNR$SampleCollectionWeek == "Week.1")

LocationWeek3Species<-LocationWeek3Species[,c(1,5,18:ncol(LocationWeek3Species))]
LocationWeek3Species<-na.omit(LocationWeek3Species)
#
LocationWeek3Species.mrpp<-mrpp(LocationWeek3Species[,3:ncol(LocationWeek3Species)], LocationWeek3Species$Location, distance = "bray")
LocationWeek3Species.mrpp

#
# 
# 
Wilcox<-pairwise.wilcox.test(LocationWeek3Species[,3], LocationWeek3Species[,2])
WilcoxTable<-as.data.frame(Wilcox$p.value)
names(WilcoxTable)<-names(LocationWeek3Species)[3]
row.names(WilcoxTable)<-"Wilcox"

library(dunn.test)
Dunn<-dunn.test(LocationWeek3Species[,5], LocationWeek3Species[,2], method = "bonferroni")


for (i in 4:length(colnames(LocationWeek3Species))){
  x<-pairwise.wilcox.test(LocationWeek3Species[,i], LocationWeek3Species[,2])
  y<-as.data.frame(x$p.value)
  names(y)<-names(LocationWeek3Species)[i]
  WilcoxTable<-cbind(WilcoxTable, y)
}

WilcoxTable<-as.data.frame(t(WilcoxTable))
names(WilcoxTable)<-c("Unadjusted_p")
WilcoxTable$FDR<-p.adjust(WilcoxTable$Unadjusted_p, method = "fdr")
WilcoxTable$Species<-row.names(WilcoxTable)
WilcoxTable<-WilcoxTable[order(WilcoxTable$Unadjusted_p),]



# Now let's figure out relative contribution of each species, fold difference between groups, and direction of difference

SpeciesCounts<-as.data.frame(colSums(LocationWeek3Species[3:length(colnames(LocationWeek3Species))]))
names(SpeciesCounts)<-"SpeciesTotal"
SpeciesCountsBigTable<-as.data.frame(colSums(LocationWeek3Species[3:length(colnames(LocationWeek3Species))]))
names(SpeciesCountsBigTable)<-"SpeciesTotal"
TotalCountsAll<-sum(SpeciesCountsBigTable$SpeciesTotal)
SpeciesCounts$Species<-row.names(SpeciesCounts)

for (i in 1:length(rownames(SpeciesCounts))){
  SpeciesCounts$Fraction[i]<-SpeciesCounts[i,1] / TotalCountsAll
}


SpeciesGroupMeans<-aggregate(.~LocationWeek3Species$Location, mean, data=LocationWeek3Species[3:ncol(LocationWeek3Species)])
row.names(SpeciesGroupMeans)<-SpeciesGroupMeans$`LocationWeek3Species$Location`
SpeciesGroupMeans<-as.data.frame(t(SpeciesGroupMeans[,2:ncol(SpeciesGroupMeans)]))
SpeciesGroupMeans$Species<-row.names(SpeciesGroupMeans)
OverallMeans<-as.data.frame(colMeans(LocationWeek3Species[3:ncol(LocationWeek3Species)]))
OverallMeans$Species<-row.names(OverallMeans)
names(OverallMeans)<-c("OverallMean", "Species")
SpeciesGroupMeans<-merge(SpeciesGroupMeans, OverallMeans, by = "Species", all.x = TRUE)
SpeciesGroupMeans$Cincinnati<-as.numeric(SpeciesGroupMeans$Cincinnati)
SpeciesGroupMeans$Hangzhou<-as.numeric(SpeciesGroupMeans$Hangzhou)
# SpeciesGroupMeans$FoldChange<-(SpeciesGroupMeans$Hangzhou -SpeciesGroupMeans$Cincinnati) / SpeciesGroupMeans$OverallMean
# SpeciesGroupMeans$Log2<-log2(SpeciesGroupMeans$Hangzhou/SpeciesGroupMeans$OverallMean)
SpeciesGroupMeans$FC<-foldchange(SpeciesGroupMeans$Hangzhou, SpeciesGroupMeans$Cincinnati)
SpeciesGroupMeans$logratio<-foldchange2logratio(SpeciesGroupMeans$FC)

SpeciesGroupTable<-merge(SpeciesGroupMeans, SpeciesCounts, by = "Species", all.x = TRUE)
SpeciesGroupTable<-merge(SpeciesGroupTable, WilcoxTable, by = "Species", all.x = TRUE)
SpeciesGroupTable<-SpeciesGroupTable[order(SpeciesGroupTable$Unadjusted_p),]
SpeciesGroupTable<-SpeciesGroupTable[,c(1,8,4, 2,3,6,9,10)]
names(SpeciesGroupTable)<-c("Species", "Abundance", "OverallMean", "Cincinnati.Mean", "Hangzhou.Mean", "LogFC", "p_Unadjusted", "FDR")
SpeciesGroupTable$Abundance<-SpeciesGroupTable$Abundance * 100
SpeciesGroupTable$LogFC[is.na(SpeciesGroupTable$LogFC)]<-0



SpeciesGroupTable<-SpeciesGroupTable[order(SpeciesGroupTable$LogFC, decreasing = TRUE),]

SigSpeciesGroupTable<-subset(SpeciesGroupTable, SpeciesGroupTable$FDR < 0.1)

SigSpeciesGroupTable$LogFC[SigSpeciesGroupTable$LogFC == Inf]<-0
SigSpeciesGroupTable$LogFC[SigSpeciesGroupTable$LogFC == -Inf]<-0


SigSpeciesGroupTable$HigherIn<-ifelse(SigSpeciesGroupTable$LogFC > 0, "Hangzhou", "Cincinnati")
SigSpeciesGroupTable$HigherIn<-factor(SigSpeciesGroupTable$HigherIn, levels = c("Hangzhou", "Cincinnati"))

# SigSpeciesGroupTable$LogFC[SigSpeciesGroupTable$LogFC == Inf]<-10
# SigSpeciesGroupTable$LogFC[SigSpeciesGroupTable$LogFC == -Inf]<--10


FilteredSigSpeciesGroupTable<-subset(SigSpeciesGroupTable, SigSpeciesGroupTable$LogFC >= 8 | SigSpeciesGroupTable$LogFC <= -3)
LocationSpecies_LogFC<- ggplot(FilteredSigSpeciesGroupTable, 
                               aes(x = reorder(Species, LogFC), y = LogFC, fill = HigherIn)) +
  geom_bar(colour="black", stat="identity") +
  #guides(fill=FALSE) + scale_fill_hue(l=40) +
  guides(fill=guide_legend(override.aes=list(colour=NA))) +
  coord_flip() +
  labs( x = NULL) +
  ylab(expression(atop("Log Fold Change"))) +
  guides(fill=guide_legend(title="More Abundant In:")) +
  #geom_hline(yintercept = c(-0.2, 0.2), color = "firebrick", lty = 2) +
  theme(axis.text.y = element_text(size= 12, color="black")) +
  theme(legend.title = element_text(size=12)) +
  theme(legend.text = element_text(size = 11)) #+ ylim(-1.8, 1.5)
#labs( caption = "*FDR < 0.05 and SDA Effect Size > 0.5 \n All Samples")

LocationSpeciesWeek3_LogFC<-LocationSpecies_LogFC + scale_fill_manual(values=col)

ggsave(filename = "LocationSpecies_LogFC_Week3.pdf", plot = LocationSpecies_LogFC, width = 10,
       height = 14, limitsize = FALSE)



SpeciesVolcFC_Week3<-  EnhancedVolcano(as.data.frame(SigSpeciesGroupTable),
                                       lab = SigSpeciesGroupTable$Species,
                                       x = 'LogFC',
                                       y = 'FDR',
                                       pCutoff = 10e-3,
                                       FCcutoff = 2,
                                       pointSize = 2,
                                       labSize = 3.5, 
                                       title = NULL,
                                       subtitle = NULL)

SpeciesVolcFC_Week3<-SpeciesVolcFC + ylim(0,10)

ggsave(filename = "SpeciesVolcanoLogFCPlotWeek1.pdf", plot = SpeciesVolcFC, limitsize = FALSE)

save.image(file = "NICUData20250209")


library(sda)

metadata <- LocationWeek3Species[, 1:2]
cts <- as.matrix(LocationWeek3Species[, -(1:2)])
rownames(cts) <- metadata$Sample

cts_l2 <- glog2(cts)
grps <- LocationWeek3Species$Location


res <- compute_ef_Location(d = cts_l2, g = LocationWeek3Species$Location, min_shrink = 0.3)

res <- arrange(left_join(res, SigSpeciesGroupTable), desc(abs(ef_shrunk))) %>%
  mutate(Species = as_factor(Species))

res<-subset(res, ! is.na(res$FDR))
resLocation<-res

resLocation$HigherIn<-ifelse(res$ef_shrunk > 0, "Hangzhou", "Cincinnati")
resLocation$HigherIn<-factor(resLocation$HigherIn, levels = c("Hangzhou", "Cincinnati"))


write.table(resLocation, file = "SignificantLocationTableNew.csv", sep = ",", row.names = FALSE, col.names = TRUE,  quote = FALSE)
write.table(SpeciesGroupTable, file = "LocationSpeciesTable.csv", sep = ",", row.names = FALSE, col.names = TRUE,  quote = FALSE)


Location_Effect<- ggplot(filter(resLocation, abs(resLocation$ef_shrunk) > 1),
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

Location_Effect<-Location_Effect + scale_fill_manual(values = tableau.colors[c(1,2)])


ggsave(filename = "Location_Effect_Week1.pdf", plot = Location_Effect, width = 10,
       height = 12, limitsize = FALSE)


# PCA of Location


col=tableau.colors[c(2,1)]
newPCA<-PCA(cts_l2, scale.unit = FALSE, ncp = 5, graph = FALSE)
pcaPlot<-fviz_pca_ind(newPCA,
                      geom.ind = "point",
                      pointsize = 1.5,
                      point.alph = 0.1,
                      title = "Unsupervised Principal Coordinate Analysis",
                      subtitle = "Week 1 Stool Samples Grouped by Location",
                      col.ind = grps,
                      addEllipses = FALSE,
                      ellipse.alpha = 0.2,
                      ellipse.type = "confidence",
                      ellipse.level = 0.95,
                      legend.title = "Location",
                      legend.size = 11,
                      mean.point = FALSE,
                      palette = col,
                      axes.linetype = "blank"
)
LocationPCA<- pcaPlot + scale_fill_manual(values = tableau.colors[c(2,1)])

pdf("LocationPCAWeek1.pdf")
print(LocationPCA)
dev.off()
# ######

GBSLocation<-ggplot(LocationWeek1Species,
                    aes(x=Location, y=as.numeric(Streptococcus.agalactiae ), fill=Location)) + geom_boxplot(lwd=1,aes(color=factor(Location),fill = NA), outlier.size = 3)  +
  stat_summary(fun=mean, geom="point", shape=5, size=8) + scale_fill_manual(values = col) + scale_colour_manual(values = col ) +
  geom_point(size=4,aes(color = factor(Location))) + xlab(NULL) + #facet_grid(rows = . ~ Perinatal_Antimicrobial_Use) +
  ylab("Streptococcus.agalactiae \n") + paramsBoxWide() + scale_y_log10() + xlab(NULL) + scale_color_tableau() +
  theme(axis.title.y = element_text(face = "italic"))
ggsave(GBSLocation, file = "GBSLocation.pdf", width = 5, height = 8, limitsize = FALSE)
pairwise.wilcox.test(LocationWeek1Species$Streptococcus.agalactiae, LocationWeek1Species$Location)


GBSLocationWeek<-ggplot(subset(NICUSpeciesNR, NICUSpeciesNR$SampleType == "Stool"),
                        aes(x=Location, y=as.numeric(Streptococcus.agalactiae ), fill=Location)) + geom_boxplot(lwd=1,aes(color=factor(Location),fill = NA), outlier.size = 3)  +
  stat_summary(fun=mean, geom="point", shape=5, size=8) + scale_fill_manual(values = col) + scale_colour_manual(values = col ) +
  geom_point(size=4,aes(color = factor(Location))) + xlab(NULL) + #facet_grid(rows = . ~ Perinatal_Antimicrobial_Use) +
  ylab("Streptococcus.agalactiae \n") + paramsBoxWide() + scale_y_log10() + xlab(NULL) + facet_grid(rows = SampleType ~ SampleCollectionWeek) + scale_color_tableau() +
  theme(axis.title.y = element_text(face = "italic"))
ggsave(GBSLocationWeek, file = "GBSLocationWeek.pdf", width = 8, height = 8, limitsize = FALSE)
pairwise.wilcox.test(NICUSpeciesNR$Streptococcus.agalactiae, NICUSpeciesNR$Location) # p = 0.012

LmonoLocation<-ggplot(LocationWeek1Species,
                      aes(x=Location, y=as.numeric(Listeria.monocytogenes ), fill=Location)) + geom_boxplot(lwd=1,aes(color=factor(Location),fill = NA), outlier.size = 3)  +
  stat_summary(fun=mean, geom="point", shape=5, size=8) + scale_fill_manual(values = col) + scale_colour_manual(values = col ) +
  geom_point(size=4,aes(color = factor(Location))) + xlab(NULL) + #facet_grid(rows = . ~ Perinatal_Antimicrobial_Use) +
  ylab("Listeria.monocytogenes \n") + paramsBoxWide() + scale_y_log10() + xlab(NULL) + scale_color_tableau() +
  theme(axis.title.y = element_text(face = "italic"))
ggsave(LmonoLocation, file = "ListeriaLocation.pdf", width = 5, height = 8, limitsize = FALSE)


LmonoLocationWeek<-ggplot(subset(NICUSpeciesNR, NICUSpeciesNR$SampleType == "Stool"),
                          aes(x=Location, y=as.numeric(Listeria.monocytogenes ), fill=Location)) + geom_boxplot(lwd=1,aes(color=factor(Location),fill = NA), outlier.size = 3)  +
  stat_summary(fun=mean, geom="point", shape=5, size=8) + scale_fill_manual(values = col) + scale_colour_manual(values = col ) +
  geom_point(size=4,aes(color = factor(Location))) + xlab(NULL) + #facet_grid(rows = . ~ Perinatal_Antimicrobial_Use) +
  ylab("Listeria.monocytogenes \n") + paramsBoxWide() + scale_y_log10() + xlab(NULL) + facet_grid(rows = . ~ SampleCollectionWeek) + scale_color_tableau() +
  theme(axis.title.y = element_text(face = "italic"))
ggsave(LmonoLocationWeek, file = "LmonoLocationWeek.pdf", width = 8, height = 8, limitsize = FALSE)

BbreveLocationWeek<-ggplot(subset(NICUSpeciesNR, NICUSpeciesNR$SampleType == "Stool"),
                           aes(x=Location, y=as.numeric(Bifidobacterium.breve ), fill=Location)) + geom_boxplot(lwd=1,aes(color=factor(Location),fill = NA), outlier.size = 3)  +
  stat_summary(fun=mean, geom="point", shape=5, size=8) + scale_fill_manual(values = col) + scale_colour_manual(values = col ) +
  geom_point(size=4,aes(color = factor(Location))) + xlab(NULL) + #facet_grid(rows = . ~ Perinatal_Antimicrobial_Use) +
  ylab("Bifidobacterium.breve \n") + paramsBoxWide() + scale_y_log10() + xlab(NULL) + facet_grid(rows = . ~ SampleCollectionWeek) + scale_color_tableau() +
  theme(axis.title.y = element_text(face = "italic"))
ggsave(BbreveLocationWeek, file = "BbreveLocationWeek.pdf", width = 8, height = 8, limitsize = FALSE)

VatypicaLocationWeek<-ggplot(subset(NICUSpeciesNR, NICUSpeciesNR$SampleType == "Stool"),
                             aes(x=Location, y=as.numeric(Veillonella.atypica ), fill=Location)) + geom_boxplot(lwd=1,aes(color=factor(Location),fill = NA), outlier.size = 3)  +
  stat_summary(fun=mean, geom="point", shape=5, size=8) + scale_fill_manual(values = col) + scale_colour_manual(values = col ) +
  geom_point(size=4,aes(color = factor(Location))) + xlab(NULL) + #facet_grid(rows = . ~ Perinatal_Antimicrobial_Use) +
  ylab("Veillonella.atypica \n") + paramsBoxWide() + scale_y_log10() + xlab(NULL) + facet_grid(rows = . ~ SampleCollectionWeek) + scale_color_tableau() +
  theme(axis.title.y = element_text(face = "italic"))
ggsave(VatypicaLocationWeek, file = "VatypicaLocationWeek.pdf", width = 8, height = 8, limitsize = FALSE)


AbaumanniiLocationWeek<-ggplot(subset(NICUSpeciesNR, NICUSpeciesNR$SampleType == "Stool"),
                               aes(x=Location, y=as.numeric(Acinetobacter.baumannii ), fill=Location)) + geom_boxplot(lwd=1,aes(color=factor(Location),fill = NA), outlier.size = 3)  +
  stat_summary(fun=mean, geom="point", shape=5, size=8) + scale_fill_manual(values = col) + scale_colour_manual(values = col ) +
  geom_point(size=4,aes(color = factor(Location))) + xlab(NULL) + #facet_grid(rows = . ~ Perinatal_Antimicrobial_Use) +
  ylab("Acinetobacter.baumannii \n") + paramsBoxWide() + scale_y_log10() + xlab(NULL) + facet_grid(rows = . ~ SampleCollectionWeek) + scale_color_tableau() +
  theme(axis.title.y = element_text(face = "italic"))
ggsave(AbaumanniiLocationWeek, file = "AbaumanniiLocationWeek.pdf", width = 8, height = 8, limitsize = FALSE)

# save.image(file = "NICUData20220516")


# Do PCA for US samples: SampleType and Week

## Next do axilla, groin, and stool. (US only) to see what changes with antibiotics or from week 1 to week 3

AxillaGroinStoolSpecies<-subset(NICUSpeciesNR, NICUSpeciesNR$SampleType %in% c("Axilla", "Groin", "Stool"))
AxillaGroinStoolSpecies$TypeWeekLocation<-paste(AxillaGroinStoolSpecies$TypeWeek, AxillaGroinStoolSpecies$Location, sep = "-")
AxillaGroinStoolSpecies<-AxillaGroinStoolSpecies[, c(1:17, ncol(AxillaGroinStoolSpecies), 18:ncol(AxillaGroinStoolSpecies)-1)]
AxillaGroinStoolSpecies$AnyMilk.1<-NULL

# AxillaGroinStoolSpecies<-subset(NICUSpeciesNR, NICUSpeciesNR$PostNatalAntibioticsNew %in% c("No.Infant.Abx", "Low.Infant.Abx", "High.Infant.Abx"))

library(factoextra)

AbxSampleInd<-AxillaGroinStoolSpecies$Sample[AxillaGroinStoolSpecies$PostNatalAntibiotics %in% c("Infant.Abx")]
NoAbxSampleInd<-AxillaGroinStoolSpecies$Sample[AxillaGroinStoolSpecies$PostNatalAntibiotics == "No.Infant.Abx"]


metadata <- AxillaGroinStoolSpecies[, 1:18]
cts <- as.matrix(AxillaGroinStoolSpecies[, -(1:18)])
rownames(cts) <- metadata$Sample
cts_l2 <- glog2(cts)
grps <- AxillaGroinStoolSpecies$TypeWeekLocation

col=rep(tableau.colors,2)
newPCA<-PCA(cts_l2, scale.unit = FALSE, ncp = 5, graph = FALSE)
pcaPlot<-fviz_pca_ind(newPCA,
                      geom.ind = "point", 
                      pointsize = 1.5,
                      point.alph = 0.1,
                      title = "Unsupervised Principal Coordinate Analysis",
                      subtitle = "Samples colored by Sample Type and Infant Age",
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
SampleTypeWeekPCA<- pcaPlot #+ xlim(-55,95) + ylim(-45,60)

pdf("SampleTypeWeekPCA_AllSamples.pdf") 
print(SampleTypeWeekPCA)
dev.off() 

col=rep(tableau.colors,2)
pcaPlot<-fviz_pca_ind(newPCA,
                      geom.ind = "point", 
                      pointsize = 1.5,
                      point.alph = 0.1,
                      title = "Unsupervised Principal Coordinate Analysis",
                      subtitle = "Samples Colored by Sample Type and Infant Age ; No Antibiotics",
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
                      select.ind = list(name = NoAbxSampleInd),
)


SampleTypeWeekPCA<- pcaPlot #+ xlim(-55,95) + ylim(-45,60)

pdf("SampleTypeWeekPCA_NOABX.pdf") 
print(SampleTypeWeekPCA)
dev.off() 

# ABX samples
pcaPlot<-fviz_pca_ind(newPCA,
                      geom.ind = "point", 
                      pointsize = 1.5,
                      point.alph = 0.1,
                      title = "Unsupervised Principal Coordinate Analysis",
                      subtitle = "Samples colored by Sample Type and Infant Age : Any Antibiotics",
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
                      select.ind = list(name = AbxSampleInd)
)
SampleTypeWeekPCA<- pcaPlot #+ xlim(-55,95) + ylim(-45,60)

pdf("SampleTypeWeekPCA_ABX.pdf") 
print(SampleTypeWeekPCA)
dev.off() 


## Next do axilla, groin, and stool.  to see what changes with antibiotics or from week 1 to week 3

library(FactoMineR)

Week1AxillaGroinStoolSpecies<-subset(NICUSpeciesNR, NICUSpeciesNR$SampleType %in% c("Axilla", "Groin", "Stool") &
                                       NICUSpeciesNR$SampleCollectionWeek == "Week.1")

Week1AxillaGroinStoolSpecies$TypeAbxLocation<-paste(Week1AxillaGroinStoolSpecies$TypeAbx, Week1AxillaGroinStoolSpecies$Location, sep = "-")
Week1AxillaGroinStoolSpecies<-Week1AxillaGroinStoolSpecies[, c(1:17, ncol(Week1AxillaGroinStoolSpecies), 18:ncol(Week1AxillaGroinStoolSpecies)-1)]
Week1AxillaGroinStoolSpecies$AnyMilk.1<-NULL

AbxSampleInd<-Week1AxillaGroinStoolSpecies$Sample[Week1AxillaGroinStoolSpecies$PostNatalAntibiotics == "Infant.Abx"]
NoAbxSampleInd<-Week1AxillaGroinStoolSpecies$Sample[Week1AxillaGroinStoolSpecies$PostNatalAntibiotics == "No.Infant.Abx"]

metadata <- Week1AxillaGroinStoolSpecies[, 1:18]
cts <- as.matrix(Week1AxillaGroinStoolSpecies[, -(1:18)])
rownames(cts) <- metadata$Sample
cts_l2 <- glog2(cts)
grps <- Week1AxillaGroinStoolSpecies$TypeAbxLocation

col=rep(tableau.colors,2)
newPCA<-PCA(cts_l2, scale.unit = FALSE, ncp = 5, graph = FALSE)
pcaPlot<-fviz_pca_ind(newPCA,
                      geom.ind = "point", 
                      pointsize = 1.5,
                      point.alph = 0.1,
                      title = "Unsupervised Principal Coordinate Analysis",
                      subtitle = "Week 1 Samples colored by Sample Type and Abx Exposure",
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
                      
)
SampleTypeAbxWk1PCA<- pcaPlot #+ xlim(-55,80) + ylim(-45,95)  

pdf("SampleTypeAbxWk1PCA_AllSamples.pdf") 
print(SampleTypeAbxWk1PCA)
dev.off() 


# This should be same figure as above but no antibiotics 

## ok there's an error because a sample with antibiotics shows up.
col<-tableau.colors
pcaPlot<-fviz_pca_ind(newPCA,
                      geom.ind = "point", 
                      pointsize = 1.5,
                      point.alph = 0.1,
                      title = "Unsupervised Principal Coordinate Analysis",
                      subtitle = "Week 1 Samples colored by Sample Type : No Antibiotic Exposure",
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
                      select.ind = list(name = NoAbxSampleInd)
)
SampleTypeAbxWk1PCA<- pcaPlot #+ xlim(-55,80) + ylim(-45,95) 

pdf("SampleTypeWk1PCA_NOAbx.pdf") 
print(SampleTypeAbxWk1PCA)
dev.off() 

col<-tableau.colors
pcaPlot<-fviz_pca_ind(newPCA,
                      geom.ind = "point", 
                      pointsize = 1.5,
                      point.alph = 0.1,
                      title = "Unsupervised Principal Coordinate Analysis",
                      subtitle = "Week 1 Samples colored by Sample Type : Antibiotic Exposed Infants",
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
                      select.ind = list(name = AbxSampleInd)
)
SampleTypeYesAbxWk1PCA<- pcaPlot + xlim(-55,80) + ylim(-45,95) 

pdf("SampleTypeWk1PCA_AbxExposed.pdf") 
print(SampleTypeYesAbxWk1PCA)
dev.off() 


#HERE 4/19/22 8:29AM

## Week 3 PCA ##
Week3AxillaGroinStoolSpecies<-subset(NICUSpeciesNR, NICUSpeciesNR$SampleType %in% c("Axilla", "Groin", "Stool") &
                                       NICUSpeciesNR$SampleCollectionWeek == "Week.3")

Week3AxillaGroinStoolSpecies$TypeAbxLocation<-paste(Week3AxillaGroinStoolSpecies$TypeAbx, Week3AxillaGroinStoolSpecies$Location, sep = "-")
Week3AxillaGroinStoolSpecies<-Week3AxillaGroinStoolSpecies[, c(1:17, ncol(Week3AxillaGroinStoolSpecies), 18:ncol(Week3AxillaGroinStoolSpecies)-1)]
Week3AxillaGroinStoolSpecies$AnyMilk.1<-NULL

AbxSampleInd<-Week3AxillaGroinStoolSpecies$Sample[Week3AxillaGroinStoolSpecies$PostNatalAntibiotics == "Infant.Abx"]
NoAbxSampleInd<-Week3AxillaGroinStoolSpecies$Sample[Week3AxillaGroinStoolSpecies$PostNatalAntibiotics == "No.Infant.Abx"]


metadata <- Week3AxillaGroinStoolSpecies[, 1:18]
cts <- as.matrix(Week3AxillaGroinStoolSpecies[, -(1:18)])
rownames(cts) <- metadata$Sample
cts_l2 <- glog2(cts)
grps <- Week3AxillaGroinStoolSpecies$TypeAbxLocation

col=rep(tableau.colors,2)
newPCA<-PCA(cts_l2, scale.unit = FALSE, ncp = 5, graph = FALSE)
pcaPlot<-fviz_pca_ind(newPCA,
                      geom.ind = "point", 
                      pointsize = 1.5,
                      point.alph = 0.1,
                      title = "Unsupervised Principal Coordinate Analysis",
                      subtitle = "Week 3 Samples colored by Sample Type and Abx Exposure",
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
SampleTypeAbxWk3PCA<- pcaPlot  + xlim(-50,80) + ylim(-75,55)

pdf("SampleTypeAbxWk3PCA_AllSamples.pdf") 
print(SampleTypeAbxWk3PCA)
dev.off() 

col=tableau.colors
newPCA<-PCA(cts_l2, scale.unit = FALSE, ncp = 5, graph = FALSE)
pcaPlot<-fviz_pca_ind(newPCA,
                      geom.ind = "point", 
                      pointsize = 1.5,
                      point.alph = 0.1,
                      title = "Unsupervised Principal Coordinate Analysis",
                      subtitle = "Week 3 Samples colored by Sample Type : Antibiotic Exposed Infants",
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
                      select.ind = list(name = AbxSampleInd)
)
SampleTypeAbxWk3PCA<- pcaPlot + xlim(-50,80) + ylim(-75,55)

pdf("SampleTypeAbxWk3PCA_AbxExposed.pdf") 
print(SampleTypeAbxWk3PCA)
dev.off() 

col=tableau.colors
newPCA<-PCA(cts_l2, scale.unit = FALSE, ncp = 5, graph = FALSE)
pcaPlot<-fviz_pca_ind(newPCA,
                      geom.ind = "point", 
                      pointsize = 1.5,
                      point.alph = 0.1,
                      title = "Unsupervised Principal Coordinate Analysis",
                      subtitle = "Week 3 Samples colored by Sample Type : No Antibiotic Exposure",
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
                      select.ind = list(name = NoAbxSampleInd)
)
SampleTypeAbxWk3PCA<- pcaPlot + xlim(-50,80) + ylim(-75,55)


pdf("SampleTypeAbxWk3PCA_NoAbxExposure.pdf") 
print(SampleTypeAbxWk3PCA)
dev.off() 



# OK look at week 1 stools to see what's different:
# This is both Hangzhou and Cincinnati

Week1Stools<-subset(NICUSpeciesNR, NICUSpeciesNR$SampleType == "Stool" & NICUSpeciesNR$SampleCollectionWeek =="Week.1")
Week1Stools<-Week1Stools[,c(1,10,18:ncol(Week1Stools))]
Week1Stools<-na.omit(Week1Stools)


# p = 0.001 if No.Abx versus Abx
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
WilcoxTable$Species<-row.names(WilcoxTable)
WilcoxTable<-WilcoxTable[order(WilcoxTable$Unadjusted_p),]



# Now let's figure out relative contribution of each species, fold difference between groups, and direction of difference


SpeciesCounts<-as.data.frame(colSums(Week1Stools[3:length(colnames(Week1Stools))]))
names(SpeciesCounts)<-"SpeciesTotal"
SpeciesCountsBigTable<-as.data.frame(colSums(Week1Stools[3:length(colnames(Week1Stools))]))
names(SpeciesCountsBigTable)<-"SpeciesTotal"
TotalCountsAll<-sum(SpeciesCountsBigTable$SpeciesTotal)
SpeciesCounts$Species<-row.names(SpeciesCounts)

for (i in 1:length(rownames(SpeciesCounts))){
  SpeciesCounts$Fraction[i]<-SpeciesCounts[i,1] / TotalCountsAll
}


SpeciesGroupMeans<-aggregate(.~Week1Stools$PostNatalAntibiotics, mean, data=Week1Stools[3:ncol(Week1Stools)])
row.names(SpeciesGroupMeans)<-SpeciesGroupMeans$`Week1Stools$PostNatalAntibiotics`
SpeciesGroupMeans<-as.data.frame(t(SpeciesGroupMeans[,2:ncol(SpeciesGroupMeans)]))
SpeciesGroupMeans$Species<-row.names(SpeciesGroupMeans)
OverallMeans<-as.data.frame(colMeans(Week1Stools[3:ncol(Week1Stools)]))
OverallMeans$Species<-row.names(OverallMeans)
names(OverallMeans)<-c("OverallMean", "Species")
SpeciesGroupMeans<-merge(SpeciesGroupMeans, OverallMeans, by = "Species", all.x = TRUE)
SpeciesGroupMeans$No.Infant.Abx<-as.numeric(SpeciesGroupMeans$No.Infant.Abx)
SpeciesGroupMeans$Infant.Abx<-as.numeric(SpeciesGroupMeans$Infant.Abx)

# SpeciesGroupMeans$FoldChange<-(SpeciesGroupMeans$Infant.Abx -SpeciesGroupMeans$No.Infant.Abx) / SpeciesGroupMeans$OverallMean
# SpeciesGroupMeans$Log2<-log2(SpeciesGroupMeans$Infant.Abx/SpeciesGroupMeans$OverallMean)

SpeciesGroupMeans$FC<-foldchange(SpeciesGroupMeans$Infant.Abx, SpeciesGroupMeans$No.Infant.Abx)
SpeciesGroupMeans$logratio<-foldchange2logratio(SpeciesGroupMeans$FC)


SpeciesGroupTable<-merge(SpeciesGroupMeans, SpeciesCounts, by = "Species", all.x = TRUE)
SpeciesGroupTable<-merge(SpeciesGroupTable, WilcoxTable, by = "Species", all.x = TRUE)
SpeciesGroupTable<-SpeciesGroupTable[order(SpeciesGroupTable$Unadjusted_p),]
SpeciesGroupTable<-SpeciesGroupTable[,c(1,8,4, 2,3,5,9,10)]
names(SpeciesGroupTable)<-c("Species", "Abundance", "OverallMean", "No.Infant.Abx.Mean", "Infant.Abx.Mean", "Fold.Change", "p_Unadjusted", "FDR")
SpeciesGroupTable$Abundance<-SpeciesGroupTable$Abundance * 100
SpeciesGroupTable$Fold.Change[is.na(SpeciesGroupTable$Fold.Change)]<-0
SpeciesGroupTable<-SpeciesGroupTable[order(SpeciesGroupTable$Fold.Change, decreasing = TRUE),]

SigSpeciesGroupTable<-subset(SpeciesGroupTable, SpeciesGroupTable$p_Unadjusted < 0.2)

# Ok try out new effect size calculation from Klaus
# there's probably a way to shorten this part but we'll use his approach and first split out class and counts:

metadata <- Week1Stools[, 1:2]
cts <- as.matrix(Week1Stools[, -(1:2)])
rownames(cts) <- metadata$Sample

cts_l2 <- glog2(cts)
grps <- Week1Stools$PostNatalAntibiotics


res <- compute_ef_Abx(d = cts_l2, g = Week1Stools$PostNatalAntibiotics, min_shrink = 0.3)

res <- arrange(left_join(res, SigSpeciesGroupTable), desc(abs(ef_shrunk))) %>%
  mutate(Species = as_factor(Species))

res<-subset(res, ! is.na(res$FDR)) 
resWk1Abx<-res

resWk1Abx$HigherIn<-ifelse(res$ef_shrunk > 0, "Infant.Abx", "No.Infant.Abx")
resWk1Abx$HigherIn<-factor(resWk1Abx$HigherIn, levels = c("Infant.Abx", "No.Infant.Abx"))


write.table(resWk1Abx, file = "SignificantWk1AbxTableNew.csv", sep = ",", row.names = FALSE, col.names = TRUE,  quote = FALSE)
write.table(SpeciesGroupTable, file = "Wk1SpeciesTable.csv", sep = ",", row.names = FALSE, col.names = TRUE,  quote = FALSE)


Week1Abx_Effect<- ggplot(filter(resWk1Abx, abs(resWk1Abx$ef_shrunk) > 0.5),
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

Week1Abx_Effect<-Week1Abx_Effect + scale_fill_manual(values = tableau.colors[c(1,2)])
ggsave(filename = "Week1Abx_Effect.pdf", plot = Week1Abx_Effect, width = 10, 
       height = 8, limitsize = FALSE)




Week1Abx_LogFC<- ggplot(filter(SigSpeciesGroupTable, SigSpeciesGroupTable$LogFC >= 7 | SigSpeciesGroupTable$LogFC <= -7), 
                        aes(x = reorder(Species, LogFC), y = LogFC, fill = HigherIn)) +
  geom_bar(colour="black", stat="identity") +
  #guides(fill=FALSE) + scale_fill_hue(l=40) +
  guides(fill=guide_legend(override.aes=list(colour=NA))) +
  coord_flip() +
  labs( x = NULL) +
  ylab(expression(atop("Log Fold Change"))) +
  guides(fill=guide_legend(title="More Abundant In:")) +
  #geom_hline(yintercept = c(-0.2, 0.2), color = "firebrick", lty = 2) +
  theme(axis.text.y = element_text(size= 12, color="black")) +
  theme(legend.title = element_text(size=12)) +
  theme(legend.text = element_text(size = 11)) #+ ylim(-1.8, 1.5)
#labs( caption = "*FDR < 0.05 and SDA Effect Size > 0.5 \n All Samples")

LocationSpecies_LogFC<-LocationSpecies_LogFC + scale_fill_manual(values=col)

ggsave(filename = "LocationSpecies_LogFC.pdf", plot = LocationSpecies_LogFC, width = 10,
       height = 14, limitsize = FALSE)



SpeciesVolcFC<-  EnhancedVolcano(as.data.frame(SigSpeciesGroupTable),
                                 lab = SigSpeciesGroupTable$Species,
                                 x = 'LogFC',
                                 y = 'FDR',
                                 pCutoff = 10e-3,
                                 FCcutoff = 2,
                                 pointSize = 2,
                                 labSize = 3.5, 
                                 title = NULL,
                                 subtitle = NULL)

SpeciesVolcFC<-SpeciesVolcFC + ylim(0,10)

ggsave(filename = "SpeciesVolcanoLogFCPlot.pdf", plot = SpeciesVolcFC, limitsize = FALSE)






# OK look at week 3 stools to see what's different:
# this is both Cincinnati and Hangzhou. It will be important to separate those out because the antibiotic use is much different

Week3Stools<-subset(NICUSpeciesNR, NICUSpeciesNR$SampleType == "Stool" & NICUSpeciesNR$SampleCollectionWeek =="Week.3")
# Week3Stools<-Week3Stools[,c(1,9,18:ncol(Week3Stools))] ## 9 is incorrect here. That's SampleCollectionWeek
Week3Stools<-Week3Stools[,c(1,10,18:ncol(Week3Stools))] # also, there must have been an NA in the table, so did na.omit
Week3Stools<-na.omit(Week3Stools)


# Pick up here 20220517

# p = 0.04 if No.Abx versus Abx
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
WilcoxTable$Species<-row.names(WilcoxTable)
WilcoxTable<-WilcoxTable[order(WilcoxTable$Unadjusted_p),]



# Now let's figure out relative contribution of each species, fold difference between groups, and direction of difference

SpeciesCounts<-as.data.frame(colSums(Week3Stools[3:length(colnames(Week3Stools))]))
names(SpeciesCounts)<-"SpeciesTotal"
SpeciesCountsBigTable<-as.data.frame(colSums(Week3Stools[3:length(colnames(Week3Stools))]))
names(SpeciesCountsBigTable)<-"SpeciesTotal"
TotalCountsAll<-sum(SpeciesCountsBigTable$SpeciesTotal)
SpeciesCounts$Species<-row.names(SpeciesCounts)

for (i in 1:length(rownames(SpeciesCounts))){
  SpeciesCounts$Fraction[i]<-SpeciesCounts[i,1] / TotalCountsAll
}



SpeciesGroupMeans<-aggregate(.~Week3Stools$PostNatalAntibiotics, mean, data=Week3Stools[3:ncol(Week3Stools)])
row.names(SpeciesGroupMeans)<-SpeciesGroupMeans$`Week3Stools$PostNatalAntibiotics`
SpeciesGroupMeans<-as.data.frame(t(SpeciesGroupMeans[,2:ncol(SpeciesGroupMeans)]))
SpeciesGroupMeans$Species<-row.names(SpeciesGroupMeans)
OverallMeans<-as.data.frame(colMeans(Week3Stools[3:ncol(Week3Stools)]))
OverallMeans$Species<-row.names(OverallMeans)
names(OverallMeans)<-c("OverallMean", "Species")
SpeciesGroupMeans<-merge(SpeciesGroupMeans, OverallMeans, by = "Species", all.x = TRUE)
SpeciesGroupMeans$No.Infant.Abx<-as.numeric(SpeciesGroupMeans$No.Infant.Abx)
SpeciesGroupMeans$Infant.Abx<-as.numeric(SpeciesGroupMeans$Infant.Abx)
SpeciesGroupMeans$FoldChange<-(SpeciesGroupMeans$Infant.Abx -SpeciesGroupMeans$No.Infant.Abx) / SpeciesGroupMeans$OverallMean
SpeciesGroupMeans$Log2<-log2(SpeciesGroupMeans$Infant.Abx/SpeciesGroupMeans$OverallMean)
SpeciesGroupTable<-merge(SpeciesGroupMeans, SpeciesCounts, by = "Species", all.x = TRUE)
SpeciesGroupTable<-merge(SpeciesGroupTable, WilcoxTable, by = "Species", all.x = TRUE)
SpeciesGroupTable<-SpeciesGroupTable[order(SpeciesGroupTable$Unadjusted_p),]
SpeciesGroupTable<-SpeciesGroupTable[,c(1,8,4, 2,3,5,9,10)]
names(SpeciesGroupTable)<-c("Species", "Abundance", "OverallMean", "No.Infant.Abx.Mean", "Infant.Abx.Mean", "Fold.Change", "p_Unadjusted", "FDR")
SpeciesGroupTable$Abundance<-SpeciesGroupTable$Abundance * 100
SpeciesGroupTable$Fold.Change[is.na(SpeciesGroupTable$Fold.Change)]<-0
SpeciesGroupTable<-SpeciesGroupTable[order(SpeciesGroupTable$Fold.Change, decreasing = TRUE),]

SigSpeciesGroupTable<-subset(SpeciesGroupTable, SpeciesGroupTable$p_Unadjusted < 0.05)

# Ok try out new effect size calculation from Klaus
# there's probably a way to shorten this part but we'll use his approach and first split out class and counts:

metadata <- Week3Stools[, 1:2]
cts <- as.matrix(Week3Stools[, -(1:2)])
rownames(cts) <- metadata$Sample

cts_l2 <- glog2(cts)
grps <- Week3Stools$PostNatalAntibiotics


res <- compute_ef_Abx(d = cts_l2, g = Week3Stools$PostNatalAntibiotics, min_shrink = 0.3)

res <- arrange(left_join(res, SigSpeciesGroupTable), desc(abs(ef_shrunk))) %>%
  mutate(Species = as_factor(Species))

res<-subset(res, ! is.na(res$FDR)) 
resWk3Abx<-res

resWk3Abx$HigherIn<-ifelse(res$ef_shrunk > 0, "Infant.Abx", "No.Infant.Abx")
resWk3Abx$HigherIn<-factor(resWk3Abx$HigherIn, levels = c("Infant.Abx", "No.Infant.Abx"))


write.table(resWk3Abx, file = "SignificantWk3AbxTableNew.csv", sep = ",", row.names = FALSE, col.names = TRUE,  quote = FALSE)
write.table(SpeciesGroupTable, file = "Wk3SpeciesTable.csv", sep = ",", row.names = FALSE, col.names = TRUE,  quote = FALSE)


Week3Abx_Effect<- ggplot(filter(resWk3Abx, resWk3Abx$ef_shrunk > 0.2 | resWk3Abx$ef_shrunk < -0.5),
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

Week3Abx_Effect<-Week3Abx_Effect + scale_fill_manual(values = tableau.colors[c(1,2)])
ggsave(filename = "Week3Abx_Effect.pdf", plot = Week3Abx_Effect, width = 10, 
       height = 8, limitsize = FALSE)


# Make a dummy table so can set 0 counts to 1

DummySpeciesNR<-NICUSpeciesNR

# Try the new GestationalAge numeric column with a species

ggplot(NICUSpeciesNR, aes(x=GestationTime, y=log(Bifidobacterium.longum))) + 
  geom_point()+
  geom_smooth(method=lm)


ggplot(NICUSpeciesNR, aes(x=GestationTime, y=log(Escherichia.coli
))) + 
  geom_point()+
  geom_smooth(method=lm)

DummySpeciesNR$Bifidobacterium.longum[DummySpeciesNR$Bifidobacterium.longum==0]<-1
BifidoAbx<-ggplot(subset(DummySpeciesNR, DummySpeciesNR$SampleType %in% c("Axilla", "Groin", "Stool")),
                  aes(x=PostNatalAntibiotics, y=as.numeric(Bifidobacterium.longum ), fill=PostNatalAntibiotics)) + geom_boxplot(lwd=1,aes(color=factor(PostNatalAntibiotics),fill = NA), outlier.size = 3)  +
  stat_summary(fun=mean, geom="point", shape=5, size=8) + scale_fill_manual(values = col) + scale_colour_manual(values = col ) +
  geom_point(size=4,aes(color = factor(PostNatalAntibiotics))) + xlab(NULL) + facet_grid(rows = SampleType ~ SampleCollectionWeek + Location) +
  ylab("Bifidobacterium.longum \n") + paramsBoxWide() + scale_y_log10() + xlab(NULL)  + scale_color_tableau() +
  theme(axis.title.y = element_text(face = "italic"))
ggsave(BifidoAbx, file = "BifidoAbx.pdf", width = 8, height = 8, limitsize = FALSE)

DummySpeciesNR$Streptococcus.sp.HMSC065E03[DummySpeciesNR$Streptococcus.sp.HMSC065E03==0]<-1
StrepHMSCAbx<-ggplot(subset(DummySpeciesNR, DummySpeciesNR$SampleType %in% c("Axilla", "Groin", "Stool")),
                     aes(x=PostNatalAntibiotics, y=as.numeric(Streptococcus.sp.HMSC065E03 ), fill=PostNatalAntibiotics)) + geom_boxplot(lwd=1,aes(color=factor(PostNatalAntibiotics),fill = NA), outlier.size = 3)  +
  stat_summary(fun=mean, geom="point", shape=5, size=8) + scale_fill_manual(values = col) + scale_colour_manual(values = col ) +
  geom_point(size=4,aes(color = factor(PostNatalAntibiotics))) + xlab(NULL) + facet_grid(rows = SampleType ~ SampleCollectionWeek) +
  ylab("Streptococcus.sp.HMSC065E03 \n") + paramsBoxWide() + scale_y_log10() + xlab(NULL)  + scale_color_tableau() +
  theme(axis.title.y = element_text(face = "italic"))
ggsave(StrepHMSCAbx, file = "StrepHMSCAbx.pdf", width = 8, height = 8, limitsize = FALSE)

DummySpeciesNR$Streptococcus.infantis[DummySpeciesNR$Streptococcus.infantis==0]<-1
StrepInfAbx<-ggplot(subset(DummySpeciesNR, DummySpeciesNR$SampleType %in% c("Axilla", "Groin", "Stool")),
                    aes(x=PostNatalAntibiotics, y=as.numeric(Streptococcus.infantis ), fill=PostNatalAntibiotics)) + geom_boxplot(lwd=1,aes(color=factor(PostNatalAntibiotics),fill = NA), outlier.size = 3)  +
  stat_summary(fun=mean, geom="point", shape=5, size=8) + scale_fill_manual(values = col) + scale_colour_manual(values = col ) +
  geom_point(size=4,aes(color = factor(PostNatalAntibiotics))) + xlab(NULL) + facet_grid(rows = SampleType ~ SampleCollectionWeek) +
  ylab("Streptococcus.infantis \n") + paramsBoxWide() + scale_y_log10() + xlab(NULL)  + scale_color_tableau() +
  theme(axis.title.y = element_text(face = "italic"))
ggsave(StrepInfAbx, file = "StrepInfantisAbx.pdf", width = 8, height = 8, limitsize = FALSE)

DummySpeciesNR$Staphylococcus.aureus[DummySpeciesNR$Staphylococcus.aureus==0]<-1
StaphAbx<-ggplot(subset(DummySpeciesNR, DummySpeciesNR$SampleType %in% c("Axilla", "Groin", "Stool")),
                 aes(x=PostNatalAntibiotics, y=as.numeric(Staphylococcus.aureus ), fill=PostNatalAntibiotics)) + geom_boxplot(lwd=1,aes(color=factor(PostNatalAntibiotics),fill = NA), outlier.size = 3)  +
  stat_summary(fun=mean, geom="point", shape=5, size=8) + scale_fill_manual(values = col) + scale_colour_manual(values = col ) +
  geom_point(size=4,aes(color = factor(PostNatalAntibiotics))) + xlab(NULL) + #facet_grid(rows = . ~ Perinatal_Antimicrobial_Use) +
  ylab("Staphylococcus.aureus \n") + paramsBoxWide() + scale_y_log10() + xlab(NULL) + facet_grid(rows = SampleType ~ SampleCollectionWeek) + scale_color_tableau() +
  theme(axis.title.y = element_text(face = "italic"))
ggsave(StaphAbx, file = "StaphAbx.pdf", width = 12, height = 10, limitsize = FALSE)

DummySpeciesNR$Klebsiella.pneumoniae[DummySpeciesNR$Klebsiella.pneumoniae==0]<-1
KpneumoAbx<-ggplot(subset(DummySpeciesNR, DummySpeciesNR$SampleType %in% c("Axilla", "Groin", "Stool")),
                   aes(x=PostNatalAntibiotics, y=as.numeric(Klebsiella.pneumoniae ), fill=PostNatalAntibiotics)) + geom_boxplot(lwd=1,aes(color=factor(PostNatalAntibiotics),fill = NA), outlier.size = 3)  +
  stat_summary(fun=mean, geom="point", shape=5, size=8) + scale_fill_manual(values = col) + scale_colour_manual(values = col ) +
  geom_point(size=4,aes(color = factor(PostNatalAntibiotics))) + xlab(NULL) + #facet_grid(rows = . ~ Perinatal_Antimicrobial_Use) +
  ylab("Klebsiella.pneumoniae \n") + paramsBoxWide() + scale_y_log10() + xlab(NULL) + facet_grid(rows = SampleType ~ SampleCollectionWeek) + scale_color_tableau() +
  theme(axis.title.y = element_text(face = "italic"))
ggsave(KpneumoAbx, file = "KpneumoAbx.pdf", width = 12, height = 10, limitsize = FALSE)

DummySpeciesNR$Staphylococcus.epidermidis[DummySpeciesNR$Staphylococcus.epidermidis==0]<-1
StaphEpiAbx<-ggplot(subset(DummySpeciesNR, DummySpeciesNR$SampleType %in% c("Axilla", "Groin", "Stool")),
                    aes(x=PostNatalAntibiotics, y=as.numeric(Staphylococcus.epidermidis ), fill=PostNatalAntibiotics)) + geom_boxplot(lwd=1,aes(color=factor(PostNatalAntibiotics),fill = NA), outlier.size = 3)  +
  stat_summary(fun=mean, geom="point", shape=5, size=8) + scale_fill_manual(values = col) + scale_colour_manual(values = col ) +
  geom_point(size=4,aes(color = factor(PostNatalAntibiotics))) + xlab(NULL) + #facet_grid(rows = . ~ Perinatal_Antimicrobial_Use) +
  ylab("Staphylococcus.epidermidis \n") + paramsBoxWide() + scale_y_log10() + xlab(NULL) + facet_grid(rows = SampleType ~ SampleCollectionWeek) + scale_color_tableau() +
  theme(axis.title.y = element_text(face = "italic"))
ggsave(StaphEpiAbx, file = "StaphEpiAbx.pdf", width = 12, height = 10, limitsize = FALSE)

DummySpeciesNR$Enterococcus.faecalis[DummySpeciesNR$Enterococcus.faecalis==0]<-1
EfaecalisAbx<-ggplot(subset(DummySpeciesNR, DummySpeciesNR$SampleType %in% c("Axilla", "Groin", "Stool")),
                     aes(x=PostNatalAntibiotics, y=as.numeric(Enterococcus.faecalis ), fill=PostNatalAntibiotics)) + geom_boxplot(lwd=1,aes(color=factor(PostNatalAntibiotics),fill = NA), outlier.size = 3)  +
  stat_summary(fun=mean, geom="point", shape=5, size=8) + scale_fill_manual(values = col) + scale_colour_manual(values = col ) +
  geom_point(size=4,aes(color = factor(PostNatalAntibiotics))) + xlab(NULL) + #facet_grid(rows = . ~ Perinatal_Antimicrobial_Use) +
  ylab("Enterococcus.faecalis \n") + paramsBoxWide() + scale_y_log10() + xlab(NULL) + facet_grid(rows = SampleType ~ SampleCollectionWeek) + scale_color_tableau() +
  theme(axis.title.y = element_text(face = "italic"))
ggsave(EfaecalisAbx, file = "EfaecalisAbx.pdf", width = 12, height = 10, limitsize = FALSE)

DummySpeciesNR$Enterococcus.faecium[DummySpeciesNR$Enterococcus.faecium==0]<-1
EfaeciumAbx<-ggplot(subset(DummySpeciesNR, DummySpeciesNR$SampleType %in% c("Axilla", "Groin", "Stool")),
                    aes(x=PostNatalAntibiotics, y=as.numeric(Enterococcus.faecium), fill=PostNatalAntibiotics)) + geom_boxplot(lwd=1,aes(color=factor(PostNatalAntibiotics),fill = NA), outlier.size = 3)  +
  stat_summary(fun=mean, geom="point", shape=5, size=8) + scale_fill_manual(values = col) + scale_colour_manual(values = col ) +
  geom_point(size=4,aes(color = factor(PostNatalAntibiotics))) + xlab(NULL) + #facet_grid(rows = . ~ Perinatal_Antimicrobial_Use) +
  ylab("Enterococcus.faecium\n") + paramsBoxWide() + scale_y_log10() + xlab(NULL) + facet_grid(rows = SampleType ~ SampleCollectionWeek) + scale_color_tableau() +
  theme(axis.title.y = element_text(face = "italic"))
ggsave(EfaeciumAbx, file = "EfaeciumAbx.pdf", width = 12, height = 10, limitsize = FALSE)

DummySpeciesNR$Bifidobacterium.breve[DummySpeciesNR$Bifidobacterium.breve==0]<-1
BBreveAbx<-ggplot(subset(DummySpeciesNR, DummySpeciesNR$SampleType %in% c("Axilla", "Groin", "Stool")),
                  aes(x=PostNatalAntibiotics, y=as.numeric(Bifidobacterium.breve), fill=PostNatalAntibiotics)) + geom_boxplot(lwd=1,aes(color=factor(PostNatalAntibiotics),fill = NA), outlier.size = 3)  +
  stat_summary(fun=mean, geom="point", shape=5, size=8) + scale_fill_manual(values = col) + scale_colour_manual(values = col ) +
  geom_point(size=4,aes(color = factor(PostNatalAntibiotics))) + xlab(NULL) + #facet_grid(rows = . ~ Perinatal_Antimicrobial_Use) +
  ylab("Bifidobacterium.breve\n") + paramsBoxWide() + scale_y_log10() + xlab(NULL) + facet_grid(rows = SampleType ~ SampleCollectionWeek) + scale_color_tableau() +
  theme(axis.title.y = element_text(face = "italic"))
ggsave(BBreveAbx, file = "BBreveAbx.pdf", width = 12, height = 10, limitsize = FALSE)

DummySpeciesNR$Escherichia.coli[DummySpeciesNR$Escherichia.coli==0]<-1
EcoliAbx<-ggplot(subset(DummySpeciesNR, DummySpeciesNR$SampleType %in% c("Axilla", "Groin", "Stool")),
                 aes(x=PostNatalAntibiotics, y=as.numeric(Escherichia.coli), fill=PostNatalAntibiotics)) + geom_boxplot(lwd=1,aes(color=factor(PostNatalAntibiotics),fill = NA), outlier.size = 3)  +
  stat_summary(fun=mean, geom="point", shape=5, size=8) + scale_fill_manual(values = col) + scale_colour_manual(values = col ) +
  geom_point(size=4,aes(color = factor(PostNatalAntibiotics))) + xlab(NULL) + #facet_grid(rows = . ~ Perinatal_Antimicrobial_Use) +
  ylab("Escherichia.coli\n") + paramsBoxWide() + scale_y_log10() + xlab(NULL) + facet_grid(rows = SampleType ~ SampleCollectionWeek) + scale_color_tableau() +
  theme(axis.title.y = element_text(face = "italic"))
ggsave(EcoliAbx, file = "EcoliAbx.pdf", width = 12, height = 10, limitsize = FALSE)

DummySpeciesNR$Bacteroides.vulgatus[DummySpeciesNR$Bacteroides.vulgatus==0]<-1
BvulgatusAbx<-ggplot(subset(DummySpeciesNR, DummySpeciesNR$SampleType %in% c("Axilla", "Groin", "Stool")),
                     aes(x=PostNatalAntibiotics, y=as.numeric(Bacteroides.vulgatus), fill=PostNatalAntibiotics)) + geom_boxplot(lwd=1,aes(color=factor(PostNatalAntibiotics),fill = NA), outlier.size = 3)  +
  stat_summary(fun=mean, geom="point", shape=5, size=8) + scale_fill_manual(values = col) + scale_colour_manual(values = col ) +
  geom_point(size=4,aes(color = factor(PostNatalAntibiotics))) + xlab(NULL) + #facet_grid(rows = . ~ Perinatal_Antimicrobial_Use) +
  ylab("Bacteroides.vulgatus\n") + paramsBoxWide() + scale_y_log10() + xlab(NULL) + facet_grid(rows = SampleType ~ SampleCollectionWeek) + scale_color_tableau() +
  theme(axis.title.y = element_text(face = "italic"))
ggsave(BvulgatusAbx, file = "BvulgatusAbx.pdf", width = 12, height = 10, limitsize = FALSE)

DummySpeciesNR$Veillonella.atypica[DummySpeciesNR$Veillonella.atypica==0]<-1
VatypicaAbx<-ggplot(subset(DummySpeciesNR, DummySpeciesNR$SampleType %in% c("Axilla", "Groin", "Stool")),
                    aes(x=PostNatalAntibiotics, y=as.numeric(Veillonella.atypica), fill=PostNatalAntibiotics)) + geom_boxplot(lwd=1,aes(color=factor(PostNatalAntibiotics),fill = NA), outlier.size = 3)  +
  stat_summary(fun=mean, geom="point", shape=5, size=8) + scale_fill_manual(values = col) + scale_colour_manual(values = col ) +
  geom_point(size=4,aes(color = factor(PostNatalAntibiotics))) + xlab(NULL) + #facet_grid(rows = . ~ Perinatal_Antimicrobial_Use) +
  ylab("Veillonella.atypica\n") + paramsBoxWide() + scale_y_log10() + xlab(NULL) + facet_grid(rows = SampleType ~ SampleCollectionWeek) + scale_color_tableau() +
  theme(axis.title.y = element_text(face = "italic"))
ggsave(VatypicaAbx, file = "VatypicaAbx.pdf", width = 12, height = 10, limitsize = FALSE)


# DummySpeciesNR$Pseudomonas.aeruginosa[DummySpeciesNR$Pseudomonas.aeruginosa==0]<-1
PseudoAbx<-ggplot(subset(DummySpeciesNR, DummySpeciesNR$SampleType %in% c("Axilla", "Groin", "Stool")),
                  aes(x=PostNatalAntibiotics, y=as.numeric(Pseudomonas.aeruginosa), fill=PostNatalAntibiotics)) + geom_boxplot(lwd=1,aes(color=factor(PostNatalAntibiotics),fill = NA), outlier.size = 3)  +
  stat_summary(fun=mean, geom="point", shape=5, size=8) + scale_fill_manual(values = col) + scale_colour_manual(values = col ) +
  geom_point(size=4,aes(color = factor(PostNatalAntibiotics))) + xlab(NULL) + #facet_grid(rows = . ~ Perinatal_Antimicrobial_Use) +
  ylab("Pseudomonas.aeruginosa\n") + paramsBoxWide() + scale_y_log10() + xlab(NULL) + facet_grid(rows = SampleType ~ SampleCollectionWeek) + scale_color_tableau() +
  theme(axis.title.y = element_text(face = "italic"))
ggsave(PseudoAbx, file = "PseudoAbx.pdf", width = 12, height = 10, limitsize = FALSE)

col<-tableau.colors[c(3,4)]
PseudoAbx2<-ggplot(subset(DummySpeciesNR, DummySpeciesNR$SampleType %in% c("Axilla", "Groin", "Stool")),
                   aes(x=SampleCollectionWeek, y=as.numeric(Pseudomonas.aeruginosa), fill=SampleCollectionWeek)) + geom_boxplot(lwd=1,aes(color=factor(SampleCollectionWeek),fill = NA), outlier.size = 3)  +
  stat_summary(fun=mean, geom="point", shape=5, size=8) + scale_fill_manual(values = col) + scale_colour_manual(values = col ) +
  geom_point(size=4,aes(color = factor(SampleCollectionWeek))) + xlab(NULL) + #facet_grid(rows = . ~ Perinatal_Antimicrobial_Use) +
  ylab("Pseudomonas.aeruginosa\n") + paramsBoxWide() + scale_y_log10() + xlab(NULL) + facet_grid(rows = PostNatalAntibiotics ~ SampleType) + 
  theme(axis.title.y = element_text(face = "italic")) + scale_fill_manual(values = tableau.colors[c(3,4)])
ggsave(PseudoAbx2, file = "PseudoAbx2.pdf", width = 12, height = 10, limitsize = FALSE)

DummySpeciesNR$Escherichia.fergusonii[DummySpeciesNR$Escherichia.fergusonii==0]<-1
EfergAbx<-ggplot(subset(DummySpeciesNR, DummySpeciesNR$SampleType %in% c("Axilla", "Groin", "Stool")),
                 aes(x=PostNatalAntibiotics, y=as.numeric(Escherichia.fergusonii), fill=PostNatalAntibiotics)) + geom_boxplot(lwd=1,aes(color=factor(PostNatalAntibiotics),fill = NA), outlier.size = 3)  +
  stat_summary(fun=mean, geom="point", shape=5, size=8) + scale_fill_manual(values = col) + scale_colour_manual(values = col ) +
  geom_point(size=4,aes(color = factor(PostNatalAntibiotics))) + xlab(NULL) + #facet_grid(rows = . ~ Perinatal_Antimicrobial_Use) +
  ylab("Escherichia.fergusonii\n") + paramsBoxWide() + scale_y_log10() + xlab(NULL) + facet_grid(rows = SampleType ~ SampleCollectionWeek) + scale_color_tableau() +
  theme(axis.title.y = element_text(face = "italic"))
ggsave(EfergAbx, file = "EfergAbx.pdf", width = 12, height = 10, limitsize = FALSE)

DummySpeciesNR$Streptococcus.mitis[DummySpeciesNR$Streptococcus.mitis==0]<-1
SmitisAbx<-ggplot(subset(DummySpeciesNR, DummySpeciesNR$SampleType %in% c("Axilla", "Groin", "Stool")),
                  aes(x=PostNatalAntibiotics, y=as.numeric(Streptococcus.mitis), fill=PostNatalAntibiotics)) + geom_boxplot(lwd=1,aes(color=factor(PostNatalAntibiotics),fill = NA), outlier.size = 3)  +
  stat_summary(fun=mean, geom="point", shape=5, size=8) + scale_fill_manual(values = col) + scale_colour_manual(values = col ) +
  geom_point(size=4,aes(color = factor(PostNatalAntibiotics))) + xlab(NULL) + #facet_grid(rows = . ~ Perinatal_Antimicrobial_Use) +
  ylab("Streptococcus.mitis\n") + paramsBoxWide() + scale_y_log10() + xlab(NULL) + facet_grid(rows = SampleType ~ SampleCollectionWeek) + scale_color_tableau() +
  theme(axis.title.y = element_text(face = "italic"))
ggsave(SmitisAbx, file = "SmitisAbxAbx.pdf", width = 12, height = 10, limitsize = FALSE)


# Location graphs

EcoliLocation<-ggplot(subset(DummySpeciesNR, DummySpeciesNR$SampleType %in% c("Axilla", "Groin", "Stool")),
                      aes(x=Location, y=as.numeric(Escherichia.coli), fill=Location)) + geom_boxplot(lwd=1,aes(color=factor(Location),fill = NA), outlier.size = 3)  +
  stat_summary(fun=mean, geom="point", shape=5, size=8) + scale_fill_manual(values = col) + scale_colour_manual(values = col ) +
  geom_point(size=4,aes(color = factor(Location))) + xlab(NULL) + #facet_grid(rows = . ~ Perinatal_Antimicrobial_Use) +
  ylab("Escherichia.coli\n") + paramsBoxWide()  + xlab(NULL) + facet_grid(rows = SampleType ~ SampleCollectionWeek) + scale_color_tableau() +
  theme(axis.title.y = element_text(face = "italic")) + scale_y_log10()
ggsave(EcoliLocation, file = "EcoliLocation.pdf", width = 12, height = 10, limitsize = FALSE)

EcoliTable<-aggregate(Escherichia.coli ~ PatientID + Location + SampleCollectionWeek, data = NICUSpeciesNR, sum) # for SSW

DummySpeciesNR<-NICUSpeciesNR
StaphLocation<-ggplot(subset(DummySpeciesNR, DummySpeciesNR$SampleType %in% c("Axilla", "Groin", "Stool")),
                      aes(x=Location, y=as.numeric(Staphylococcus.aureus ), fill=Location)) + geom_boxplot(lwd=1,aes(color=factor(Location),fill = NA), outlier.size = 3)  +
  stat_summary(fun=mean, geom="point", shape=5, size=8) + scale_fill_manual(values = col) + scale_colour_manual(values = col ) +
  geom_point(size=4,aes(color = factor(Location))) + xlab(NULL) + #facet_grid(rows = . ~ Perinatal_Antimicrobial_Use) +
  ylab("Staphylococcus.aureus \n") + paramsBoxWide() + scale_y_log10() + xlab(NULL) + facet_grid(rows = . ~ SampleCollectionWeek) + scale_color_tableau() +
  theme(axis.title.y = element_text(face = "italic"))
ggsave(StaphLocation, file = "StaphLocation.pdf", width = 8, height = 6, limitsize = FALSE)

StaphLocationSampleType<-ggplot(subset(DummySpeciesNR, DummySpeciesNR$SampleType %in% c("Axilla", "Groin", "Stool")),
                                aes(x=Location, y=as.numeric(Staphylococcus.aureus ), fill=Location)) + geom_boxplot(lwd=1,aes(color=factor(Location),fill = NA), outlier.size = 3)  +
  stat_summary(fun=mean, geom="point", shape=5, size=8) + scale_fill_manual(values = col) + scale_colour_manual(values = col ) +
  geom_point(size=4,aes(color = factor(Location))) + xlab(NULL) + #facet_grid(rows = . ~ Perinatal_Antimicrobial_Use) +
  ylab("Staphylococcus.aureus \n") + paramsBoxWide() + scale_y_log10() + xlab(NULL) + facet_grid(rows = SampleType ~ SampleCollectionWeek) + scale_color_tableau() +
  theme(axis.title.y = element_text(face = "italic"))
ggsave(StaphLocationSampleType, file = "StaphLocationSampleType.pdf", width = 8, height = 12, limitsize = FALSE)


StaphLocation<-ggplot(subset(DummySpeciesNR, DummySpeciesNR$SampleType %in% c("Axilla", "Groin", "Stool")),
                      aes(x=Location, y=as.numeric(Staphylococcus.aureus ), fill=Location)) + geom_boxplot(lwd=1,aes(color=factor(Location),fill = NA), outlier.size = 3)  +
  stat_summary(fun=mean, geom="point", shape=5, size=8) + scale_fill_manual(values = col) + scale_colour_manual(values = col ) +
  geom_point(size=4,aes(color = factor(Location))) + xlab(NULL) + #facet_grid(rows = . ~ Perinatal_Antimicrobial_Use) +
  ylab("Staphylococcus.aureus \n") + paramsBoxWide() + scale_y_log10() + xlab(NULL) + facet_grid(rows = . ~ SampleCollectionWeek) + scale_color_tableau() +
  theme(axis.title.y = element_text(face = "italic"))
ggsave(StaphLocation, file = "StaphLocation.pdf", width = 8, height = 6, limitsize = FALSE)

StaphEpiLocationSampleType<-ggplot(subset(DummySpeciesNR, DummySpeciesNR$SampleType %in% c("Axilla", "Groin", "Stool")),
                                   aes(x=Location, y=as.numeric(Staphylococcus.epidermidis ), fill=Location)) + geom_boxplot(lwd=1,aes(color=factor(Location),fill = NA), outlier.size = 3)  +
  stat_summary(fun=mean, geom="point", shape=5, size=8) + scale_fill_manual(values = col) + scale_colour_manual(values = col ) +
  geom_point(size=4,aes(color = factor(Location))) + xlab(NULL) + #facet_grid(rows = . ~ Perinatal_Antimicrobial_Use) +
  ylab("Staphylococcus.epidermidis \n") + paramsBoxWide() + scale_y_log10() + xlab(NULL) + facet_grid(rows = SampleType ~ SampleCollectionWeek) + scale_color_tableau() +
  theme(axis.title.y = element_text(face = "italic"))
ggsave(StaphEpiLocationSampleType, file = "StaphEpiLocationSampleType.pdf", width = 8, height = 12, limitsize = FALSE)


StaphEpiLocation<-ggplot(subset(DummySpeciesNR, DummySpeciesNR$SampleType %in% c("Axilla", "Groin", "Stool")),
                         aes(x=Location, y=as.numeric(Staphylococcus.epidermidis ), fill=Location)) + geom_boxplot(lwd=1,aes(color=factor(Location),fill = NA), outlier.size = 3)  +
  stat_summary(fun=mean, geom="point", shape=5, size=8) + scale_fill_manual(values = col) + scale_colour_manual(values = col ) +
  geom_point(size=4,aes(color = factor(Location))) + xlab(NULL) + #facet_grid(rows = . ~ Perinatal_Antimicrobial_Use) +
  ylab("Staphylococcus.epidermidis \n") + paramsBoxWide() + scale_y_log10() + xlab(NULL) + facet_grid(rows = . ~ SampleCollectionWeek) + scale_color_tableau() +
  theme(axis.title.y = element_text(face = "italic"))
ggsave(StaphEpiLocation, file = "StaphEpiLocation.pdf", width = 8, height = 6, limitsize = FALSE)


StaphEpiLocationSampleType<-ggplot(subset(DummySpeciesNR, DummySpeciesNR$SampleType %in% c("Axilla", "Groin", "Stool")),
                                   aes(x=Location, y=as.numeric(Staphylococcus.epidermidis ), fill=Location)) + geom_boxplot(lwd=1,aes(color=factor(Location),fill = NA), outlier.size = 3)  +
  stat_summary(fun=mean, geom="point", shape=5, size=8) + scale_fill_manual(values = col) + scale_colour_manual(values = col ) +
  geom_point(size=4,aes(color = factor(Location))) + xlab(NULL) + #facet_grid(rows = . ~ Perinatal_Antimicrobial_Use) +
  ylab("Staphylococcus.epidermidis \n") + paramsBoxWide() + scale_y_log10() + xlab(NULL) + facet_grid(rows = SampleType ~ SampleCollectionWeek) + scale_color_tableau() +
  theme(axis.title.y = element_text(face = "italic"))
ggsave(StaphEpiLocationSampleType, file = "StaphEpiLocationSampleType.pdf", width = 8, height = 12, limitsize = FALSE)


KlebPneumoniaeLocation<-ggplot(subset(DummySpeciesNR, DummySpeciesNR$SampleType %in% c("Axilla", "Groin", "Stool")),
                               aes(x=Location, y=as.numeric(Klebsiella.pneumoniae ), fill=Location)) + geom_boxplot(lwd=1,aes(color=factor(Location),fill = NA), outlier.size = 3)  +
  stat_summary(fun=mean, geom="point", shape=5, size=8) + scale_fill_manual(values = col) + scale_colour_manual(values = col ) +
  geom_point(size=4,aes(color = factor(Location))) + xlab(NULL) + #facet_grid(rows = . ~ Perinatal_Antimicrobial_Use) +
  ylab("Klebsiella.pneumoniae \n") + paramsBoxWide() + scale_y_log10() + xlab(NULL) + facet_grid(rows = . ~ SampleCollectionWeek) + scale_color_tableau() +
  theme(axis.title.y = element_text(face = "italic"))
ggsave(KlebPneumoniaeLocation, file = "KlebPneumoniaeLocation.pdf", width = 8, height = 6, limitsize = FALSE)

KlebLocationSampleType<-ggplot(subset(DummySpeciesNR, DummySpeciesNR$SampleType %in% c("Axilla", "Groin", "Stool")),
                               aes(x=Location, y=as.numeric(Klebsiella.pneumoniae ), fill=Location)) + geom_boxplot(lwd=1,aes(color=factor(Location),fill = NA), outlier.size = 3)  +
  stat_summary(fun=mean, geom="point", shape=5, size=8) + scale_fill_manual(values = col) + scale_colour_manual(values = col ) +
  geom_point(size=4,aes(color = factor(Location))) + xlab(NULL) + #facet_grid(rows = . ~ Perinatal_Antimicrobial_Use) +
  ylab("Klebsiella.pneumoniae \n") + paramsBoxWide() + scale_y_log10() + xlab(NULL) + facet_grid(rows = SampleType ~ SampleCollectionWeek) + scale_color_tableau() +
  theme(axis.title.y = element_text(face = "italic"))
ggsave(KlebLocationSampleType, file = "KlebLocationSampleType.pdf", width = 8, height = 12, limitsize = FALSE)

KlebOxLocationSampleType<-ggplot(subset(DummySpeciesNR, DummySpeciesNR$SampleType %in% c("Axilla", "Groin", "Stool")),
                                 aes(x=Location, y=as.numeric(Klebsiella.oxytoca ), fill=Location)) + geom_boxplot(lwd=1,aes(color=factor(Location),fill = NA), outlier.size = 3)  +
  stat_summary(fun=mean, geom="point", shape=5, size=8) + scale_fill_manual(values = col) + scale_colour_manual(values = col ) +
  geom_point(size=4,aes(color = factor(Location))) + xlab(NULL) + #facet_grid(rows = . ~ Perinatal_Antimicrobial_Use) +
  ylab("Klebsiella oxytoca \n") + paramsBoxWide() + scale_y_log10() + xlab(NULL) + facet_grid(rows = SampleType ~ SampleCollectionWeek) + scale_color_tableau() +
  theme(axis.title.y = element_text(face = "italic"))
ggsave(KlebOxLocationSampleType, file = "KlebOxLocationSampleType.pdf", width = 8, height = 12, limitsize = FALSE)


EcloacaeLocationSampleType<-ggplot(subset(DummySpeciesNR, DummySpeciesNR$SampleType %in% c("Axilla", "Groin", "Stool")),
                                   aes(x=Location, y=as.numeric(Enterobacter.cloacae ), fill=Location)) + geom_boxplot(lwd=1,aes(color=factor(Location),fill = NA), outlier.size = 3)  +
  stat_summary(fun=mean, geom="point", shape=5, size=8) + scale_fill_manual(values = col) + scale_colour_manual(values = col ) +
  geom_point(size=4,aes(color = factor(Location))) + xlab(NULL) + #facet_grid(rows = . ~ Perinatal_Antimicrobial_Use) +
  ylab("Enterobacter.cloacae \n") + paramsBoxWide() + scale_y_log10() + xlab(NULL) + facet_grid(rows = SampleType ~ SampleCollectionWeek) + scale_color_tableau() +
  theme(axis.title.y = element_text(face = "italic"))
ggsave(EcloacaeLocationSampleType, file = "EcloacaeLocationSampleType.pdf", width = 8, height = 12, limitsize = FALSE)


AcinetoLocationSampleType<-ggplot(subset(DummySpeciesNR, DummySpeciesNR$SampleType %in% c("Axilla", "Groin", "Stool")),
                                  aes(x=Location, y=as.numeric(Acinetobacter.baumannii), fill=Location)) + geom_boxplot(lwd=1,aes(color=factor(Location),fill = NA), outlier.size = 3)  +
  stat_summary(fun=mean, geom="point", shape=5, size=8) + scale_fill_manual(values = col) + scale_colour_manual(values = col ) +
  geom_point(size=4,aes(color = factor(Location))) + xlab(NULL) + #facet_grid(rows = . ~ Perinatal_Antimicrobial_Use) +
  ylab("Acinetobacter.baumannii\n") + paramsBoxWide() + scale_y_log10() + xlab(NULL) + facet_grid(rows = SampleType ~ SampleCollectionWeek) + scale_color_tableau() +
  theme(axis.title.y = element_text(face = "italic"))
ggsave(AcinetoLocationSampleType, file = "AcinetoLocationSampleType.pdf", width = 8, height = 12, limitsize = FALSE)

GBSLocationSampleType<-ggplot(subset(DummySpeciesNR, DummySpeciesNR$SampleType %in% c("Axilla", "Groin", "Stool")),
                              aes(x=Location, y=as.numeric(Streptococcus.agalactiae), fill=Location)) + geom_boxplot(lwd=1,aes(color=factor(Location),fill = NA), outlier.size = 3)  +
  stat_summary(fun=mean, geom="point", shape=5, size=8) + scale_fill_manual(values = col) + scale_colour_manual(values = col ) +
  geom_point(size=4,aes(color = factor(Location))) + xlab(NULL) + #facet_grid(rows = . ~ Perinatal_Antimicrobial_Use) +
  ylab("Streptococcus.agalactiae\n") + paramsBoxWide() + xlab(NULL) + facet_grid(rows = SampleType ~ SampleCollectionWeek) + scale_color_tableau() +
  theme(axis.title.y = element_text(face = "italic")) + scale_y_log10()
ggsave(GBSLocationSampleType, file = "GBSLocationSampleType.pdf", width = 8, height = 12, limitsize = FALSE)

EfaecalisLocationSampleType<-ggplot(subset(DummySpeciesNR, DummySpeciesNR$SampleType %in% c("Axilla", "Groin", "Stool")),
                                    aes(x=Location, y=as.numeric(Enterococcus.faecalis), fill=Location)) + geom_boxplot(lwd=1,aes(color=factor(Location),fill = NA), outlier.size = 3)  +
  stat_summary(fun=mean, geom="point", shape=5, size=8) + scale_fill_manual(values = col) + scale_colour_manual(values = col ) +
  geom_point(size=4,aes(color = factor(Location))) + xlab(NULL) + #facet_grid(rows = . ~ Perinatal_Antimicrobial_Use) +
  ylab("Enterococcus faecalis\n") + paramsBoxWide() + scale_y_log10() + xlab(NULL) + facet_grid(rows = SampleType ~ SampleCollectionWeek) + scale_color_tableau() +
  theme(axis.title.y = element_text(face = "italic"))
ggsave(EfaecalisLocationSampleType, file = "EfaecalisLocationSampleType.pdf", width = 8, height = 12, limitsize = FALSE)

EfaeciumLocationSampleType<-ggplot(subset(DummySpeciesNR, DummySpeciesNR$SampleType %in% c("Axilla", "Groin", "Stool")),
                                   aes(x=Location, y=as.numeric(Enterococcus.faecium), fill=Location)) + geom_boxplot(lwd=1,aes(color=factor(Location),fill = NA), outlier.size = 3)  +
  stat_summary(fun=mean, geom="point", shape=5, size=8) + scale_fill_manual(values = col) + scale_colour_manual(values = col ) +
  geom_point(size=4,aes(color = factor(Location))) + xlab(NULL) + #facet_grid(rows = . ~ Perinatal_Antimicrobial_Use) +
  ylab("Enterococcus faecium\n") + paramsBoxWide() + scale_y_log10() + xlab(NULL) + facet_grid(rows = SampleType ~ SampleCollectionWeek) + scale_color_tableau() +
  theme(axis.title.y = element_text(face = "italic"))
ggsave(EfaeciumLocationSampleType, file = "EfaeciumLocationSampleType.pdf", width = 8, height = 12, limitsize = FALSE)



SerratiamLocationSampleType<-ggplot(subset(DummySpeciesNR, DummySpeciesNR$SampleType %in% c("Axilla", "Groin", "Stool")),
                                    aes(x=Location, y=as.numeric(Serratia.marcescens), fill=Location)) + geom_boxplot(lwd=1,aes(color=factor(Location),fill = NA), outlier.size = 3)  +
  stat_summary(fun=mean, geom="point", shape=5, size=8) + scale_fill_manual(values = col) + scale_colour_manual(values = col ) +
  geom_point(size=4,aes(color = factor(Location))) + xlab(NULL) + #facet_grid(rows = . ~ Perinatal_Antimicrobial_Use) +
  ylab("Serratia marcescens\n") + paramsBoxWide() + scale_y_log10() + xlab(NULL) + facet_grid(rows = SampleType ~ SampleCollectionWeek) + scale_color_tableau() +
  theme(axis.title.y = element_text(face = "italic"))
ggsave(SerratiamLocationSampleType, file = "SerratiamLocationSampleType.pdf", width = 8, height = 12, limitsize = FALSE)

SpyogenesLocationSampleType<-ggplot(subset(DummySpeciesNR, DummySpeciesNR$SampleType %in% c("Axilla", "Groin", "Stool")),
                                    aes(x=Location, y=as.numeric(Streptococcus.pyogenes), fill=Location)) + geom_boxplot(lwd=1,aes(color=factor(Location),fill = NA), outlier.size = 3)  +
  stat_summary(fun=mean, geom="point", shape=5, size=8) + scale_fill_manual(values = col) + scale_colour_manual(values = col ) +
  geom_point(size=4,aes(color = factor(Location))) + xlab(NULL) + #facet_grid(rows = . ~ Perinatal_Antimicrobial_Use) +
  ylab("Streptococcus.pyogenes\n") + paramsBoxWide()  + xlab(NULL) + facet_grid(rows = SampleType ~ SampleCollectionWeek) + scale_color_tableau() +
  theme(axis.title.y = element_text(face = "italic")) + scale_y_log10()
ggsave(SpyogenesLocationSampleType, file = "SpyogenesLocationSampleType.pdf", width = 8, height = 12, limitsize = FALSE)

CparapsilosisLocationSampleType<-ggplot(subset(DummySpeciesNR, DummySpeciesNR$SampleType %in% c("Axilla", "Groin", "Stool")),
                                        aes(x=Location, y=as.numeric(Candida.parapsilosis), fill=Location)) + geom_boxplot(lwd=1,aes(color=factor(Location),fill = NA), outlier.size = 3)  +
  stat_summary(fun=mean, geom="point", shape=5, size=8) + scale_fill_manual(values = col) + scale_colour_manual(values = col ) +
  geom_point(size=4,aes(color = factor(Location))) + xlab(NULL) + #facet_grid(rows = . ~ Perinatal_Antimicrobial_Use) +
  ylab("Candida.parapsilosis\n") + paramsBoxWide()  + xlab(NULL) + facet_grid(rows = SampleType ~ SampleCollectionWeek) + scale_color_tableau() +
  theme(axis.title.y = element_text(face = "italic")) + scale_y_log10()
ggsave(CparapsilosisLocationSampleType, file = "CparapsilosisLocationSampleType.pdf", width = 8, height = 12, limitsize = FALSE)

CalbicansLocationSampleType<-ggplot(subset(DummySpeciesNR, DummySpeciesNR$SampleType %in% c("Axilla", "Groin", "Stool")),
                                    aes(x=Location, y=as.numeric(Candida.albicans), fill=Location)) + geom_boxplot(lwd=1,aes(color=factor(Location),fill = NA), outlier.size = 3)  +
  stat_summary(fun=mean, geom="point", shape=5, size=8) + scale_fill_manual(values = col) + scale_colour_manual(values = col ) +
  geom_point(size=4,aes(color = factor(Location))) + xlab(NULL) + #facet_grid(rows = . ~ Perinatal_Antimicrobial_Use) +
  ylab("Candida.albicans\n") + paramsBoxWide()  + xlab(NULL) + facet_grid(rows = SampleType ~ SampleCollectionWeek) + scale_color_tableau() +
  theme(axis.title.y = element_text(face = "italic")) + scale_y_log10()
ggsave(CalbicansLocationSampleType, file = "CalbicansLocationSampleType.pdf", width = 8, height = 12, limitsize = FALSE)




# Calculate statistical significance of differnece in each species by sample type and week (courtesy of ChatGPT)
species_list <- c("Staphylococcus.aureus", "Enterococcus.faecalis", "Enterococcus.faecium", "Serratia.marcescens","Streptococcus.agalactiae", "Enterobacter.cloacae",
                  "Acinetobacter.baumannii", "Klebsiella.pneumoniae", "Klebsiella.oxytoca", "Staphylococcus.epidermidis", "Staphylococcus.aureus", "Escherichia.coli", 
                  "Candida.albicans", "Candida.parapsilosis", "Staphylococcus.capitis","Streptococcus.pyogenes" )


# Loop through each species
for (species in species_list) {
  
  # Initialize an empty list to store results for the current species
  species_results <- list()
  
  # Loop through unique SampleType and SampleCollectionWeek combinations
  for (sample_type in unique(DummySpeciesNR$SampleType)) {
    for (week in unique(DummySpeciesNR$SampleCollectionWeek)) {
      
      # Subset data for the given SampleType and SampleCollectionWeek
      subset_data <- DummySpeciesNR %>%
        filter(SampleType == sample_type, SampleCollectionWeek == week)
      
      # Ensure there are at least two unique locations to compare
      if (length(unique(subset_data$Location)) > 1) {
        
        # Kruskal-Wallis test
        kw_test <- kruskal.test(subset_data[[species]] ~ subset_data$Location)
        
        # Wilcoxon rank-sum test (if only two locations exist)
        wilcox_test <- NULL
        if (length(unique(subset_data$Location)) == 2) {
          wilcox_test <- wilcox.test(subset_data[[species]] ~ subset_data$Location)
        }
        
        # Store results with updated significance notation
        species_results[[paste(sample_type, week, sep = "_")]] <- list(
          "SampleType_Week" = paste(sample_type, week, sep = "_"),
          "Kruskal-Wallis p-value" = kw_test$p.value,
          "KW_Significance" = ifelse(kw_test$p.value < 0.001, "***",
                                     ifelse(kw_test$p.value < 0.01, "**",
                                            ifelse(kw_test$p.value < 0.05, "*", ""))),
          "Wilcoxon p-value" = if (!is.null(wilcox_test)) wilcox_test$p.value else NA,
          "Wilcox_Significance" = if (!is.null(wilcox_test)) {
            ifelse(wilcox_test$p.value < 0.001, "***",
                   ifelse(wilcox_test$p.value < 0.01, "**",
                          ifelse(wilcox_test$p.value < 0.05, "*", "")))
          } else NA
        )
      }
    }
  }
  
  # Convert results to a dataframe
  species_df <- do.call(rbind, lapply(names(species_results), function(name) {
    as.data.frame(species_results[[name]])
  }))
  
  # Assign the dataframe to a variable named after the species
  assign(species, species_df, envir = .GlobalEnv)
  
  # Print the results for review
  print(paste("Results for:", species))
  print(species_df)
  
  # Save results to a CSV file with flagged significant values
  write.csv(species_df, paste0(species, "_results.csv"), row.names = FALSE)
}

## Graph some genera 

EneterococcusLocationSampleType<-ggplot(subset(DummyGenusNR, DummyGenusNR$SampleType %in% c("Axilla", "Groin", "Stool")),
                                        aes(x=Location, y=as.numeric(Enterococcus), fill=Location)) + geom_boxplot(lwd=1,aes(color=factor(Location),fill = NA), outlier.size = 3)  +
  stat_summary(fun=mean, geom="point", shape=5, size=8) + scale_fill_manual(values = col) + scale_colour_manual(values = col ) +
  geom_point(size=4,aes(color = factor(Location))) + xlab(NULL) + #facet_grid(rows = . ~ Perinatal_Antimicrobial_Use) +
  ylab("Enterococcus\n") + paramsBoxWide() + scale_y_log10() + xlab(NULL) + facet_grid(rows = SampleType ~ SampleCollectionWeek) + scale_color_tableau() +
  theme(axis.title.y = element_text(face = "italic"))
ggsave(EneterococcusLocationSampleType, file = "EneterococcusLocationSampleType.pdf", width = 8, height = 12, limitsize = FALSE)


StaphylococcusLocationSampleType<-ggplot(subset(DummyGenusNR, DummyGenusNR$SampleType %in% c("Axilla", "Groin", "Stool")),
                                         aes(x=Location, y=as.numeric(Staphylococcus), fill=Location)) + geom_boxplot(lwd=1,aes(color=factor(Location),fill = NA), outlier.size = 3)  +
  stat_summary(fun=mean, geom="point", shape=5, size=8) + scale_fill_manual(values = col) + scale_colour_manual(values = col ) +
  geom_point(size=4,aes(color = factor(Location))) + xlab(NULL) + #facet_grid(rows = . ~ Perinatal_Antimicrobial_Use) +
  ylab("Staphylococcus\n") + paramsBoxWide() + scale_y_log10() + xlab(NULL) + facet_grid(rows = SampleType ~ SampleCollectionWeek) + scale_color_tableau() +
  theme(axis.title.y = element_text(face = "italic"))
ggsave(StaphylococcusLocationSampleType, file = "StaphylococcusLocationSampleType.pdf", width = 8, height = 12, limitsize = FALSE)


KlebsiellaLocationSampleType<-ggplot(subset(DummyGenusNR, DummyGenusNR$SampleType %in% c("Axilla", "Groin", "Stool")),
                                     aes(x=Location, y=as.numeric(Klebsiella), fill=Location)) + geom_boxplot(lwd=1,aes(color=factor(Location),fill = NA), outlier.size = 3)  +
  stat_summary(fun=mean, geom="point", shape=5, size=8) + scale_fill_manual(values = col) + scale_colour_manual(values = col ) +
  geom_point(size=4,aes(color = factor(Location))) + xlab(NULL) + #facet_grid(rows = . ~ Perinatal_Antimicrobial_Use) +
  ylab("Klebsiella\n") + paramsBoxWide() + scale_y_log10() + xlab(NULL) + facet_grid(rows = SampleType ~ SampleCollectionWeek) + scale_color_tableau() +
  theme(axis.title.y = element_text(face = "italic"))
ggsave(KlebsiellaLocationSampleType, file = "KlebsiellaLocationSampleType.pdf", width = 8, height = 12, limitsize = FALSE)





# Calculate statistical significance of differnece in each genus by sample type and week (courtesy of ChatGPT)
genus_list <- c("Staphylococcus", "Enterococcus", "Serratia","Streptococcus", "Enterobacter",
                "Acinetobacter", "Klebsiella", "Staphylococcus", "Escherichia")


# Loop through each genus
for (genus in genus_list) {
  
  # Initialize an empty list to store results for the current genus
  genus_results <- list()
  
  # Loop through unique SampleType and SampleCollectionWeek combinations
  for (sample_type in unique(DummyGenusNR$SampleType)) {
    for (week in unique(DummyGenusNR$SampleCollectionWeek)) {
      
      # Subset data for the given SampleType and SampleCollectionWeek
      subset_data <- DummyGenusNR %>%
        filter(SampleType == sample_type, SampleCollectionWeek == week)
      
      # Ensure there are at least two unique locations to compare
      if (length(unique(subset_data$Location)) > 1) {
        
        # Kruskal-Wallis test
        kw_test <- kruskal.test(subset_data[[genus]] ~ subset_data$Location)
        
        # Wilcoxon rank-sum test (if only two locations exist)
        wilcox_test <- NULL
        if (length(unique(subset_data$Location)) == 2) {
          wilcox_test <- wilcox.test(subset_data[[genus]] ~ subset_data$Location)
        }
        
        # Store results with updated significance notation
        genus_results[[paste(sample_type, week, sep = "_")]] <- list(
          "SampleType_Week" = paste(sample_type, week, sep = "_"),
          "Kruskal-Wallis p-value" = kw_test$p.value,
          "KW_Significance" = ifelse(kw_test$p.value < 0.001, "***",
                                     ifelse(kw_test$p.value < 0.01, "**",
                                            ifelse(kw_test$p.value < 0.05, "*", ""))),
          "Wilcoxon p-value" = if (!is.null(wilcox_test)) wilcox_test$p.value else NA,
          "Wilcox_Significance" = if (!is.null(wilcox_test)) {
            ifelse(wilcox_test$p.value < 0.001, "***",
                   ifelse(wilcox_test$p.value < 0.01, "**",
                          ifelse(wilcox_test$p.value < 0.05, "*", "")))
          } else NA
        )
      }
    }
  }
  
  # Convert results to a dataframe
  genus_df <- do.call(rbind, lapply(names(genus_results), function(name) {
    as.data.frame(genus_results[[name]])
  }))
  
  # Assign the dataframe to a variable named after the genus
  assign(genus, genus_df, envir = .GlobalEnv)
  
  # Print the results for review
  print(paste("Results for:", genus))
  print(genus_df)
  
  # Save results to a CSV file with flagged significant values
  write.csv(genus_df, paste0(genus, "_results.csv"), row.names = FALSE)
}


### Ok check out maturation of microbiome from antibiotic negative babies at different gestational age

AxillaGroinStoolNoAbx<-subset(NICUSpeciesNR, NICUSpeciesNR$SampleType %in% c("Axilla") &
                                NICUSpeciesNR$PostNatalAntibiotics == "No.Infant.Abx")
AxillaGroinStoolNoAbx$GestationTypeWeek<-paste(AxillaGroinStoolNoAbx$GestationCohort, AxillaGroinStoolNoAbx$TypeWeek, sep = "-")
GestationTypeWeekCol<-grep("GestationTypeWeek", colnames(AxillaGroinStoolNoAbx))
AxillaGroinStoolNoAbx$GestationTypeWeek<-factor(AxillaGroinStoolNoAbx$GestationTypeWeek, levels = 
                                                  c("28-32 Weeks-Axilla-Week.1", "33-36 Weeks-Axilla-Week.1", "28-32 Weeks-Axilla-Week.3", "33-36 Weeks-Axilla-Week.3"))

AxillaGroinStoolNoAbx<-AxillaGroinStoolNoAbx[,c(1:17,GestationTypeWeekCol,18:(GestationTypeWeekCol-1))]

metadata <- AxillaGroinStoolNoAbx[, 1:18]
cts <- as.matrix(AxillaGroinStoolNoAbx[, -(1:18)])
rownames(cts) <- metadata$SampleID
cts_l2 <- glog2(cts)
grps <- AxillaGroinStoolNoAbx$GestationTypeWeek

col=tableau.colors
newPCA<-PCA(cts_l2, scale.unit = FALSE, ncp = 5, graph = FALSE)
pcaPlot<-fviz_pca_ind(newPCA,
                      geom.ind = "point", 
                      pointsize = 1.5,
                      point.alph = 0.1,
                      title = "Unsupervised Principal Coordinate Analysis",
                      subtitle = "Gestational Age and Sample Type : No Antibiotics Samples",
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
SampleTypeGestation<- pcaPlot 

pdf("GestationWeekAxilla.pdf") 
print(SampleTypeGestation)
dev.off() 


AxillaGroinStoolNoAbx<-subset(NICUSpeciesNR, NICUSpeciesNR$SampleType %in% c("Stool") &
                                NICUSpeciesNR$PostNatalAntibiotics == "No.Infant.Abx")
AxillaGroinStoolNoAbx$GestationTypeWeek<-paste(AxillaGroinStoolNoAbx$GestationCohort, AxillaGroinStoolNoAbx$TypeWeek, sep = "-")
GestationTypeWeekCol<-grep("GestationTypeWeek", colnames(AxillaGroinStoolNoAbx))
AxillaGroinStoolNoAbx$GestationTypeWeek<-factor(AxillaGroinStoolNoAbx$GestationTypeWeek, levels = 
                                                  c("28-32 Weeks-Stool-Week.1", "33-36 Weeks-Stool-Week.1", "28-32 Weeks-Stool-Week.3", "33-36 Weeks-Stool-Week.3"))

AxillaGroinStoolNoAbx<-AxillaGroinStoolNoAbx[,c(1:17,GestationTypeWeekCol,18:(GestationTypeWeekCol-1))]

metadata <- AxillaGroinStoolNoAbx[, 1:18]
cts <- as.matrix(AxillaGroinStoolNoAbx[, -(1:18)])
rownames(cts) <- metadata$SampleID
cts_l2 <- glog2(cts)
grps <- AxillaGroinStoolNoAbx$GestationTypeWeek

col=rep(tableau.colors,"blue", "green")
newPCA<-PCA(cts_l2, scale.unit = FALSE, ncp = 5, graph = FALSE)
pcaPlot<-fviz_pca_ind(newPCA,
                      geom.ind = "point", 
                      pointsize = 1.5,
                      point.alph = 0.1,
                      title = "Unsupervised Principal Coordinate Analysis",
                      subtitle = "Gestational Age and Sample Type : No Antibiotics Samples",
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
SampleTypeGestation<- pcaPlot 

pdf("GestationWeekStool.pdf") 
print(SampleTypeGestation)
dev.off() 

save.image(file = "NICUData20250209")



