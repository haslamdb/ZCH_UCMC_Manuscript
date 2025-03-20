# bsi_microbiome_comparison.R
# -----------------------------------------------------------
# Purpose:
#   1) Read sample metadata (AllSampleKey)
#   2) Import species-level data from Kraken2 (_species_abundance.txt)
#   3) Merge with BSI data, remove contaminants
#   4) Filter, rarefy, compute Bray-Curtis, do PCA
#   5) Conduct basic BSI vs. Microbiome comparisons (O:E ratio)
#   6) Generate some plots (boxplots, PCA)

# NOTE:
#   - The script uses BSI data from "BSIData20250206.csv" plus
#     multiple setwd calls for different local environments. 
#     Comment out the ones you don’t need.
#   - Repetitive code (for instance, re-labelling species) can be
#     consolidated into a utility function in a separate file if desired.
#   - If you have the previously reorganized “Reorganized R Code,” you can
#     integrate these steps or remove duplication.

# -----------------------------------------------------------
# SECTION 0: SETUP
# -----------------------------------------------------------

## (Optional) setwd calls for your environment:
# setwd("C:/Users/dbhas/OneDrive/Documents/Alignments/KrakenAlignments/Kraken2")
# setwd("C:/Users/HASI9S/OneDrive/Documents/Alignments/KrakenAlignments/Kraken2")
# setwd("P:/Documents/Alignments/KrakenAlignments/Kraken2")
# setwd("~/pCloudDrive/Documents/Alignments/KrakenAlignments/Kraken2")

## Libraries you might need (assumes loaded already, else do library(...)):
# library(tidyverse)
# library(vegan)
# library(phyloseq) ... etc.

# -----------------------------------------------------------
# SECTION 1: READ SAMPLE METADATA
# -----------------------------------------------------------

NewSampleKey <- read.csv("AllNICUSampleKeyRevised20250206_for_HangzhouCincinnatiSamples.csv",
                         header = TRUE, stringsAsFactors = FALSE)
NewSampleKey$PostNatalAntibiotics <- factor(
  NewSampleKey$PostNatalAntibiotics,
  levels = c("No.Infant.Abx", "Infant.Abx")
)
NewSampleKey$PostNatalAntibioticsNew <- factor(
  NewSampleKey$PostNatalAntibioticsNew,
  levels = c("No.Infant.Abx", "Low.Infant.Abx", "High.Infant.Abx")
)
NewSampleKey$Sample <- NewSampleKey$SampleID

# Another file that might track antibiotic usage
NewSampleKeyAbx <- read.csv("SampleKeyAbx20231106.csv", header = FALSE)

AllSampleKey <- NewSampleKey

AllSampleKey$GestationCohort <- factor(AllSampleKey$GestationCohort,
                                       levels = c("23-27 Weeks","28-32 Weeks","33-36 Weeks","37-42 Weeks"))
AllSampleKey$SampleCollectionWeek <- factor(AllSampleKey$SampleCollectionWeek,
                                            levels = c("Week.1","Week.3"))
AllSampleKey$SampleType <- factor(AllSampleKey$SampleType,
                                  levels = c("Nares","Axilla","Groin","Stool"))
AllSampleKey$MaternalAntibiotics <- factor(AllSampleKey$MaternalAntibiotics,
                                           levels = c("No.Mat.Abx","Mat.Abx"))

## Create convenience columns
AllSampleKey$TypeWeek <- paste(AllSampleKey$SampleType,
                               AllSampleKey$SampleCollectionWeek, sep = "-")
AllSampleKey$TypeAbx <- paste(AllSampleKey$SampleType,
                              AllSampleKey$PostNatalAntibiotics, sep = "-")
AllSampleKey$TypeAbx <- factor(
  AllSampleKey$TypeAbx,
  levels = c("Axilla-No.Infant.Abx","Axilla-Infant.Abx",
             "Groin-No.Infant.Abx","Groin-Infant.Abx",
             "Stool-No.Infant.Abx","Stool-Infant.Abx",
             "Nares-No.Infant.Abx","Nares-Infant.Abx")
)

# Restrict to Stool, Groin, Axilla (remove duplicates)
AllSampleKey <- subset(
  AllSampleKey,
  !duplicated(AllSampleKey$SampleID) & AllSampleKey$SampleType %in% c("Stool","Groin","Axilla")
)

# -----------------------------------------------------------
# SECTION 2: ADDITIONAL HUMAN-REACTIVE SPECIES
# -----------------------------------------------------------
HumanReactiveSpecies <- read.csv("../HumanReactiveKraken2.csv")
HumanReactiveSpecies$SubSpecies <- gsub("_", ".", HumanReactiveSpecies$SubSpecies)
HumanReactiveSpecies$SubSpecies <- gsub("-", ".", HumanReactiveSpecies$SubSpecies)
HumanReactiveSpecies$SubSpecies <- gsub("/",".",HumanReactiveSpecies$SubSpecies,fixed=TRUE)
HumanReactiveSpecies$SubSpecies <- gsub("X,","",HumanReactiveSpecies$SubSpecies,fixed=TRUE)
HumanReactiveSpecies$SubSpecies <- gsub("\\[|\\]","",HumanReactiveSpecies$SubSpecies)
HumanReactiveSpecies$SubSpecies <- gsub("\\.\\.",".",HumanReactiveSpecies$SubSpecies)
HumanReactiveSpecies$SubSpecies <- gsub(" ",".",HumanReactiveSpecies$SubSpecies,fixed=TRUE)

# -----------------------------------------------------------
# SECTION 3: READ SPECIES COUNTS FROM KRAKEN2
# -----------------------------------------------------------

AllK2Files <- list.files()
SpeciesFileList3 <- grep("_species_abundance.txt", AllK2Files)
K2SpeciesFiles <- AllK2Files[SpeciesFileList3]
K2FileList <- gsub("_species_abundance.txt","",K2SpeciesFiles)

NICUFiles <- as.character(AllSampleKey$SampleID)
K2NICUFiles <- K2FileList[K2FileList %in% NICUFiles]
ALLNICUFiles <- K2NICUFiles

MissingNICU <- subset(NICUFiles, !NICUFiles %in% c(K2NICUFiles))
# Optionally write out missing sample IDs
# write.csv(MissingNICU, file="MissingNICUFiles20210730.txt")

# Loop to read each <sample>_species_abundance.txt
for(f in seq_along(K2NICUFiles)){
  fnr <- paste0(K2NICUFiles[f], "_species_abundance.txt")
  filename <- K2NICUFiles[f]
  # assign to environment - you might prefer storing in a list
  assign(filename,
         read.csv(fnr, sep="\t", header=TRUE, stringsAsFactors=FALSE))
}

# -----------------------------------------------------------
# SECTION 4: MERGE ALL SPECIES TABLES
# -----------------------------------------------------------
NewSpeciesTable <- read.csv(paste0(K2NICUFiles[1], "_species_abundance.txt"),
                            sep="\t", header=TRUE)
names(NewSpeciesTable)[1] <- "Species"
NewSpeciesTable$Species <- gsub(" ", ".", NewSpeciesTable$Species,fixed=TRUE)
NewSpeciesTable$Species <- gsub("\\.\\.",".",NewSpeciesTable$Species)
NewSpeciesTable$Species <- gsub("_",".",NewSpeciesTable$Species)
NewSpeciesTable$Species <- gsub("-",".",NewSpeciesTable$Species)
NewSpeciesTable$Species <- gsub("/",".",NewSpeciesTable$Species,fixed=TRUE)
NewSpeciesTable$Species <- gsub("X,","",NewSpeciesTable$Species,fixed=TRUE)
NewSpeciesTable$Species <- gsub("\\[|\\]","",NewSpeciesTable$Species)
NewSpeciesTable <- as.data.frame(NewSpeciesTable[, c(1,6)])
names(NewSpeciesTable) <- c("Species", K2NICUFiles[1])
NewSpeciesTable <- subset(NewSpeciesTable, !duplicated(NewSpeciesTable$Species))

for(i in 2:length(K2NICUFiles)){
  x <- get(K2NICUFiles[i])  # from assign above
  x <- as.data.frame(x[, c(1,6)])
  names(x) <- c("Species", paste(K2NICUFiles[i]))
  x$Species <- gsub(" ", ".", x$Species, fixed=TRUE)
  x$Species <- gsub("\\.\\.",".", x$Species)
  x$Species <- gsub("_",".",x$Species)
  x$Species <- gsub("-",".",x$Species)
  x$Species <- gsub("/",".",x$Species,fixed=TRUE)
  x$Species <- gsub("X,","",x$Species,fixed=TRUE)
  x$Species <- gsub("\\[|\\]","",x$Species)
  x <- subset(x, !duplicated(x$Species))
  NewSpeciesTable <- merge(x, NewSpeciesTable, by="Species", all=TRUE)
  NewSpeciesTable[is.na(NewSpeciesTable)] <- 0
}

row.names(NewSpeciesTable) <- NewSpeciesTable$Species

# -----------------------------------------------------------
# SECTION 5: MERGE BSI DATA, FILTER, CREATE RAREFIED TABLE
# -----------------------------------------------------------
BSIRaw <- read.csv("BSIData20250206.csv", header=TRUE)
BSI550kx <- BSIRaw[, c(2,3)] * 500000
BSI5k <- as.data.frame(cbind(BSIRaw$Species, BSI550kx))
names(BSI5k) <- c("Species","ZCH.BSI","UC.BSI")

SpeciesBSITable <- merge(NewSpeciesTable, BSI5k, by="Species", all=TRUE)
row.names(SpeciesBSITable) <- SpeciesBSITable$Species
SpeciesBSITable$Species <- NULL

save.image(file="NICUData20250209")

# Transpose for further filtering
SpeciesBSITableTrans <- as.data.frame(t(SpeciesBSITable))

# Remove Salinibacter, Homo, human-reactive from the transposed data
SaliniColBSI <- grep("Salinibacter.ruber", names(SpeciesBSITableTrans))
if(length(SaliniColBSI)>0) {
  SpeciesBSITableTrans <- SpeciesBSITableTrans[,-SaliniColBSI]
}
HumanColBSI <- grep("Homo.sapiens", names(SpeciesBSITableTrans))
if(length(HumanColBSI)>0) {
  SpeciesBSITableTrans <- SpeciesBSITableTrans[,-HumanColBSI]
}
HumanReactiveColsBSI <- which(colnames(SpeciesBSITableTrans) %in% HumanReactiveSpecies$SubSpecies)
if(length(HumanReactiveColsBSI)>0) {
  SpeciesBSITableTrans <- SpeciesBSITableTrans[,-HumanReactiveColsBSI]
}

# Filter: only keep species present in >= 10% of samples
TenPercentCutoff <- floor(nrow(SpeciesBSITableTrans)/20)
NonZeroCounts <- sapply(SpeciesBSITableTrans, function(col) sum(col>0))
TenPercentNotZero <- which(NonZeroCounts >= TenPercentCutoff)
SpeciesBSITableTransReduced <- SpeciesBSITableTrans[,TenPercentNotZero]

# Remove samples with row sums < 10000
LowSamples <- which(rowSums(SpeciesBSITableTransReduced) <= 10000)
SpeciesBSITableTransReduced <- SpeciesBSITableTransReduced[-LowSamples,]
SpeciesBSITableTransReduced[is.na(SpeciesBSITableTransReduced)] <- 0

# Optional noise removal and rarefaction
SpeciesBSITableTransReducedNR <- as.data.frame(
  t(noise.removal(t(SpeciesBSITableTransReduced), 0.001))
)
minCount <- min(rowSums(SpeciesBSITableTransReducedNR))
SpeciesBSITableTransReducedNR <- data.frame(rrarefy(SpeciesBSITableTransReducedNR, 10000))
SpeciesBSITableTransReducedNR$SampleID <- row.names(SpeciesBSITableTransReducedNR)

# Merge with metadata
BSIMetaRows <- which(SpeciesBSITableTransReducedNR$SampleID %in% AllSampleKey$SampleID)
BSISampleKey <- AllSampleKey[BSIMetaRows,]

SpeciesBSITableNR <- merge(AllSampleKey, SpeciesBSITableTransReducedNR,
                           by="SampleID", all.y=TRUE)

# If SampleID has “BSI” in it, mark TypeWeek as BSI
BSIRows <- grep("BSI", SpeciesBSITableNR$SampleID)
SpeciesBSITableNR$TypeWeek[BSIRows] <- "BSI"

# Adjust location for ZCH vs UC
ZCHBSIRow <- grep("ZCH.BSI", SpeciesBSITableNR$SampleID)
UCBSIRow  <- grep("UC.BSI",  SpeciesBSITableNR$SampleID)
SpeciesBSITableNR$Location[ZCHBSIRow] <- "Hangzhou"
SpeciesBSITableNR$Location[UCBSIRow]  <- "Cincinnati"
SpeciesBSITableNR$TypeWeekLocation <- paste(SpeciesBSITableNR$TypeWeek,
                                            SpeciesBSITableNR$Location, sep="-")

# -----------------------------------------------------------
# SECTION 6: BRAY-CURTIS DIST & HEATMAP
# -----------------------------------------------------------
# Example aggregator for centroids
Centroids <- aggregate(
  SpeciesBSITableNR[, c(18:ncol(SpeciesBSITableNR))],
  by=list(SpeciesBSITableNR$TypeWeekLocation),
  FUN=mean
)
row.names(Centroids) <- Centroids$Group.1
Centroids$Group.1 <- NULL
CentroidsTrans <- t(Centroids)
NACol <- grep("NA", colnames(CentroidsTrans))
# if needed: CentroidsTrans <- CentroidsTrans[,-NACol]

CentroidsTrans <- as.matrix(CentroidsTrans)
CentroidsTrans <- as.data.frame(round(CentroidsTrans))
df_imputed <- CentroidsTrans
df_imputed[is.na(df_imputed)] <- 0

df_transposed <- t(df_imputed)
dist_bc <- vegan::vegdist(df_transposed, method="bray")
dist_matrix <- as.matrix(dist_bc)
dist_matrix[dist_matrix==0] <- NA

write.csv(dist_matrix, file="BSIBrayCurtisDist.csv", row.names=TRUE)

my_colors_reversed <- rev(diverging.colors(100))
heatmap <- pheatmap(dist_matrix,
                    color=my_colors_reversed,
                    na_col="white")

pdf("BrayCurtisBSIHeatmap.pdf")
print(heatmap)
dev.off()

# -----------------------------------------------------------
# SECTION 7: OBSERVED VS EXPECTED RATIO (O:E) FOR BSI
# -----------------------------------------------------------
# Example approach for Cincinnati vs Hangzhou comparisons
BSIRowsFlag <- (SpeciesBSITableNR$TypeWeek=="BSI")
SpeciesBSITableNRNoBSI <- SpeciesBSITableNR[!BSIRowsFlag, ]

LocationCentroids <- aggregate(
  SpeciesBSITableNRNoBSI[, c(18:(ncol(SpeciesBSITableNRNoBSI)-1))],
  by=list(SpeciesBSITableNRNoBSI$Location), FUN=mean
)
row.names(LocationCentroids) <- LocationCentroids$Group.1
LocationCentroids$Group.1 <- NULL
LocationCentroidsTrans <- as.data.frame(round(t(LocationCentroids)))
LocationCentroidsTrans$Species <- row.names(LocationCentroidsTrans)

BSIColsCentroids <- grep("BSI", colnames(CentroidsTrans))
BSICentroids <- as.data.frame(CentroidsTrans[, BSIColsCentroids])
BSICentroids$Species <- row.names(BSICentroids)
LocationBSICentroids <- merge(LocationCentroidsTrans, BSICentroids,
                              by="Species", all=TRUE)
LocationBSICentroids <- na.omit(LocationBSICentroids)
row.names(LocationBSICentroids) <- LocationBSICentroids$Species

# Same logic for e.g. Cincinnati vs Hangzhou
# Proportional test snippet below

# 1) Cincinnati
CincinnatiCentroids <- LocationBSICentroids[, grep("Cincinnati", colnames(LocationBSICentroids))]
names(CincinnatiCentroids) <- c("Microbiome.Cincinnati", "BSI.Cincinnati")
BSI_total <- sum(CincinnatiCentroids$BSI.Cincinnati)
Micro_total <- sum(CincinnatiCentroids$Microbiome.Cincinnati)

pvals <- numeric(nrow(CincinnatiCentroids))
obs_exp_ratio <- numeric(nrow(CincinnatiCentroids))
cohens_h <- numeric(nrow(CincinnatiCentroids))

for(i in seq_len(nrow(CincinnatiCentroids))){
  k_i <- CincinnatiCentroids$BSI.Cincinnati[i]
  p_i <- CincinnatiCentroids$Microbiome.Cincinnati[i]/Micro_total
  if(p_i>0 && BSI_total>0) {
    testres <- binom.test(k_i, BSI_total, p_i)
    pvals[i] <- testres$p.value
    p_obs <- k_i/BSI_total
    obs_exp_ratio[i] <- if(p_i==0) Inf else p_obs/p_i
    p_exp <- p_i
    cohens_h[i] <- 2*(asin(sqrt(p_obs)) - asin(sqrt(p_exp)))
  } else {
    pvals[i] <- NA
    obs_exp_ratio[i] <- NA
    cohens_h[i] <- NA
  }
}
pvals_adj <- p.adjust(pvals, method="BH")
CincinnatiCentroids$pval <- pvals
CincinnatiCentroids$pval_adj <- pvals_adj
CincinnatiCentroids$ObsExp_ratio <- obs_exp_ratio
CincinnatiCentroids$Cohens_h <- cohens_h
write.csv(CincinnatiCentroids, file="CincinnatiObservedExpected20250209.csv")

# 2) Hangzhou
HangzhouCentroids <- LocationBSICentroids[, grep("Hangzhou", colnames(LocationBSICentroids))]
names(HangzhouCentroids) <- c("Microbiome.Hangzhou","BSI.Hangzhou")
BSI_total <- sum(HangzhouCentroids$BSI.Hangzhou)
Micro_total <- sum(HangzhouCentroids$Microbiome.Hangzhou)

pvals <- numeric(nrow(HangzhouCentroids))
obs_exp_ratio <- numeric(nrow(HangzhouCentroids))
cohens_h <- numeric(nrow(HangzhouCentroids))

for(i in seq_len(nrow(HangzhouCentroids))){
  k_i <- HangzhouCentroids$BSI.Hangzhou[i]
  p_i <- HangzhouCentroids$Microbiome.Hangzhou[i]/Micro_total
  if(p_i>0 && BSI_total>0){
    testres <- binom.test(k_i, BSI_total, p_i)
    pvals[i] <- testres$p.value
    p_obs <- k_i/BSI_total
    obs_exp_ratio[i] <- if(p_i==0) Inf else p_obs/p_i
    p_exp <- p_i
    cohens_h[i] <- 2*(asin(sqrt(p_obs)) - asin(sqrt(p_exp)))
  } else {
    pvals[i] <- NA
    obs_exp_ratio[i] <- NA
    cohens_h[i] <- NA
  }
}
pvals_adj <- p.adjust(pvals, method="BH")
HangzhouCentroids$pval <- pvals
HangzhouCentroids$pval_adj <- pvals_adj
HangzhouCentroids$ObsExp_ratio <- obs_exp_ratio
HangzhouCentroids$Cohens_h <- cohens_h
write.csv(HangzhouCentroids, file="HangzhouObservedExpected20250209.csv")

save.image(file="NICUData20250209")

# -----------------------------------------------------------
# SECTION 8: PCA & PLOTS BY LOCATION
# -----------------------------------------------------------
# Breaking out data by location, removing BSI, making PCA, etc.

CincinnatiSpeciesBSITable <- subset(SpeciesBSITableNR, Location=="Cincinnati")
CincinnatiSpeciesBSITable <- CincinnatiSpeciesBSITable[
  , c(1, ncol(CincinnatiSpeciesBSITable), 2:(ncol(CincinnatiSpeciesBSITable)-1))
]
CincinnatiSpeciesBSITable$TypeWeekLocation <- factor(
  CincinnatiSpeciesBSITable$TypeWeekLocation,
  levels = c("BSI-Cincinnati","Axilla-Week.1-Cincinnati","Axilla-Week.3-Cincinnati",
             "Groin-Week.1-Cincinnati","Groin-Week.3-Cincinnati",
             "Stool-Week.1-Cincinnati","Stool-Week.3-Cincinnati")
)
CincinnatiSpeciesBSITableNoBSI <- subset(
  CincinnatiSpeciesBSITable,
  TypeWeekLocation != "BSI-Cincinnati"
)
CincinnatiSpeciesBSITableNoBSI$TypeWeekLocation <- factor(
  CincinnatiSpeciesBSITableNoBSI$TypeWeekLocation,
  levels = c("Axilla-Week.1-Cincinnati","Axilla-Week.3-Cincinnati",
             "Groin-Week.1-Cincinnati","Groin-Week.3-Cincinnati",
             "Stool-Week.1-Cincinnati","Stool-Week.3-Cincinnati")
)

metadata <- CincinnatiSpeciesBSITable[, 1:19]
cts <- as.matrix(CincinnatiSpeciesBSITable[, -(1:19)])
rownames(cts) <- metadata$Sample
cts_l2 <- glog2(cts)
grps <- CincinnatiSpeciesBSITable$TypeWeekLocation

col <- c("black", new.col(7)[2:7])
newPCA <- PCA(cts_l2, scale.unit=FALSE, ncp=5, graph=FALSE)
pcaPlot <- fviz_pca_ind(
  newPCA,
  geom.ind="point",
  pointsize=1.5,
  point.alpha=0.1,
  title="Unsupervised Principal Coordinate Analysis",
  subtitle="Samples colored by Body Site and Sample Collection Week",
  col.ind=grps,
  addEllipses=FALSE,
  ellipse.alpha=0.2,
  ellipse.type="confidence",
  ellipse.level=0.95,
  legend.title="Location",
  legend.size=11,
  mean.point=TRUE,
  palette=col,
  axes.linetype="blank"
) + theme(legend.key=element_blank()) +
  scale_shape_manual(values = c(17, rep(16, length(unique(grps))-1)))

CincinnatiBSIPCA <- pcaPlot + scale_fill_manual(values=col)
pdf("CincinnatiBSIPCA.pdf")
print(CincinnatiBSIPCA)
dev.off()

# Repeating for NoBSI
# ...
# The same approach done for Hangzhou
# ...
# Then combining them in one big PCA

# -----------------------------------------------------------
# SECTION 9: QUICK BOX PLOTS FOR BSI-DIFFERENTIAL ORGS
# -----------------------------------------------------------

BSIDifferentialOrgs <- c(
  "Staphylococcus.aureus","Staphylococcus.capitis","Staphylococcus.epidermidis",
  "Streptococcus.pyogenes","Klebsiella.pneumoniae","Klebsiella.oxytoca",
  "Escherichia.coli","Enterococcus.faecium","Enterococcus.faecalis","Candida.parapsilosis",
  "Serratia.marcescens"
)

col <- tableau.colors  # or define your color palette

## If you have a combined dataframe e.g. 'DummySpeciesNR'
## loop for boxplot vs. location
for(org in BSIDifferentialOrgs){
  temp_df <- subset(DummySpeciesNR, SampleType %in% c("Axilla","Groin","Stool"))
  p <- ggplot(temp_df,
              aes_string(
                x="Location",
                y=paste0("as.numeric(", org, ")"),
                color="Location"
              )) +
    geom_boxplot(fill=NA, lwd=1, outlier.size=3) +
    stat_summary(fun=mean, geom="point", shape=5, size=8) +
    scale_colour_manual(values=col) +
    geom_point(size=4) +
    xlab(NULL) +
    ylab(paste0(org, " \n")) +
    stat_compare_means(
      aes(label=paste0("p = ", after_stat(p.format))),
      method="wilcox.test", label.x.npc="center", label.y.npc="top"
    ) +
    paramsBoxWide() +
    scale_y_log10() +
    facet_grid(rows=vars(SampleType), cols=vars(SampleCollectionWeek)) +
    theme(strip.background=element_rect(fill=NA,color=NA),
          strip.text=element_text(size=14, face="bold"),
          axis.title.y=element_text(face="italic"))
  
  out_fn <- paste0(org, "_SampleTypeLocation.pdf")
  ggsave(out_fn, plot=p, width=8, height=12, limitsize=FALSE)
}

## Alternative simpler facet:
for(org in BSIDifferentialOrgs){
  temp_df <- subset(DummySpeciesNR, SampleType %in% c("Axilla","Groin","Stool"))
  p <- ggplot(temp_df,
              aes_string(
                x="Location",
                y=paste0("as.numeric(", org, ")"),
                color="Location"
              )) +
    geom_boxplot(fill=NA, lwd=1, outlier.size=3) +
    stat_summary(fun=mean, geom="point", shape=5, size=8) +
    scale_colour_manual(values=col) +
    geom_point(size=4) +
    xlab(NULL) +
    ylab(paste0(org, " \n")) +
    stat_compare_means(
      aes(label=paste0("p = ", after_stat(p.format))),
      method="wilcox.test", label.x.npc="center", label.y.npc="top"
    ) +
    paramsBoxWide() +
    scale_y_log10() +
    theme(axis.title.y=element_text(face="italic"))
  
  out_fn <- paste0(org, "_Location.pdf")
  ggsave(out_fn, plot=p, width=5, height=11, limitsize=FALSE)
}

# Done. Save final environment:
# save.image(file="NICUData20250209.RData")
