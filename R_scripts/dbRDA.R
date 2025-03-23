##################################################
# 1) Load packages
##################################################

# install.packages("vegan")
library(vegan)

##################################################
# 2) Get data and metadata
##################################################
setwd("C:/users/dbhas//Documents/Code/ZCH_UCMC_Manuscript/")

# Create a community (abundance) matrix
community_mat <- read.csv("data/NICUSpeciesReduced.csv")

# Create a metadata data.frame
metadata <- read.csv("metadata/AllNICUSampleKey20250206.csv")
metadata <- subset(metadata, metadata$SampleID %in% community_mat$SampleID)

community_mat <- subset(community_mat, community_mat$SampleID %in% metadata$SampleID)

##################################################
# 3) Prepare the data for dbRDA
##################################################
# Create a community data matrix (species matrix) without the SampleID column
# First, save the SampleID for later reference
sample_ids <- community_mat$SampleID

# Remove the SampleID column from the community matrix
species_data <- community_mat[, !colnames(community_mat) %in% c("SampleID", "SubjectID", "Subject")]

# Make sure the metadata and community matrix are in the same order
metadata <- metadata[match(sample_ids, metadata$SampleID), ]

##################################################
# 4) Run partial dbRDA with capscale to control for Subject
##################################################
# Using Condition(SubjectID) to control for repeated measures
dbrda_result <- capscale(species_data ~ SampleType + Location + GestationCohort + 
                         SampleCollectionWeek + MaternalAntibiotics + PostNatalAbxCohort + 
                         BSI_30D + NEC_30D + AnyMilk + PICC + UVC + Delivery +
                         Condition(Subject), 
                         data = metadata,
                         distance = "bray")

##################################################
# 5) Review the results
##################################################
# Print a summary of the dbRDA
summary(dbrda_result)

# Test the significance of the model (or each term) via permutation
# For overall model test
anova(dbrda_result, permutations = 999)

# For individual terms test
anova(dbrda_result, by="terms", permutations=999)

# Inspect the proportion of variance explained by constraints vs. unconstrained
dbrda_result

##################################################
# 6) Basic visualization
##################################################
# Ordination biplot with sites (samples) and biplot arrows for the constraints
plot(dbrda_result, 
     display = c("sites", "bp", "cn"),  # 'bp' = biplot arrows, 'cn' = centroids
     main = "Distance-based Redundancy Analysis (dbRDA)")

# Extract coordinates for custom plotting:
dbRDA_sites <- scores(dbrda_result, display = "sites")  # sample scores
dbRDA_species <- scores(dbrda_result, display = "species")  # taxa scores
dbRDA_biplot <- scores(dbrda_result, display = "bp")    # biplot arrows

# Optional: Create a better visualization with ggplot2
library(ggplot2)
library(ggrepel)

# # Convert site scores to a data frame and add metadata
sites_df <- as.data.frame(dbRDA_sites)
sites_df$SampleID <- sample_ids
sites_df$SubjectID <- metadata$SubjectID
sites_df$SampleType <- metadata$SampleType  # Add other relevant variables

# # Create the plot
ggplot(sites_df, aes(x = CAP1, y = CAP2, color = SampleType)) +
  geom_point(size = 3, alpha = 0.7) +
  geom_text_repel(aes(label = SampleID), size = 3) +
  theme_minimal() +
  labs(title = "dbRDA of Microbiome Data (Controlled for Subject)",
       x = "dbRDA1", y = "dbRDA2")

save.image(file = "dbRDA.RData")
