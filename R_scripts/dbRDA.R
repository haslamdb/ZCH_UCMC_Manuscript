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
rownames(metadata) <- rownames(community_mat)

##################################################
# 3) Run dbRDA with capscale
##################################################
# capscale(formula, data, distance) is typically used like:
#    capscale(community_data ~ predictors, data = metadata, distance = "...")
# If you want partial dbRDA (controlling for certain variables),
# use: capscale(community_data ~ Var1 + Condition(Var2), ...)

# For full dbRDA on Var1, Var2, Var3:
dbrda_result <- capscale(community_mat ~ Var1 + Var2 + Var3,
                         data = metadata,
                         distance = "bray")

##################################################
# 4) Review the results
##################################################
# Print a summary of the dbRDA
summary(dbrda_result)

# Test the significance of the model (or each term) via permutation
anova(dbrda_result, permutations = 999)       # tests the overall dbRDA
anova(dbrda_result, by="terms", permutations=999)  # tests each predictor

# Inspect the proportion of variance explained by constraints vs. unconstrained
dbrda_result

##################################################
# 5) Basic visualization
##################################################
# Ordination biplot with sites (samples) and biplot arrows for the constraints
plot(dbrda_result, 
     display = c("sites", "bp", "cn"),  # 'bp' = biplot arrows, 'cn' = centroids
     main = "Distance-based Redundancy Analysis (dbRDA)")

# Alternatively, you can extract coordinates for custom plotting:
dbRDA_sites <- scores(dbrda_result, display = "sites")  # sample scores
dbRDA_species <- scores(dbrda_result, display = "species")  # taxa scores
dbRDA_biplot <- scores(dbrda_result, display = "bp")    # biplot arrows

# Now you can use ggplot2 or base R graphics to plot 'dbRDA_sites' etc. 
