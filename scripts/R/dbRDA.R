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

# Add these before running the capscale function

# Option 1: Hellinger transformation
# Good for abundance data, downweights rare species
species_data_hellinger <- decostand(species_data, method = "hellinger")

# Option 2: CLR (Centered Log-Ratio) transformation
# For compositional data, handles zero values with pseudocounts
# First add a small pseudocount to zeros
species_data_clr <- decostand(species_data + 0.5, method = "clr")

# Option 3: Wisconsin double standardization
# Commonly used for community data
species_data_wisconsin <- wisconsin(species_data)

# Option 4: Log transformation
# Reduces the impact of extremely abundant taxa
species_data_log <- log1p(species_data)  # log(x+1) to handle zeros





##################################################
# 4) Run partial dbRDA with capscale to control for Subject
##################################################
# Using Condition(SubjectID) to control for repeated measures
dbrda_result <- capscale(species_data_hellinger ~ SampleType + Location + GestationCohort + 
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

# Get a summary of the analysis including residuals information
summary_result <- summary(dbrda_result)
print(summary_result)

# Extract the proportion of variation explained vs. residual
constrained <- sum(dbrda_result$CCA$eig)
unconstrained <- sum(dbrda_result$CA$eig)
total <- constrained + unconstrained
print(paste("Proportion constrained:", round(constrained/total, 4)))
print(paste("Proportion residual (unexplained):", round(unconstrained/total, 4)))

# Extract residuals distance of each sample from the model
residuals_dist <- sqrt(rowSums(scores(dbrda_result, display="lc", choices=1:ncol(scores(dbrda_result, display="lc")))^2))
sample_residuals <- data.frame(SampleID = sample_ids, Residual = residuals_dist)
sample_residuals <- sample_residuals[order(-sample_residuals$Residual), ]
head(sample_residuals, 10)  # Show 10 samples with largest residuals

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
sites_df$SampleType <- metadata$SampleType 
sites_df$SampleCollectionWeek <- metadata$SampleCollectionWeek
sites_df$Location <- metadata$Location


# # Create the plot
ggplot(sites_df, aes(x = CAP1, y = CAP2, color = SampleType, shape = SampleCollectionWeek)) +
  geom_point(size = 3, alpha = 0.7) +
  geom_text_repel(aes(label = SampleID), size = 3) +
  theme_minimal() +
  labs(title = "dbRDA of Microbiome Data (Controlled for Subject)",
       x = "dbRDA1", y = "dbRDA2")


# Plot variance partitioning
# Plot eigenvalues

# 1. Extract the eigenvalues and calculate proportion of variance explained
eigenvals <- eigenvals(dbrda_result)
explainedvar <- eigenvals / sum(eigenvals)
barplot(explainedvar[1:10], names.arg = paste0("Axis ", 1:10),
        main = "Proportion of Variance Explained by Each Axis")

# 2. Run the permutation test for individual terms
set.seed(123)  # For reproducibility
perm_test <- anova(dbrda_result, by="terms", permutations=999)
print(perm_test)

# 3. Calculate the proportion of constrained variance explained by each variable

# Look at the structure of perm_test to see the actual column names
print(names(perm_test))

# Check the full output to understand the structure
print(perm_test)

# Modified approach - using the correct column names
# Usually the R² or SumOfSqs column represents the variance explained
importance <- data.frame(
  Variable = rownames(perm_test)[-nrow(perm_test)],  # Exclude the last row (residual)
  F_value = perm_test$F[-nrow(perm_test)],
  P_value = perm_test$`Pr(>F)`[-nrow(perm_test)]
)

# If SumOfSqs exists, use it as a proxy for variance explained
if("SumOfSqs" %in% names(perm_test)) {
  importance$Variance_Explained <- perm_test$SumOfSqs[-nrow(perm_test)] / 
    sum(perm_test$SumOfSqs[-nrow(perm_test)])
} else if("R2" %in% names(perm_test)) {
  importance$Variance_Explained <- perm_test$R2[-nrow(perm_test)] / 
    sum(perm_test$R2[-nrow(perm_test)])
} else {
  # If neither exists, we'll just use F-values as a proxy
  importance$Variance_Explained <- perm_test$F[-nrow(perm_test)] / 
    sum(perm_test$F[-nrow(perm_test)])
  print("Note: Using F-values as proxy for variance explained")
}

# Sort by variance explained (descending)
importance <- importance[order(-importance$Variance_Explained), ]
print(importance)

write.csv(importance, file = "dbRDA_Importance.csv")

# Create a bar plot of variance explained by each variable
# Sort data by Variance_Explained in descending order
# Sort data by Variance_Explained in descending order (highest values first)
# Sort data by Variance_Explained in descending order
importance_sorted <- importance[order(-importance$Variance_Explained),]

# Create PDF device
pdf("results/dbRDA_variance_explained_plot.pdf", width = 10, height = 8)

# Create a horizontal line + dot plot
# Use reverse ordering for y-axis to get highest values at top
plot(importance_sorted$Variance_Explained, 
     seq(nrow(importance_sorted), 1, -1),  # Reverse y-axis ordering
     xlim = c(0, max(importance_sorted$Variance_Explained) * 1.1),
     yaxt = "n",  # Hide y-axis, we'll create our own
     ylab = "",
     xlab = "Proportion of Variance Explained",
     main = "Proportion of Constrained Variance Explained by Each Variable",
     type = "n")  # Don't plot points yet

# Add horizontal lines from y-axis to each point
segments(0, seq(nrow(importance_sorted), 1, -1), 
         importance_sorted$Variance_Explained, seq(nrow(importance_sorted), 1, -1),
         col = "skyblue", lwd = 2)

# Add points for each value
points(importance_sorted$Variance_Explained, seq(nrow(importance_sorted), 1, -1),
       pch = 16, col = "blue", cex = 1.5)

# Create y-axis with variable names - now with reversed ordering
axis(2, at = seq(nrow(importance_sorted), 1, -1), 
     labels = importance_sorted$Variable, 
     las = 2, cex.axis = 0.7)

# Add significance indicators - using the reversed y positions
significant <- importance_sorted$P_value < 0.05
significant_positions <- which(significant)
if(length(significant_positions) > 0) {
  text(importance_sorted$Variance_Explained[significant] + max(importance_sorted$Variance_Explained) * 0.02, 
       seq(nrow(importance_sorted), 1, -1)[significant_positions], 
       labels = "*", 
       col = "red", 
       cex = 1.5)
}

# Add grid lines for better readability
grid(nx = NULL, ny = NULL, lty = 2, col = "lightgray")

# Close the PDF device
dev.off()


# Sort by variance explained (descending)
importance <- importance[order(-importance$Variance_Explained), ]
print(importance)

# 4. Create a bar plot of variance explained by each variable
barplot(importance$Variance_Explained, 
        names.arg = importance$Variable, 
        las = 2,  # Rotate labels for better readability
        cex.names = 0.7,  # Smaller font for variable names
        main = "Proportion of Constrained Variance Explained by Each Variable")

# 5. For categorical variables with multiple levels, look at scores of centroids
centroids <- scores(dbrda_result, display = "cn")
if(!is.null(centroids)) {
  print("Centroid distances from origin:")
  # Calculate distance of each centroid from the origin
  centroid_distances <- sqrt(rowSums(centroids^2))
  centroid_importance <- data.frame(
    Level = names(centroid_distances),
    Distance = centroid_distances
  )
  # Sort by distance (descending)
  centroid_importance <- centroid_importance[order(-centroid_importance$Distance), ]
  print(centroid_importance)
}



save.image(file = "dbRDA.RData")
