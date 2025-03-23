# Random forest with feature selection
# assess importance of various clinical features on microbiome differences
# in this case, we're retaining body site and postnatal age as features, but this could be modified
# in this case we have not scaled the count data. Consider doing so

import pandas as pd
import numpy as np
import scipy.spatial
from sklearn.ensemble import RandomForestRegressor
from sklearn.preprocessing import StandardScaler, LabelEncoder
from sklearn.model_selection import train_test_split
import matplotlib.pyplot as plt
import seaborn as sns
import microbiome_transform as mt

# Load microbiome abundance data (rows: samples, columns: species/genus)
def select_best_transformation(df):
    """
    Automatically selects the best microbiome transformation based on data properties.
    
    Parameters:
    - df: Pandas DataFrame (raw microbiome count data)

    Returns:
    - Transformed DataFrame
    """
    zero_fraction = (df == 0).sum().sum() / df.size  # % of zero values
    total_reads_variation = df.sum(axis=1).std() / df.sum(axis=1).mean()  # Variability in sequencing depth

    if zero_fraction > 0.30:  
        print("⚠️ High zero fraction detected! Using CLR Transformation.")
        return mt.clr_transform(df)
    elif total_reads_variation > 0.5:  
        print("⚠️ Large sequencing depth variation detected! Applying Rarefaction.")
        return mt.rarefaction(df)
    else:
        print("✅ Data appears normalized, applying Total Sum Scaling (TSS).")
        return mt.tss_transform(df)

# Load microbiome count data
microbiome_df = pd.read_csv("../data/NICUSpeciesReduced.csv", index_col=0)

# Auto-select and apply transformation
microbiome_transformed = select_best_transformation(microbiome_df)

# Load clinical metadata
metadata_df = pd.read_csv("../metadata/AllNICUSampleKey20250206.csv", index_col=0)
subject_id_col = "Subject"  

# Check for missing values in metadata
print(metadata_df.isnull().sum())


# transform using CLR, VST, or other method from microbiome_transform module
microbiome_clr = mt.clr_transform(microbiome_df)
microbiome_tss = mt.tss_transform(microbiome_df)
microbiome_vst = mt.vst_transform_r(microbiome_df) # this will fail if zero counts in all features
microbiome_vst = mt.debug_vst_transform(microbiome_df) # has to be run from command prompt, not powershell

# Merge microbiome data and metadata on Sample ID
# here we're using clr transformed data
data = microbiome_clr.merge(metadata_df, left_index=True, right_index=True)

# Calculate pairwise Bray-Curtis distance matrix
bray_curtis_dist = scipy.spatial.distance.pdist(microbiome_clr, metric='braycurtis')
bray_curtis_dist = scipy.spatial.distance.squareform(bray_curtis_dist)  # Convert to square matrix

# Convert distance matrix into a vector (response variable)
y = bray_curtis_dist[np.triu_indices_from(bray_curtis_dist, k=1)]  # Take upper triangle (pairwise distances)

# Prepare metadata for regression (drop body site and postnatal age to remove dominance)
# predictors = metadata_df.drop(columns=["body_site", "postnatal_age"])  


# Encode categorical variables
categorical_features = ["SampleType", "Location", "GestationCohort", "SampleCollectionWeek", 
                        "MaternalAntibiotics", "PostNatalAbxCohort", "BSI_30D", "NEC_30D", "AnyMilk", "PICC", "UVC"]

predictors = metadata_df[categorical_features]

for col in categorical_features:
    predictors[col] = LabelEncoder().fit_transform(predictors[col])

# Scale numerical variables
scaler = StandardScaler()
X_scaled = scaler.fit_transform(predictors)

# Convert X into pairwise differences
X_diff = np.abs(X_scaled[:, None, :] - X_scaled[None, :, :])  # Pairwise absolute differences
X_diff = X_diff[np.triu_indices_from(bray_curtis_dist, k=1)]  # Take upper triangle

# Split data into training and test sets
X_train, X_test, y_train, y_test = train_test_split(X_diff, y, test_size=0.2, random_state=42)

# Train Random Forest Regressor
rf = RandomForestRegressor(n_estimators=100, random_state=42)
rf.fit(X_train, y_train)

# Feature importance analysis
feature_importances = pd.DataFrame({"Feature": predictors.columns, "Importance": rf.feature_importances_})
feature_importances = feature_importances.sort_values(by="Importance", ascending=False)

# Plot feature importances
plt.figure(figsize=(10, 6))
sns.barplot(x="Importance", y="Feature", data=feature_importances[:15], palette="viridis")
plt.xlabel("Feature Importance Score")
plt.ylabel("Features")
plt.title("Top Features Driving Microbiome Differences")
plt.show()

# Create line+dot plot (PyCaret style) for feature importance
plt.figure(figsize=(12, 8))
top_n = min(15, len(feature_importances))
top_features = feature_importances.head(top_n)

top_features = top_features.iloc[::-1]


# Create horizontal line + dot plot
plt.hlines(y=range(top_n), 
           xmin=0, 
           xmax=top_features["Importance"].values,  
           color="skyblue", 
           alpha=0.7, 
           linewidth=2)

plt.plot(top_features["Importance"].values,  
         range(top_n), 
         "o", 
         markersize=10, 
         color="blue", 
         alpha=0.8)

# Add feature names
plt.yticks(range(top_n), top_features["Feature"].values)
plt.xlabel("Feature Importance Score")
plt.title("Clinical Variables Associated with Microbiome Composition Differences")  
# Add values next to dots
for i, importance in enumerate(top_features["Importance"].values):  
    plt.text(importance + 0.001, i, f"{importance:.4f}", va='center')
    
plt.tight_layout()
plt.savefig("../results/feature_importance_plot.pdf", bbox_inches="tight")  
plt.show()
plt.close()

#