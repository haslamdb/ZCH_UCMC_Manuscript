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

# Load microbiome abundance data (rows: samples, columns: species/genus)
microbiome_df = pd.read_csv("microbiome_abundance.csv", index_col=0)

# Load metadata
metadata_df = pd.read_csv("metadata.csv", index_col=0)

# Merge microbiome data and metadata on Sample ID
data = microbiome_df.merge(metadata_df, left_index=True, right_index=True)

# Calculate pairwise Bray-Curtis distance matrix
bray_curtis_dist = scipy.spatial.distance.pdist(microbiome_df, metric='braycurtis')
bray_curtis_dist = scipy.spatial.distance.squareform(bray_curtis_dist)  # Convert to square matrix

# Convert distance matrix into a vector (response variable)
y = bray_curtis_dist[np.triu_indices_from(bray_curtis_dist, k=1)]  # Take upper triangle (pairwise distances)

# Prepare metadata for regression (drop body site and postnatal age to remove dominance)
# predictors = metadata_df.drop(columns=["body_site", "postnatal_age"])  
predictors = metadata_df # no dropped columns

# Encode categorical variables
categorical_features = ["feeding", "antibiotics", "PICC_placement"]
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
