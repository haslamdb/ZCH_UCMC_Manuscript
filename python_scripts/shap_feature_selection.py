# SHAP (SHapley Additive exPlanations) for Feature Importance in Microbiome Differences


import pandas as pd
import numpy as np
import scipy.spatial
import shap
import matplotlib.pyplot as plt
import seaborn as sns
from sklearn.ensemble import RandomForestRegressor
from sklearn.preprocessing import StandardScaler, LabelEncoder
from sklearn.model_selection import train_test_split

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
predictors = metadata_df.drop(columns=["body_site", "postnatal_age"])  

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

# Compute SHAP values
explainer = shap.Explainer(rf, X_train)
shap_values = explainer(X_test)

# SHAP Summary Plot
plt.figure(figsize=(10, 6))
shap.summary_plot(shap_values, X_test, feature_names=predictors.columns)

# SHAP Feature Importance Bar Plot
shap_values_mean = np.abs(shap_values.values).mean(axis=0)
shap_importance = pd.DataFrame({"Feature": predictors.columns, "SHAP Importance": shap_values_mean})
shap_importance = shap_importance.sort_values(by="SHAP Importance", ascending=False)

plt.figure(figsize=(10, 6))
sns.barplot(x="SHAP Importance", y="Feature", data=shap_importance[:15], palette="viridis")
plt.xlabel("SHAP Importance Score")
plt.ylabel("Features")
plt.title("Top Features Driving Microbiome Differences (SHAP)")
plt.show()

# SHAP Dependence Plot for a Key Feature
shap.dependence_plot("antibiotics", shap_values, X_test, feature_names=predictors.columns)
