## SHaP to identify clinical features associated with specific microbe abundance

import pandas as pd
import numpy as np
import shap
import seaborn as sns
import matplotlib.pyplot as plt
from sklearn.ensemble import RandomForestRegressor
from sklearn.preprocessing import StandardScaler, LabelEncoder
from sklearn.model_selection import train_test_split

# Load microbiome abundance data (rows: samples, columns: species/genus)
microbiome_df = pd.read_csv("microbiome_abundance.csv", index_col=0)

# Load metadata
metadata_df = pd.read_csv("metadata.csv", index_col=0)

# Merge microbiome data and metadata on Sample ID
data = microbiome_df.merge(metadata_df, left_index=True, right_index=True)

# Select a specific microbe abundance as the outcome variable
target_microbe = "Staphylococcus.aureus"  # Replace with any species of interest
y = data[target_microbe]

# Select clinical metadata as predictors
X = metadata_df.copy()

# Encode categorical variables (e.g., feeding, antibiotic use)
categorical_features = ["Location", "SampleType", "SampleCollectionWeek", "GestationCohort", "PostNatalAbxCohort", 
                        "BSI_30D", "NEC_30D", "AnyMilk", "PICC", "UVC"]
for col in categorical_features:
    X[col] = LabelEncoder().fit_transform(X[col])

# Standardize numerical variables
scaler = StandardScaler()
X_scaled = scaler.fit_transform(X)

# Split data into training and test sets
X_train, X_test, y_train, y_test = train_test_split(X_scaled, y, test_size=0.2, random_state=42)

# Train a Random Forest Regressor
rf = RandomForestRegressor(n_estimators=100, random_state=42)
rf.fit(X_train, y_train)

# Compute SHAP values
explainer = shap.Explainer(rf, X_train)
shap_values = explainer(X_test)

# SHAP Summary Plot: Identifying key clinical variables associated with Staphylococcus aureus
shap.summary_plot(shap_values, X_test, feature_names=X.columns)

# SHAP Feature Importance Bar Plot
shap_values_mean = np.abs(shap_values.values).mean(axis=0)
shap_importance = pd.DataFrame({"Feature": X.columns, "SHAP Importance": shap_values_mean})
shap_importance = shap_importance.sort_values(by="SHAP Importance", ascending=False)

plt.figure(figsize=(10, 6))
sns.barplot(x="SHAP Importance", y="Feature", data=shap_importance[:15], palette="viridis")
plt.xlabel("SHAP Importance Score")
plt.ylabel("Features")
plt.title(f"Clinical Variables Associated with {target_microbe} Abundance")
plt.show()

# Correlation Heatmap: Checking raw correlations
correlation_matrix = pd.concat([y, X], axis=1).corr()
plt.figure(figsize=(8, 6))
sns.heatmap(correlation_matrix[[target_microbe]].sort_values(by=target_microbe, ascending=False), 
            annot=True, cmap="coolwarm", vmin=-1, vmax=1)
plt.title(f"Correlation of Clinical Variables with {target_microbe}")
plt.show()
