# SHAP (SHapley Additive exPlanations) for Feature Importance

import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import seaborn as sns
from sklearn.ensemble import RandomForestClassifier
from sklearn.inspection import permutation_importance
from sklearn.preprocessing import StandardScaler, LabelEncoder
from sklearn.model_selection import train_test_split

print("Starting feature importance analysis...")

# Load the metadata
metadata_df = pd.read_csv("metadata/AllNICUSampleKeyRevised20250320_for_HanzhouCincinnatiSamples.csv")
print(f"Loaded metadata with shape: {metadata_df.shape}")

# Define categorical features
categorical_features = [
    "SampleType", "Location", "GestationCohort", "SampleCollectionWeek", 
    "PostNatalAbxCohort", "TypeWeek", "TypeAbx", "AnyMilk", 
    "PICC", "UVC"
]

# Convert categorical variables to numeric
for col in categorical_features:
    if col in metadata_df.columns:
        metadata_df[col] = metadata_df[col].astype(str)
        le = LabelEncoder()
        metadata_df[col] = le.fit_transform(metadata_df[col])

# Convert the GestationalAge to numeric (weeks)
metadata_df['GestationalAge_Numeric'] = metadata_df['GestationalAge'].str.extract(r'(\d+)').astype(float)
print(f"Converted GestationalAge to numeric")

# We'll analyze what factors are associated with sample type 
print(f"SampleType values: {metadata_df['SampleType'].value_counts()}")

# Select features for the model, excluding SampleType itself
feature_cols = [col for col in categorical_features if col in metadata_df.columns and col != 'SampleType']
feature_cols.append('GestationalAge_Numeric')
print(f"Using features: {feature_cols}")

# Select rows with complete data
predictors = metadata_df[feature_cols].copy()
target = metadata_df['SampleType']  # Predict the sample type

# Handle missing values
for col in predictors.columns:
    if predictors[col].isnull().sum() > 0:
        print(f"Column {col} has {predictors[col].isnull().sum()} missing values")
        if pd.api.types.is_numeric_dtype(predictors[col]):
            predictors[col] = predictors[col].fillna(predictors[col].median())
        else:
            predictors[col] = predictors[col].fillna(predictors[col].mode()[0])

print(f"Final data shape: {predictors.shape} with target shape: {target.shape}")

# Split the data
X_train, X_test, y_train, y_test = train_test_split(
    predictors, target, test_size=0.3, random_state=42, stratify=target
)

# Train a RandomForest model
print("Training Random Forest model...")
model = RandomForestClassifier(n_estimators=100, random_state=42)
model.fit(X_train, y_train)

# Get feature importance from the model directly
feature_importance = pd.DataFrame({
    'Feature': predictors.columns,
    'Importance': model.feature_importances_
})
feature_importance = feature_importance.sort_values('Importance', ascending=False)

print("Model's built-in feature importance:")
print(feature_importance)

# Calculate permutation importance
print("Calculating permutation importance...")
perm_importance = permutation_importance(
    model, X_test, y_test, n_repeats=10, random_state=42
)

# Create permutation importance dataframe
perm_importance_df = pd.DataFrame({
    'Feature': predictors.columns,
    'Importance': perm_importance.importances_mean
})
perm_importance_df = perm_importance_df.sort_values('Importance', ascending=False)

print("Permutation feature importance:")
print(perm_importance_df)

# Save results to CSV
feature_importance.to_csv("results/feature_importance_sampletype.csv", index=False)
perm_importance_df.to_csv("results/permutation_importance_sampletype.csv", index=False)
print("Saved feature importance results to CSV files")

# Generate visualizations
plt.figure(figsize=(10, 6))
sns.barplot(x='Importance', y='Feature', data=feature_importance)
plt.title('Feature Importance for Sample Type Prediction')
plt.tight_layout()
plt.savefig('results/feature_importance_plot.pdf')

plt.figure(figsize=(10, 6))
sns.barplot(x='Importance', y='Feature', data=perm_importance_df)
plt.title('Permutation Feature Importance for Sample Type Prediction')
plt.tight_layout()
plt.savefig('results/permutation_importance_plot.pdf')

print("Generated and saved feature importance plots")
print("Feature importance analysis completed successfully.")
