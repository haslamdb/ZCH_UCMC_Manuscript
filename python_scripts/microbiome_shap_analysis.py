# microbiome_shap_analysis.py

import pandas as pd
import numpy as np
import shap
import microbiome_transform as mt  # Import the transformation module
import seaborn as sns
import matplotlib.pyplot as plt
from sklearn.ensemble import RandomForestRegressor
from sklearn.preprocessing import StandardScaler, LabelEncoder
from sklearn.model_selection import train_test_split

def select_best_transformation(df):
    """
    Automatically selects the best microbiome transformation based on data properties.
    
    Parameters:
    - df: Pandas DataFrame (raw microbiome count data)

    Returns:
    - Transformed DataFrame
    """
    zero_fraction = (df == 0).sum().sum() / df.size  # % of zero values in the dataset
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
microbiome_df=pd.read_csv("../data/NICUSpeciesReduced.csv", index_col=0)


# Auto-select and apply transformation
microbiome_transformed = select_best_transformation(microbiome_df)
microbiome_transformed.to_csv("microbiome_abundance_transformed.csv")

# Load clinical metadata
metadata_df = pd.read_csv("../metadata/AllNICUSampleKey20250206.csv", index_col=0)

# Merge metadata with transformed microbiome data
data = microbiome_transformed.merge(metadata_df, left_index=True, right_index=True)

# Define the target microbe (modify this based on your research question)
target_microbe = "Staphylococcus.aureus"  # Change species as needed
y = data[target_microbe]
X = metadata_df.copy()

# Encode categorical variables
categorical_features = ["SampleID", "Subject","SampleType","Location",  "GestationCohort",   "SampleCollectionWeek", 
                         "MaternalAntibiotics","PostNatalAbxCohort", "BSI_30D", "NEC_30D", "AnyMilk", "PICC", "UVC"]
# keep only the columns listed in categorical_features
X = X[categorical_features]
y = data[target_microbe]




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

# SHAP Summary Plot
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
