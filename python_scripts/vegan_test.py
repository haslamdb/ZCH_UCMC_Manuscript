import pandas as pd
import numpy as np
import os
from vegan import adonis2
from scipy.spatial.distance import squareform
from skbio.stats.distance import DistanceMatrix

print("Testing vegan integration...")

# Check for available data files
data_files = []
for root, dirs, files in os.walk('.'):
    if '.git' in dirs:
        dirs.remove('.git')  # Skip git directory
    for file in files:
        if file.endswith('.csv'):
            data_files.append(os.path.join(root, file))

print(f"Found {len(data_files)} CSV files:")
for file in data_files[:10]:  # Show just the first 10
    print(f"  - {file}")

# Try to load the microbiome normalized data
try:
    microbiome_df = pd.read_csv("results/variance_analysis/microbiome_normalized_presence_absence.csv", index_col=0)
    print(f"Loaded microbiome data with shape: {microbiome_df.shape}")
except FileNotFoundError:
    print("Microbiome normalized file not found. Looking for alternative...")
    
    # Look for any microbiome file
    microbiome_files = [f for f in data_files if 'microbiome' in f.lower()]
    if microbiome_files:
        print(f"Found alternative microbiome files: {microbiome_files}")
        microbiome_df = pd.read_csv(microbiome_files[0], index_col=0)
        print(f"Loaded alternative microbiome data with shape: {microbiome_df.shape}")
    else:
        print("No microbiome data files found.")
        exit(1)

# Try to test adonis2 with a distance matrix
print("\nTesting adonis2 with a small distance matrix...")

# Create a small test dataset
n_samples = 20
np.random.seed(42)
group_a = np.random.normal(0, 1, (10, 5))
group_b = np.random.normal(2, 1, (10, 5))
data = np.vstack([group_a, group_b])
sample_ids = [f'sample_{i}' for i in range(n_samples)]

# Create distance matrix
from scipy.spatial.distance import pdist, squareform
dist = pdist(data, metric='euclidean')
dist_matrix = squareform(dist)

# Create metadata
metadata = pd.DataFrame({
    'group': ['A'] * 10 + ['B'] * 10,
    'value': np.random.normal(0, 1, n_samples)
}, index=sample_ids)

# Convert to distance matrix format
skdist = DistanceMatrix(dist_matrix, ids=sample_ids)

try:
    print("\nPreparing distance matrix...")
    from rpy2.robjects.numpy2ri import numpy2rpy
    from rpy2.robjects import r
    import rpy2.robjects as ro
    
    # Convert distance matrix to R dist object
    dist_vector = ro.FloatVector(skdist.condensed_form())
    r_dist = r['as.dist'](r['matrix'](dist_vector, 
                                   nrow=int(n_samples*(n_samples-1)/2)))
    
    print("Running adonis2 directly through R...")
    # Run adonis2 in R directly
    from rpy2.robjects.packages import importr
    vegan = importr('vegan')
    
    # Convert metadata to R
    from rpy2.robjects import pandas2ri
    from rpy2.robjects.conversion import localconverter
    with localconverter(ro.default_converter + pandas2ri.converter):
        r_metadata = pandas2ri.py2rpy(metadata)
    
    # Run adonis2 through R
    result = vegan.adonis2(r_dist, r_metadata, permutations=99)
    
    # Convert result back to pandas
    with localconverter(ro.default_converter + pandas2ri.converter):
        pd_result = pandas2ri.rpy2py(result)
    
    print("Success! Result:")
    print(pd_result)
    
except Exception as e:
    print(f"Error in direct R test: {e}")

print("\nVegan integration test complete.")