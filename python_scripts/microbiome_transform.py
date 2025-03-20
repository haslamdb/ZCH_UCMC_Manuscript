# Module for microbiome count data transformations
import numpy as np
import pandas as pd
from scipy.stats import gmean
from sklearn.utils import resample
import subprocess

def clr_transform(df, pseudocount=1):
    """
    Perform Centered Log-Ratio (CLR) transformation on microbiome data.
    
    Parameters:
    - df: Pandas DataFrame (samples as rows, species as columns)
    - pseudocount: Small value added to avoid log(0) issues (default: 1)

    Returns:
    - Transformed DataFrame
    """
    df = df + pseudocount  # Add pseudocount to avoid log(0)
    clr_df = np.log(df) - np.log(gmean(df, axis=1, keepdims=True))
    return pd.DataFrame(clr_df, index=df.index, columns=df.columns)

def tss_transform(df):
    """
    Perform Total Sum Scaling (TSS) normalization (relative abundance).
    
    Parameters:
    - df: Pandas DataFrame (samples as rows, species as columns)
    
    Returns:
    - Transformed DataFrame
    """
    return df.div(df.sum(axis=1), axis=0)

def vst_transform_r(input_csv, output_csv="microbiome_abundance_vst.csv"):
    """
    Perform Variance-Stabilizing Transformation (VST) using DESeq2 in R.
    
    Parameters:
    - input_csv: Path to input CSV file (raw microbiome counts)
    - output_csv: Path to output CSV file (VST-transformed data)

    Returns:
    - None (saves transformed data as CSV)
    """
    r_script = f"""
    library(DESeq2)
    df <- read.csv("{input_csv}", row.names=1)
    dds <- DESeqDataSetFromMatrix(countData = df, colData = NULL, design = ~ 1)
    vsd <- varianceStabilizingTransformation(dds)
    write.csv(assay(vsd), "{output_csv}")
    """
    with open("vst_transform.R", "w") as f:
        f.write(r_script)
    
    subprocess.run(["Rscript", "vst_transform.R"], check=True)
    print(f"VST-transformed data saved to {output_csv}")

def log_transform(df, pseudocount=1):
    """
    Perform log normalization on microbiome count data.
    
    Parameters:
    - df: Pandas DataFrame (samples as rows, species as columns)
    - pseudocount: Small value added to avoid log(0) issues (default: 1)

    Returns:
    - Log-transformed DataFrame
    """
    return np.log1p(df + pseudocount)  # log1p avoids log(0) issues

def rarefaction(df, min_depth=None):
    """
    Perform rarefaction (subsampling) to normalize sequencing depth across samples.
    
    Parameters:
    - df: Pandas DataFrame (samples as rows, species as columns)
    - min_depth: Minimum sequencing depth to subsample to (default: min(sample read counts))

    Returns:
    - Rarefied DataFrame (subsampled counts)
    """
    if min_depth is None:
        min_depth = df.sum(axis=1).min()  # Set to minimum total read count across samples

    rarefied_df = df.apply(lambda row: resample(row, replace=False, n_samples=min_depth) if row.sum() >= min_depth else row, axis=1)
    return rarefied_df.fillna(0).astype(int)
