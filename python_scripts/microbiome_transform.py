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

def vst_transform_r(df, r_path=None, output_csv="microbiome_abundance_vst.csv"):
    """
    Perform Variance-Stabilizing Transformation (VST) using DESeq2 in R.
    
    Parameters:
    - df: Pandas DataFrame with microbiome counts
    - r_path: Path to Rscript executable (e.g., "C:/Program Files/R/R-4.2.2/bin/Rscript.exe")
    - output_csv: Path to output CSV file (VST-transformed data)

    Returns:
    - Transformed DataFrame
    """
    import os
    import subprocess
    import tempfile
    
    # Save input DataFrame to a temporary CSV
    temp_dir = tempfile.mkdtemp()
    input_csv = os.path.join(temp_dir, "temp_input.csv")
    df.to_csv(input_csv)
    
    # Create R script with full path
    r_script_path = os.path.join(temp_dir, "vst_transform.R")
    r_script = f"""
    library(DESeq2)
    df <- read.csv("{input_csv}", row.names=1)
    # Convert to integer matrix as required by DESeq2
    df <- round(as.matrix(df))
    # Create DESeq2 dataset
    dds <- DESeqDataSetFromMatrix(countData = df, colData = data.frame(row.names=colnames(df)), design = ~ 1)
    # Apply VST transformation
    vsd <- varianceStabilizingTransformation(dds)
    # Save results
    write.csv(assay(vsd), "{output_csv}")
    """
    
    with open(r_script_path, "w") as f:
        f.write(r_script)
    
    # Execute R script
    try:
        r_command = "Rscript" if r_path is None else r_path
        cmd = [r_command, r_script_path]
        
        # Run the subprocess with full error output
        process = subprocess.run(
            cmd, 
            check=False,  # Don't raise exception yet
            capture_output=True,
            text=True
        )
        
        if process.returncode != 0:
            print(f"R Error Output:\n{process.stderr}")
            raise RuntimeError(f"R script failed with exit code {process.returncode}")
            
        print(f"VST-transformed data saved to {output_csv}")
        
        # Read and return the transformed data
        return pd.read_csv(output_csv, index_col=0)
        
    except FileNotFoundError:
        print("\nERROR: Rscript executable not found! Please provide the full path to Rscript.exe")
        print("Example usage:")
        print('vst_transform_r(df, r_path="C:/Program Files/R/R-4.2.2/bin/Rscript.exe")')
        print("\nAlternative Python-only method:")
        print("Consider using the Python package 'skbio' for similar transformations:")
        print("from skbio.stats.composition import clr")
        print("transformed_data = clr(df + 1)  # CLR transform with pseudocount")
        return None
    finally:
        # Clean up temporary files
        try:
            os.remove(input_csv)
            os.remove(r_script_path)
            os.rmdir(temp_dir)
        except:
            pass  # Ignore cleanup errors

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
