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

def vst_transform_r(df_or_path, output_csv=None):
    """
    Perform Variance-Stabilizing Transformation (VST) using DESeq2 in R.
    
    Parameters:
    - df_or_path: Pandas DataFrame or path to input CSV
    - output_csv: Path to output CSV file (optional)
    
    Returns:
    - DataFrame with VST-transformed data or None if error
    """
    import os
    import pandas as pd
    import tempfile
    import subprocess
    
    # Create temporary directory and files with normalized paths
    temp_dir = tempfile.mkdtemp()
    
    # Handle input - either DataFrame or path
    is_dataframe = isinstance(df_or_path, pd.DataFrame)
    if is_dataframe:
        # If DataFrame, save to temp CSV
        temp_input = os.path.join(temp_dir, "temp_input.csv")
        df_or_path.to_csv(temp_input)
        input_path = temp_input
    else:
        # If string path, use as is
        input_path = df_or_path
    
    # Set output path
    if output_csv is None:
        temp_output = os.path.join(temp_dir, "temp_output.csv")
    else:
        temp_output = output_csv
    
    # Normalize paths to use forward slashes for R
    input_path_r = input_path.replace("\\", "/")
    output_path_r = temp_output.replace("\\", "/")
    
    # Create R script file
    r_script_path = os.path.join(temp_dir, "vst_transform.R")
    
    r_script = f"""
    # Load required libraries
    library(DESeq2)
    
    # Read the CSV file
    df <- try(read.csv("{input_path_r}", row.names=1), silent=TRUE)
    if(inherits(df, "try-error")) {{
        cat("ERROR: Failed to read input CSV file\\n")
        quit(status=1)
    }}
    
    # Convert to integer matrix for DESeq2
    df <- round(as.matrix(df))
    
    # Check for all-zero rows or columns
    if(any(rowSums(df) == 0)) {{
        cat("WARNING: Some rows have all zeros - removing them\\n")
        df <- df[rowSums(df) > 0, ]
    }}
    
    if(any(colSums(df) == 0)) {{
        cat("WARNING: Some columns have all zeros - removing them\\n")
        df <- df[, colSums(df) > 0]
    }}
    
    # Create DESeq2 dataset with minimal design
    sample_info <- data.frame(row.names=colnames(df), condition=rep(1, ncol(df)))
    
    # Create DESeq2 dataset
    dds <- try(DESeqDataSetFromMatrix(
        countData = df,
        colData = sample_info,
        design = ~ 1
    ), silent=TRUE)
    
    if(inherits(dds, "try-error")) {{
        cat("ERROR: Failed to create DESeq2 dataset\\n")
        cat("Possible issues:\\n")
        cat("- Data may not be counts (should be integers)\\n")
        cat("- Some rows or columns might be all zeros\\n")
        quit(status=1)
    }}
    
    # Apply VST transformation
    vsd <- try(varianceStabilizingTransformation(dds), silent=TRUE)
    
    if(inherits(vsd, "try-error")) {{
        cat("ERROR: VST transformation failed\\n")
        quit(status=1)
    }}
    
    # Write transformed data to output CSV
    transformed_data <- assay(vsd)
    try(write.csv(transformed_data, "{output_path_r}"), silent=TRUE)
    
    if(!file.exists("{output_path_r}")) {{
        cat("ERROR: Failed to write output file\\n")
        quit(status=1)
    }}
    
    cat("SUCCESS: VST transformation completed\\n")
    """
    
    # Write the R script to file with UTF-8 encoding
    with open(r_script_path, "w", encoding="utf-8") as f:
        f.write(r_script)
    
    # Execute R script
    try:
        # Run the R script
        result = subprocess.run(
            ["Rscript", r_script_path],
            check=False,
            capture_output=True,
            text=True
        )
        
        # Display R output
        if result.stdout:
            print("R Output:")
            print(result.stdout)
            
        # Handle errors
        if result.returncode != 0:
            print(f"R Error Output:")
            print(result.stderr)
            print("VST transformation failed. See R error output above.")
            return None
            
        # Read the transformed data
        if os.path.exists(temp_output):
            transformed_df = pd.read_csv(temp_output, index_col=0)
            print("VST transformation completed successfully.")
            return transformed_df
        else:
            print("Error: Output file not created.")
            return None
            
    except FileNotFoundError:
        print("\nERROR: Rscript executable not found.")
        print("Make sure R and DESeq2 are installed, and Rscript is in your PATH.")
        print("\nTry these steps:")
        print("1. Install R from https://cran.r-project.org/")
        print("2. In R, install DESeq2: BiocManager::install('DESeq2')")
        print("3. Add R bin directory to your PATH environment variable")
        return None
    finally:
        # Clean up temporary files
        try:
            if is_dataframe and os.path.exists(temp_input):
                os.remove(temp_input)
            if output_csv is None and os.path.exists(temp_output):
                os.remove(temp_output)
            if os.path.exists(r_script_path):
                os.remove(r_script_path)
            os.rmdir(temp_dir)
        except:
            pass  # Ignore cleanup errors

# run this if vst_transform_r fails
# run in command prompt rather than windows powershell (which can't find R script)
def debug_vst_transform(df, debug_mode=True):
    """
    Debug-friendly VST transformation with alternative methods when DESeq2 fails.
    
    Parameters:
    - df: Pandas DataFrame with microbiome counts
    - debug_mode: Print detailed debugging information
    
    Returns:
    - DataFrame with transformed data or None if all methods fail
    """
    import os
    import pandas as pd
    import numpy as np
    import tempfile
    import subprocess
    import time
    
    # Create a unique name for temporary files
    timestamp = int(time.time())
    temp_dir = tempfile.mkdtemp()
    input_csv = os.path.join(temp_dir, f"input_{timestamp}.csv")
    output_csv = os.path.join(temp_dir, f"output_{timestamp}.csv")
    
    # Print data summary if in debug mode
    if debug_mode:
        print("\n--- Data Summary ---")
        print(f"DataFrame shape: {df.shape}")
        print(f"Non-zero values: {(df > 0).sum().sum()} / {df.size} ({(df > 0).sum().sum() / df.size:.2%})")
        print(f"Zero values: {(df == 0).sum().sum()} / {df.size} ({(df == 0).sum().sum() / df.size:.2%})")
        print(f"Min value: {df.values.min()}")
        print(f"Max value: {df.values.max()}")
        print(f"Mean value: {df.values.mean():.2f}")
        print(f"Samples with all zeros: {(df.sum(axis=1) == 0).sum()}")
        print(f"Features with all zeros: {(df.sum(axis=0) == 0).sum()}")
        print(f"Data type: {df.dtypes.iloc[0]}")
    
    # Save to CSV file
    df.to_csv(input_csv)
    
    # Normalize paths for R
    input_r_path = input_csv.replace("\\", "/")
    output_r_path = output_csv.replace("\\", "/")
    
    # Create a more robust R script with extensive error checking
    r_script_path = os.path.join(temp_dir, f"vst_script_{timestamp}.R")
    
    r_script = f"""
    # Set error handling to get more information
    options(error = function() {{ 
      traceback(3)
      if(!interactive()) quit(status = 1, save = "no")
    }})
    
    # Load DESeq2
    suppressPackageStartupMessages(library(DESeq2))
    suppressPackageStartupMessages(library(SummarizedExperiment))
    
    # Read input data
    cat("Reading input data from:", "{input_r_path}", "\\n")
    counts <- read.csv("{input_r_path}", row.names=1, check.names=FALSE)
    
    # Print data summary
    cat("\\nData summary:\\n")
    cat("Dimensions:", nrow(counts), "x", ncol(counts), "\\n")
    cat("Non-zero count:", sum(counts > 0), "/", nrow(counts) * ncol(counts), "\\n")
    cat("Min value:", min(counts), "\\n")
    cat("Max value:", max(counts), "\\n")
    
    # Check for and remove zero-sum rows
    zero_rows <- rowSums(counts) == 0
    if(sum(zero_rows) > 0) {{
        cat("Removing", sum(zero_rows), "rows with all zeros\\n")
        counts <- counts[!zero_rows, ]
    }}
    
    # Check for and remove zero-sum columns
    zero_cols <- colSums(counts) == 0
    if(sum(zero_cols) > 0) {{
        cat("Removing", sum(zero_cols), "columns with all zeros\\n")
        counts <- counts[, !zero_cols]
    }}
    
    # Enforce integer mode
    counts_int <- round(as.matrix(counts))
    
    # Create a simple design with one condition
    colData <- data.frame(row.names=colnames(counts_int), condition=factor(rep("A", ncol(counts_int))))
    
    # Alternative 1: Try standard DESeq2 approach
    tryCatch({{
        cat("\\nMethod 1: Standard VST\\n")
        # Create DESeq dataset
        dds <- DESeqDataSetFromMatrix(countData = counts_int, 
                                     colData = colData,
                                     design = ~ 1)
        
        # Try variance stabilizing transformation
        vsd <- varianceStabilizingTransformation(dds)
        
        # Get transformed values
        transformed <- assay(vsd)
        
        # Write to output
        write.csv(transformed, "{output_r_path}")
        cat("Success: Standard VST transformation completed\\n")
    }}, error = function(e) {{
        cat("Method 1 failed:", conditionMessage(e), "\\n")
        
        # Alternative 2: Try with fitType="local"
        tryCatch({{
            cat("\\nMethod 2: VST with fitType='local'\\n")
            dds <- DESeqDataSetFromMatrix(countData = counts_int, 
                                         colData = colData,
                                         design = ~ 1)
            
            # Try variance stabilizing transformation with local fit
            vsd <- varianceStabilizingTransformation(dds, fitType="local")
            transformed <- assay(vsd)
            write.csv(transformed, "{output_r_path}")
            cat("Success: VST transformation with local fit completed\\n")
        }}, error = function(e) {{
            cat("Method 2 failed:", conditionMessage(e), "\\n")
            
            # Alternative 3: Filter low-expression genes more aggressively
            tryCatch({{
                cat("\\nMethod 3: More aggressive filtering + VST\\n")
                
                # More aggressive filtering (remove rows with low counts)
                keep <- rowSums(counts_int >= 10) >= 3
                if(sum(keep) < 2) {{
                    # If too strict, try lower threshold
                    keep <- rowSums(counts_int >= 5) >= 2
                }}
                
                cat("Keeping", sum(keep), "rows after filtering\\n")
                if(sum(keep) >= 2) {{
                    # Only proceed if we have at least 2 rows
                    counts_filtered <- counts_int[keep, ]
                    
                    dds <- DESeqDataSetFromMatrix(countData = counts_filtered, 
                                                 colData = colData,
                                                 design = ~ 1)
                    
                    vsd <- varianceStabilizingTransformation(dds)
                    transformed <- assay(vsd)
                    write.csv(transformed, "{output_r_path}")
                    cat("Success: VST transformation with filtering completed\\n")
                }} else {{
                    cat("Too few rows after filtering, trying alternative transformation\\n")
                    # Fall through to alternative transformation
                }}
            }}, error = function(e) {{
                cat("Method 3 failed:", conditionMessage(e), "\\n")
                
                # Alternative 4: Use rlog instead of VST
                tryCatch({{
                    cat("\\nMethod 4: Using rlog instead of VST\\n")
                    dds <- DESeqDataSetFromMatrix(countData = counts_int, 
                                                 colData = colData,
                                                 design = ~ 1)
                    
                    # Try rlog transformation instead
                    rld <- rlog(dds, blind=TRUE)
                    transformed <- assay(rld)
                    write.csv(transformed, "{output_r_path}")
                    cat("Success: rlog transformation completed\\n")
                }}, error = function(e) {{
                    cat("Method 4 failed:", conditionMessage(e), "\\n")
                    
                    # Alternative 5: Fallback to simpler log2 transformation
                    tryCatch({{
                        cat("\\nMethod 5: Simple log2 transformation\\n")
                        # Simple log2 transformation with pseudocount
                        log2_counts <- log2(counts + 1)
                        write.csv(log2_counts, "{output_r_path}")
                        cat("Success: Simple log2 transformation completed\\n")
                    }}, error = function(e) {{
                        cat("All transformation methods failed.\\n")
                        cat("Final error:", conditionMessage(e), "\\n")
                    }})
                }})
            }})
        }})
    }})
    
    # Final check - did we create output?
    if(file.exists("{output_r_path}")) {{
        cat("\\nTransformation complete - output saved to:", "{output_r_path}", "\\n")
    }} else {{
        cat("\\nFailed to create output file.\\n")
    }}
    """
    
    # Write the R script
    with open(r_script_path, "w", encoding="utf-8") as f:
        f.write(r_script)
    
    # Run the R script
    try:
        if debug_mode:
            print("\n--- Running R Script ---")
        
        result = subprocess.run(
            ["Rscript", r_script_path],
            check=False,
            capture_output=True,
            text=True
        )
        
        # Print R output for debugging
        if result.stdout:
            print("\n--- R Output ---")
            print(result.stdout)
        
        if result.stderr and debug_mode:
            print("\n--- R Error Output ---")
            print(result.stderr)
        
        # Check if output file was created
        if os.path.exists(output_csv):
            print("\n--- Transformation successful! ---")
            transformed_df = pd.read_csv(output_csv, index_col=0)
            
            if debug_mode:
                print(f"Transformed data shape: {transformed_df.shape}")
                print(f"Transformed min: {transformed_df.values.min():.4f}")
                print(f"Transformed max: {transformed_df.values.max():.4f}")
                print(f"Transformed mean: {transformed_df.values.mean():.4f}")
            
            return transformed_df
        else:
            print("\n--- Transformation failed: No output file created ---")
            
            # Try Python-based fallback method
            print("\n--- Trying Python fallback method ---")
            print("Using log(x+1) transformation as fallback")
            
            # Remove zero rows and columns to match R processing
            df_clean = df.loc[df.sum(axis=1) > 0, df.sum(axis=0) > 0]
            transformed = np.log1p(df_clean)
            
            print(f"Fallback transformation shape: {transformed.shape}")
            print(f"Fallback min: {transformed.values.min():.4f}")
            print(f"Fallback max: {transformed.values.max():.4f}")
            print(f"Fallback mean: {transformed.values.mean():.4f}")
            
            return pd.DataFrame(transformed, index=df_clean.index, columns=df_clean.columns)
            
    except FileNotFoundError:
        print("\nERROR: Rscript not found. Make sure R is installed and in your PATH.")
        
        # Python fallback if R isn't available
        print("Using Python log(x+1) transformation as fallback")
        return pd.DataFrame(
            np.log1p(df.loc[df.sum(axis=1) > 0, df.sum(axis=0) > 0]),
            index=df.loc[df.sum(axis=1) > 0].index,
            columns=df.loc[:, df.sum(axis=0) > 0].columns
        )
        
    finally:
        # Clean up temporary files
        try:
            if os.path.exists(input_csv):
                os.remove(input_csv)
            if os.path.exists(output_csv):
                os.remove(output_csv)
            if os.path.exists(r_script_path):
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
