# Example of migrating an existing script to use the unified configuration
# This shows "before" and "after" versions of code sections

# --------------- BEFORE MIGRATION ---------------
# Hard-coded parameters scattered throughout the script:

# Hard-coded paths
DATA_DIR = "../data"
RESULTS_DIR = "../results"
KRAKEN_DIR = "../KrakenAlignments/Kraken2"

# Hard-coded file names
microbiome_df = pd.read_csv("../data/NICUSpeciesReduced.csv", index_col=0)
metadata_df = pd.read_csv("../metadata/AllNICUSampleKey20250206.csv", index_col=0)

# Hard-coded analysis parameters
RAREFACTION_DEPTH = 10000
NOISE_THRESHOLD = 0.001
MIN_PREVALENCE = 0.05

# Hard-coded target organisms
key_organisms = ["Klebsiella.pneumoniae", "Staphylococcus.aureus", "Escherichia.coli", 
                "Klebsi