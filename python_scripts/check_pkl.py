#!/usr/bin/env python3

import pickle
import os

# Load the file
with open('results/variance_analysis/variance_analysis_results.pkl', 'rb') as f:
    data = pickle.load(f)

# Print the type and structure
print(f"Type: {type(data)}")

if isinstance(data, dict):
    print("\nKeys:")
    for key in data.keys():
        print(f"- {key}")
        if isinstance(data[key], dict):
            print(f"  Subkeys in {key}:")
            for subkey in data[key].keys():
                print(f"  - {subkey}")
else:
    print("Not a dictionary")