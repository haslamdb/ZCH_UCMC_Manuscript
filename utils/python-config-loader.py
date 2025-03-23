#!/usr/bin/env python
# utils/config.py - Configuration loader for Python scripts

import os
import yaml
import re

def resolve_path(path, variables):
    """Resolve variable references in paths like ${variable_name}"""
    if not isinstance(path, str):
        return path
        
    # Replace ${variable} with its value
    pattern = r'\${([^}]*)}'
    
    def replace(match):
        var_name = match.group(1)
        if var_name in variables:
            return variables[var_name]
        return match.group(0)
    
    return re.sub(pattern, replace, path)

def expand_paths(config):
    """Recursively expand all path variables in the config"""
    if isinstance(config, dict):
        # First, collect all potential path variables
        path_vars = {}
        
        # Start with paths section
        if 'paths' in config:
            for key, value in config['paths'].items():
                path_vars[key] = resolve_path(value, path_vars)
            
        # Then process all other sections recursively
        for key, value in config.items():
            if isinstance(value, (dict, list)):
                config[key] = expand_paths(value)
            elif isinstance(value, str):
                config[key] = resolve_path(value, path_vars)
                
    elif isinstance(config, list):
        for i, item in enumerate(config):
            if isinstance(item, (dict, list, str)):
                config[i] = expand_paths(item)
                
    return config

def load_config(config_path="config.yaml"):
    """
    Load configuration from YAML file and resolve path variables
    
    Parameters:
    -----------
    config_path : str
        Path to the YAML configuration file
        
    Returns:
    --------
    dict
        Configuration dictionary with resolved paths
    """
    # Get the directory of the script being run
    script_dir = os.path.dirname(os.path.abspath(__file__))
    
    # If config_path is not absolute, make it relative to the script directory
    if not os.path.isabs(config_path):
        # Go one level up from utils directory
        project_dir = os.path.dirname(script_dir)
        config_path = os.path.join(project_dir, config_path)
    
    # Load the YAML file
    with open(config_path, 'r') as file:
        config = yaml.safe_load(file)
    
    # Expand path variables
    config = expand_paths(config)
    
    # Create directories if they don't exist
    if 'paths' in config:
        for path_name, path_value in config['paths'].items():
            if path_name.endswith('_dir') and path_value != ".":
                os.makedirs(path_value, exist_ok=True)
    
    return config

if __name__ == "__main__":
    # Example usage
    config = load_config()
    print("Loaded configuration:")
    print(f"Project directory: {config['paths']['project_dir']}")
    print(f"Results directory: {config['paths']['results_dir']}")
    print(f"Key organisms: {', '.join(config['key_organisms'])}")
