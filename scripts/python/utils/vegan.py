import rpy2.robjects as robjects
from rpy2.robjects.packages import importr
from rpy2.robjects import pandas2ri
from rpy2.robjects.conversion import localconverter
import pandas as pd
import numpy as np

# Activate pandas conversion
pandas2ri.activate()

# Import the vegan package
vegan = importr('vegan')

# Custom wrapper for rda function
def rda(formula, data, scale=False):
    """
    Wrapper around vegan.rda to handle pandas DataFrame conversion.
    """
    # Convert formula to R format
    r_formula = robjects.Formula(formula)
    
    # Convert data to R dataframe
    with localconverter(robjects.default_converter + pandas2ri.converter):
        r_data = pandas2ri.py2rpy(data)
    
    # Call R's rda function
    result = vegan.rda(r_formula, data=r_data, scale=scale)
    
    return result

# Custom wrapper for varpart function
def varpart(y, *args):
    """
    Wrapper around vegan.varpart to handle pandas DataFrame conversion.
    """
    # Convert y to R dataframe
    with localconverter(robjects.default_converter + pandas2ri.converter):
        r_y = pandas2ri.py2rpy(y)
    
    # Convert each X matrix to R
    r_args = []
    for arg in args:
        with localconverter(robjects.default_converter + pandas2ri.converter):
            r_arg = pandas2ri.py2rpy(arg)
        r_args.append(r_arg)
    
    # Call varpart with appropriate number of arguments
    if len(r_args) == 2:
        result = vegan.varpart(r_y, r_args[0], r_args[1])
    elif len(r_args) == 3:
        result = vegan.varpart(r_y, r_args[0], r_args[1], r_args[2])
    elif len(r_args) == 4:
        result = vegan.varpart(r_y, r_args[0], r_args[1], r_args[2], r_args[3])
    else:
        raise ValueError("varpart supports 2-4 explanatory tables")
    
    return result

# Custom wrapper for adonis2 function
def adonis2(formula, data, permutations=999, method="bray", by="terms", parallel=1):
    """
    Wrapper around vegan.adonis2 to handle pandas DataFrame conversion.
    """
    # Convert formula to R format
    r_formula = robjects.Formula(formula)
    
    # Convert data to R dataframe
    with localconverter(robjects.default_converter + pandas2ri.converter):
        r_data = pandas2ri.py2rpy(data)
    
    # Call adonis2
    result = vegan.adonis2(r_formula, data=r_data, permutations=permutations, 
                          method=method, by=by, parallel=parallel)
    
    # Convert result to pandas DataFrame
    with localconverter(robjects.default_converter + pandas2ri.converter):
        pd_result = pandas2ri.rpy2py(result)
    
    return pd_result