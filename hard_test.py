"""
Rosetta Hard Test - Implementation of deseq2 wrapper via rpy2
Contributor: Catherine Chi Chung (GSoC 2026)
Environment: R 4.4.2, Python 3.10.12

The script validates the Tier 1 API prototype and makes sure robust
data exchange between Pandas and Bioconductor S4 objects
"""

import os
import pandas as pd
import numpy as np
import pytest

# Important: Patch for rpy2 compatibility with Pandas 2.0+
pd.DataFrame.iteritems = pd.DataFrame.items

# Hardening R environment pathing
os.environ['R_HOME'] = '/usr/lib/R'

import rpy2.robjects as robjects
from rpy2.robjects.packages import importr
from rpy2.robjects import pandas2ri
from rpy2.robjects.conversion import localconverter

def deseq2(counts, metadata, design):
    """
    A hardened Python wrapper for the DESeq2 pipeline
    Implements pre-check validation to prevent R-level segmentation faults
    """
    # Pre-validation: Catch common user errors before hitting R C-API
    if not isinstance(counts, pd.DataFrame):
        raise TypeError("counts must be a pandas DataFrame")
    
    # 1. Non-negative integer check
    if (counts < 0).any().any():
        raise ValueError("Counts cannot contain negative values")
    if not np.array_equal(counts, counts.astype(int)):
        raise ValueError("Counts must be integers")

    # 2. Index consistency check
    if not all(counts.columns == metadata.index):
        raise ValueError("Sample names (counts columns) must match metadata index")
    
    # 3. Design formula validation
    design_factor = design.replace('~', '').strip()
    if design_factor not in metadata.columns:
        raise ValueError(f"Design factor '{design_factor}' not found in metadata columns")

    # Initialize Bioconductor packages
    deseq2_pkg = importr('DESeq2')
    stats = importr('stats')
    
    with localconverter(robjects.default_converter + pandas2ri.converter):
        r_counts = robjects.conversion.py2rpy(counts)
        r_meta = robjects.conversion.py2rpy(metadata)

        # Build DESeqDataSet
        dds = deseq2_pkg.DESeqDataSetFromMatrix(
            countData=r_counts, 
            colData=r_meta, 
            design=robjects.Formula(design)
        )
        
        # Inject dds into R global env for complex S4 manipulation
        robjects.globalenv['dds'] = dds

        robjects.r('''
            dds <- estimateSizeFactors(dds)
            dds <- estimateDispersionsGeneEst(dds)
                   
            # WORKAROUND: Explicit R-side type coercion to 'as.vector'
            # Resolves a 'Signature Mismatch' error in R's S4 method dispatch for dispersions()<-
            dispersions(dds) <- as.vector(mcols(dds)$dispGeneEst)
                   
            dds <- nbinomWaldTest(dds)
            res_obj <- results(dds)
        ''')
        
        # Safe data frame conversion back to Python
        res_r = robjects.r('as.data.frame(res_obj)')
        return robjects.conversion.rpy2py(res_r)

# --- Test Funcs ---

def test_success_workflow():
    """Validates the standard DESeq2 workflow with simulated counts"""
    counts = pd.DataFrame({
        's1': [10, 20, 30], 's2': [11, 21, 31],
        's3': [100, 200, 300], 's4': [105, 205, 305]
    }, index=['G1', 'G2', 'G3'])
    meta = pd.DataFrame({'cond': ['A', 'A', 'B', 'B']}, index=['s1', 's2', 's3', 's4'])
    
    result = deseq2(counts, meta, "~ cond")
    assert isinstance(result, pd.DataFrame)
    assert 'log2FoldChange' in result.columns
    assert 'padj' in result.columns

def test_validation_negative_counts():
    """Ensures the wrapper catches invalid data before R crashes"""
    counts = pd.DataFrame({'s1': [-5, 10], 's2': [5, 10]}, index=['G1', 'G2'])
    meta = pd.DataFrame({'cond': ['A', 'B']}, index=['s1', 's2'])
    with pytest.raises(ValueError, match="Counts cannot contain negative values"):
        deseq2(counts, meta, "~ cond")

def test_validation_mismatched_samples():
    """Checks for sample-metadata index alignment"""
    counts = pd.DataFrame({'s1': [10, 20], 's2': [11, 21]}, index=['G1', 'G2'])
    meta = pd.DataFrame({'cond': ['A', 'B']}, index=['sample_X', 'sample_Y'])
    with pytest.raises(ValueError, match="Sample names.*must match"):
        deseq2(counts, meta, "~ cond")

def test_validation_float_counts():
    """Make sure the script can stop non integer input, to avoid wrong calculate in R"""
    counts = pd.DataFrame({'s1': [10.5, 20], 's2': [11, 21]}, index=['G1', 'G2'])
    meta = pd.DataFrame({'cond': ['A', 'B']}, index=['s1', 's2'])
    with pytest.raises(ValueError, match="Counts must be integers"):
        deseq2(counts, meta, "~ cond")

def test_validation_missing_design_column():
    """Test for not existing column name"""
    counts = pd.DataFrame({'s1': [10, 20], 's2': [11, 21]}, index=['G1', 'G2'])
    meta = pd.DataFrame({'cond': ['A', 'B']}, index=['s1', 's2'])
    # provide a not existing column name 'treatment'
    with pytest.raises(ValueError, match="not found in metadata columns"):
        deseq2(counts, meta, "~ treatment")