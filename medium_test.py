import os
import sys

# Set R environment variables
os.environ['R_HOME'] = '/usr/lib/R'
os.environ['LD_LIBRARY_PATH'] = '/usr/local/lib/R/lib:' + os.environ.get('LD_LIBRARY_PATH', '')

try:
    import rpy2.robjects as robjects
    from rpy2.robjects.packages import importr
    from rpy2.robjects import pandas2ri
except ImportError as e:
    print(f"Failed to import rpy2: {e}")
    sys.exit(1)

import pandas as pd

# 1. Import R packages
base = importr('base')
deseq2 = importr('DESeq2')
airway = importr('airway')

print("--- Running DESeq2 analysis in Python via rpy2 ---")

# 2. Load airway data
robjects.r('data("airway")')
airway_obj = robjects.r['airway']

# 3. DESeq2 Pipeline
# formula: ~ dex
dds = deseq2.DESeqDataSet(airway_obj, design=robjects.Formula('~ dex'))
dds = deseq2.DESeq(dds)
res = deseq2.results(dds)

# 4. Data Transfer: Force R-side conversion to DF before pulling to Pandas
robjects.globalenv['res_r'] = res
res_df_r = robjects.r('as.data.frame(res_r)')

# Activate pandas conversion globally
pandas2ri.activate()
res_pandas = pandas2ri.rpy2py(res_df_r)

# Filter out NAs in padj
res_pandas = res_pandas.dropna(subset=['padj'])
res_pandas = res_pandas.sort_values('padj')

print("\n--- Python (rpy2) Success: Top 5 Differentially Expressed Genes ---")
print(top5_py := res_pandas.head(5))