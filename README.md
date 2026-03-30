# Rosetta-Test-Implementation

This repository contains my implementation of the **Easy**, **Medium**, and **Hard** qualification tests for the **Rosetta** project. 

## Overview
The goal of these tests is to demonstrate the ability to bridge the gap between **Python (Pandas)** and **R (Bioconductor/DESeq2)** using `rpy2`, while ensuring a robust and reproducible development environment.

## Environment & Reproducibility
To resolve system-level dependencies and shared library linking issues (`libR.so`), I developed a **Dockerized workflow** based on Ubuntu 22.04.

### How to Run
1. **Build the Docker Image:**
   ```bash
   docker build -t rosetta-test .
2. **Execute Hard Test (Unit Testing):**
   ```bash
   docker run -it -v $(pwd):/app rosetta-test pytest hard_test.py
3. **Execute Medium Test (Biological Validation):**
   ```bash
   docker run -it -v $(pwd):/app rosetta-test python3 medium_test.py

## Technical Challenges Resolved
During the implementation, I identified and resolved several critical hurdles:

Pandas 2.x Compatibility: Implemented a monkey patch for pd.DataFrame.iteritems to support modern Python data stacks where this method is deprecated.

S4 Object Interoperability: Managed complex DESeq2 RS4 classes by using robjects.globalenv and explicit type coercion (as.vector) to satisfy R method signatures for dispersions<-.

Edge Case Fitting: Customized the DESeq2 workflow (estimateDispersionsGeneEst) to handle low-variance simulated datasets, ensuring 100% test coverage even in extreme edge cases.

## Results
Easy Test: Basic R environment verification (Success).

Medium Test: Successfully processed the airway dataset and matched Bioconductor's statistical output in Python.

Hard Test: Implemented defensive programming with input validation and achieved 2 passed in pytest.
