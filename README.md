## Overview
This repository contains the data and R scripts used for the analyses presented in the following study:

"Leaf ionomic data improve genomic prediction for soybean lines and their progeny under contrasting water conditions"

This study investigated the usefulness of leaf ionomic information for improving the prediction of phenotypic and genotypic values in soybean.

## Scripts
The R scripts are numbered using a two-part numbering system (X.Y), where the first number indicates the analytical group and the second number indicates the execution order within that group.

The scripts should be executed sequentially in numerical order, starting from 1.1 and proceeding through 4.9.

The scripts are organized into the following four groups:

1. Principal Component Analysis and Correlation Analyses

2. Prediction Accuracy Within Each Year–Treatment Combination

3. Prediction Accuracy Across Years and Water Treatments

4. Prediction Accuracy for the Progeny Population

## Requirements and Environment Setup
The analyses were performed using R (version 4.3.3).

The required R packages and their versions are recorded in renv.lock, which was generated using the renv package.

To reproduce the computational environment, install renv and run the following commands in the R console from the project directory:

install.packages("renv")
renv::restore()

This will install the package versions recorded in renv.lock into the project-specific library.

The .Rprofile file and the renv/ directory are required for activating the project-specific package environment.
