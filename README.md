# Estimating the diagonal of the inverse
 
This repository contains the numerical experiments for the `Estimating the diagonal of the inverse` project, carried out in the context of the MATH-453 Computational Linear Algebra course at EPFL.

## Repository description
- `code` - Implementation of the experiments
  - `utils` - Subfolder with utilities functions
  - `comparison.m` - Main function for comparison of the three estimators
  - `lanczos_mc.m` - Main function for Lanczos MC estimator
  - `lanczos.m` - Main function for Lanczos estimator
  - `mc.m` - Main function for MC estimator
  - `main.m` - Main script to run all the experiments
- `figures` - Plots of the obtained results
- `matrices` - Folder containing the matrices used in the experiments in `.mat` format

## Reproducibility of the results
We provide a unique Matlab script `main.m` to reproduce the results shown in the report. We fixed the seed once and for all to produce the plots in the report, hence running such script will reproduce exactly those plots, in the same order as they appear in the document. Because we tested multiple matrices in our reports, a variable `matname` is present in `main.m`. To ease your job, this is by default set to `nos3`.

## Authors
- Federico Betti
- Theophile Boinnard
