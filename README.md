## DARth ABC: Analysis of the Dated Archaeological Record through Approximate Bayesian Computation

This repository contains a set of scripts in R used in the work "Analysis of the Abundance of Radiocarbon Samples without the Sum of Probability Distributions" by Miguel de Navascués, Concetta Burgarella and Mattias Jakobsson ([doi:10.5281/zenodo.13381596](https://zenodo.org/doi/10.5281/zenodo.13381596)).

### Scripts

- *DARthABC.R* contains the functions to simulate and calculate summary statistics for ABC
- *MethodEvaluation.R* is a script that analyses a toy model with different models (probability distribution model vs. poisson model) and summary statistics (based on SPD vs. based on CRA dates)


### Requirements

R packages:
- rcarbon
- abcrf
- weights
- doSNOW
- doParallel
- doRNG
- hexbin
- latticeExtra