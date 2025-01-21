## DARth ABC: Analysis of the Dated Archaeological Record through Approximate Bayesian Computation

This repository contains a set of scripts in R that implement the method described in "Analysis of the Abundance of Radiocarbon Samples without the Sum of Probability Distributions" by Miguel de Navascués, Concetta Burgarella and Mattias Jakobsson (2024) [doi:10.5281/zenodo.13381596](https://doi.org/10.5281/zenodo.13381596).

### Vignette

A vignette is included with a simple example. After reading the publication, use this as a starting point for undersanding how to apply the functions in your own analysis.

### Scripts

#### Main file: DARthABC.R

File that contains the functions needed for the analysis of the abundance of radiocarbon data.

#### Files to reproduce analysis presented in Navascués *et al.* (2024) [doi:10.5281/zenodo.13381596](https://doi.org/10.5281/zenodo.13381596)  

- *MethodEvaluation.R* analyses a toy model with different models (probability distribution model *vs.* Poisson (data as counts) model) and summary statistics (based on SPD *vs.* based on CRA dates).
- *MethodEvaluationPlots.R* plots the results from the analysis performed in *MethodEvaluation.R*.
- *Bevan_all.R* analyses the data from Bevan *et al.* (2017) [doi:10.1073/pnas.1709190114](https://doi.org/10.1073/pnas.1709190114). Requires to download the data set published in Bevan (2017) [doi:10.14324/000.ds.10025178](https://doi.org/10.14324/000.ds.10025178). 
- *Bevan_cereals.R* analyses a subset of those data that correponds to samples identified as barley or wheat.

#### Other scrpits

Other files correspond to work in progress and will be described in the future.

### Requirements (R packages)

These scripts have been tested in R 4.1.2. The following packages are required (with corresponding version tested):

DARthABC:

- rcarbon (1.5.1)
- weights (4.1.2)
- extraDistr (1.9.1)

Scripts that apply the method to data:

- abcrf (1.9)
- doSNOW (1.0.20)
- doParallel (1.0.17)
- doRNG (1.8.6)

