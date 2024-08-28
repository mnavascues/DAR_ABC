## DARth ABC: Analysis of the Dated Archaeological Record through Approximate Bayesian Computation

This repository contains a set of scripts in R that implement the method described in "Analysis of the Abundance of Radiocarbon Samples without the Sum of Probability Distributions" by Miguel de Navascués, Concetta Burgarella and Mattias Jakobsson (2024) [doi:10.5281/zenodo.13381596](https://doi.org/10.5281/zenodo.13381596).

### Scripts

#### Application file: Example.R

File that contains an example application. Start here (after reading the publication) if you intend to apply the method to your data.

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

DARthABC:

- rcarbon
- weights
- extraDistr

Scripts that apply the moethod to data:

- abcrf
- doSNOW
- doParallel
- doRNG