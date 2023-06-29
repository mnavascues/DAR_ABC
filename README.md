## DARth ABC

This repository contains a set of scripts in R used in the work "A Model-Based Framework for the Analysis of the Dated Archaeological Record" by Miguel de Navascués and Concetta Burgarella.

### Scrpts

- *DARthABC.R* contains the functions to simulate and calculate summary statistics for ABC
- *MethodEvaluation.R* is a script that analyses a toy model with different models (probability distribution model vs. poisson model) and summary statistics (based on SPD vs. based on CRA dates)
- *Malmo.R* is a script for the analysis of archaeological data from Malmö area (ref to be added)


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