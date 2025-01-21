################################################################################
# Application example
################################################################################

# setup - load data and functions
##################################

# Load DARth ABC functions
source("scripts/DARthABC.R") # loads rcarbon package

# Load data from rcarbon package for an example (substitute this with the loading of you own data)
data(euroevol)
# get data in format expected by functions in DARthABC.R
target_data = list(CRA = euroevol$C14Age, error = euroevol$C14SD)

# Define time range of interest (YBP) (substitute here for the range of your interest)
time_range = c(8000, 4000)

# Calculates summary statistics on CRAs with function from DARthABC.R
target_sumstat = get_sumstats(target_data, time_range)

# Generate reference table for exponential model for ABC analysis
##################################################################
# for other models explore DARthABC.R to find the appropiate functions

# size of the reference table. as a rule of thumb 10000 simulations should be enough
num_of_sims = 10000 

# load packages necessary for running simulations in parallel (or rewrite the script to run them sequentially):
require(doSNOW)
require(doParallel)
require(doRNG)

# set up the cluster,
# change number of cores accordingly to your computing capacity (set to 5 for safety)
# change configuration of cluster according to your operating system if necessary
ncores = 5
cl <- makeCluster(ncores, type = "FORK")  
registerDoParallel(cl)  
registerDoRNG(seed = 33333)

# foreach loop to run the simulations
reftable <- foreach(sim = seq_len(num_of_sims), .combine = rbind) %dopar% {

  # sample parameters from prior (change priors to those appropriate for your data if necessary)
  rate = runif(1, -0.001, 0.001) # exponential rate from an uniform distribution
  lambda_0 = 10^runif(1, log10(0.001), log10(1)) # initial lambda value from a log-uniform distribution
  params = cbind(lambda_0, rate)
  
  # calculates parameter lambda for the whole time range
  lambda_t = get_exponential_lambda_t(lambda_0, rate, time_range)
  
  # simulates the true ages of the data set
  dates = sim_dates_lambda(lambda_t, time_range)
  
  # simulates the radiocarbon ages (CRA) of the data set
  dates_CRA = sim_CRA(dates, errors = target_data$error)
  
  # calculates the summary statistics for the simulated data
  sumstats = get_sumstats(dates_CRA, time_range)
  
  cbind(params, sumstats)
}
stopCluster(cl) 

# Approximate Bayesian Computation
###################################

# load abcrf pacjages for ABC via random forest
require(abcrf)

# separate reference table in parameters and summary statistics
rate = reftable$rate
lambda_0 = reftable$lambda_0
sumstats = reftable[names(target_sumstat)]

# train random forest to predict rate parameter 
RFmodel = regAbcrf(rate~., data.frame(rate,sumstats),
                   ntree = 1000, paral = TRUE)

# calculate posterior probability distribution fro observed (target) data
results_p = predict(RFmodel,
                    obs=target_sumstat,
                    training=data.frame(rate,sumstats),
                    paral=T)

cat("Point estimate of growth rate is", results_p$med, "with 95% credibility interval (",results_p$quantile[1],",",results_p$quantile[2],")")

# repeat ABCRF steps for lambda_0