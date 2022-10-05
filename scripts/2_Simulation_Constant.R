library(rcarbon)
library(weights)
library(doSNOW)
library(doParallel)
library(doRNG)
source("scripts/sim14c.R")

set.seed(24)
ncores = 30

load(file = "results/Bevan_time_range_BP.rda")
load(file = "results/Bevan_num_of_sims.rda")
load(file = "results/Bevan_lambda_prior.rda")

# make a reference table constant model
if (!file.exists("results/Bevan_constant_model_reftable.rda") ){
  # setup parallel computing
  cl <- makeCluster(ncores, type="FORK")  
  registerDoParallel(cl)  
  registerDoRNG(seed = 1234567)
  reftable <- foreach(sim=seq_len(num_of_sims), .combine=rbind) %dopar% {
    gc()
    lambda = exp(runif(1,log(lambda_min),log(lambda_max)))
    demography = rep(lambda, length(t_))
    ss = sim_all(demography, t_, SPD=F, errors=dates$Error, runm=100, window=100)
    params = as.data.frame(lambda)
    cbind(params,ss)
  }
  stopCluster(cl) 
  save(reftable,file="results/Bevan_constant_model_reftable.rda")
}


