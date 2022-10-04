library(rcarbon)
library(weights)
library(doSNOW)
library(doParallel)
library(doRNG)
source("scripts/sim14c.R")

set.seed(24)
ncores = 30

load(file = "results/Bevan_dates.rda")
load(file = "results/Bevan_time_range_BP.rda")
load(file = "results/Bevan_num_of_sims.rda")
load(file = "results/Bevan_lambda_prior.rda")

# make a reference table logistic expansion model
if (!file.exists("results/Bevan_logistic_model_reftable.rda") ){
  # setup parallel computing
  cl <- makeCluster(ncores, type="FORK")  
  registerDoParallel(cl)  
  registerDoRNG(seed = 1234567)
  reftable <- foreach(sim=seq_len(num_of_sims), .combine=rbind) %dopar% {
    gc()
    demograhy = get_logistic_model(lambda_min,lambda_max,time_range_BP,expansion=T)
    ss = sim_all(demograhy$lambda_t, t_, SPD=F, errors=dates$Error, runm=100, window=100)
    params = as.data.frame(t(c(demograhy$lambda_0,
                               demograhy$lambda_f,
                               demograhy$lambda_K,
                               demograhy$rate)))
    names(params) = c("lambda_0","lambda_f","K","r")
    cbind(params,ss)
  }
  stopCluster(cl) 
  save(reftable,file="results/Bevan_logistic_model_reftable.rda")
}
