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

load(file="results/Bevan_m_hat.rda")


num_of_periods = round(m_hat)
skyline_years = get_time_of_change(num_of_periods, time_range_BP, intervals="regular")

# make reference tables for piecewise exponential models
if (!file.exists("results/Bevan_piecewise_fixed_m_model_reftable.rda") ){
  # setup parallel computing
  cl <- makeCluster(ncores, type="FORK")  
  registerDoParallel(cl)  
  registerDoRNG(seed = sample(1:1000000,1) )
  reftable <- foreach(sim=seq_len(num_of_sims), .combine=rbind) %dopar% {
    gc()
    demography = get_piecewise_exponential_model(num_of_periods, 
                                                 lambda_min, 
                                                 lambda_max, 
                                                 time_range_BP, 
                                                 intervals = "regular", 
                                                 skyline = F)
    ss = sim_all(demography$lambda_t, 
                 t_, 
                 SPD = F, 
                 errors = dates$Error, 
                 runm = 100, 
                 window = 100)
    params = as.data.frame(t(c(demography$lambda_skyline,
                               demography$rate_skyline,
                               num_of_periods)))
    names(params) = c(paste0("lambda", skyline_years),
                      paste0("rate", skyline_years[1:(length(skyline_years)-1)]+(skyline_years[1]-skyline_years[2])/2),
                      "num_of_periods")
    cbind(params,ss)
  }
  save(reftable, file="results/Bevan_piecewise_fixed_m_model_reftable.rda")
  stopCluster(cl) 
}




