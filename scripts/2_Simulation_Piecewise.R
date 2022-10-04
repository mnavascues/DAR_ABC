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

# max_num_of_periods = 1000
# skyline_years = seq(time_range_BP[1], time_range_BP[2], by = -100)

num_of_periods = round((time_range_BP[1] - time_range_BP[2]) / 100 / 4)
skyline_years = get_time_of_change(num_of_periods, time_range_BP, 
                                   intervals = "regular")


# make reference tables for piecewise exponential models
if (!file.exists("results/Bevan_piecewise_model_reftable.rda") ){
  # setup parallel computing
  cl <- makeCluster(ncores, type="FORK")  
  registerDoParallel(cl)  
  registerDoRNG(seed = sample(1:1000000,1) )
  reftable <- foreach(sim=seq_len(num_of_sims), .combine=rbind) %dopar% {
    gc()
    # num_of_periods = round(exp(runif(1,log(2),log(max_num_of_periods))))
    demograhy = get_piecewise_exponential_model(num_of_periods,
                                                lambda_min,
                                                lambda_max,
                                                time_range_BP,
                                                intervals = "regular",
                                                skyline = F)
    ss = sim_all(demograhy$lambda_t, 
                 t_, 
                 SPD = F, 
                 errors = dates$Error, 
                 runm = 100, 
                 window = 100)
    params = as.data.frame(t(c(demograhy$lambda_skyline,
                               demograhy$rate_skyline,
                               num_of_periods)))
    names(params) = c(paste0("lambda", skyline_years),
                      paste0("rate", skyline_years[1:(length(skyline_years)-1)]+(skyline_years[1]-skyline_years[2])/2),
                      "num_of_periods")
    #names(params) = c(paste0("lambda", skyline_years),
    #                  paste0("rate", skyline_years),
    #                  "num_of_periods")
    cbind(params,ss)
  }
  save(reftable, file="results/Bevan_piecewise_model_reftable.rda")
  stopCluster(cl) 
}




