library(rcarbon)
library(weights)
library(doSNOW)
library(doParallel)
library(doRNG)
source("scripts/sim14c.R")

set.seed(24)
ncores = 30

load(file = "results/dates.rda")
load(file = "results/time_range_BP.rda")
load(file = "results/num_of_sims.rda")
load(file = "results/lambda_prior.rda")


gisp2 = read.table(file = "data/climate/gisp2_temp_accum_alley2000.txt", skip = 75, nrows = 1707 - 75)
names(gisp2) = c("YBP","Temperature")
gisp2$YBP = gisp2$YBP * 1000
dates_4_interpolation = seq(time_range_BP[1], time_range_BP[2] + 1,by = -1)
temperature = with(gisp2, data.frame(approx(YBP,Temperature,xout=dates_4_interpolation)))



# make reference tables for dynamic models
if (!file.exists("results/dynamic_model_reftable.rda") ){
  # setup parallel computing
  cl <- makeCluster(ncores, type="FORK")  
  registerDoParallel(cl)  
  registerDoRNG(seed = sample(1:1000000,1) )
  reftable <- foreach(sim=seq_len(num_of_sims), .combine=rbind) %dopar% {
    gc()
    demography = get_dynamic_model(lambda_min=0.001, lambda_max=1,
                                   b0_min=0.00001, b0_max=0.01,
                                   b1_min=0.00001, b1_max=0.01,
                                   b2_min=0.00001, b2_max=0.01,
                                   d=temperature$y, time_range=time_range_BP)
    #plot(dates_4_interpolation,demography$lambda_t, type="l",  xlim=time_range_BP)
    ss = sim_all(demography$lambda_t, 
                 t_, 
                 SPD = F, 
                 errors = dates$Error, 
                 runm = 100, 
                 window = 100)
    params = as.data.frame(t(c(demography$lambda_0,
                               demography$b0,
                               demography$b1,
                               demography$b2,
                               demography$d)))
    names(params) = c("lambda_0","b0","b1","b2","d")
    cbind(params,ss)
  }
  save(reftable, file="results/dynamic_model_reftable.rda")
  stopCluster(cl) 
}




