library(rcarbon)
library(weights)
library(doSNOW)
library(doParallel)
library(doRNG)
source("scripts/sim14c.R")

set.seed(24)
ncores = 25

load(file = "results/Cereals_dates.rda")
load(file = "results/Cereals_time_range_BP.rda")
load(file = "results/num_of_sims.rda")
load(file = "results/Cereals_lambda_prior.rda")

num_of_periods = round((time_range_BP[1] - time_range_BP[2]) / 100 / 4)
skyline_years = get_time_of_change(num_of_periods, time_range_BP, 
                                   intervals = "regular")

# make reference tables for piecewise exponential models
if (!file.exists("results/Cereals_correlated_model_reftable.rda") ){
  # setup parallel computing
  cl <- makeCluster(ncores, type="FORK")  
  registerDoParallel(cl)  
  registerDoRNG(seed = sample(1:1000000,1) )
  reftable <- foreach(sim=seq_len(num_of_sims), .combine=rbind) %dopar% {
    gc()
 
    na_in_demography = T
    alpha = beta = 1
    pi = NA
    demography = get_piecewise_exponential_model_2_categories(num_of_periods,
                                                              lambda_min,
                                                              lambda_max,
                                                              model = "correlated",
                                                              alpha, beta,
                                                              pi = NULL,
                                                              time_range_BP)
    
    #plot(demography$lambda_t_B, type = "l", lwd = 2, ylim = c(0, lambda_max))
    #lines(demography$lambda_t_A, lwd = 2, col = "red")
    
    na_in_demography = any(c(is.na(demography$lambda_t_A),
                             is.na(demography$lambda_t_B)))
    if (na_in_demography){
      ss = rep(NA,598)
      ds = seq(time_range_BP[1],time_range_BP[2]-1,-100)[-1]+50
      probs = seq(0, 1, 0.02)
      names(ss) = c("count_A",  paste0("hist", ds, "_A"),
                    paste0("delta_h", ds[-1], "_A"), paste0("quantile", probs, "_A"),
                    "mean_A", "sd_A", "count_B", paste0("hist", ds, "_B"),
                    paste0("delta_h", ds[-1], "_B"), paste0("quantile", probs, "_B"),
                    "mean_B", "sd_B", "count", paste0("hist", ds),
                    paste0("delta_h", ds[-1]), paste0("quantile", probs),
                    "mean", "sd", "PcorH", "KcorH", "ScorH",
                    "PcovH","KcovH","ScovH","PcorD","KcorD","ScorD",
                    "PcovD","KcovD","ScovD","PcorQ","KcorQ","ScorQ",
                    "PcovQ","KcovQ","ScovQ","ratioABtot",
                    paste0("ratioABhist", ds,))
      ss=t(ss)
    }else{
      ss = sim_2_categories(demography$lambda_t_A, demography$lambda_t_B, t_,
                            calCurves = "intcal20", errors = cereals_dates$Error)    
    }
    #ds = seq(time_range_BP[1],time_range_BP[2]-1,-100)[-1]+50
    #plot(ds, ss[paste0("hist", ds, "_A")], type="l",lwd=2)
    #lines(ds, ss[paste0("hist", ds, "_B")], col="red",lwd=2)
    #ss$PcorH
    params = as.data.frame(t(c(alpha, beta,
                               demography$pi,
                               demography$lambda_values_A,
                               demography$growth_rates_A,
                               demography$lambda_values_B,
                               demography$growth_rates_B)))
    names(params) = c("alpha", "beta","pi",
                  paste0("lambda", skyline_years, "_A"),
                  paste0("rate", skyline_years[1:(length(skyline_years)-1)]+(skyline_years[1]-skyline_years[2])/2, "_A"),
                  paste0("lambda", skyline_years, "_B"),
                  paste0("rate", skyline_years[1:(length(skyline_years)-1)]+(skyline_years[1]-skyline_years[2])/2, "_B"))
    cbind(params,ss)
    
  }
  reftable = reftable[!is.na(reftable$count),]
  reftable = reftable[reftable$count!=1,]
  
  save(reftable, file="results/Cereals_correlated_model_reftable.rda")
  stopCluster(cl) 
}
# load( file="results/Cereals_correlated_model_reftable.rda")
















