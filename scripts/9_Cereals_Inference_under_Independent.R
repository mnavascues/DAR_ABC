library(abcrf)
source("../../Misc_R_tools/color.R")
source("scripts/sim14c.R")

load(file = "results/num_of_sims.rda")

# load target (i.e. observed) summary statistics
load(file = "results/Cereals_sumstats.rda")

# lead reference tables
load(file = "results/Cereals_independent_model_reftable.rda")
reftable = reftable[!is.na(reftable$count_A),]
reftable = reftable[reftable$count_A!=1,]
reftable = reftable[!is.na(reftable$count_B),]
reftable = reftable[reftable$count_B!=1,]




load(file = "results/Cereals_time_range_BP.rda")
num_of_periods = round((time_range_BP[1] - time_range_BP[2]) / 100 / 4)
skyline_years = get_time_of_change(num_of_periods, time_range_BP, 
                                   intervals = "regular")

if ( !file.exists("results/Cereals_independent_posterior.rda") ){

  sumstats = rbind(reftable[names(cereals_sumstats_2_categories)])
  
  # A
  lambda_error_A = rep(NA,length(skyline_years))
  lambda_hat_A = rep(NA,length(skyline_years))
  lambda_95low_A = rep(NA,length(skyline_years))
  lambda_95upp_A = rep(NA,length(skyline_years))
  for (i in seq_along(skyline_years)){
    param_name = paste0("lambda",skyline_years[i],"_A")
    param_index = which(names(reftable)==param_name)
    param = log10(reftable[param_index])
    names(param) = "param"
    RF_lambda = regAbcrf(param~., data.frame(param,sumstats),
                         ntree = 1000, paral = TRUE)
    
    posterior_lambda = predict(RF_lambda, cereals_sumstats_2_categories,
                               training = data.frame(param,sumstats),
                               paral = TRUE, rf.weights = FALSE) 
    
    lambda_error_A[i] = RF_lambda$model.rf$prediction.error
    lambda_hat_A[i] = 10^(posterior_lambda$med[1])
    lambda_95low_A[i] = 10^(posterior_lambda$quantiles[1])
    lambda_95upp_A[i] = 10^(posterior_lambda$quantiles[2])
  }
  
  # B
  lambda_error_B = rep(NA,length(skyline_years))
  lambda_hat_B = rep(NA,length(skyline_years))
  lambda_95low_B = rep(NA,length(skyline_years))
  lambda_95upp_B = rep(NA,length(skyline_years))
  for (i in seq_along(skyline_years)){
    param_name = paste0("lambda",skyline_years[i],"_B")
    param_index = which(names(reftable)==param_name)
    param = log10(reftable[param_index])
    names(param) = "param"
    RF_lambda = regAbcrf(param~., data.frame(param,sumstats),
                         ntree = 1000, paral = TRUE)
    
    posterior_lambda = predict(RF_lambda, cereals_sumstats_2_categories,
                               training = data.frame(param,sumstats),
                               paral = TRUE, rf.weights = FALSE) 
    
    lambda_error_B[i] = RF_lambda$model.rf$prediction.error
    lambda_hat_B[i] = 10^(posterior_lambda$med[1])
    lambda_95low_B[i] = 10^(posterior_lambda$quantiles[1])
    lambda_95upp_B[i] = 10^(posterior_lambda$quantiles[2])
  }
  
  
  save(skyline_years,
       lambda_error_A, lambda_hat_A, lambda_95low_A, lambda_95upp_A,
       lambda_error_B, lambda_hat_B, lambda_95low_B, lambda_95upp_B,  
       file="results/Cereals_independent_posterior.rda")
  
  
}
load(file = "results/Cereals_independent_posterior.rda")






load(file = "results/Cereals_spd.rda")
load(file = "results/Cereals_time_range_BP.rda")

pdf(file="results/Cereals_independent_model_result.pdf", width=10, height=5)
par(mar=c(4.5, 4.5, 1, 1) + 0.1)
#plot(triticum_spd$grid$calBP, triticum_spd$grid$PrDens, xlim = time_range_BP, ylim = c(0.0001, 1), log = "y", type="l", xlab="Years cal BP", ylab=expression(lambda), col="grey", lwd = 2)
plot(skyline_years, lambda_hat_A, xlim = time_range_BP, ylim = c(0.0001, 1), log = "y",
      type="l", xlab="Years cal BP", ylab=expression(lambda), col = PCI_blue, lwd = 2)
lines(skyline_years, lambda_95low_A, lty = 2, lwd = 2, col = PCI_blue)
lines(skyline_years, lambda_95upp_A, lty = 2, lwd = 2, col = PCI_blue)

#lines(hordeum_spd$grid$calBP, hordeum_spd$grid$PrDens, col="orange", lwd = 2)
lines(skyline_years, lambda_hat_B, lwd = 2)
lines(skyline_years, lambda_95low_B, lty = 2, lwd = 2)
lines(skyline_years, lambda_95upp_B, lty = 2, lwd = 2)


dev.off()
