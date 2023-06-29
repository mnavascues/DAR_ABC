library(abcrf)
source("../../Misc_R_tools/color.R")
source("scripts/sim14c.R")

load(file = "results/num_of_sims.rda")

# load target (i.e. observed) summary statistics
load(file = "results/Cereals_sumstats.rda")

# lead reference tables
load(file = "results/Cereals_interdependent_model_reftable.rda")
reftable = reftable[!is.na(reftable$count_A),]
reftable = reftable[reftable$count_A!=1,]
reftable = reftable[!is.na(reftable$count_B),]
reftable = reftable[reftable$count_B!=1,]




load(file = "results/Cereals_time_range_BP.rda")
num_of_periods = round((time_range_BP[1] - time_range_BP[2]) / 100 / 4)
skyline_years = get_time_of_change(num_of_periods, time_range_BP, 
                                   intervals = "regular")

if ( !file.exists("results/Cereals_interdependent_parameter_estimates.rda") ){

  sumstats = rbind(reftable[names(cereals_sumstats_2_categories)])
  
  # total
  lambda_error = rep(NA,length(skyline_years))
  lambda_hat = rep(NA,length(skyline_years))
  lambda_95low = rep(NA,length(skyline_years))
  lambda_95upp = rep(NA,length(skyline_years))
  for (i in seq_along(skyline_years)){
    results_file = paste0("results/Cereals_interdependent_posterior_log_lambda_",
                          skyline_years[i], ".rda")
    if ( !file.exists(results_file) ){
      param_name_A = paste0("lambda",skyline_years[i],"_A")
      param_name_B = paste0("lambda",skyline_years[i],"_B")
      param_index_A = which(names(reftable)==param_name_A)
      param_index_B = which(names(reftable)==param_name_B)
      
      param = log10(reftable[param_index_A]+reftable[param_index_B])
      names(param) = "param"
      RF_log_lambda = regAbcrf(param~., data.frame(param,sumstats),
                           ntree = 1000, paral = TRUE)
      
      posterior_log_lambda = predict(RF_log_lambda, cereals_sumstats_2_categories,
                                 training = data.frame(param,sumstats),
                                 paral = TRUE, rf.weights = TRUE) 
      save(RF_log_lambda, posterior_log_lambda,
           file = results_file)
    }else{load(file = results_file)}
    gc()
    lambda_error[i] = RF_log_lambda$model.rf$prediction.error
    lambda_hat[i] = 10^(posterior_log_lambda$med[1])
    lambda_95low[i] = 10^(posterior_log_lambda$quantiles[1])
    lambda_95upp[i] = 10^(posterior_log_lambda$quantiles[2])
  }
  
  # A
  lambda_A_error = rep(NA,length(skyline_years))
  lambda_A_hat = rep(NA,length(skyline_years))
  lambda_A_95low = rep(NA,length(skyline_years))
  lambda_A_95upp = rep(NA,length(skyline_years))
  for (i in seq_along(skyline_years)){
    results_file = paste0("results/Cereals_interdependent_posterior_log_lambda_A_",
                          skyline_years[i], ".rda")
    if ( !file.exists(results_file) ){
      param_name_A = paste0("lambda",skyline_years[i],"_A")
      param_index_A = which(names(reftable)==param_name_A)
      
      param = log10(reftable[param_index_A])
      names(param) = "param"
      RF_log_lambda = regAbcrf(param~., data.frame(param,sumstats),
                           ntree = 1000, paral = TRUE)
      
      posterior_log_lambda = predict(RF_log_lambda, cereals_sumstats_2_categories,
                                 training = data.frame(param,sumstats),
                                 paral = TRUE, rf.weights = TRUE) 
      
      save(RF_log_lambda, posterior_log_lambda,
           file = results_file)
    }else{load(file = results_file)}
    
    gc()
    
    lambda_A_error[i] = RF_log_lambda$model.rf$prediction.error
    lambda_A_hat[i] = 10^(posterior_log_lambda$med[1])
    lambda_A_95low[i] = 10^(posterior_log_lambda$quantiles[1])
    lambda_A_95upp[i] = 10^(posterior_log_lambda$quantiles[2])
  }


  # B
  lambda_B_error = rep(NA,length(skyline_years))
  lambda_B_hat = rep(NA,length(skyline_years))
  lambda_B_95low = rep(NA,length(skyline_years))
  lambda_B_95upp = rep(NA,length(skyline_years))
  for (i in seq_along(skyline_years)){
    results_file = paste0("results/Cereals_interdependent_posterior_log_lambda_B_",
                          skyline_years[i], ".rda")
    if ( !file.exists(results_file) ){
      param_name_B = paste0("lambda",skyline_years[i],"_B")
      param_index_B = which(names(reftable)==param_name_B)
      
      param = log10(reftable[param_index_B])
      names(param) = "param"
      RF_log_lambda = regAbcrf(param~., data.frame(param,sumstats),
                           ntree = 1000, paral = TRUE)
      
      posterior_log_lambda = predict(RF_log_lambda, cereals_sumstats_2_categories,
                                 training = data.frame(param,sumstats),
                                 paral = TRUE, rf.weights = TRUE) 
      
      save(RF_log_lambda, posterior_log_lambda,
           file = results_file)
    }else{load(file = results_file)}
    gc()    
    lambda_B_error[i] = RF_log_lambda$model.rf$prediction.error
    lambda_B_hat[i] = 10^(posterior_log_lambda$med[1])
    lambda_B_95low[i] = 10^(posterior_log_lambda$quantiles[1])
    lambda_B_95upp[i] = 10^(posterior_log_lambda$quantiles[2])
  }
  
  
  # pi
  pi_error = rep(NA,length(skyline_years))
  pi_hat = rep(NA,length(skyline_years))
  pi_95low = rep(NA,length(skyline_years))
  pi_95upp = rep(NA,length(skyline_years))
  for (i in seq_along(skyline_years)){
    results_file = paste0("results/Cereals_interdependent_posterior_pi_",
                          skyline_years[i], ".rda")
    if ( !file.exists(results_file) ){
      param_name_A = paste0("lambda",skyline_years[i],"_A")
      param_name_B = paste0("lambda",skyline_years[i],"_B")
      param_index_A = which(names(reftable)==param_name_A)
      param_index_B = which(names(reftable)==param_name_B)
      
      p = reftable[param_index_A]/(reftable[param_index_A]+reftable[param_index_B])
      param = log(p/(1-p))
      names(param) = "param"
      RF_logit_pi = regAbcrf(param~., data.frame(param,sumstats),
                             ntree = 1000, paral = TRUE)
      
      posterior_logit_pi = predict(RF_logit_pi, cereals_sumstats_2_categories,
                                   training = data.frame(param,sumstats),
                                   paral = TRUE, rf.weights = TRUE) 
      
      save(RF_logit_pi, posterior_logit_pi,
           file = results_file)
    }else{load(file = results_file)}
    gc()    
    pi_error[i] = RF_logit_pi$model.rf$prediction.error
    pi_hat[i] = exp(posterior_logit_pi$med[1])/(1+exp(posterior_logit_pi$med[1]))  
    pi_95low[i] = exp(posterior_logit_pi$quantiles[1])/(1+exp(posterior_logit_pi$quantiles[1]))  
    pi_95upp[i] = exp(posterior_logit_pi$quantiles[2])/(1+exp(posterior_logit_pi$quantiles[2]))   
    print(pi_hat[i])
    print(pi_95low[i] )
    print(pi_95upp[i])
  }
  
  
  save(skyline_years,
       lambda_error, lambda_hat, lambda_95low, lambda_95upp,
       lambda_A_error, lambda_A_hat, lambda_A_95low, lambda_A_95upp,
       lambda_B_error, lambda_B_hat, lambda_B_95low, lambda_B_95upp,
       pi_error, pi_hat, pi_95low, pi_95upp,  
       file="results/Cereals_interdependent_parameter_estimates.rda")
  
  
}else{load(file = "results/Cereals_interdependent_parameter_estimates.rda")}






load(file = "results/Cereals_spd.rda")
load(file = "results/Cereals_time_range_BP.rda")

pdf(file="results/Cereals_interdependent_model_result_lambda.pdf", width=10, height=5)
par(mar=c(4.5, 4.5, 1, 1) + 0.1)
#plot(triticum_spd$grid$calBP, triticum_spd$grid$PrDens, xlim = time_range_BP, ylim = c(0.0001, 1), log = "y", type="l", xlab="Years cal BP", ylab=expression(lambda), col="grey", lwd = 2)
plot(hordeum_spd$grid$calBP, hordeum_spd$grid$PrDens, xlim = time_range_BP, ylim = c(0.0001, 1), log = "y",
     type="l", xlab="Years cal BP", ylab=expression(lambda), col="grey", lwd = 2)
lines(triticum_spd$grid$calBP, triticum_spd$grid$PrDens, lty = 2, lwd = 2, col = "grey")
lines(skyline_years, lambda_hat, col = PCI_blue, lwd = 2)
lines(skyline_years, lambda_95low, lty = 2, lwd = 2, col = PCI_blue)
lines(skyline_years, lambda_95upp, lty = 2, lwd = 2, col = PCI_blue)
text(time_range_BP[1],1,"a",cex=2)
dev.off()

pdf(file="results/Cereals_interdependent_model_result_pi.pdf", width=10, height=5)
par(mar=c(4.5, 4.5, 1, 1) + 0.1)
#plot(triticum_spd$grid$calBP, triticum_spd$grid$PrDens, xlim = time_range_BP, ylim = c(0.0001, 1), log = "y", type="l", xlab="Years cal BP", ylab=expression(lambda), col="grey", lwd = 2)
plot(skyline_years, pi_hat, xlim = time_range_BP, ylim = c(0, 1), 
     type="l", xlab="Years cal BP", ylab=expression(pi), col = PCI_blue, lwd = 2)
lines(skyline_years, pi_95low, lty = 2, lwd = 2, col = PCI_blue)
lines(skyline_years, pi_95upp, lty = 2, lwd = 2, col = PCI_blue)
text(time_range_BP[1],1,"b",cex=2)
dev.off()
