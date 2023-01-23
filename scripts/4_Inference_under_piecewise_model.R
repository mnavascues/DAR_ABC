library(abcrf)
source("../../Misc_R_tools/color.R")
source("scripts/sim14c.R")

load(file = "results/num_of_sims.rda")

# load target (i.e. observed) summary statistics
load(file = "results/sumstats.rda")

# lead reference tables for piecewise model
load(file = "results/piecewise_model_reftable.rda")
sumstats = reftable[names(all_sumstats_c14)]

load(file = "results/time_range_BP.rda")

num_of_periods = round((time_range_BP[1] - time_range_BP[2]) / 100 / 4)
skyline_years = get_time_of_change(num_of_periods, time_range_BP, intervals="regular")


if ( !file.exists(file="results/parameter_estimates_piecewise.rda") ){
  lambda_hat      = rep(NA,length(skyline_years))
  lambda_95CI_low = rep(NA,length(skyline_years))
  lambda_95CI_upp = rep(NA,length(skyline_years))
  lambda_error    = rep(NA,length(skyline_years))
  for (i in seq_len(length(skyline_years))){
    results_file = paste0("results/piecewise_posterior_lambda_",
                          skyline_years[i], ".rda")
    if ( !file.exists(results_file) ){
      param_name = paste0("lambda",skyline_years[i])
      param_index = which(names(reftable)==param_name)
      param = log10(reftable[param_index])
      names(param) = "param"
      RF_log10lambda = regAbcrf(param~., data.frame(param,sumstats),
                                ntree = 1000, paral = TRUE)
      posterior_log10lambda = predict(RF_log10lambda, all_sumstats_c14,
                                      training = data.frame(param,sumstats),
                                      paral = TRUE, rf.weights = TRUE) 
      
      save(RF_log10lambda, posterior_log10lambda, file = results_file)
    }else{load(file = results_file)}
    lambda_error[i] = RF_log10lambda$model.rf$prediction.error
    rm(RF_log10lambda);gc()
    lambda_hat[i] = 10^(posterior_log10lambda$med[1])
    lambda_95CI_low[i] = 10^(posterior_log10lambda$quantiles[1])
    lambda_95CI_upp[i] = 10^(posterior_log10lambda$quantiles[2])
  }
  
  rate_param_names = paste0("rate", skyline_years[1:(length(skyline_years)-1)]+(skyline_years[1]-skyline_years[2])/2)
  rate_error    = rep(NA,length(rate_param_names))
  rate_hat      = rep(NA,length(rate_param_names))
  rate_95CI_low = rep(NA,length(rate_param_names))
  rate_95CI_upp = rep(NA,length(rate_param_names))
  for (i in seq_along(rate_param_names)){
    results_file = paste0("results/piecewise_posterior_",
                          rate_param_names[i], ".rda")
    if ( !file.exists(results_file) ){
      param_name = rate_param_names[i]
      param_index = which(names(reftable)==param_name)
      param = reftable[param_index]
      names(param) = "param"
      RF_rate = regAbcrf(param~., data.frame(param,sumstats),
                                ntree = 1000, paral = TRUE)
      posterior_rate = predict(RF_rate, all_sumstats_c14,
                                      training = data.frame(param,sumstats),
                                      paral = TRUE, rf.weights = TRUE) 
      
      save(RF_rate, posterior_rate, file = results_file)
    }else{load(file = results_file)}
    rate_error[i] = RF_rate$model.rf$prediction.error
    rm(RF_rate);gc()
    rate_hat[i] = (posterior_rate$med[1])
    rate_95CI_low[i] = (posterior_rate$quantiles[1])
    rate_95CI_upp[i] = (posterior_rate$quantiles[2])
  }
  
  
  
  save(lambda_hat,
       lambda_95CI_low, lambda_95CI_upp,
       lambda_error,
       rate_hat,
       rate_95CI_low, rate_95CI_upp,
       rate_error,
       file="results/parameter_estimates_piecewise.rda")
}else{load(file="results/parameter_estimates_piecewise.rda")}
########## OLD FROM HERE
load(file = "results/spd.rda")
load(file = "results/time_range_BP.rda")

pdf(file="results/piecewise_model_result_lambda.pdf", width=10, height=5)
par(mar=c(4.5, 4.5, 1, 1) + 0.1)
plot(allspd$grid$calBP, allspd$grid$PrDens, xlim = time_range_BP, ylim = c(0.04, 6), log = "y",
     type="l", xlab="Years cal BP", ylab=expression(lambda), col="grey", lwd = 2)
lines(skyline_years, lambda_hat, col = PCI_blue, lwd = 2)
lines(skyline_years, lambda_95CI_low, lty = 2, lwd = 2, col = PCI_blue)
lines(skyline_years, lambda_95CI_upp, lty = 2, lwd = 2, col = PCI_blue)
dev.off()







#skyline_years_midpoint = skyline_years[1:(length(skyline_years)-1)]+(skyline_years[1]-skyline_years[2])/2


pdf(file="results/piecewise_model_result_rate.pdf", width=10, height=5)
step_wise_years = c(skyline_years[1],rep(skyline_years[2:(length(skyline_years)-1)],each=2),skyline_years[length(skyline_years)])
par(mar=c(4.5, 4.5, 1, 1) + 0.1)
plot(step_wise_years, rep(rate_hat,each=2),
     xlab = "Years cal BP", ylab=expression(italic(r)),
     lwd = 2, type = "l", xlim = time_range_BP,
     ylim = c(-0.006, 0.006),
     #ylim=c(min(rate_95CI_low),max(rate_95CI_upp)),
     col=PCI_blue)
lines(step_wise_years,rep(rate_95CI_low,each=2),lty=3, lwd=2,col=PCI_blue)
lines(step_wise_years,rep(rate_95CI_upp,each=2),lty=3, lwd=2,col=PCI_blue)
abline(h=0,col="gray")
points( skyline_years[which(rate_95CI_low>0)]-196.5, rep(0,sum(rate_95CI_low>0)),
        pch="*", cex=2)
points( skyline_years[which(rate_95CI_upp<0)]-196.5, rep(0,sum(rate_95CI_upp<0)),
        pch="*", cex=2)
dev.off()







for (i in seq_len(length(skyline_years))){
  results_file = paste0("results/piecewise_posterior_lambda_", skyline_years[i], ".rda")
  load(results_file)
  param_name = paste0("lambda",skyline_years[i])
  param_index = which(names(reftable)==param_name)
  
  pdf_file_name = paste0("results/posterior_lambda_",skyline_years[i],"_piecewise.pdf")
  pdf(file=pdf_file_name, width=4, height=4)
  par(mar=c(4.5, 4.5, 1, 1) + 0.1)
  breaks= seq(-3,1.1,0.02)
  
  hist(log10(t(reftable[param_index])),
       breaks = breaks,
       main = "",
       xlab = bquote(log[10]*lambda[.(skyline_years[i])]),
       ylim = c(0,10),
       col = adjustcolor( "gray", alpha.f = 0.6), freq = F)
  wtd.hist(log10(t(reftable[param_index])),
           breaks = breaks,
           col=PCI_t_blue,
           weight = posterior_log10lambda$weights,
           add=T, freq=F)
  box()
  dev.off()
}


for (i in seq_along(rate_param_names)){
  results_file = paste0("results/piecewise_posterior_", rate_param_names[i], ".rda")
  load(results_file)
  param_name = rate_param_names[i]
  param_index = which(names(reftable)==param_name)
  
  pdf_file_name = paste0("results/posterior_r_",skyline_years[i],"_",skyline_years[i+1],"_piecewise.pdf")
  pdf(file=pdf_file_name, width=4, height=4)
  par(mar=c(4.5, 4.5, 1, 1) + 0.1)
  breaks= seq(-0.006,0.006,0.0002)
  
  hist(t(reftable[param_index]),
       breaks = breaks,
       main = "",
       xlab = bquote(italic(r)["("*.(skyline_years[i])-.(skyline_years[i+1])*")"]),
       ylim = c(0,900),
       col = adjustcolor( "gray", alpha.f = 0.6), freq = F)
  wtd.hist(t(reftable[param_index]),
           breaks = breaks,
           col=PCI_t_blue,
           weight = posterior_rate$weights,
           add=T, freq=F)
  box()
  dev.off()
}


