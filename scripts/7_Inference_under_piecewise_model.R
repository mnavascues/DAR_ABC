library(abcrf)
source("../../Misc_R_tools/color.R")
source("scripts/sim14c.R")

load(file = "results/Bevan_num_of_sims.rda")

# load target (i.e. observed) summary statistics
load(file = "results/Bevan_sumstats.rda")

# lead reference tables for piecewise model
load(file = "results/Bevan_piecewise_model_reftable.rda")

load(file = "results/Bevan_time_range_BP.rda")
#load(file="results/Bevan_m_hat.rda")
num_of_periods = round((time_range_BP[1] - time_range_BP[2]) / 100 / 4)
skyline_years = get_time_of_change(num_of_periods, time_range_BP, intervals="regular")

if ( !file.exists("results/Bevan_piecewise_posterior.rda") ){

  sumstats = reftable[names(all_sumstats_c14)]
  
  lambda_error = rep(NA,length(skyline_years))
  lambda_hat = rep(NA,length(skyline_years))
  lambda_95low = rep(NA,length(skyline_years))
  lambda_95upp = rep(NA,length(skyline_years))
  
  for (i in seq_along(skyline_years)){
    param_name = paste0("lambda",skyline_years[i])
    param_index = which(names(reftable)==param_name)
    param = log10(reftable[param_index])
    names(param) = "param"
    RF_lambda = regAbcrf(param~., data.frame(param,sumstats),
                         ntree = 1000, paral = TRUE)
    
    posterior_lambda = predict(RF_lambda, all_sumstats_c14,
                               training = data.frame(param,sumstats),
                               paral = TRUE, rf.weights = FALSE) 
    
    lambda_error[i] = RF_lambda$model.rf$prediction.error
    lambda_hat[i] = 10^(posterior_lambda$med[1])
    lambda_95low[i] = 10^(posterior_lambda$quantiles[1])
    lambda_95upp[i] = 10^(posterior_lambda$quantiles[2])
  }
  save(skyline_years, lambda_error, lambda_hat, lambda_95low, lambda_95upp, 
       file="results/Bevan_piecewise_posterior.rda")
  
  
}
load(file = "results/Bevan_piecewise_posterior.rda")

load(file = "results/Bevan_spd.rda")
load(file = "results/Bevan_time_range_BP.rda")

pdf(file="results/Bevan_piecewise_model_result.pdf", width=10, height=5)
par(mar=c(4.5, 4.5, 1, 1) + 0.1)
plot(allspd$grid$calBP, allspd$grid$PrDens, xlim = time_range_BP, ylim = c(0.01, 10), log = "y",
     type="l", xlab="Years cal BP", ylab=expression(lambda), col="grey", lwd = 2)
lines(skyline_years, lambda_hat, col = PCI_blue, lwd = 2)
lines(skyline_years, lambda_95low, lty = 2, lwd = 2, col = PCI_blue)
lines(skyline_years, lambda_95upp, lty = 2, lwd = 2, col = PCI_blue)
dev.off()

pdf(file="results/Bevan_piecewise_model_result_error.pdf", width=10, height=5)
par(mar=c(4.5, 4.5, 1, 1) + 0.1)
plot(skyline_years, lambda_error,
     type = "l", xlim = time_range_BP, lwd = 2,
     xlab = "Years cal BP",
     ylab = expression("OOB error "*log[10](lambda)), col = PCI_blue)
dev.off()






rate_param_names = paste0("rate", skyline_years[1:(length(skyline_years)-1)]+(skyline_years[1]-skyline_years[2])/2)
skyline_years_midpoint = skyline_years[1:(length(skyline_years)-1)]+(skyline_years[1]-skyline_years[2])/2

if ( !file.exists("results/Bevan_piecewise_posterior_r.rda") ){
  sumstats = reftable[names(all_sumstats_c14)]
  rate_error = rep(NA,length(rate_param_names))
  rate_hat = rep(NA,length(rate_param_names))
  rate_95low = rep(NA,length(rate_param_names))
  rate_95upp = rep(NA,length(rate_param_names))
  for (i in seq_along(rate_param_names)){
    param_name = rate_param_names[i]
    param_index = which(names(reftable)==param_name)
    param = reftable[param_index]
    names(param) = "param"
    RF_rate = regAbcrf(param~., data.frame(param,sumstats),
                       ntree = 1000, paral = TRUE, ncores = 24)
    
    posterior_rate = predict(RF_rate, all_sumstats_c14,
                             training = data.frame(param,sumstats),
                             paral = TRUE, ncores = 24, rf.weights = FALSE) 
    
    rate_error[i] = RF_rate$model.rf$prediction.error
    rate_hat[i] = posterior_rate$med[1]
    rate_95low[i] = posterior_rate$quantiles[1]
    rate_95upp[i] = posterior_rate$quantiles[2]
  }
  save(skyline_years_midpoint, rate_error, rate_hat, rate_95low, rate_95upp, 
       file="results/Bevan_piecewise_posterior_r.rda")
}
load(file="results/Bevan_piecewise_posterior_r.rda")

pdf(file="results/Bevan_piecewise_model_result_rate_error.pdf", width=10, height=5)
par(mar=c(4.5, 4.5, 1, 1) + 0.1)
plot(skyline_years_midpoint,rate_error,
     xlab = "Years cal BP", ylab=expression("OOB error ("*italic(r)*")"),
     type="l",xlim=time_range_BP,lwd=2, col=PCI_blue)
dev.off()



pdf(file="results/Bevan_piecewise_model_result_rate.pdf", width=10, height=5)
step_wise_years = c(skyline_years[1],rep(skyline_years[2:(length(skyline_years)-1)],each=2),skyline_years[length(skyline_years)])
par(mar=c(4.5, 4.5, 1, 1) + 0.1)
plot(step_wise_years, rep(rate_hat,each=2),
     xlab = "Years cal BP", ylab=expression(italic(r)),
     lwd=2, type="l", xlim = time_range_BP, ylim=c(-0.018,0.018), col=PCI_blue)
lines(step_wise_years,rep(rate_95low,each=2),lty=3, lwd=2,col=PCI_blue)
lines(step_wise_years,rep(rate_95upp,each=2),lty=3, lwd=2,col=PCI_blue)
abline(h=0,col="gray")
points( skyline_years[which(rate_95low>0)]-70, rep(-0.0185,sum(rate_95low>0)),
        pch="*", cex=2.5, col=PCI_blue )
points( skyline_years[which(rate_95upp<0)]-70, rep(-0.0185,sum(rate_95upp<0)),
        pch="*", cex=2.5, col=PCI_blue )
dev.off()

