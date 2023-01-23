library(abcrf)
source("../../Misc_R_tools/color.R")
source("scripts/sim14c.R")

load(file = "results/num_of_sims.rda")

# load target (i.e. observed) summary statistics
load(file = "results/Cereals_sumstats.rda")

# lead reference tables
load(file = "results/Cereals_proportional_model_reftable.rda")
reftable = reftable[!is.na(reftable$count_A),]
reftable = reftable[reftable$count_A!=1,]
reftable = reftable[!is.na(reftable$count_B),]
reftable = reftable[reftable$count_B!=1,]




load(file = "results/Cereals_time_range_BP.rda")
num_of_periods = round((time_range_BP[1] - time_range_BP[2]) / 100 / 4)
skyline_years = get_time_of_change(num_of_periods, time_range_BP, 
                                   intervals = "regular")

if ( !file.exists("results/Cereals_proportional_posterior.rda") ){

  sumstats = rbind(reftable[names(cereals_sumstats_2_categories)])
  
  
  # reminder:
  # logit :  log(x/(1-x))
  # inverse logit : exp(x)/(1+exp(x))
  # pi
  logit_pi = log(reftable["pi"]/(1-reftable["pi"]))
  names(logit_pi) = "logit_pi"
  RF_lambda = regAbcrf(logit_pi~., data.frame(logit_pi,sumstats),
                       ntree = 1000, paral = TRUE)
  posterior_logit_pi = predict(RF_lambda, cereals_sumstats_2_categories,
                             training = data.frame(logit_pi,sumstats),
                             paral = TRUE, rf.weights = FALSE)  
  # lambda
  lambda_error = rep(NA,length(skyline_years))
  lambda_hat = rep(NA,length(skyline_years))
  lambda_95low = rep(NA,length(skyline_years))
  lambda_95upp = rep(NA,length(skyline_years))
  for (i in seq_along(skyline_years)){
    param_name_A = paste0("lambda",skyline_years[i],"_A")
    param_name_B = paste0("lambda",skyline_years[i],"_B")
    param_index_A = which(names(reftable)==param_name_A)
    param_index_B = which(names(reftable)==param_name_B)
    
    param = log10(reftable[param_index_A]+reftable[param_index_B])
    names(param) = "param"
    RF_lambda = regAbcrf(param~., data.frame(param,sumstats),
                         ntree = 1000, paral = TRUE)
    
    posterior_lambda = predict(RF_lambda, cereals_sumstats_2_categories,
                               training = data.frame(param,sumstats),
                               paral = TRUE, rf.weights = FALSE) 
    
    lambda_error[i] = RF_lambda$model.rf$prediction.error
    lambda_hat[i] = 10^(posterior_lambda$med[1])
    lambda_95low[i] = 10^(posterior_lambda$quantiles[1])
    lambda_95upp[i] = 10^(posterior_lambda$quantiles[2])
  }
  

  
  save(skyline_years, posterior_logit_pi,
       lambda_error, lambda_hat, lambda_95low, lambda_95upp,
       file="results/Cereals_proportional_posterior.rda")
  
  
}
load(file = "results/Cereals_proportional_posterior.rda")


posterior_pi_median = exp(posterior_logit_pi$med)/(1+exp(posterior_logit_pi$med))
posterior_pi_quantiles = exp(posterior_logit_pi$quantiles)/(1+exp(posterior_logit_pi$quantiles))


load(file = "results/Cereals_spd.rda")
load(file = "results/Cereals_time_range_BP.rda")

pdf(file="results/Cereals_proportional_model_result.pdf", width=10, height=5)
par(mar=c(4.5, 4.5, 1, 1) + 0.1)
#plot(triticum_spd$grid$calBP, triticum_spd$grid$PrDens, xlim = time_range_BP, ylim = c(0.0001, 1), log = "y", type="l", xlab="Years cal BP", ylab=expression(lambda), col="grey", lwd = 2)
plot(skyline_years, lambda_hat, xlim = time_range_BP, ylim = c(0.0001, 1), log = "y",
      type="l", xlab="Years cal BP", ylab=expression(lambda), col = PCI_blue, lwd = 2)
lines(skyline_years, lambda_95low, lty = 2, lwd = 2, col = PCI_blue)
lines(skyline_years, lambda_95upp, lty = 2, lwd = 2, col = PCI_blue)

text(3000,0.0006,expression(hat(pi)*"=0.40 (0.22, 0.50)"), cex=1.5)

dev.off()
