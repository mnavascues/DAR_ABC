library(abcrf)
source("../../Misc_R_tools/color.R")
source("scripts/sim14c.R")

load(file = "results/Bevan_num_of_sims.rda")

# load target (i.e. observed) summary statistics
load(file = "results/Bevan_sumstats.rda")

# lead reference tables for constant model
load(file = "results/Bevan_constant_model_reftable.rda")

if ( !file.exists("results/Bevan_constant_posterior.rda") ){

  sumstats = reftable[names(all_sumstats_c14)]
  param = log10(reftable$lambda)
  # hist(log10(param))
  
  RF_log10lambda = regAbcrf(param~., data.frame(param,sumstats),
                               ntree = 5000, paral = TRUE)
  
  posterior_log10lambda = predict(RF_log10lambda, all_sumstats_c14,
                                  training = data.frame(param,sumstats),
                                  paral = TRUE, rf.weights = TRUE) 
  
  save(RF_log10lambda, posterior_log10lambda, file="results/Bevan_constant_posterior.rda")
}
load(file="results/Bevan_constant_posterior.rda")
(posterior_log10lambda)

lambda_hat = 10^posterior_log10lambda$med[1]
lambda_CI = 10^posterior_log10lambda$quantiles
save(lambda_hat, lambda_CI, file="results/Bevan_lambda_hat_constant.rda")

load(file="results/Bevan_lambda_hat_constant.rda")

pdf(file="results/Bevan_posterior_lambda_constant.pdf", width=10, height=5)
breaks= seq(-3,1.1,0.1)
hist(log10(reftable$lambda),
     breaks = breaks,
     main="",#expression("Posterior probabilty number of periods ("*italic(m)*")"),
     xlab=expression(log[10]*lambda),
     ylim=c(0,3),
     col=adjustcolor( "gray", alpha.f = 0.6),freq=F)
wtd.hist(log10(reftable$lambda),
         breaks = breaks,
         col=PCI_t_blue,
         weight = posterior_log10lambda$weights,
         add=T, freq=F)
box()
dev.off()

load(file = "results/Bevan_spd.rda")
load(file = "results/Bevan_time_range_BP.rda")

pdf(file="results/Bevan_constant_model_result.pdf", width=10, height=5)
par(mar=c(4.5, 4.5, 1, 1) + 0.1)
plot(allspd$grid$calBP, allspd$grid$PrDens, xlim = time_range_BP, ylim = c(0.01,10), log = "y",
     type="l", xlab="Years cal BP", ylab=expression(lambda), col="grey", lwd=2)
abline(h=lambda_hat, col=PCI_blue, lwd=2)
abline(h=lambda_CI[1], col=PCI_blue, lwd=2, lty=2)
abline(h=lambda_CI[2], col=PCI_blue, lwd=2, lty=2)
dev.off()
