library(abcrf)
source("../../Misc_R_tools/color.R")
source("scripts/sim14c.R")

load(file = "results/Bevan_num_of_sims.rda")

# load target (i.e. observed) summary statistics
load(file = "results/Bevan_sumstats.rda")

# lead reference tables for exponential model
load(file = "results/Bevan_exponential_model_reftable.rda")

if ( !file.exists("results/Bevan_exponential_posterior.rda") ){

  sumstats = reftable[names(all_sumstats_c14)]
  
  # params: lambda_0"     "lambda_f"     "r"
  param = log10(reftable$lambda_0)
  RF_log10lambda_0 = regAbcrf(param~., data.frame(param,sumstats),
                              ntree = 5000, paral = TRUE)
  posterior_log10lambda_0 = predict(RF_log10lambda_0, all_sumstats_c14,
                                    training = data.frame(param,sumstats),
                                    paral = TRUE, rf.weights = TRUE) 
  param = log10(reftable$lambda_f)
  RF_log10lambda_f = regAbcrf(param~., data.frame(param,sumstats),
                              ntree = 5000, paral = TRUE)
  posterior_log10lambda_f = predict(RF_log10lambda_f, all_sumstats_c14,
                                    training = data.frame(param,sumstats),
                                    paral = TRUE, rf.weights = TRUE) 
  
  param = reftable$r
  RF_r = regAbcrf(param~., data.frame(param,sumstats),
                  ntree = 5000, paral = TRUE)
  posterior_r = predict(RF_r, all_sumstats_c14,
                        training = data.frame(param,sumstats),
                        paral = TRUE, rf.weights = TRUE) 
  
  save(RF_r, posterior_r,
       RF_log10lambda_0, posterior_log10lambda_0,
       RF_log10lambda_f, posterior_log10lambda_f, file="results/Bevan_exponential_posterior.rda")
}
#load(file="results/Bevan_exponential_posterior.rda")
(posterior_log10lambda_0)
(posterior_log10lambda_f)
(posterior_r)

lambda_0_hat = 10^posterior_log10lambda_0$med[1]
lambda_0_CI = 10^posterior_log10lambda_0$quantiles
lambda_f_hat = 10^posterior_log10lambda_f$med[1]
lambda_f_CI = 10^posterior_log10lambda_f$quantiles
r_hat = posterior_r$med[1]
r_CI = posterior_r$quantiles
save(r_hat, r_CI,
     lambda_0_hat, lambda_0_CI,
     lambda_f_hat, lambda_f_CI, file="results/Bevan_lambda_hat_exponential.rda")

load(file="results/Bevan_lambda_hat_exponential.rda")

pdf(file="results/Bevan_posterior_lambda_0_exponential.pdf", width=10, height=5)
breaks= seq(-3,1.1,0.1)
hist(log10(reftable$lambda_0),
     breaks = breaks,
     main="",#expression("Posterior probabilty number of periods ("*italic(m)*")"),
     xlab=expression(log[10]*lambda[0]),
     ylim=c(0,4),
     col=adjustcolor( "gray", alpha.f = 0.6),freq=F)
wtd.hist(log10(reftable$lambda_0),
         breaks = breaks,
         col=PCI_t_blue,
         weight = posterior_log10lambda_0$weights,
         add=T, freq=F)
box()
dev.off()

pdf(file="results/Bevan_posterior_lambda_f_exponential.pdf", width=10, height=5)
breaks= seq(-3,1.1,0.1)
hist(log10(reftable$lambda_f),
     breaks = breaks,
     main="",#expression("Posterior probabilty number of periods ("*italic(m)*")"),
     xlab=expression(log[10]*lambda["f"]),
     ylim=c(0,1),
     col=adjustcolor( "gray", alpha.f = 0.6),freq=F)
wtd.hist(log10(reftable$lambda_f),
         breaks = breaks,
         col=PCI_t_blue,
         weight = posterior_log10lambda_f$weights,
         add=T, freq=F)
box()
dev.off()

pdf(file="results/Bevan_posterior_r_exponential.pdf", width=10, height=5)
breaks= seq(0,0.0009,0.00002)
hist(reftable$r,
     breaks = breaks,
     main="",#expression("Posterior probabilty number of periods ("*italic(m)*")"),
     xlab=expression(italic("r")),
     ylim=c(0,10000),
     col=adjustcolor( "gray", alpha.f = 0.6),freq=F)
wtd.hist(reftable$r,
         breaks = breaks,
         col=PCI_t_blue,
         weight = posterior_r$weights,
         add=T, freq=F)
box()
dev.off()

load(file = "results/Bevan_spd.rda")
load(file = "results/Bevan_time_range_BP.rda")

pdf(file="results/Bevan_exponential_model_result.pdf", width=10, height=5)
par(mar=c(4.5, 4.5, 1, 1) + 0.1)
plot(allspd$grid$calBP, allspd$grid$PrDens, xlim = time_range_BP, ylim = c(0.01,10), log = "y",
     type="l", xlab="Years cal BP", ylab=expression(lambda), col="grey", lwd=2)

lambda = lambda_0_hat*exp(r_hat*1:(time_range_BP[1]-time_range_BP[2]))
lines(t_, lambda, col=PCI_blue, lwd=2)

lambda = lambda_0_CI[1]*exp(r_CI[1]*1:(time_range_BP[1]-time_range_BP[2]))
lines(t_, lambda, col=PCI_blue, lwd=2, lty=2)

lambda = lambda_0_CI[2]*exp(r_CI[2]*1:(time_range_BP[1]-time_range_BP[2]))
lines(t_, lambda, col=PCI_blue, lwd=2, lty=2)

dev.off()

