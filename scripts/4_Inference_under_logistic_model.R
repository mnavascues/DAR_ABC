library(abcrf)
source("../../Misc_R_tools/color.R")
source("scripts/sim14c.R")

load(file = "results/num_of_sims.rda")

# load target (i.e. observed) summary statistics
load(file = "results/sumstats.rda")

# lead reference tables for logistic model
load(file = "results/logistic_model_reftable.rda")
rows2keep = complete.cases(reftable)
reftable = reftable[rows2keep,]

sumstats = reftable[names(all_sumstats_c14)]

if ( !file.exists("results/logistic_posterior_lambda_0.rda") ){

  param = log10(reftable$lambda_0)
  RF_log10lambda_0 = regAbcrf(param~., data.frame(param,sumstats),
                              ntree = 5000, paral = TRUE)
  posterior_log10lambda_0 = predict(RF_log10lambda_0, all_sumstats_c14,
                                    training = data.frame(param,sumstats),
                                    paral = TRUE, rf.weights = TRUE) 
  save(RF_log10lambda_0, posterior_log10lambda_0,
       file="results/logistic_posterior_lambda_0.rda")
}
rm(RF_log10lambda_0);gc()
if ( !file.exists("results/logistic_posterior_K.rda") ){
  param = log10(reftable$K)
  RF_log10K = regAbcrf(param~., data.frame(param,sumstats),
                              ntree = 5000, paral = TRUE)
  posterior_log10K = predict(RF_log10K, all_sumstats_c14,
                                    training = data.frame(param,sumstats),
                                    paral = TRUE, rf.weights = TRUE) 
  save(RF_log10K, posterior_log10K,
       file="results/logistic_posterior_K.rda")
}
rm(RF_log10K);gc()
if ( !file.exists("results/logistic_posterior_r.rda") ){
  param = reftable$r
  RF_r = regAbcrf(param~., data.frame(param,sumstats),
                       ntree = 5000, paral = TRUE)
  posterior_r = predict(RF_r, all_sumstats_c14,
                             training = data.frame(param,sumstats),
                             paral = TRUE, rf.weights = TRUE) 
  save(RF_r, posterior_r,
       file="results/logistic_posterior_r.rda")
}
rm(RF_r);gc()

load("results/logistic_posterior_lambda_0.rda")
rm(RF_log10lambda_0);gc()
load("results/logistic_posterior_K.rda")
rm(RF_log10K);gc()
load("results/logistic_posterior_r.rda")
rm(RF_r);gc()
(posterior_log10lambda_0)
(posterior_log10K)
(posterior_r)

if ( !file.exists("results/parameter_estimates_logistic.rda") ){
  lambda_0_hat = 10^posterior_log10lambda_0$med[1]
  lambda_0_CI = 10^posterior_log10lambda_0$quantiles
  K_hat = 10^posterior_log10K$med[1]
  K_CI = 10^posterior_log10K$quantiles
  r_hat = posterior_r$med[1]
  r_CI = posterior_r$quantiles
  save(lambda_0_hat, lambda_0_CI,
       K_hat, K_CI,
       r_hat, r_CI, file="results/parameter_estimates_logistic.rda")
}



pdf(file="results/posterior_lambda_0_logistic.pdf", width=4, height=4)
par(mar=c(4.5, 4.5, 1, 1) + 0.1)
breaks= seq(-3,1.1,0.05)
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


pdf(file="results/posterior_K_logistic.pdf", width=4, height=4)
par(mar=c(4.5, 4.5, 1, 1) + 0.1)
breaks= seq(-3,2,0.05)
hist(log10(reftable$K),
     breaks = breaks,
     main="",#expression("Posterior probabilty number of periods ("*italic(m)*")"),
     xlab=expression(log[10]*italic(K)),
     ylim=c(0,5),
     col=adjustcolor( "gray", alpha.f = 0.6),freq=F)
wtd.hist(log10(reftable$K),
         breaks = breaks,
         col=PCI_t_blue,
         weight = posterior_log10K$weights,
         add=T, freq=F)
box()
dev.off()



pdf(file="results/posterior_r_logistic.pdf", width=4, height=4)
par(mar=c(4.5, 4.5, 1, 1) + 0.1)
breaks= seq(0,0.002,0.00003)
hist(reftable$r,
     breaks = breaks,
     main="",#expression("Posterior probabilty number of periods ("*italic(m)*")"),
     xlab=expression(italic(r)),
     ylim=c(0,3500),
     col=adjustcolor( "gray", alpha.f = 0.6),freq=F)
wtd.hist(reftable$r,
         breaks = breaks,
         col=PCI_t_blue,
         weight = posterior_r$weights,
         add=T, freq=F)
box()
dev.off()
