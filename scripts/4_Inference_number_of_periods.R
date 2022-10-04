library(abcrf)
source("../../Misc_R_tools/color.R")
source("scripts/sim14c.R")

load(file = "results/Bevan_num_of_sims.rda")

# load target (i.e. observed) summary statistics
load(file = "results/Bevan_sumstats.rda")

# lead reference tables for piecewise model
load(file = "results/Bevan_piecewise_model_reftable.rda")

if ( !file.exists("results/Bevan_num_of_periods_posterior.rda") ){

  sumstats = reftable[names(all_sumstats_c14)]
  param = log10(reftable$num_of_periods)
  # hist(log10(param))
  
  RF_log10num_of_periods = regAbcrf(param~., data.frame(param,sumstats),
                               ntree = 5000, paral = TRUE)
  
  posterior_log10num_of_periods = predict(RF_log10num_of_periods, all_sumstats_c14,
                                          training = data.frame(param,sumstats),
                                          paral = TRUE, rf.weights = TRUE) 
  
  save(RF_log10num_of_periods, posterior_log10num_of_periods, file="results/Bevan_num_of_periods_posterior.rda")
}
load(file="results/Bevan_num_of_periods_posterior.rda")
(posterior_log10num_of_periods)
m_hat = 10^posterior_log10num_of_periods$med[1]
m_CI = 10^posterior_log10num_of_periods$quantiles
save(m_hat, m_CI, file="results/Bevan_m_hat.rda")




pdf(file="results/Bevan_posterior_m.pdf", width=10, height=5)
par(mar=c(4.5, 4.5, 1, 1) + 0.1)
breaks= seq(0.3,3,0.05)
hist(log10(reftable$num_of_periods),
     breaks = breaks,
     main="",#expression("Posterior probabilty number of periods ("*italic(m)*")"),
     xlab=expression(log[10]*italic(m)),
     col=adjustcolor( "gray", alpha.f = 0.6),freq=F)
wtd.hist(log10(reftable$num_of_periods),
         breaks = breaks,
         col=PCI_t_blue,
         weight = posterior_log10num_of_periods$weights,
         add=T, freq=F)
box()
dev.off()
