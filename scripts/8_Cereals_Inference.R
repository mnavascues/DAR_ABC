library(abcrf)
source("../../Misc_R_tools/color.R")
source("scripts/sim14c.R")

load(file = "results/num_of_sims.rda")

# load target (i.e. observed) summary statistics
load(file = "results/Cereals_sumstats.rda")

# lead reference tables for piecewise model
load(file = "results/Cereals_piecewise_model_reftable.rda")




#if ( !file.exists("results/Cereals_num_of_periods_posterior.rda") ){

  sumstats = reftable[names(cereals_sumstats_c14)]
  param = reftable$alpha/(reftable$alpha+reftable$beta)
  # hist(log10(param))
  
  RF_beta_mean = regAbcrf(param~., data.frame(param,sumstats),
                           ntree = 5000, paral = TRUE)
  
  posterior_beta_mean = predict(RF_beta_mean, cereals_sumstats_c14,
                                 training = data.frame(param,sumstats),
                                 paral = TRUE, rf.weights = TRUE) 
  
  par(mar=c(4.5, 4.5, 1, 1) + 0.1)
  breaks= seq(-5,3,0.2)
  hist(reftable$alpha/(reftable$alpha+reftable$beta),
       #breaks = breaks,
       xlim=c(0,1),ylim=c(0,10),
       main="",#expression("Posterior probabilty number of periods ("*italic(m)*")"),
       xlab=expression(log[10]*alpha),
       col=adjustcolor( "gray", alpha.f = 0.6),freq=F)
  wtd.hist(reftable$alpha/(reftable$alpha+reftable$beta),
           #breaks = breaks,
           col=PCI_t_blue,
           weight = posterior_beta_mean$weights,
           add=T, freq=F)
  
  
  
  #save(RF_log10num_of_periods, posterior_log10num_of_periods, file="results/Cereals_num_of_periods_posterior.rda")


  
  
  
  
  
    
  
  
load(file="results/Cereals_num_of_periods_posterior.rda")
(posterior_log10num_of_periods)
m_hat = 10^posterior_log10num_of_periods$med[1]
m_CI = 10^posterior_log10num_of_periods$quantiles
save(m_hat, m_CI, file="results/Cereals_m_hat.rda")




pdf(file="results/Cereals_posterior_m.pdf", width=10, height=5)
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

