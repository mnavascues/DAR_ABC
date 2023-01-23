library(abcrf)
source("../../Misc_R_tools/color.R")
source("scripts/sim14c.R")

load(file = "results/num_of_sims.rda")

# load target (i.e. observed) summary statistics
load(file = "results/sumstats.rda")

# lead reference tables for dynamic model
load(file = "results/dynamic_model_reftable.rda")
rows2keep = complete.cases(reftable)
reftable = reftable[rows2keep,]



sumstats = reftable[names(all_sumstats_c14)]

if ( !file.exists("results/dynamic_posterior_lambda_0.rda") ){
  param = log10(reftable$lambda_0)
  RF_log10lambda_0 = regAbcrf(param~., data.frame(param,sumstats),
                              ntree = 5000, paral = TRUE)
  posterior_log10lambda_0 = predict(RF_log10lambda_0, all_sumstats_c14,
                                    training = data.frame(param,sumstats),
                                    paral = TRUE, rf.weights = TRUE) 
  save(RF_log10lambda_0, posterior_log10lambda_0,
       file="results/dynamic_posterior_lambda_0.rda")
}
rm(RF_log10lambda_0);gc()
if ( !file.exists("results/dynamic_posterior_beta_0.rda") ){
  param = log10(reftable$b0)
  RF_log10b0 = regAbcrf(param~., data.frame(param,sumstats),
                        ntree = 5000, paral = TRUE)
  posterior_log10b0 = predict(RF_log10b0, all_sumstats_c14,
                              training = data.frame(param,sumstats),
                              paral = TRUE, rf.weights = TRUE) 
  save(RF_log10b0, posterior_log10b0,
       file="results/dynamic_posterior_beta_0.rda")
}
rm(RF_log10b0);gc()
if ( !file.exists("results/dynamic_posterior_beta_1.rda") ){
  param = log10(-reftable$b1)
  RF_log10b1 = regAbcrf(param~., data.frame(param,sumstats),
                        ntree = 5000, paral = TRUE)
  posterior_log10b1 = predict(RF_log10b1, all_sumstats_c14,
                              training = data.frame(param,sumstats),
                              paral = TRUE, rf.weights = TRUE) 
  save(RF_log10b1, posterior_log10b1,
       file="results/dynamic_posterior_beta_1.rda")
}
rm(RF_log10b1);gc()
if ( !file.exists("results/dynamic_posterior_beta_2.rda") ){
  param = log10(-reftable$b2)
  RF_log10b2 = regAbcrf(param~., data.frame(param,sumstats),
                        ntree = 5000, paral = TRUE)
  posterior_log10b2 = predict(RF_log10b2, all_sumstats_c14,
                              training = data.frame(param,sumstats),
                              paral = TRUE, rf.weights = TRUE) 
  save(RF_log10b2, posterior_log10b2,
       file="results/dynamic_posterior_beta_2.rda")
}
rm(RF_log10b2);gc()
if ( !file.exists("results/dynamic_posterior_d.rda") ){
  param = reftable$d
  RF_d = regAbcrf(param~., data.frame(param,sumstats),
                  ntree = 5000, paral = TRUE)
  posterior_d = predict(RF_d, all_sumstats_c14,
                        training = data.frame(param,sumstats),
                        paral = TRUE, rf.weights = TRUE) 
  save(RF_d, posterior_d,
       file="results/dynamic_posterior_d.rda")
}
rm(RF_d);gc()

load("results/dynamic_posterior_lambda_0.rda")
rm(RF_log10lambda_0);gc()
load("results/dynamic_posterior_beta_0.rda")
rm(RF_log10b0);gc()
load("results/dynamic_posterior_beta_1.rda")
rm(RF_log10b1);gc()
load("results/dynamic_posterior_beta_2.rda")
rm(RF_log10b2);gc()
load("results/dynamic_posterior_d.rda")
rm(RF_d);gc()

lambda_0_hat = 10^posterior_log10lambda_0$med[1]
lambda_0_CI = 10^posterior_log10lambda_0$quantiles
beta_0_hat = 10^posterior_log10b0$med[1]
beta_0_CI = 10^posterior_log10b0$quantiles
beta_1_hat = 10^posterior_log10b1$med[1]
beta_1_CI = 10^posterior_log10b1$quantiles
beta_2_hat = 10^posterior_log10b2$med[1]
beta_2_CI = 10^posterior_log10b2$quantiles
d_hat = posterior_d$med[1]
d_CI = posterior_d$quantiles
save(lambda_0_hat, lambda_0_CI,
     beta_0_hat, beta_0_CI,
     beta_1_hat, beta_1_CI,
     beta_2_hat, beta_2_CI,
     d_hat, d_CI, file="results/parameter_estimates_dynamic.rda")

pdf(file="results/posterior_lambda_0_dynamic.pdf", width=4, height=4)
par(mar=c(4.5, 4.5, 1, 1) + 0.1)
breaks= seq(-3,0,0.05)
hist(log10(reftable$lambda_0),
     breaks = breaks,
     main="",#expression("Posterior probabilty number of periods ("*italic(m)*")"),
     xlab=expression(log[10]*lambda[0]),
     ylim=c(0,3),
     col=adjustcolor( "gray", alpha.f = 0.6),freq=F)
wtd.hist(log10(reftable$lambda_0),
         breaks = breaks,
         col=PCI_t_blue,
         weight = posterior_log10lambda_0$weights,
         add=T, freq=F)
box()
dev.off()

pdf(file="results/posterior_beta_0_dynamic.pdf", width=4, height=4)
par(mar=c(4.5, 4.5, 1, 1) + 0.1)
breaks= seq(-5,-2,0.05)
hist(log10(reftable$b0),
     breaks = breaks,
     main="",#expression("Posterior probabilty number of periods ("*italic(m)*")"),
     xlab=expression(log[10]*beta[0]),
     ylim=c(0,2),
     col=adjustcolor( "gray", alpha.f = 0.6),freq=F)
wtd.hist(log10(reftable$b0),
         breaks = breaks,
         col=PCI_t_blue,
         weight = posterior_log10b0$weights,
         add=T, freq=F)
box()
dev.off()

pdf(file="results/posterior_beta_1_dynamic.pdf", width=4, height=4)
par(mar=c(4.5, 4.5, 1, 1) + 0.1)
breaks= seq(-5,-2,0.05)
hist(log10(-reftable$b1),
     breaks = breaks,
     main="",#expression("Posterior probabilty number of periods ("*italic(m)*")"),
     xlab=expression(log[10]*(beta[1])),
     ylim=c(0,2),
     col=adjustcolor( "gray", alpha.f = 0.6),freq=F)
wtd.hist(log10(-reftable$b1),
         breaks = breaks,
         col=PCI_t_blue,
         weight = posterior_log10b1$weights,
         add=T, freq=F)
box()
dev.off()


pdf(file = "results/posterior_beta_2_dynamic.pdf", width = 4, height = 4)
par(mar=c(4.5, 4.5, 1, 1) + 0.1)
breaks = seq(-5, -2, 0.05)
hist(log10(-reftable$b2),
     breaks = breaks,
     main = "",
     xlab = expression(log[10]*(beta[2])),
     ylim = c(0,.8),
     col = adjustcolor("gray", alpha.f = 0.6), freq = F)
wtd.hist(log10(-reftable$b2),
         breaks = breaks,
         col = PCI_t_blue,
         weight = posterior_log10b2$weights,
         add = T, freq = F)
box()
dev.off()

load(file = "results/time_range_BP.rda")
gisp2 = read.table(file = "data/climate/gisp2_temp_accum_alley2000.txt", skip = 75, nrows = 1707 - 75)
names(gisp2) = c("YBP","Temperature")
gisp2$YBP = gisp2$YBP * 1000
dates_4_interpolation = seq(time_range_BP[1], time_range_BP[2] + 1,by = -1)
temperature = with(gisp2, data.frame(approx(YBP,Temperature,xout=dates_4_interpolation)))

pdf(file="results/posterior_d_dynamic.pdf", width=4, height=4)
par(mar=c(4.5, 4.5, 1, 1) + 0.1)
breaks= seq(-4,4,0.1)
hist(reftable$d-mean(temperature$y),
     breaks = breaks,
     main="",#expression("Posterior probabilty number of periods ("*italic(m)*")"),
     xlab=expression(italic(d)),
     ylim=c(0,.6), xlim=c(-2.5,2.5),
     col=adjustcolor( "gray", alpha.f = 0.6),freq=F)
wtd.hist(reftable$d-mean(temperature$y),
         breaks = breaks,
         col=PCI_t_blue,
         weight = posterior_d$weights,
         add=T, freq=F)
box()
dev.off()










