library(abcrf)
library(hexbin)
library(latticeExtra)
source("../../Misc_R_tools/color.R")
source("scripts/sim14c.R")

load(file = "results/num_of_sims.rda")

# load target (i.e. observed) summary statistics
load(file = "results/sumstats.rda")

# lead reference tables for piecewise model
load(file = "results/piecewise_model_reftable.rda")

load(file = "results/time_range_BP.rda")
num_of_periods = round((time_range_BP[1] - time_range_BP[2]) / 100 / 4)
skyline_years = get_time_of_change(num_of_periods, time_range_BP, intervals="regular")

OOB_year = 5646 #skyline_years[round(length(skyline_years)/2)]

sumstats = reftable[names(all_sumstats_c14)]
param_name = paste0("lambda",OOB_year)
param_index = which(names(reftable)==param_name)
param = log10(reftable[param_index])
names(param) = "param"

load(file = paste0("results/piecewise_posterior_lambda_",OOB_year,".rda"))

df = data.frame(lambda = param$param, lambda_hat = RF_log10lambda$model.rf$predictions)
pdf(file="results/piecewise_OOB_lambda.pdf", width=6, height=5)
hexbinplot(lambda_hat~lambda, data=df,
           ylab=expression(log[10]*hat(lambda)),
           xlab=expression(log[10]*lambda),
           trans=log, inv=exp,
           panel= function(...){
             panel.hexbinplot(...)
             panel.xyplot(x = seq(-3, 1.1,0.01), y = seq(-3, 1.1,0.01),
                          cex = 0.1, col = PCI_blue, fill = PCI_blue)
             panel.xyplot(x = -3, y = 1.1, pch="a", cex = 3, col = "black")
           })
dev.off()



#

sumstats = reftable[names(all_sumstats_c14)]
rate_param_names = paste0("rate", skyline_years[1:(length(skyline_years)-1)]+(skyline_years[1]-skyline_years[2])/2)



param_name = rate_param_names[length(rate_param_names)/2]
param_index = which(names(reftable)==param_name)
param = reftable[param_index]
names(param) = "param"

load(file=paste0("results/piecewise_posterior_",param_name,".rda"))

df = data.frame(rate = param$param, rate_hat = RF_rate$model.rf$predictions)
pdf(file="results/piecewise_OOB_rate.pdf", width=5, height=5)
hexbinplot(rate_hat~rate, data=df,
           ylab=expression(hat(italic(r))),
           xlab=expression(italic(r)),
           xlim=c(-0.006, 0.006),ylim=c(-0.006, 0.006),
           trans=log, inv=exp,
           panel= function(...){
             panel.hexbinplot(...)
             panel.xyplot(x = seq(-0.006, 0.006,0.00001), y = seq(-0.006, 0.006, 0.00001),
                          cex = 0.1, col = PCI_blue, fill = PCI_blue)
             panel.xyplot(x = -0.005, y = 0.005, pch="b", cex = 3, col = "black")
           })
dev.off()

