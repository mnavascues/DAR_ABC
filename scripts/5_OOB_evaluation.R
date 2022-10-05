library(abcrf)
library(hexbin)
library(latticeExtra)
source("../../Misc_R_tools/color.R")
source("scripts/sim14c.R")

load(file = "results/Bevan_num_of_sims.rda")

# load target (i.e. observed) summary statistics
load(file = "results/Bevan_sumstats.rda")

# lead reference tables for piecewise model
load(file = "results/Bevan_piecewise_model_reftable.rda")

load(file = "results/Bevan_time_range_BP.rda")
num_of_periods = round((time_range_BP[1] - time_range_BP[2]) / 100 / 4)
skyline_years = get_time_of_change(num_of_periods, time_range_BP, intervals="regular")

OOB_year = 5500 #skyline_years[round(length(skyline_years)/2)]

sumstats = reftable[names(all_sumstats_c14)]
param_name = paste0("lambda",OOB_year)
param_index = which(names(reftable)==param_name)
param = log10(reftable[param_index])
names(param) = "param"
if ( !file.exists("results/Bevan_piecewise_example_RF_lambda.rda") ){
  RF_lambda = regAbcrf(param~., data.frame(param,sumstats),
                       ntree = 1000, paral = TRUE)
  save(RF_lambda, 
       file="results/Bevan_piecewise_example_RF_lambda.rda")
}
load(file = "results/Bevan_piecewise_example_RF_lambda.rda")

df = data.frame(lambda = param$param, lambda_hat = RF_lambda$model.rf$predictions)
pdf(file="results/Bevan_piecewise_OOB_lambda.pdf", width=6, height=5)
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

param_name = paste0("rate",5893-(5893-5500)/2)
param_index = which(names(reftable)==param_name)
param = reftable[param_index]
names(param) = "param"
if ( !file.exists("results/Bevan_piecewise_example_RF_rate.rda") ){
  RF_rate = regAbcrf(param~., data.frame(param,sumstats),
                     ntree = 1000, paral = TRUE, ncores = 24)
    
  save(RF_rate, 
       file="results/Bevan_piecewise_example_RF_rate.rda")
}
load(file="results/Bevan_piecewise_example_RF_rate.rda")

df = data.frame(rate = param$param, rate_hat = RF_rate$model.rf$predictions)
pdf(file="results/Bevan_piecewise_OOB_rate.pdf", width=5, height=5)
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

