library(abcrf)
library(ggplot2)
library(viridis)
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

OOB_year = 5250 #skyline_years[round(length(skyline_years)/2)]

sumstats = reftable[names(all_sumstats_c14)]
param_name = paste0("lambda",OOB_year)
param_index = which(names(reftable)==param_name)
param = log10(reftable[param_index])
names(param) = "param"

load(file = paste0("results/piecewise_posterior_lambda_",OOB_year,".rda"))

df = data.frame(lambda = param$param, lambda_hat = RF_log10lambda$model.rf$predictions)



pdf(file="results/piecewise_OOB_lambda.pdf", width=6, height=5)
MSE = round(RF_log10lambda$model.rf$prediction.error, 2)
R2 = round(RF_log10lambda$model.rf$r.squared, 2)
my_breaks_1=c(500,50,5)
oobplot = ggplot(df, aes(lambda,lambda_hat))
oobplot2d = oobplot + 
  stat_bin2d(bins=100, na.rm = T, show.legend=T) + 
  scale_fill_gradientn(colours=viridis::viridis(32), trans="log", name="counts",
                       breaks = round(my_breaks_1),labels = round(my_breaks_1)
                       ) +
  xlab(expression(paste("true", " ", log[10](lambda[5250])))) +
  ylab(expression(paste(log[10](hat(lambda)[5250]), " ","(OOB prediction)")))+
  xlim(-3,1) +
  ylim(-3,1) +
  geom_abline(intercept = 0, slope = 1, color = "black",  linetype="dashed") +
  annotate("text", x=0, y=-2, label= paste("MSE =",MSE),  size=4) +
  annotate("text", x=0, y=-2.3, label= paste("R2 =",R2),  size=4) +
  annotate("text", x=-3, y=1, label= "a",  size=8) +
  theme_bw() + theme(panel.border = element_rect(colour = "black", fill=NA), 
                     panel.grid.major = element_blank(),
                     panel.grid.minor = element_blank(), axis.line = element_line(colour = "black"),
                     legend.position = "right", axis.text = element_text(size = 12),
                     axis.title=element_text(size=12))
oobplot2d
dev.off()
#

sumstats = reftable[names(all_sumstats_c14)]
rate_param_names = paste0("rate", skyline_years[1:(length(skyline_years)-1)]+(skyline_years[1]-skyline_years[2])/2)

OOB_r = 14

param_name = rate_param_names[OOB_r]
param_index = which(names(reftable)==param_name)
param = reftable[param_index]
names(param) = "param"

year_interval = c(skyline_years[OOB_r],
                  skyline_years[OOB_r+1])

load(file=paste0("results/piecewise_posterior_",param_name,".rda"))

df = data.frame(rate = param$param, rate_hat = RF_rate$model.rf$predictions)

pdf(file="results/piecewise_OOB_rate.pdf", width=6, height=5)
MSE = round(RF_rate$model.rf$prediction.error, 6)
R2 = round(RF_rate$model.rf$r.squared, 2)
my_breaks_1=c(500,50,5)
oobplot = ggplot(df, aes(rate,rate_hat))
oobplot2d = oobplot + 
  stat_bin2d(bins=100, na.rm = T, show.legend=T) + 
  scale_fill_gradientn(colours=viridis::viridis(32), trans="log", name="counts",
                       breaks = round(my_breaks_1),labels = round(my_breaks_1)
  ) +
  xlab(expression(paste("true", " ", italic(r)[4854-4458] ))) +
  ylab(expression(paste(hat(italic(r))[4854-4458], " ","(OOB prediction)")))+
  xlim(-5.814e-03,5.814e-03) +
  ylim(-5.814e-03,5.814e-03) +
  geom_abline(intercept = 0, slope = 1, color = "black",  linetype="dashed") +
  annotate("text", x=0.003, y=-0.0045, label= paste("MSE =",MSE),  size=4) +
  annotate("text", x=0.003, y=-0.0053, label= paste("R2 =",R2),  size=4) +
  annotate("text", x=-5.814e-03, y=5.814e-03, label= "b",  size=8) +
  theme_bw() + theme(panel.border = element_rect(colour = "black", fill=NA), 
                     panel.grid.major = element_blank(),
                     panel.grid.minor = element_blank(), axis.line = element_line(colour = "black"),
                     legend.position = "right", axis.text = element_text(size = 12),
                     axis.title=element_text(size=12))
oobplot2d
dev.off()

