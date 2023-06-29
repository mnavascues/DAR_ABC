library(abcrf)
library(ggplot2)
library(viridis)
source("../../Misc_R_tools/color.R")
source("scripts/sim14c.R")

load(file = "results/num_of_sims.rda")

# load target (i.e. observed) summary statistics
load(file = "results/Cereals_sumstats.rda")

# lead reference tables for interdependent model
load(file = "results/Cereals_interdependent_model_reftable.rda")
reftable = reftable[!is.na(reftable$count_A),]
reftable = reftable[reftable$count_A!=1,]
reftable = reftable[!is.na(reftable$count_B),]
reftable = reftable[reftable$count_B!=1,]

load(file = "results/Cereals_time_range_BP.rda")
num_of_periods = round((time_range_BP[1] - time_range_BP[2]) / 100 / 4)
skyline_years = get_time_of_change(num_of_periods, time_range_BP, intervals="regular")

OOB_year = 3250 #skyline_years[round(length(skyline_years)/2)]

sumstats = reftable[names(cereals_sumstats_c14)]
param_name = paste0("lambda",OOB_year)
param_index_A = which(names(reftable)==paste0(param_name,"_A"))
param_index_B = which(names(reftable)==paste0(param_name,"_B"))
param = log10(reftable[param_index_A]+reftable[param_index_B])
names(param) = "param"

load(file = paste0("results/Cereals_interdependent_posterior_log_lambda_",OOB_year,".rda"))

df = data.frame(lambda = param$param, lambda_hat = RF_log_lambda$model.rf$predictions)



pdf(file="results/cereals_OOB_lambda.pdf", width=6, height=5)
MSE = round(RF_log_lambda$model.rf$prediction.error, 2)
R2 = round(RF_log_lambda$model.rf$r.squared, 2)
my_breaks_1=c(500,50,5)
oobplot = ggplot(df, aes(lambda,lambda_hat))
oobplot2d = oobplot + 
  stat_bin2d(bins=100, na.rm = T, show.legend=T) + 
  scale_fill_gradientn(colours=viridis::viridis(32), trans="log", name="counts",
                       breaks = round(my_breaks_1),labels = round(my_breaks_1)
                       ) +
  xlab(expression(paste("true", " ", log[10](lambda[3250])))) +
  ylab(expression(paste(log[10](hat(lambda)[3250]), " ","(OOB prediction)")))+
  xlim(-3,0.3) +
  ylim(-3,0.3) +
  geom_abline(intercept = 0, slope = 1, color = "black",  linetype="dashed") +
  annotate("text", x=0, y=-2, label= paste("MSE =",MSE),  size=4) +
  annotate("text", x=0, y=-2.3, label= paste("R2 =",R2),  size=4) +
  annotate("text", x=-3, y=0.3, label= "a",  size=8) +
  theme_bw() + theme(panel.border = element_rect(colour = "black", fill=NA), 
                     panel.grid.major = element_blank(),
                     panel.grid.minor = element_blank(), axis.line = element_line(colour = "black"),
                     legend.position = "right", axis.text = element_text(size = 12),
                     axis.title=element_text(size=12))
oobplot2d
dev.off()
#


# PI

param = reftable[param_index_A] / (reftable[param_index_A] + reftable[param_index_B])
names(param) = "param"

load(file=paste0("results/Cereals_interdependent_posterior_pi_",OOB_year,".rda"))

pi_hat = exp(RF_logit_pi$model.rf$predictions) / (1 + exp(RF_logit_pi$model.rf$predictions))

df = data.frame(rate = param$param, pi_hat = pi_hat )

pdf(file="results/cereals_OOB_pi.pdf", width=6, height=5)
MSE = round(RF_logit_pi$model.rf$prediction.error, 6)
R2 = round(RF_logit_pi$model.rf$r.squared, 2)
my_breaks_1=c(500,50,5)
oobplot = ggplot(df, aes(rate,pi_hat))
oobplot2d = oobplot + 
  stat_bin2d(bins=100, na.rm = T, show.legend=T) + 
  scale_fill_gradientn(colours=viridis::viridis(32), trans="log", name="counts",
                       breaks = round(my_breaks_1),labels = round(my_breaks_1)
  ) +
  xlab(expression(paste("true", " ", pi[3250] ))) +
  ylab(expression(paste(hat(pi)[3250], " ","(OOB prediction)")))+
  xlim(0,1) +
  ylim(0,1) +
  geom_abline(intercept = 0, slope = 1, color = "black",  linetype="dashed") +
  annotate("text", x=0.6, y=0, label= paste("MSE =",MSE),  size=4) +
  annotate("text", x=0.9, y=0, label= paste("R2 =",R2),  size=4) +
  annotate("text", x=0, y=1, label= "c",  size=8) +
  theme_bw() + theme(panel.border = element_rect(colour = "black", fill=NA), 
                     panel.grid.major = element_blank(),
                     panel.grid.minor = element_blank(), axis.line = element_line(colour = "black"),
                     legend.position = "right", axis.text = element_text(size = 12),
                     axis.title=element_text(size=12))
oobplot2d
dev.off()

