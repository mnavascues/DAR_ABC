source("scripts/DARthABC.R")
require(rbenchmark)
require(abcrf)
require(ggplot2)
results_directory = "results/method_evaluation"
dir.create(results_directory)

# MODEL DESCRIPTION
time_range = c(2450, 2400)
errors = 10
lambda_0 = 1
rate = 0.04
time_range_CRA = get_time_range_CRA(time_range, error = max(errors))
lambda_t = get_exponential_lambda_t(lambda_0, rate, time_range)

plot_lambda_model_file = paste0(results_directory, "/plot_lambda_model.pdf")
if (!file.exists(plot_lambda_model_file)){
  pdf(file = plot_lambda_model_file, width = 7, height = 4)
  par(mar=c(4.5, 4.5, 1, 1) + 0.1)
  plot(seq(time_range[1], time_range[2], -1),
       lambda_t,
       xlim = time_range, ylim=c(0,10),
       ylab = expression(lambda),
       xlab = "Years BP (calibrated)",
       type="l", col="grey")
  points(seq(time_range[1], time_range[2], -1),
         lambda_t)
  text(time_range[1],10,"a",cex=1.5)
  #text(2425,8,expression(bold(lambda)),cex=3)
  text(2425,8,expression(lambda[italic(t)]*"=e"^0.04*""^italic(t)),cex=2.5)
  
  dev.off()  
}

  

plot_R_file = paste0(results_directory, "/plot_R.pdf")
if (!file.exists(plot_R_file)){
  DAR = sim_dates_lambda(lambda_t, time_range)
  pdf(file = plot_R_file, width = 7, height = 4)
  par(mar=c(4.5, 4.5, 1, 1) + 0.1)
  h = hist(DAR,
           breaks = seq(time_range[1],time_range[2],-1),
           xlim = time_range,
           xlab = "Years BP (calibrated)",
           main=""
  )
  text(time_range[1],max(h$counts), "b", cex = 1.5)
  text(2425, max(h$counts)-1, expression(bold(italic(R))), cex = 2.5)
  box()
  dev.off()
}
  
  

plot_Rprime_file = paste0(results_directory, "/plot_Rprime.pdf")
if (!file.exists(plot_Rprime_file)){
  DAR_prime = sim_CRA(DAR, calCurves="intcal20", errors)
  pdf(file = plot_Rprime_file, width = 7, height = 4)
  par(mar=c(4.5, 4.5, 1, 1) + 0.1)
  h = hist(DAR_prime$CRA,
           breaks = seq(time_range_CRA[1],time_range_CRA[2],-1),
           xlim = time_range,
           xlab = "Years BP (CRA)",
           main=""
  )
  text(time_range[1],max(h$counts),"c",cex=1.5)
  text(2425,max(h$counts)-0.5,expression(bold(italic("R'"))),cex=2.5)
  box()
  dev.off()
}


  
  


# SPD AS SUMMARY STATISTICS

time_range = c(7000, 5000)
number_of_replicates = 30 # X
errors = 30
true_lambda_0 = 0.01
true_r = 0.003
time_range_CRA = get_time_range_CRA(time_range, error = max(errors))
lambda_t = get_exponential_lambda_t(true_lambda_0, true_r, time_range)
p_t = transform_to_pdf(lambda_t)
true_p_0 = p_t[1]
true_p_f = tail(p_t, 1)
rm(p_t);gc()
true_lambda_f = tail(lambda_t, 1)
n = round(sum(lambda_t))

bench_results_file = paste0(results_directory, "/bench_results.rda")
load(file = bench_results_file)
(bench_results)

rate_min = -0.01
rate_max = 0.01
reftable_file = paste0(results_directory, "/reftable_exponential_prob_model.rda")

OOBplot_p_hist_file = paste0(results_directory, "/OOB_p_hist.pdf")
if (!file.exists(OOBplot_p_hist_file)){
  load(file = reftable_file)
  sumstats = reftable[paste0("hist",seq(time_range_CRA[1], time_range_CRA[2], -1))]
  param    = reftable$rate
  RF_model_p_hist_file = paste0(results_directory, "/RF_model_p_hist.rda")
  load(file = RF_model_p_hist_file)

  pdf(file = OOBplot_p_hist_file, width = 5, height = 4)
  df = data.frame(rate = param, rate_hat = RF_p_hist$model.rf$predictions)
  MSE = round(RF_p_hist$model.rf$prediction.error, 10)
  R2 = round(RF_p_hist$model.rf$r.squared, 3)
  my_breaks_1=c(625,125,25,5,1)
  oobplot = ggplot(df, aes(rate,rate_hat))
  oobplot2d = oobplot + 
    stat_bin2d(bins=100, na.rm = T, show.legend=T) + 
    scale_fill_gradientn(colours=viridis::viridis(32), trans="log", name="counts",
                         breaks = round(my_breaks_1), labels = round(my_breaks_1)
    ) +
    xlab(expression(paste("true", " ", italic(r)))) +
    ylab(expression(paste(hat(italic(r)), " ","(OOB prediction)")))+
    xlim(min(param),max(param)) +
    ylim(min(param),max(param)) +
    geom_abline(intercept = 0, slope = 1, color = "black",  linetype="dashed") +
    annotate("text", x=0.005, y=-0.005, label= paste("MSE =",MSE),  size=4) +
    annotate("text", x=0.005, y=-0.006, label= bquote(paste(italic("R")^"2"," = ",.(R2))),  size=4) +
    annotate("text", x=min(param), y=max(param), label= "b",  size=8) +
    theme_bw() + theme(panel.border = element_rect(colour = "black", fill=NA), 
                       panel.grid.major = element_blank(),
                       panel.grid.minor = element_blank(), axis.line = element_line(colour = "black"),
                       legend.position = "right", axis.text = element_text(size = 12),
                       axis.title=element_text(size=12))
  oobplot2d
  dev.off()
}








OOBplot_p_spd_file = paste0(results_directory, "/OOB_p_spd.pdf")
if (!file.exists(OOBplot_p_spd_file)){
  load(file = reftable_file)
  sumstats = reftable[paste0("spd",seq(time_range[1], time_range[2], -1))]
  param    = reftable$rate
  RF_model_p_spd_file = paste0(results_directory, "/RF_model_p_spd.rda")
  load(file = RF_model_p_spd_file)

  pdf(file=OOBplot_p_spd_file, width = 5, height = 4)
  df = data.frame(rate = param, rate_hat = RF_p_spd$model.rf$predictions)
  MSE = round(RF_p_spd$model.rf$prediction.error, 10)
  R2 = round(RF_p_spd$model.rf$r.squared, 3)
  my_breaks_1=c(625,125,25,5,1)
  oobplot = ggplot(df, aes(rate,rate_hat))
  oobplot2d = oobplot + 
    stat_bin2d(bins=100, na.rm = T, show.legend=T) + 
    scale_fill_gradientn(colours=viridis::viridis(32), trans="log", name="counts",
                         breaks = round(my_breaks_1), labels = round(my_breaks_1)
    ) +
    xlab(expression(paste("true", " ", italic(r)))) +
    ylab(expression(paste(hat(italic(r)), " ","(OOB prediction)")))+
    xlim(min(param),max(param)) +
    ylim(min(param),max(param)) +
    geom_abline(intercept = 0, slope = 1, color = "black",  linetype="dashed") +
    annotate("text", x=0.005, y=-0.005, label= paste("MSE =",MSE),  size=4) +
    annotate("text", x=0.005, y=-0.006, label= bquote(paste(italic("R")^"2"," = ",.(R2))),  size=4) +
    annotate("text", x=min(param), y=max(param), label= "a",  size=8) +
    theme_bw() + theme(panel.border = element_rect(colour = "black", fill=NA), 
                       panel.grid.major = element_blank(),
                       panel.grid.minor = element_blank(), axis.line = element_line(colour = "black"),
                       legend.position = "right", axis.text = element_text(size = 12),
                       axis.title=element_text(size=12))
  oobplot2d
  dev.off()
}  





OOBplot_p_hist2_file = paste0(results_directory, "/OOB_p_hist2.pdf")
if (!file.exists(OOBplot_p_hist2_file)){
  load(file = reftable_file)
  sumstats = reftable[3994:4516]
  param    = reftable$rate
  RF_model_p_hist2_file = paste0(results_directory, "/RF_model_p_hist2.rda")
  load(file = RF_model_p_hist2_file)

  pdf(file=OOBplot_p_hist2_file, width = 5, height = 4)
  df = data.frame(rate = param, rate_hat = RF_p_hist2$model.rf$predictions)
  MSE = round(RF_p_hist2$model.rf$prediction.error, 10)
  R2 = round(RF_p_hist2$model.rf$r.squared, 3)
  my_breaks_1=c(625,125,25,5,1)
  oobplot = ggplot(df, aes(rate,rate_hat))
  oobplot2d = oobplot + 
    stat_bin2d(bins=100, na.rm = T, show.legend=T) + 
    scale_fill_gradientn(colours=viridis::viridis(32), trans="log", name="counts",
                         breaks = round(my_breaks_1), labels = round(my_breaks_1)
    ) +
    xlab(expression(paste("true", " ", italic(r)))) +
    ylab(expression(paste(hat(italic(r)), " ","(OOB prediction)")))+
    xlim(min(param),max(param)) +
    ylim(min(param),max(param)) +
    geom_abline(intercept = 0, slope = 1, color = "black",  linetype="dashed") +
    annotate("text", x=0.005, y=-0.005, label= paste("MSE =",MSE),  size=4) +
    annotate("text", x=0.005, y=-0.006, label= bquote(paste(italic("R")^"2"," = ",.(R2))),  size=4) +
    annotate("text", x=min(param), y=max(param), label= "b",  size=8) +
    theme_bw() + theme(panel.border = element_rect(colour = "black", fill=NA), 
                       panel.grid.major = element_blank(),
                       panel.grid.minor = element_blank(), axis.line = element_line(colour = "black"),
                       legend.position = "right", axis.text = element_text(size = 12),
                       axis.title=element_text(size=12))
  oobplot2d
  dev.off()
}

#################################################################################################

# Comparison of p model vs. lambda model
lambda_min = 0.005
lambda_max = 5
rate_min = -0.005
rate_max = 0.005
num_of_sims = 100000
reftable_file = paste0(results_directory, "/reftable_exponential_2_exponential_models.rda")

OOBplot_lambda_rate_file = paste0(results_directory, "/OOB_lambda_rate.pdf")
if (!file.exists(OOBplot_lambda_rate_file)){
  load(file = reftable_file)
  param    = reftable$rate
  RF_model_lambda_file = paste0(results_directory, "/RF_model_lambda.rda")
  load(file = RF_model_lambda_file)
  
  pdf(file=OOBplot_lambda_rate_file, width = 5, height = 4)
  df = data.frame(rate = param, rate_hat = RF_lambda$model.rf$predictions)
  MSE = round(RF_lambda$model.rf$prediction.error, 10)
  R2 = round(RF_lambda$model.rf$r.squared, 3)
  my_breaks_1=c(625,125,25,5,1)
  oobplot = ggplot(df, aes(rate,rate_hat))
  oobplot2d = oobplot + 
    stat_bin2d(bins=100, na.rm = T, show.legend=T) + 
    scale_fill_gradientn(colours=viridis::viridis(32), trans="log", name="counts",
                         breaks = round(my_breaks_1), labels = round(my_breaks_1)
    ) +
    xlab(expression(paste("true", " ", italic(r)))) +
    ylab(expression(paste(hat(italic(r)), " ","(OOB prediction)")))+
    xlim(-0.005, 0.005) +
    ylim(-0.005, 0.005) +
    geom_abline(intercept = 0, slope = 1, color = "black",  linetype="dashed") +
    annotate("text", x=0.0025, y=-0.0025, label= paste("MSE =",MSE),  size=4) +
    annotate("text", x=0.0025, y=-0.0035, label= bquote(paste(italic("R")^"2"," = ",.(R2))),  size=4) +
    annotate("text", x=-0.005, y=0.005, label= "a",  size=8) +
    theme_bw() + theme(panel.border = element_rect(colour = "black", fill=NA), 
                       panel.grid.major = element_blank(),
                       panel.grid.minor = element_blank(), axis.line = element_line(colour = "black"),
                       legend.position = "right", axis.text = element_text(size = 12),
                       axis.title=element_text(size=12))
  oobplot2d
  dev.off()
}

OOBplot_lambda_0_file = paste0(results_directory, "/OOB_lambda_0.pdf")
if (!file.exists(OOBplot_lambda_rate_file)){
  load(file = reftable_file)
  param    = log10(reftable$lambda_0)
  RF_model_lambda_0_file = paste0(results_directory, "/RF_model_lambda_0.rda")
  load(file = RF_model_lambda_0_file)
  
  pdf(file=OOBplot_lambda_0_file, width = 5, height = 4)
  df = data.frame(rate = param, rate_hat = RF_lambda_0$model.rf$predictions)
  MSE = round(RF_lambda_0$model.rf$prediction.error, 10)
  R2 = round(RF_lambda_0$model.rf$r.squared, 3)
  my_breaks_1=c(625,125,25,5,1)
  oobplot = ggplot(df, aes(rate,rate_hat))
  oobplot2d = oobplot + 
    stat_bin2d(bins=100, na.rm = T, show.legend=T) + 
    scale_fill_gradientn(colours=viridis::viridis(32), trans="log", name="counts",
                         breaks = round(my_breaks_1), labels = round(my_breaks_1)
    ) +
    xlab(expression(paste("true", " ", log[10]*lambda[0]))) +
    ylab(expression(paste(log[10]*hat(lambda)[0], " ","(OOB prediction)")))+
    xlim(log10(0.005), log10(5)) +
    ylim(log10(0.005), log10(5)) +
    geom_abline(intercept = 0, slope = 1, color = "black",  linetype="dashed") +
    annotate("text", x=0, y=-1.5, label= paste("MSE =",MSE),  size=4) +
    annotate("text", x=0, y=-1.7, label= bquote(paste(italic("R")^"2"," = ",.(R2))),  size=4) +
    annotate("text", x=log10(0.005), y=log10(5), label= "b",  size=8) +
    theme_bw() + theme(panel.border = element_rect(colour = "black", fill=NA), 
                       panel.grid.major = element_blank(),
                       panel.grid.minor = element_blank(), axis.line = element_line(colour = "black"),
                       legend.position = "right", axis.text = element_text(size = 12),
                       axis.title=element_text(size=12))
  oobplot2d
  dev.off()
}


OOBplot_p_rate_file = paste0(results_directory, "/OOB_p_rate.pdf")
if (!file.exists(OOBplot_p_rate_file)){
  load(file = reftable_file)
  param    = reftable$rate
  RF_model_p_file = paste0(results_directory, "/RF_model_p.rda")
  load(file = RF_model_p_file)
  
  pdf(file=OOBplot_p_rate_file, width = 5, height = 4)
  df = data.frame(rate = param, rate_hat = RF_p$model.rf$predictions)
  MSE = round(RF_p$model.rf$prediction.error, 10)
  R2 = round(RF_p$model.rf$r.squared, 3)
  my_breaks_1=c(625,125,25,5,1)
  oobplot = ggplot(df, aes(rate,rate_hat))
  oobplot2d = oobplot + 
    stat_bin2d(bins=100, na.rm = T, show.legend=T) + 
    scale_fill_gradientn(colours=viridis::viridis(32), trans="log", name="counts",
                         breaks = round(my_breaks_1), labels = round(my_breaks_1)
    ) +
    xlab(expression(paste("true", " ", italic(r)))) +
    ylab(expression(paste(hat(italic(r)), " ","(OOB prediction)")))+
    xlim(-0.005, 0.005) +
    ylim(-0.005, 0.005) +
    geom_abline(intercept = 0, slope = 1, color = "black",  linetype="dashed") +
    annotate("text", x=0.0025, y=-0.0025, label= paste("MSE =",MSE),  size=4) +
    annotate("text", x=0.0025, y=-0.0035, label= bquote(paste(italic("R")^"2"," = ",.(R2))),  size=4) +
    annotate("text", x=-0.005, y=0.005, label= "c",  size=8) +
    theme_bw() + theme(panel.border = element_rect(colour = "black", fill=NA), 
                       panel.grid.major = element_blank(),
                       panel.grid.minor = element_blank(), axis.line = element_line(colour = "black"),
                       legend.position = "right", axis.text = element_text(size = 12),
                       axis.title=element_text(size=12))
  oobplot2d
  dev.off()
}


OOBplot_p_0_file = paste0(results_directory, "/OOB_p_0.pdf")
if (!file.exists(OOBplot_p_0_file)){
  load(file = reftable_file)
  param    = log10(reftable$p_0)
  RF_model_p_0_file = paste0(results_directory, "/RF_model_p_0.rda")
  load(file = RF_model_p_0_file)
  
  pdf(file=OOBplot_p_0_file, width = 5, height = 4)
  df = data.frame(rate = param, rate_hat = RF_p_0$model.rf$predictions)
  MSE = round(RF_p_0$model.rf$prediction.error, 10)
  R2 = round(RF_p_0$model.rf$r.squared, 3)
  my_breaks_1=c(625,125,25,5,1)
  oobplot = ggplot(df, aes(rate,rate_hat))
  oobplot2d = oobplot + 
    stat_bin2d(bins=100, na.rm = T, show.legend=T) + 
    scale_fill_gradientn(colours=viridis::viridis(32), trans="log", name="counts",
                         breaks = round(my_breaks_1), labels = round(my_breaks_1)
    ) +
    xlab(expression(paste("true", " ", log[10]*italic(p)[0]))) +
    ylab(expression(paste(log[10]*hat(italic(p)[0]), " ","(OOB prediction)")))+
    xlim(min(param), max(param)) +
    ylim(min(param), max(param)) +
    geom_abline(intercept = 0, slope = 1, color = "black",  linetype="dashed") +
    annotate("text", x=-3.5, y=-4.7, label= paste("MSE =",MSE),  size=4) +
    annotate("text", x=-3.5, y=-5, label= bquote(paste(italic("R")^"2"," = ",.(R2))),  size=4) +
    annotate("text", x=min(param),  y=max(param), label= "d",  size=8) +
    theme_bw() + theme(panel.border = element_rect(colour = "black", fill=NA), 
                       panel.grid.major = element_blank(),
                       panel.grid.minor = element_blank(), axis.line = element_line(colour = "black"),
                       legend.position = "right", axis.text = element_text(size = 12),
                       axis.title=element_text(size=12))
  oobplot2d
  dev.off()
}



########################################################################

















num_of_param = 5
num_of_sumstats = (ncol(reftable) - num_of_param) / 2
lambda_sumstats_names = names(reftable)[(1:num_of_sumstats) + num_of_param]
p_sumstats_names = names(reftable)[(1:num_of_sumstats) + num_of_param + num_of_sumstats]

results_RF_model_p_file = paste0(results_directory, "/results_RF_model_p.rda")
load(file=results_RF_model_p_file)
results_RF_model_lambda_file = paste0(results_directory, "/results_RF_model_lambda.rda")
load(file=results_RF_model_lambda_file)


true_r
quantiles = seq(0,1,0.005)
p_model_quantiles = rep(NA,length(quantiles))
lambda_model_quantiles = rep(NA,length(quantiles))
for (i in seq_along(quantiles)){
  p_model_quantiles[i] = sum(results_p$quantiles[,i]>true_r)/length(results_p$expectation)
  lambda_model_quantiles[i] = sum(results_lambda$quantiles[,i]>true_r)/length(results_lambda$expectation)
}

plot(quantiles,p_model_quantiles,
     type="l",
     ylab="coverage probability",
     xlab="credibility interval quantiles",
     lwd=2,col="blue")
lines(quantiles,lambda_model_quantiles,
      lwd=2,col="red", lty=2)
abline(a=0, b=1, lty=3)




#separation = 0.1
#plot(seq_len(number_of_replicates)-separation,
#     results_p$med[seq_len(number_of_replicates)],
#     ylab=expression(italic(r)),
#     xlab=expression("Pseudo-observed data sets with "*italic(n)*"=1343"),
#     ylim=c(0.0025, 0.0035))
#segments(seq_len(number_of_replicates)-separation,
#         results_p$quantiles[seq_len(number_of_replicates),6],
#         seq_len(number_of_replicates)-separation,
#         results_p$quantiles[seq_len(number_of_replicates),196])
#points(seq_len(number_of_replicates)+separation,
#       results_lambda$med[seq_len(number_of_replicates)], col="blue")
#segments(seq_len(number_of_replicates)+separation,
#         results_lambda$quantiles[seq_len(number_of_replicates),6],
#         seq_len(number_of_replicates)+separation,
#         results_lambda$quantiles[seq_len(number_of_replicates),196], col="blue")
#abline(h=true_r,col="grey")



















true_lambda_0
true_p_0

results_RF_model_p_0_file = paste0(results_directory, "/results_RF_model_p_0.rda")
load(file = results_RF_model_p_0_file)
results_RF_model_lambda_0_file = paste0(results_directory, "/results_RF_model_lambda_0.rda")
load(file = results_RF_model_lambda_0_file)

separation = 0.1
plot(seq_len(number_of_replicates)-separation,
     results_p_0$med[seq_len(number_of_replicates)]*true_lambda_0/true_p_0,
     ylab=expression(lambda[0]*" or rescaled "*italic(p)[0]),
     xlab=expression("Pseudo-observed data sets with "*italic(n)*"=1343"),
     ylim=c(0.00, 0.03))
segments(seq_len(number_of_replicates)-separation,
         results_p_0$quantiles[seq_len(number_of_replicates),1]*true_lambda_0/true_p_0,
         seq_len(number_of_replicates)-separation,
         results_p_0$quantiles[seq_len(number_of_replicates),2]*true_lambda_0/true_p_0)
points(seq_len(number_of_replicates)+separation,
       results_lambda_0$med[seq_len(number_of_replicates)], col="blue")
segments(seq_len(number_of_replicates)+separation,
         results_lambda_0$quantiles[seq_len(number_of_replicates),1],
         seq_len(number_of_replicates)+separation,
         results_lambda_0$quantiles[seq_len(number_of_replicates),2], col="blue")
abline(h=true_lambda_0,col="grey")













true_lambda_f
true_p_f

results_RF_model_p_f_file = paste0(results_directory, "/results_RF_model_p_f.rda")
load(file = results_RF_model_p_f_file)
results_RF_model_lambda_f_file = paste0(results_directory, "/results_RF_model_lambda_f.rda")
load(file = results_RF_model_lambda_f_file)

separation = 0.1
plot(seq_len(number_of_replicates)-separation,
     results_p_f$med[seq_len(number_of_replicates)]*true_lambda_f/true_p_f,
     ylab=expression(lambda[f]*" or rescaled "*italic(p)[f]),
     xlab=expression("Pseudo-observed data sets with "*italic(n)*"=1343"),
     ylim=c(3, 5))
segments(seq_len(number_of_replicates)-separation,
         results_p_f$quantiles[seq_len(number_of_replicates),1]*true_lambda_f/true_p_f,
         seq_len(number_of_replicates)-separation,
         results_p_f$quantiles[seq_len(number_of_replicates),2]*true_lambda_f/true_p_f)
points(seq_len(number_of_replicates)+separation,
       results_lambda_f$med[seq_len(number_of_replicates)], col="blue")
segments(seq_len(number_of_replicates)+separation,
         results_lambda_f$quantiles[seq_len(number_of_replicates),1],
         seq_len(number_of_replicates)+separation,
         results_lambda_f$quantiles[seq_len(number_of_replicates),2], col="blue")
abline(h=true_lambda_f,col="grey")









