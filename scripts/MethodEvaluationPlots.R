source("scripts/DARthABC.R")
require(rbenchmark)
require(abcrf)
require(ggplot2)
results_directory = "results/method_evaluation"
dir.create(results_directory)

# color
PCI_blue = rgb(44, 110.4, 148.3, 255, maxColorValue = 255)
PCI_t_blue = rgb(44, 110.4, 148.3, 100, maxColorValue = 255)


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
       type="l", col=PCI_t_blue)
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
           main="",
           col = PCI_t_blue
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
           main="",
           col = PCI_t_blue
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
  my_breaks_1=c(5^5,5^4,5^3,5^2,5^1,5^0)
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
  my_breaks_1=c(5^5,5^4,5^3,5^2,5^1,5^0)
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
  my_breaks_1=c(5^5,5^4,5^3,5^2,5^1,5^0)
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
  my_breaks_1=c(5^5,5^4,5^3,5^2,5^1,5^0)
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
  my_breaks_1=c(5^5,5^4,5^3,5^2,5^1,5^0)
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
  my_breaks_1=c(5^5,5^4,5^3,5^2,5^1,5^0)
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
  my_breaks_1=c(5^5,5^4,5^3,5^2,5^1,5^0)
  oobplot = ggplot(df, aes(rate,rate_hat))
  oobplot2d = oobplot + 
    stat_bin2d(bins=100, na.rm = T, show.legend=T) + 
    scale_fill_gradientn(colours=viridis::viridis(32), trans="log", name="counts",
                         breaks = round(my_breaks_1), labels = round(my_breaks_1)
    ) +
    xlab(expression(paste("true", " ", log[10]*pi[0]))) +
    ylab(expression(paste(log[10]*hat(pi)[0]), " ","(OOB prediction)"))+
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

# Analysis of PODs

results_RF_model_p_file = paste0(results_directory, "/results_RF_model_p.rda")
load(file=results_RF_model_p_file)
results_RF_model_lambda_file = paste0(results_directory, "/results_RF_model_lambda.rda")
load(file=results_RF_model_lambda_file)

plot_95CI_width = paste0(results_directory, "/plot_95CI_width.pdf")
if (!file.exists(plot_95CI_width)){
  pdf(file = plot_95CI_width, width = 5, height = 5)
  par(mar=c(4.5, 4.5, 1, 1) + 0.1)
  hist((results_p$quantiles[,96]-results_p$quantiles[,6])/true_r,
       breaks = seq(0.0,0.17,0.005),
       main="",
       xlab="relative 95%CI width",
       col = "white",
       xlim=c(0.05,0.15))
  hist((results_lambda$quantiles[,96]-results_lambda$quantiles[,6])/true_r,
       breaks = seq(0.0,0.17,0.005), add=T, col = PCI_t_blue)
  text(x=0.05, y=110, label="b", cex=2)
  legend("topright",
         legend = c("model of probabilities","model of counts"),
         fill=c("white",PCI_t_blue))
  box()
  dev.off()  
}
bias_p = mean(true_r-results_p$med)
MSE_p = mean((true_r-results_p$med)^2)
bias_lambda = mean(true_r-results_lambda$med)
MSE_lambda = mean((true_r-results_lambda$med)^2)

theoretical_quantiles = seq(0,1,0.005)
obs_quantile_lambda = obs_quantile_p = rep(NA,length(theoretical_quantiles))
for (i in seq_along(theoretical_quantiles)){

  obs_quantile_lambda[i] =sum(results_lambda$quantiles[,i] > true_r)/length(results_lambda$quantiles[,i])
  obs_quantile_p[i] =sum(results_p$quantiles[,i] > true_r)/length(results_lambda$quantiles[,i])
  
}

plot_CI_coverage = paste0(results_directory, "/plot_CI_coverage.pdf")
if (!file.exists(plot_CI_coverage)){
  pdf(file = plot_CI_coverage, width = 5, height = 5)
  par(mar=c(4.5, 4.5, 1, 1) + 0.1)
  plot(theoretical_quantiles, obs_quantile_lambda, type="l", lty=2, lwd=2, col=PCI_blue,
       xlab = "nominal CI coverage", ylab="actual CI coverage")
  lines(theoretical_quantiles, obs_quantile_p, lwd=2)
  abline(a=0,b=1,lty=3, col="grey")
  text(x=0, y=1, label="a", cex=2)
  legend("bottomright",
         legend = c("model of probabilities","model of counts"),
         col=c("black",PCI_blue),
         lwd=2,
         lty=c(1,2))
  dev.off()  
}
