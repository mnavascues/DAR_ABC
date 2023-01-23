source("scripts/sim14c.R")
source("../../Misc_R_tools/color.R")
load(file = "results/spd.rda")
load(file = "results/time_range_BP.rda")


pdf(file="results/SPD_plus_lambda_models.pdf", width=10, height=5)
par(mar=c(4.5, 4.5, 1, 1) + 0.1)
plot(allspd$grid$calBP, allspd$grid$PrDens, xlim = time_range_BP, ylim = c(0.09,5), log = "y",
     type="l", xlab="Years cal BP", ylab=expression(lambda), col="grey", lwd=2)

time_seq = seq(time_range_BP[1], time_range_BP[2] + 1,by = -1)

# DYNAMIC MODEL
load(file="results/parameter_estimates_dynamic.rda")
gisp2 = read.table(file = "data/climate/gisp2_temp_accum_alley2000.txt", skip = 75, nrows = 1707 - 75)
names(gisp2) = c("YBP","Temperature")
gisp2$YBP = gisp2$YBP * 1000
dates_4_interpolation = seq(time_range_BP[1], time_range_BP[2] + 1,by = -1)
temperature = with(gisp2, data.frame(approx(YBP,Temperature,xout=dates_4_interpolation)))
D = temperature$y - d_hat
lambda_t = array(NA, length(time_seq))
for (i in seq_along(time_seq)){
  if (i == 1) {
    lambda_t[i] = lambda_0_hat
  }else{
    lambda_t[i] = lambda_t[i-1] * exp(beta_0_hat - beta_1_hat * lambda_t[i-1]^2 - beta_2_hat * D[i]^2) 
  }
}
lines(time_seq, lambda_t, col=cbPalette1[2], lwd=2)

# LOGISTIC MODEL
load(file="results/parameter_estimates_logistic.rda")
lambda_t = K_hat*lambda_0_hat / (lambda_0_hat + (K_hat - lambda_0_hat)*exp(-r_hat*seq_along(time_seq)) )
lines(time_seq, lambda_t, col=cbPalette1[3], lwd=2)

# PIECEWISE EXPONENTIAL MODEL
load(file="results/parameter_estimates_piecewise.rda")
num_of_periods = round((time_range_BP[1] - time_range_BP[2]) / 100 / 4)
skyline_years = get_time_of_change(num_of_periods, time_range_BP, intervals="regular")
lines(skyline_years, lambda_hat, col = cbPalette1[4], lwd = 2)

legend(4000,0.5,
       legend=c("SPD","Dynamic","Logistic","Piecewise exponential"),
       col=c("grey",cbPalette1[2:4]),
       lwd=2)
dev.off()

