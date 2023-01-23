library(rcarbon)
# Climate data:
# Cite as: Alley, R.B. (2004-03-01): NOAA/WDS Paleoclimatology - GISP2 - 
# Temperature Reconstruction and Accumulation Data. [indicate subset used]. 
# NOAA National Centers for Environmental Information. 
# https://www.ncei.noaa.gov/access/paleo-search/study/2475. Accessed 28 April 2022.

# Column 1: Age (thousand years before present) 
# Column 2: Temperature in central Greenland (degrees C) 
gisp2 = read.table(file = "data/climate/gisp2_temp_accum_alley2000.txt",
                   skip = 75, nrows = 1707 - 75)
names(gisp2) = c("YBP","Temperature")
gisp2$YBP = gisp2$YBP * 1000
oldest_date = max(gisp2$YBP)

startBP = 10000
endBP = 500
dates_4_interpolation = seq(startBP, endBP, -1)



temperature = with(gisp2, data.frame(approx(YBP,Temperature,xout=dates_4_interpolation)))

plot(temperature$x,temperature$y-mean(temperature$y),type="l", xlim=c(startBP,endBP), ylim=c(-5,10))
lines(temperature$x,(temperature$y-mean(temperature$y))^2,col="gray")
abline(h=0,lty=3)


Brown_lambda = function(lambda0, b0, b1, b2, temp_series, ref_temp, time_range, time_step=-1){
  D=temp_series-ref_temp
  time_seq = seq(time_range[1],time_range[2],time_step)
  lambda = array(NA, length(time_seq))
  for (i in seq_along(time_seq)){
    if (i == 1) {
      lambda[i] = lambda0
    }else{
      #lambda[i] = lambda[i-1] * exp(b0) * exp(b1)^(lambda[i-1]^2) * exp(b2)^(D[i]^2) 
      lambda[i] = lambda[i-1] * exp(b0 + b1 * lambda[i-1]^2 + b2 * D[i]^2) 
    }
  }
  return(lambda)  
}

popsize = Brown_lambda(600,0.003,-3.333E-12,-0.003, temperature$y, mean(temperature$y), time_range=c(startBP,endBP))
plot(dates_4_interpolation,popsize, type="l",  xlim=c(startBP,endBP))


popsize = Brown_lambda(0.001,0.003,-1E-4,-0.003, temperature$y, mean(temperature$y), time_range=c(startBP,endBP))
plot(dates_4_interpolation,popsize, type="l",  xlim=c(startBP,endBP))

{
popsize = Brown_lambda(lambda0 = exp(runif(1,log(0.001),log(1))),
                       b0 = exp(runif(1,log(0.00001),log(0.01))),
                       b1 = -exp(runif(1,log(0.00001),log(0.01))),
                       b2 = -exp(runif(1,log(0.00001),log(0.01))),
                       temperature$y, rnorm(1,mean(temperature$y),sd(temperature$y)),# mean(D$y), #
                       time_range=c(startBP,endBP))
plot(dates_4_interpolation,popsize, type="l",  xlim=c(startBP,endBP))
}

popsize = Brown_lambda(lambda0 = lambda_0_hat,
                       b0 = beta_0_hat,
                       b1 = -beta_1_hat,
                       b2 = -beta_2_hat,
                       temperature$y, d_hat,# mean(D$y), #
                       time_range=c(startBP,endBP))
plot(dates_4_interpolation,popsize, type="l",  xlim=c(startBP,endBP))

