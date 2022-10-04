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
dates_4_interpolation = seq(10000, 500, -5)

D = with(gisp2, data.frame(approx(YBP,Temperature,xout=dates_4_interpolation)))

plot(D$x,D$y-mean(D$y),type="l", xlim=c(10000,500), ylim=c(-4,4))
lines(D$x,(D$y-mean(D$y))^2,col="gray")
abline(h=0,lty=3)

Brown_N = function(N0, b0, b1, b2, s, D, time_range, time_step=-5){
  time_seq = seq(time_range[1],time_range[2],time_step)
  N = array(NA, length(time_seq))
  for (i in seq_along(time_seq)){
    if (i == 1) {
      N[i] = N0
    }else{
      N[i] = N[i-1]*exp(5*b0)*exp(5*b1)^(N[i-1]^2)*exp(5*b2)^(D[i]^2)*exp(5*rnorm(1,0,s))
    }
  }
  return(N)  
}

popsize = Brown_N(600,0.003,-3.333E-12,-0.003,0.0005, D$y-mean(D$y), time_range=c(10000,500))
plot(dates_4_interpolation,popsize, type="l",  xlim=c(10000,500))

for (i in 1:10){
  theta0 = exp(runif(2,log(0.01),log(0.5)))
  b0 = exp(runif(2,log(0.001),log(0.01)))
  b1 = -exp(runif(2,log(1E-7),log(1E-1)))
  b2 = -exp(runif(2,log(0.001),log(0.01)))
  s = exp(runif(2,log(0.0001),log(0.001)))
  

  
  popsize = Brown_N(theta0,b0,b1,b2,s, D$y-mean(D$y), time_range=c(10000,500))
  plot(dates_4_interpolation,popsize, type="l",  xlim=c(10000,500))
  
}




