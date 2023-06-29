require(rcarbon)
require(weights)
library(extraDistr)

logit <- function(p) {
  log(p / (1 - p))
}
inv_logit <- function(x) {
  exp(x) / (1 + exp(x))
}


# get exponential model from prior
sample_exponential_parameters_from_priors = function(lambda_min, lambda_max, time_range, force=NULL){
  lambda = exp(runif(2,log(lambda_min),log(lambda_max)))
  if (!is.null(force)){
    if (force=="expansion") lambda = sort(lambda)
    if (force=="decline") lambda = sort(lambda, decreasing = T)
  }
  r = (log(lambda[2]) - log(lambda[1])) / (time_range[1] - time_range[2])
  return(data.frame(lambda_0 = lambda[1],
                    lambda_f = lambda[2],
                    rate     = r))
}
# get lambda values for each year in the interval (lambda_t)
get_exponential_lambda_t = function(lambda_0, r, time_range){
  t = (time_range[1] - time_range[2])
  lambda_t = c(lambda_0, lambda_0 * exp(r * seq_len(t)))
  return(lambda_t)
}
# get exponential model from prior
sample_exponential_piecewise_parameters_from_priors = function(lambda_min, lambda_max,
                                                               time_range,
                                                               num_of_periods,
                                                               intervals="dirichlet"){

  lambda_values = sample_lambda_values_piecewise_model(num_of_periods+1, lambda_min, lambda_max)
  time_of_change = get_time_of_change(num_of_periods, time_range, intervals = intervals)
  growth_rates = get_growth_rates(num_of_periods, log(lambda_values), time_of_change)
  names(lambda_values) = paste0("lambda_",seq_len(num_of_periods+1))
  names(time_of_change) = c("t_0",paste0("t_",seq_len(num_of_periods-1)),"t_f")
  names(growth_rates) = paste0("r_",seq_len(num_of_periods))
  return(data.frame(t(c(lambda_values,time_of_change,growth_rates))))
}
get_piecewise_exponential_lambda_t = function(lambda_values, time_of_change, growth_rates){
  if ((length(lambda_values)==length(time_of_change)) & length(lambda_values)==(1+length(growth_rates))){
    num_of_periods = length(growth_rates)
    lambda_t = numeric()
    for (i in 1:num_of_periods){
      lambda_t = c(lambda_t,
                   lambda_values[i] * exp(growth_rates[i]*seq_len(as.integer(abs(time_of_change[i+1]-time_of_change[i]))) ) )
    }
    lambda_t = c(lambda_t, lambda_values[num_of_periods+1])
    return(lambda_t)   
  }else{
    stop("invalid input size,must be: length(lambda_values)==length(time_of_change)==1+length(growth_rates)")
  }
}
sample_lambda_values_piecewise_model = function(n, lambda_min, lambda_max){
  lambda_values = array(NA, n)
  lambda_values[1] = exp(runif(1, log(lambda_min), log(lambda_max)))
  for (i in 2:(n)){
    alpha = runif(1, log(0.1), log(10))
    lambda_values[i] = exp(max(min(log(lambda_values[i-1]) + alpha, log(lambda_max)), log(lambda_min)))
  }
  return(lambda_values)
}
get_time_of_change = function(num_of_periods, time_range, intervals = "dirichlet"){
  if (intervals == "dirichlet"){
    time_of_change = round(time_range[1] - cumsum(extraDistr::rdirichlet(1, rep(1, num_of_periods))) * (time_range[1] - time_range[2]))
    time_of_change = c(time_range[1], time_of_change)
    while (length(time_of_change) != length(unique(time_of_change))){
      time_of_change = round(time_range[1] - cumsum(extraDistr::rdirichlet(1, rep(1, num_of_periods))) * (time_range[1] - time_range[2]))
      time_of_change = c(time_range[1], time_of_change)
    }
  }else if (intervals == "regular"){
    time_of_change = round(seq(time_range[1], time_range[2], - (time_range[1] - time_range[2]) / num_of_periods) )
  }else{
    stop("invalid intervals value. Values accepted: 'dirichlet' or 'regular'")
  }
  return(time_of_change)
}
get_growth_rates = function(num_of_periods, lambda_values, time_of_change){
  growth_rates = rep(NA, num_of_periods)
  for (i in seq_len(num_of_periods)){
    growth_rates[i] = (lambda_values[i + 1] - lambda_values[i]) / abs(time_of_change[i + 1] - time_of_change[i])
  }
  return(growth_rates)
}
# normalize lambda_t so that sum(p_t)=1
transform_to_pdf = function(lambda_t){
  p_t = lambda_t/sum(lambda_t)
  return(p_t)
}
# simulate number of samples and their age in the archaeological record giving a model with lambda_t
sim_dates_lambda = function(lambda_t, time_range){
  t_values = seq(time_range[1], time_range[2], -1)
  if (length(lambda_t)!=length(t_values)) {
    stop("lambda_t and t_values have to be of the same length")
  }
  n_t   = rpois(length(lambda_t), lambda_t)
  dates = rep(t_values, n_t)
  return(dates)
}
# simulate age of a fixed number of samples in the archaeological record giving a model with p_t
sim_dates_from_pdf = function(n, p_t, time_range){
  t_values = seq(time_range[1], time_range[2], -1)
  if (length(p_t)!=length(t_values)) {
    stop("p_t and t_values have to be of the same length")
  }
  dates = sample(x = t_values, size = n, replace = T, prob = p_t)
  return(dates)
}
# simulate Calibrated Radiocarbon Age (14C measure)
sim_CRA = function(dates, calCurves="intcal20", errors){
  if (length(errors)==1){
    errors = rep(errors, times = length(dates))
  }else{
    errors = sample(errors, size = length(dates), replace = TRUE)
  }
  x = rcarbon::uncalibrate(dates, CRAerrors = errors, calCurves = calCurves)
  return(data.frame(CRA   = x$rCRA,
                    error = x$rError))
}
# get the range of CRA expected for 99.73% (3*sigma) of ages from time range in calibrated years
get_time_range_CRA =  function(time_range, calCurves = "intcal20", error = 0){
  time_range_CRA = rcarbon::uncalibrate(time_range, calCurves = calCurves)
  time_range_CRA = c(time_range_CRA$ccCRA[1] + 3 * (time_range_CRA$ccError[1] + error),
                     time_range_CRA$ccCRA[2] - 3 * (time_range_CRA$ccError[2] + error))
  return(time_range_CRA)
}
# calculate summary statistics as SPD
get_sumstats_spd = function(dates_CRA, time_range, calCurves = "intcal20", bins = NA, runm = 100){
  if (length(dates)>0){
    dates_calibrated = rcarbon::calibrate(x         = dates_CRA$CRA,
                                          errors    = dates_CRA$error,
                                          calCurves = calCurves)
    spd = rcarbon::spd(dates_calibrated, timeRange = time_range, bins = bins, runm = runm)
    stats = as.data.frame(t(spd$grid$PrDens))
    names(stats) = paste0("spd",spd$grid$calBP)
  }else{
    stats = as.data.frame(t(rep(0, time_range[1] - time_range[2] + 1)))
    names(stats) = paste0("spd",seq(time_range[1], time_range[2], -1))
  }
  return(as.data.frame(stats))
}
# calculate summary statistics
get_sumstats_hist = function(dates_CRA, time_range_CRA, w = NULL){
  x = dates_CRA$CRA
  dates_on_range = intersect(which(x<time_range_CRA[1]), which(x>time_range_CRA[2]))
  x = x[dates_on_range]
  if (is.null(w)){
    n = length(x)
  } else{
    w = w[dates_on_range]
    n = sum(w)
  }
  if (n>0){
    breaks = c(seq(time_range_CRA[1],time_range_CRA[2],-1) + 0.5, time_range_CRA[2]-0.5)
    stats = as.data.frame(t(rev(weights::wtd.hist(x, breaks=breaks, plot=F, weight = w)$counts)))
    names(stats) = paste0("hist",seq(time_range_CRA[1], time_range_CRA[2], -1))
  }else{
    stats = as.data.frame(t(rep(0, time_range_CRA[1] - time_range_CRA[2] + 1)))
    names(stats) = paste0("hist",seq(time_range_CRA[1], time_range_CRA[2], -1))
  }
  return(stats)
}

# calculate summary statistics
get_sumstats = function(dates_CRA, time_range_CRA, window = c(50,100), w = NULL){
  x = dates_CRA$CRA
  dates_on_range = intersect(which(x<time_range_CRA[1]), which(x>time_range_CRA[2]))
  x = x[dates_on_range]
  if (is.null(w)){
    n = length(x)
  } else{
    w = w[dates_on_range]
    n = sum(w)
  }
  stats = data.frame(n=n)
  for (i in seq_along(window)){
    window_size = window[i]
    breaks = seq(time_range_CRA[1],time_range_CRA[2]-window_size,-window_size)
    dates = breaks[-1]+window_size/2
    num_of_bins = length(breaks)-1
    if (n>0){
      h = rev(weights::wtd.hist(x, breaks=breaks, plot=F, weight = w)$counts)
    }else{
      h = rep(0,num_of_bins)
    }
    dh = diff(h)
    names(h) = paste0("h",dates,"_w",window[i])
    names(dh) = paste0("dh",dates[-1],"_w",window[i])
    stats = cbind(stats,as.data.frame(t(h)),as.data.frame(t(dh)))
  }
  return(stats)
}
# calculate summary statistics for 2 DAR
get_sumstats_correlation = function(ssA, ssB, time_range_CRA, window = c(25, 50, 100)){
  stats = data.frame(ratioDARs=ssA$n/ssB$n)
  for (i in seq_along(window)){
    window_size = window[i]
    breaks = seq(time_range_CRA[1],time_range_CRA[2]-window_size,-window_size)
    dates = breaks[-1]+window_size/2
    num_of_bins = length(breaks)-1
    HA = as.vector(t(ssA[paste0("h",dates,"_w",window[i])]))
    HB = as.vector(t(ssB[paste0("h",dates,"_w",window[i])]))
    DHA = as.vector(t(ssA[paste0("dh",dates[-1],"_w",window[i])]))
    DHB = as.vector(t(ssB[paste0("dh",dates[-1],"_w",window[i])]))

    corH = cor(HA, HB, method = "pearson")
    covH = cov(HA, HB, method = "pearson")
    corH_lag1 = cor(HA[-1], HB[-length(HB)], method = "pearson")
    covH_lag1 = cov(HA[-1], HB[-length(HB)], method = "pearson")
    corH_lag2 = cor(HB[-1], HA[-length(HA)], method = "pearson")
    covH_lag2 = cov(HB[-1], HA[-length(HA)], method = "pearson")
    
    corDH = cor(DHA, DHB, method = "pearson")
    covDH = cov(DHA, DHB, method = "pearson")
    corDH_lag1 = cor(DHA[-1], DHB[-length(DHB)], method = "pearson")
    covDH_lag1 = cov(DHA[-1], DHB[-length(DHB)], method = "pearson")
    corDH_lag2 = cor(DHB[-1], DHA[-length(DHA)], method = "pearson")
    covDH_lag2 = cov(DHB[-1], DHA[-length(DHA)], method = "pearson")
    
    temp = data.frame(corH, covH, corH_lag1, covH_lag1, corH_lag2, covH_lag2,
             corDH, covDH, corDH_lag1, covDH_lag1, corDH_lag2, covDH_lag2)
    names(temp) = paste0(c("corH", "covH", "corH_lag1", "covH_lag1", "corH_lag2", "covH_lag2",
                           "corDH", "covDH", "corDH_lag1", "covDH_lag1", "corDH_lag2", "covDH_lag2"),
                         "_w", window[i])
    stats = cbind(stats,temp)
  } 
  return(stats)
}



interpret_K = function(k){
  interpretation = "Strength of evidence is:"
  if (k < 10^0) {          interpretation = paste(interpretation, "negative")
  } else if (k < 10^0.5) { interpretation = paste(interpretation, "barely worth mentioning")
  } else if (k < 10^1)   { interpretation = paste(interpretation, "substantial")
  } else if (k < 10^1.5) { interpretation = paste(interpretation, "strong")
  } else if (k < 10^2)   { interpretation = paste(interpretation, "very strong")
  } else if (k >= 10^2)  { interpretation = paste(interpretation, "decisive")
  }
  print(interpretation)
}


