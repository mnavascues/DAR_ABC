require(rcarbon)
require(weights)
library(extraDistr)

#' Sample from a prior distribution the values of the exponential model parameters
#'
#'  This function is used when an exponential model is fitted to the observed
#'  radiocarbon data. The parameters sampled from the prior distribution are the 
#'  lambda_0 at time t_0 and lambda_f at time t_f (time_range). The prior is a 
#'  log-uniform distribution between lambda_min and lambda_max. Given lambda_0,
#'  lambda_f and time_range it calculate the growth rate (r) of the exponential
#'  model.
#'
#' @param lambda_min Numeric. Lower limit for prior
#' @param lambda_max Numeric. Upper limit for prior
#' @param time_range Numeric. A vector of two values, c(t_0, t_f), expressed in 
#' calibrated years before present (i.e. before 1950). 
#' @param force String. Values "expansion" or "decline" force the dynamics of 
#' lambda to be an expansion or a decline
#'
#' @return data.frame with values of lambda_0, lambda_f and rate
#'
#' @examples
#' # General example
#' sample_exponential_parameters_from_priors(0.01, 10, c(8000, 4000))
#' 
#' # Example of forced expansion
#' sample_exponential_parameters_from_priors(0.01, 10, c(8000, 4000), force="expansion")
#' 
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

#' Calculate lambda for each year according to the exponential model
#' 
#'  This function calculates the values of lambda for each year in the given 
#'  time interval according to an exponential model given by the initial value 
#'  of lambda and the rate of exponential change
#'  
#' @param lambda_0 Numeric. initial value of lambda
#' @param r Numeric. rate of exponential change
#' @param time_range Numeric. A vector of two values, c(t_0, t_f), expressed in 
#' calibrated years before present (i.e. before 1950). 
#' 
#' @return Numeric. vector of size t_0-t_f+1 with the value of lambda for each 
#' year between t_0 and t_f (including t_0 and t_f)
#' 
#' @examples
#' lambda_t = get_exponential_lambda_t(0.01, 0.01, c(8000, 4000))
#' plot(lambda_t, type="l")
#' 
#' # Example with parameters taken from a prior
#' time_range = c(8000, 4000)
#' exponential_parameters = sample_exponential_parameters_from_priors(0.01, 10, time_range)
#' lambda_t = get_exponential_lambda_t(exponential_parameters$lambda_0,
#'                                     exponential_parameters$rate, time_range)
#' plot(lambda_t, type="l")
#' 
get_exponential_lambda_t = function(lambda_0, r, time_range){
  t = (time_range[1] - time_range[2])
  lambda_t = c(lambda_0, lambda_0 * exp(r * seq_len(t)))
  return(lambda_t)
}

#' Sample from a prior distribution the values of the logistic model parameters
#'
#'  This function is used when an logistic model is fitted to the observed
#'  radiocarbon data. The parameters sampled from the prior distribution are the 
#'  lambda_0 at time t_0 and lambda_f at time t_f (time_range), and the carrying
#'  capacity, K. The prior is a log-uniform distribution between lambda_min and 
#'  lambda_max for parameters lambda_0 and lambda_f. The prior for K is 
#'  determined by the values lambda_0 and lambda_f. Given lambda_0, lambda_f, K
#'  and time_range it calculate the growth rate (r) of the logistic model.
#'
#' @param lambda_min Numeric. Lower limit for prior
#' @param lambda_max Numeric. Upper limit for prior
#' @param time_range Numeric. A vector of two values, c(t_0, t_f), expressed in 
#' calibrated years before present (i.e. before 1950). 
#' @param force String. Values "expansion" or "decline" force the dynamics of 
#' lambda to be an expansion or a decline
#'
#' @return data.frame with values of lambda_0, lambda_f, K and rate
#'
#' @examples
#' # General example
#' sample_logistic_parameters_from_priors(0.01, 10, c(8000, 4000))
#' 
#' # Example of forced expansion
#' sample_logistic_parameters_from_priors(0.01, 10, c(8000, 4000), force="expansion")
#' 
sample_logistic_parameters_from_priors = function(lambda_min, lambda_max, time_range, force=NULL){
  lambda = exp(runif(2,log(lambda_min),log(lambda_max)))
  if (!is.null(force)){
    if (force=="expansion") lambda = sort(lambda)
    if (force=="decline") lambda = sort(lambda, decreasing = T)
  }
  K = max(lambda) + exp(runif(1,log(lambda_min),log(lambda_max)))
  while (K*lambda[1]/lambda[2]-lambda[1] <= 0 | K-lambda[1] <= 0){
    lambda = exp(runif(2,log(lambda_min),log(lambda_max)))
    if (expansion) lambda = sort(lambda)
    K = max(lambda) + exp(runif(1,log(lambda_min),log(lambda_max)))
  }
  r = (log(K-lambda[1]) - log(K*lambda[1]/lambda[2]-lambda[1])) / (time_range[1] - time_range[2])
  
  
  
  r = (log(lambda[2]) - log(lambda[1])) / (time_range[1] - time_range[2])
  return(data.frame(lambda_0 = lambda[1],
                    lambda_f = lambda[2],
                    K        = K,
                    rate     = r))
}

#' Calculate lambda for each year according to the logistic model
#' 
#'  This function calculates the values of lambda for each year in the given 
#'  time interval according to an logistic model given by the initial value 
#'  of lambda and the rate of logistic change
#'  
#' @param lambda_0 Numeric. initial value of lambda
#' @param K Numeric. carrying capacity for lambda
#' @param r Numeric. rate of logistic change
#' @param time_range Numeric. A vector of two values, c(t_0, t_f), expressed in 
#' calibrated years before present (i.e. before 1950). 
#' 
#' @return Numeric. vector of size t_0-t_f+1 with the value of lambda for each 
#' year between t_0 and t_f (including t_0 and t_f)
#' 
#' @examples
#' lambda_t = get_logistic_lambda_t(0.01, 10, 0.01, c(8000, 4000))
#' plot(lambda_t, type="l")
#' 
#' # Example with parameters taken from a prior
#' time_range = c(8000, 4000)
#' logistic_parameters = sample_logistic_parameters_from_priors(0.01, 10, time_range)
#' lambda_t = get_logistic_lambda_t(logistic_parameters$lambda_0,
#'                                  logistic_parameters$K,
#'                                  logistic_parameters$rate, time_range)
#' plot(lambda_t, type="l")
#' 
get_logistic_lambda_t = function(lambda_0, K, r, time_range){
  t = (time_range[1] - time_range[2])
  lambda_t =  c(lambda_0, K*lambda_0 / (lambda_0 + (K - lambda_0)*exp(-r*seq_len(t)) ))
  return(lambda_t)
}

#' Sample from a prior distribution the values of the piecewise exponential model parameters
#'
#'  This function is used when an piecewise exponential model is fitted to the 
#'  observed radiocarbon data. The parameters sampled from the prior 
#'  distribution are the number of periods and the lambda values at the start 
#'  and end of each period. The prior is a log-uniform distribution between 
#'  lambda_min and lambda_max. The correlation of lambda values between 
#'  consecutive periods is controlled by parameter f. Given lambda values at the
#'  start and end of each period, it calculate the growth rate (r) of the 
#'  exponential model for each period.
#'
#' @param lambda_min Numeric. Lower limit for prior
#' @param lambda_max Numeric. Upper limit for prior
#' @param time_range Numeric. A vector of two values, c(t_0, t_f), expressed in 
#' calibrated years before present (i.e. before 1950).
#' @param num_of_periods Integer. number of period in which the time interval
#' will be divided
#' @param intervals String or vector of values. Determines in which way the time
#' interval will be divided in the number of periods. The default ("dirichlet")
#' will take periods of random length following a Dirichlet(1) distribution.
#' Option "regular" will make equal (or approximately equal) periods. Finally, 
#' a numeric vector can be provided for user defined period lengths.  
#' @param f Numeric >=1. parameter to control the correlation between values of
#' lambda at consecutive time intervals. Low values make strong correlation
#' (1 makes a constant model).
#' 
#' @return data.frame with values of lambda at the times defining the periods
#' (num_of_periods+1 values), the times defining the periods (num_of_periods+1
#' values) and the exponential change rate within each interval (num_of_periods
#' values)
#'
#' @examples
#' # General example
#' sample_exponential_piecewise_parameters_from_priors(0.01, 10, c(8000, 4000), 4)
#' 
#' # Example of regular divided periods
#' sample_exponential_piecewise_parameters_from_priors(0.01, 10, c(8000, 4000), 4, intervals="regular")
#' 
#' # Example of user defined periods
#' sample_exponential_piecewise_parameters_from_priors(0.01, 10, c(8000, 4000), 4, intervals=c(8000,6000,5500,5000,4000))
#' 
sample_exponential_piecewise_parameters_from_priors = function(lambda_min, lambda_max,
                                                               time_range,
                                                               num_of_periods,
                                                               intervals="dirichlet",
                                                               f=10){

  lambda_values = sample_lambda_values_piecewise_model(num_of_periods+1, lambda_min, lambda_max, f)
  if (length(intervals)==1){
      time_of_change = get_time_of_change(num_of_periods, time_range, intervals = intervals)
  }else if (all(is.numeric(intervals),
                length(intervals)==num_of_periods+1,
                intervals[1]==time_range[1],
                intervals[length(intervals)]==time_range[length(time_range)])){
    time_of_change = intervals
  }else{
    stop("invalid intervals value. Values accepted: 'dirichlet, 'regular' or vector of times of change")
  }
  growth_rates = get_growth_rates(num_of_periods, log(lambda_values), time_of_change)
  names(lambda_values) = paste0("lambda_",seq_len(num_of_periods+1))
  names(time_of_change) = c("t_0",paste0("t_",seq_len(num_of_periods-1)),"t_f")
  names(growth_rates) = paste0("r_",seq_len(num_of_periods))
  return(data.frame(t(c(lambda_values,time_of_change,growth_rates))))
}

#' Calculate lambda for each year according to the piecewise exponential model
#' 
#'  This function calculates the values of lambda for each year in the given 
#'  time interval according to an piecewise exponential model given by the 
#'  values of lambda at the start time of each periods and the provided rate of
#'  change
#'  
#' @param lambda_values Numeric vector. value of lambda at each of the times 
#' defining the periods
#' @param time_of_change Numeric vector. times defining the periods expressed in 
#' calibrated years before present (i.e. before 1950).
#' @param growth_rates Numeric vector. exponential rate change within each period
#'  
#' @return Numeric. vector of size t_0-t_f+1 with the value of lambda for each 
#' year between first year of the first period and the last year of the last 
#' period.
#' 
#' @examples
#' lambda_t = get_piecewise_exponential_lambda_t(c(0.3, 0.2, 0.5),
#'                                               c(8000,5000,4000),
#'                                               c(-6e-5,0.001))
#' plot(lambda_t, type="l")
#' 
#' # Example with parameters taken from a prior
#' time_range = c(8000, 4000)
#' num_of_periods = 4
#' piecewise_exponential_parameters = sample_exponential_piecewise_parameters_from_priors(0.01, 10, time_range, num_of_periods)
#' lambda_t = get_piecewise_exponential_lambda_t(as.numeric(piecewise_exponential_parameters[seq_len(num_of_periods+1)]),
#'                                               as.integer(piecewise_exponential_parameters[1+num_of_periods+seq_len(num_of_periods+1)]),
#'                                               as.numeric(piecewise_exponential_parameters[2+2*num_of_periods+seq_len(num_of_periods)]))
#' plot(lambda_t, type="l")
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

#' Use a correlated prior distribution to sample several values of lambda
#' 
#'  The first lambda values is taken from a log-uniform prior and consecutive
#'  values are taken conditional to the last lambda value obtained. The prior is
#'  a log-uniform distribution between lambda_min and lambda_max. The 
#'  correlation of lambda values between consecutive periods is controlled by 
#'  parameter f.
#'
#' @param n Numeric. Number of values of lambda to be sampled
#' @param lambda_min Numeric. Lower limit for prior
#' @param lambda_max Numeric. Upper limit for prior
#' @param f Numeric >=1. parameter to control the correlation between 
#' consecutive values of lambda. Low values make strong correlation
#' (1 makes a constant value of lambda).
#' 
#' @return numeric vector with n values of lambda
#'
#' @examples
#' sample_lambda_values_piecewise_model(3,0.01,10)
#' 
#' # Example with strong correlation
#' sample_lambda_values_piecewise_model(3,0.01,100,f=1.2)
#' 
#' # Example with little correlation
#' sample_lambda_values_piecewise_model(3,0.01,100,f=10000)
#' 
sample_lambda_values_piecewise_model = function(n, lambda_min, lambda_max, f=10){
  lambda_values = array(NA, n)
  lambda_values[1] = exp(runif(1, log(lambda_min), log(lambda_max)))
  for (i in 2:(n)){
    alpha = runif(1, log(1/f), log(f))
    lambda_values[i] = exp(max(min(log(lambda_values[i-1]) + alpha, log(lambda_max)), log(lambda_min)))
  }
  return(lambda_values)
}

#' Divides a time interval in a given number of sub-intervals
#' 
#'  The interval is divided randomly according to a Dirichlet distribution D(1) 
#'  or it equal (or nearly-equal) sized intervals.
#'
#' @param num_of_periods Integer. Number of sub-intervals
#' @param time_range Numeric. A vector of two values, c(t_0, t_f)
#' @param intervals Can take only two values: "dirichlet" or "regular"
#' 
#' @return numeric vector the times at that define the start and end of the 
#' sub-intervals
#'
#' @examples
#' # Random sized sub-intervals
#' get_time_of_change(10, c(8000, 4000))
#' 
#' # Equally sized sub-intervals
#' get_time_of_change(10, c(8000, 4000), intervals = "regular")
#' 
#' # Nearly-equally sized sub-intervals
#' get_time_of_change(10, c(8000, 4002), intervals = "regular")
#' 
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

#' Calculates exponential change rate given initial and end values
#' 
#'  Given the initial and end values, and the starting and end times of the time
#'  intervals, it calculates the corresponding exponential rates of change.
#'
#' @param num_of_periods Integer. Number of consecutive periods in which
#' calculate the exponential rate
#' @param lambda_values Numeric. A vector of two or more values (with size
#' num_of_periods + 1) with the initial and final values (consecutive intervals 
#' have the same final and initial values)
#' @param time_of_change Numeric. Vector with the values of starting and ending
#' times for intervals (size num_of_periods + 1)
#' 
#' @return numeric vector (size num_of_periods) with the rates
#'
#' @examples
#' get_growth_rates(2, c(0.1, 0.001, 0.05), c(8000, 7000, 6000))
#' get_growth_rates(3, c(0.1, 0.001, 0.05, 0.04), c(8000, 7000, 6000, 5000))
#' 
get_growth_rates = function(num_of_periods, lambda_values, time_of_change){
  growth_rates = rep(NA, num_of_periods)
  for (i in seq_len(num_of_periods)){
    growth_rates[i] = (lambda_values[i + 1] - lambda_values[i]) / abs(time_of_change[i + 1] - time_of_change[i])
  }
  return(growth_rates)
}

#' Normalization to unity
#' 
#'  Normalizes input vector so that the resulting vector represents a 
#'  probability distribution.
#'
#' @param lambda_t Numeric. Vector to be normalized
#' 
#' @return normalized vector
#'
#' @examples
#' transform_to_pdf(c(10,20,30,25,15,5))
#' 
transform_to_pdf = function(lambda_t){
  p_t = lambda_t/sum(lambda_t)
  return(p_t)
}

#' Simulation of age of samples given the lambda values of each year
#' 
#'  Simulates, using a Poisson distribution, the number of samples for each year 
#'  in the interval. It outputs the vector with the true ("calibrated") ages of 
#'  each individual simulated sample
#'
#' @param lambda_t Numeric. vector with the lambda values for each year of the
#' period
#' @param time_range Numeric. A vector of two values, c(t_0, t_f), expressed in 
#' calibrated years before present (i.e. before 1950).
#' 
#' @return vector of dates of the simulated sample
#'
#' @examples
#' sim_dates_lambda(1:10, c(8010,8001))
#' 
sim_dates_lambda = function(lambda_t, time_range){
  t_values = seq(time_range[1], time_range[2], -1)
  if (length(lambda_t)!=length(t_values)) {
    stop("lambda_t and t_values have to be of the same length")
  }
  n_t   = rpois(length(lambda_t), lambda_t)
  dates = rep(t_values, n_t)
  return(dates)
}

#' Simulation of age of a fixed number of samples given by a probability distribution
#' 
#'  Simulates, using a probability distribution, the ages of a fixed number of 
#'  samples. It outputs the vector with the true ("calibrated") ages of 
#'  each individual simulated sample
#'
#' @param n Integer. number of samples to simulate
#' @param p_t Numeric. vector of probabilities (its sums should be 1)
#' @param time_range Numeric. A vector of two values, c(t_0, t_f), expressed in 
#' calibrated years before present (i.e. before 1950).
#' 
#' @return vector of dates of the simulated sample
#'
#' @examples
#' sim_dates_from_pdf(20, transform_to_pdf(1:10), c(8010,8001))
#' 
sim_dates_from_pdf = function(n, p_t, time_range){
  t_values = seq(time_range[1], time_range[2], -1)
  if (length(p_t)!=length(t_values)) {
    stop("p_t and t_values have to be of the same length")
  }
  dates = sample(x = t_values, size = n, replace = T, prob = p_t)
  return(dates)
}

#' Simulation of Conventional Radiocarbon Age (14C measure)
#' 
#'  Simulates Conventional Radiocarbon Age (14C measure) for a vector of true
#'  ages (given as calibrated BP ages) using normal distributions based on the
#'  calibration curve.
#'
#' @param dates Integer. vector of ages measured as "calibrated" ages BP (i.e. 
#' before 1950)
#' @param calCurves calibration curve identifier from rcarbon package
#' @param errors Integer. a vector of integer values from which errors will be 
#' taken randomly to set the errors of CRA
#' 
#' @return a data frame with two columns, the CRA and the error
#'
#' @examples
#' sim_CRA(rep(8000, 6), errors = 30)
#' sim_CRA(rep(8000, 6), errors = c(30,45,80))
#' sim_CRA(c(8000, 7500, 7000, 6500), errors = rpois(10,50))
#' 
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

#' Calculates the range of CRA values expected given a calibration curve
#' 
#'  Given a time range in "calibrated" years BP, it calculate the range of CRA
#'  values that is expected to contain 99.73% (3*sigma) of the CRA values
#'  obtained from the true ages in the given interval.
#'
#' @param time_range Numeric. A vector of two values, c(t_0, t_f), expressed in 
#' calibrated years before present (i.e. before 1950).
#' @param calCurves calibration curve identifier from rcarbon package
#' @param error Integer. maximum error value considered in CRA
#' 
#' @return a vector of two values of CRA
#'
#' @examples
#' get_time_range_CRA(c(8000, 4000))
#' get_time_range_CRA(c(8000, 4000), calCurves = "shcal20")
#' get_time_range_CRA(c(8000, 4000), error = 50)
#' 
get_time_range_CRA = function(time_range, calCurves = "intcal20", error = 0){
  time_range_CRA = rcarbon::uncalibrate(time_range, calCurves = calCurves)
  time_range_CRA = c(time_range_CRA$ccCRA[1] + 3 * (time_range_CRA$ccError[1] + error),
                     time_range_CRA$ccCRA[2] - 3 * (time_range_CRA$ccError[2] + error))
  return(time_range_CRA)
}

#' Calculate summary statistics from the SPD from a set of CRA
#' 
#'  Obtains the SPD using rcarbon functions and takes the value of the SPD for 
#'  each year in the given interval.
#'
#' @param dates_CRA data frame containing two columns named "CRA" and "error"
#' with the conventional radiocarbon ages of the samples
#' @param time_range Numeric. A vector of two values, c(t_0, t_f), expressed in 
#' calibrated years before present (i.e. before 1950).
#' @param calCurves calibration curve identifier from rcarbon package
#' @param bins parameter for rcarbon::spd() A vector containing the bin names 
#' associated with each radiocarbon date. If set to NA, binning is not carried 
#' out.
#' @param runm parameter for rcarbon::spd() A number indicating the window size 
#' of the moving average to smooth the SPD. If set to NA no moving average is 
#' applied. Default is 100
#' 
#' @return a data frame with the values of the SPD for evary year in the time
#' interval indicated, with columns named as "spd" plus the year in cal BP (e.g.
#' "spd8000").
#'
#' @examples
#' get_sumstats_spd( data.frame(CRA=rep(7000,10), error=rpois(10,50)), c(8000,7800))
#' 
get_sumstats_spd = function(dates_CRA, time_range, calCurves = "intcal20", bins = NA, runm = 100){
  if (length(dates_CRA)>0){
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

#' Calculates summary statistics from a set of CRA (histogram version, DEPRECATED)
#' 
#'  Deprecated function to calculate summary statistics from a set of CRA. 
#'  Superseded by function get_sumstats()
#'
#' @param dates_CRA data frame containing in column "CRA" the values of 
#' conventional radiocarbon ages in a sample
#' @param time_range_CRA Numeric. the time range in CRA over which the summary 
#' statistics will be calculated
#' @param w Numeric. weight given to each sample
#' 
#' @return data frame with one row of summary statistics (in columns)
#' 
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

#' Calculates summary statistics from a set of CRA
#' 
#'  A set of conventional radiocarbon ages is summarize into a set of statistics
#'   that are informative about the abundance of samples through time. These
#'   statistics include: the total number of samples within the time range
#'   provided (n), the number of samples in the time sub-intervals defined by 
#'   dividing the total range into windows of sizes specified by parameter
#'   "window" (H) and the difference in the number of samples between 
#'   consecutive intervals (Delta H). These can be calculated with samples 
#'   having different weights, for instance, to compensate biases due to 
#'   oversampling of certain localities.
#'
#' @param dates_CRA data frame containing in column "CRA" the values of 
#' conventional radiocarbon ages in a sample
#' @param time_range_CRA Numeric. the time range in CRA over which the summary 
#' statistics will be calculated
#' @param window Numeric. sizes of the time windows over which summary statistics
#' will be calculated
#' @param w Numeric. weight given to each sample
#' 
#' @return data frame with one row of summary statistics (in columns). First
#' column is the total number of samples (n) then subsequent columns contain the
#' number of samples in periods ("windows") of the sizes provided with @param window
#' (denoted "w" plus window/period size), two stats are calculated par window,
#' the number of samples (h) and the difference in the number of samples with 
#' the consecutive window (dh), the windows are identified with the middle year 
#' of the window 
#'
#' @examples
#' get_sumstats( data.frame(CRA = sample(5000:7000, 500), error = rpois(500, 50)), c(7000, 5000))
#' 
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

#' Calculates summary statistics from two sets of CRA
#' 
#'  The summary statistics of two sets of CRA (as calculated by function 
#'    get_sumstats() using the same @param time_range_CRA and @param window ) that
#'    measure the correlation and covariance between the two sets of summary
#'    statistics.
#'
#' @param ssA data frame of summary statistics obtained with get_sumstats()
#' @param ssB data frame of summary statistics obtained with get_sumstats()
#' @param time_range_CRA Numeric. the time range in CRA over which the summary 
#' statistics will be calculated
#' @param window Numeric. sizes of the time windows over which summary statistics
#' will be calculated
#' 
#' @return data frame with one row and summary statistics in columns, including
#' the ratio between the number of samples of the first and second sets of CRAs,
#' and the correlation and covariances between the h and dh statistics of the
#' two sets
#'
#' @examples
#' set_a = data.frame(CRA = sample(5000:7000, 500), error = rpois(500, 50))
#' set_b = data.frame(CRA = sample(6000:7000, 400), error = rpois(400, 50))
#' time_range = c(7000, 5000)
#' window = 100
#' sumstats_a = get_sumstats(set_a, time_range, window)
#' sumstats_b = get_sumstats(set_b, time_range, window)
#' 
#' get_sumstats_correlation(sumstats_a, sumstats_b, time_range, window)
#' 
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

#' Sample from a prior distribution the values of the piecewise exponential model for two categories
#' 
#'  This function is used for fitting a piecewise exponential model for two 
#'  categories of samples. The parameters of two piecewise exponential models 
#'  can have three degrees of dependency: "independent", "interdependent" or 
#'  "parallel". The parameters sampled from the prior distribution are the 
#'  lambda values at the start and end of each period, for each category.
#'  
#'  For the independent model the prior is a log-uniform distribution between 
#'  lambda_min and lambda_max, the ratio between categories (pi) is calculated 
#'  from the lambda values. 
#'  
#'  For the parallel and interdependent models first the ratio between 
#'  categories (pi) is sampled from a beta distribution (with parameters alpha 
#'  and beta). For the parallel model a single value is sampled, for the 
#'  interdependent model a value for each period is sampled. Then lambda values 
#'  for the sum of the two categories are sampled from a log-uniform 
#'  distribution between lambda_min and lambda_max, and the lambda values for 
#'  each categories are obtained using pi.
#'  
#'  Given lambda values at the start and end of each period, it calculate the 
#'  growth rate (r) of the exponential model for each period and each category.
#'
#' @param lambda_min Numeric. Lower limit for prior
#' @param lambda_max Numeric. Upper limit for prior
#' @param time_range Numeric. A vector of two values, c(t_0, t_f), expressed in 
#' calibrated years before present (i.e. before 1950).
#' @param num_of_periods Integer. number of period in which the time interval
#' will be divided
#' @param intervals String or vector of values. Determines in which way the time
#' interval will be divided in the number of periods. The default ("dirichlet")
#' will take periods of random length following a Dirichlet(1) distribution.
#' Option "regular" will make equal (or approximately equal) periods. Finally, 
#' a numeric vector can be provided for user defined period lengths.  
#' @param scenario takes values "independent", "interdependent", "parallel"
#' @param alpha numeric first parameter for the beta prior for pi
#' @param beta numeric second parameter for the beta prior for pi
#' 
#' @return data.frame with values of lambda at the times defining the periods
#' (num_of_periods+1 values), the times defining the periods (num_of_periods+1
#' values), the exponential change rate within each interval (num_of_periods
#' values) and the ratio between lambdas of the two categories for each period
#' (num_of_periods values)
#'
#' @examples
#' sample_exponential_piecewise_parameters_2_categories(0.01, 10, c(8000, 4000), 4)
#' 
#' # Missing values for alpha and beta, independent model is used instead
#' sample_exponential_piecewise_parameters_2_categories(0.01, 10, c(8000, 4000), 4, scenario = "interdependent")
#' 
#' sample_exponential_piecewise_parameters_2_categories(0.01, 10, c(8000, 4000), 4, scenario = "interdependent", alpha = 1, beta = 1)
#' 
#' sample_exponential_piecewise_parameters_2_categories(0.01, 10, c(8000, 4000), 4, scenario = "parallel", alpha = 1, beta = 1)
#' 
sample_exponential_piecewise_parameters_2_categories = function(lambda_min, lambda_max,
                                                                time_range,
                                                                num_of_periods,
                                                                intervals = "dirichlet",
                                                                scenario = "independent", # "interdependent" "parallel"
                                                                alpha = NULL, beta = NULL){ 
  
  time_of_change = get_time_of_change(num_of_periods, time_range, intervals = intervals)
  
  if ( (scenario == "parallel" | scenario == "interdependent") & !is.null(alpha) & !is.null(beta) ){
    lambda_values = sample_lambda_values_piecewise_model(num_of_periods + 1, lambda_min, lambda_max)
    if (scenario == "parallel"){
      pi = rep(rbeta(1, alpha, beta), num_of_periods + 1) 
    }else{
      pi = rbeta(num_of_periods + 1, alpha, beta)
    }
    lambda_values_A = pi * lambda_values
    lambda_values_B = (1 - pi) * lambda_values
  } else {
    if ((scenario == "parallel" | scenario == "interdependent") & ( is.null(alpha) | is.null(beta) )){
      cat("Value for alpha or beta is NULL, running model 'independent' instead'\n")
    } else if (scenario != "independent"){
      cat("Undefined model, running model 'independent'\n")
    }
    lambda_values_A = sample_lambda_values_piecewise_model(num_of_periods+1, lambda_min, lambda_max)
    lambda_values_B = sample_lambda_values_piecewise_model(num_of_periods+1, lambda_min, lambda_max)
    pi = lambda_values_A/(lambda_values_A+lambda_values_B)
  }
  growth_rates_A = get_growth_rates(num_of_periods,log(lambda_values_A),time_of_change)
  growth_rates_B = get_growth_rates(num_of_periods,log(lambda_values_B),time_of_change)
  
  return(list(time_of_change=time_of_change,
              lambda_values_A=lambda_values_A,
              lambda_values_B=lambda_values_B,
              pi=pi,
              growth_rates_A=growth_rates_A,
              growth_rates_B=growth_rates_B))
}

#' Sample from a prior distribution the values of the piecewise exponential model for two categories (DEPRECATED)
get_piecewise_exponential_model_2_categories = function(num_of_periods,
                                                        lambda_min, lambda_max,
                                                        model = "independent", # "interdependent" "parallel"
                                                        alpha = NULL, beta = NULL,
                                                        time_range){
  time_of_change = get_time_of_change(num_of_periods, time_range, intervals = "regular")
  
  if (model == "parallel" & !is.null(pi) ){
    lambda_values = get_lambda_values_piecewise_model(num_of_periods + 1, lambda_min, lambda_max)
    growth_rates = get_growth_rates(num_of_periods,log(lambda_values),time_of_change)
    
    lambda_values_A = pi * lambda_values
    growth_rates_A = get_growth_rates(num_of_periods,log(lambda_values_A),time_of_change)
    lambda_values_B = (1 - pi) * lambda_values
    growth_rates_B = get_growth_rates(num_of_periods,log(lambda_values_B),time_of_change)
    
  } else if (model == "interdependent" & !is.null(alpha) & !is.null(beta) ){
    lambda_values = get_lambda_values_piecewise_model(num_of_periods + 1, lambda_min, lambda_max)
    growth_rates = get_growth_rates(num_of_periods,log(lambda_values),time_of_change)
    
    pi = rbeta(num_of_periods + 1,alpha,beta)
    
    lambda_values_A = pi * lambda_values
    growth_rates_A = get_growth_rates(num_of_periods,log(lambda_values_A),time_of_change)
    lambda_values_B = (1 - pi) * lambda_values
    growth_rates_B = get_growth_rates(num_of_periods,log(lambda_values_B),time_of_change)
    
  } else {
    if (model == "parallel" & is.null(pi)){
      cat("Value for pi is NULL, running model 'independent' instead of 'parallel'\n")
    } else if (model == "interdependent" & ( is.null(alpha) | is.null(beta) )){
      cat("Value for alpha or beta is NULL, running model 'independent' instead of 'interdependent'\n")
    } else if (model != "independent"){
      cat("Undefined model, running model 'independent'\n")
    }
    
    lambda_values_A = sample_lambda_values_piecewise_model(num_of_periods+1,lambda_min,lambda_max)
    growth_rates_A = get_growth_rates(num_of_periods,log(lambda_values_A),time_of_change)
    
    lambda_values_B = sample_lambda_values_piecewise_model(num_of_periods+1,lambda_min,lambda_max)
    growth_rates_B = get_growth_rates(num_of_periods,log(lambda_values_B),time_of_change)
    
    pi = lambda_values_A/(lambda_values_A+lambda_values_B)
    
  }

  lambda_t_A = lambda_t_B = numeric()
  r_t_A = r_t_B = numeric()
  for (i in 1:num_of_periods){
    lambda_t_A = c(lambda_t_A,
                   lambda_values_A[i] * exp(growth_rates_A[i]*seq_len(abs(time_of_change[i+1]-time_of_change[i])) ) )
    r_t_A = c(r_t_A, rep(growth_rates_A[i],abs(time_of_change[i+1]-time_of_change[i])))
    lambda_t_B = c(lambda_t_B,
                   lambda_values_B[i] * exp(growth_rates_B[i]*seq_len(abs(time_of_change[i+1]-time_of_change[i])) ) )
    r_t_B = c(r_t_B, rep(growth_rates_B[i],abs(time_of_change[i+1]-time_of_change[i])))
  }
  
  return(list(lambda_t_A      = lambda_t_A,
              lambda_t_B      = lambda_t_B,
              alpha           = alpha,
              beta            = beta,
              pi              = if (length(pi)>1) {mean(pi)} else {pi},
              lambda_values_A = lambda_values_A,
              growth_rates_A  = growth_rates_A,
              lambda_values_B = lambda_values_B,
              growth_rates_B  = growth_rates_B))
}

#' Calculates the logit of a number 
#'
#' @param p Numeric. number between 0 and 1
#'
#' @return Numeric. logit of p
#'
#' @examples
#' # Logit of 0.5
#' logit(0.5)
#' 
#' # Logit of 10; will produce a warning and NaN result
#' logit(10)
#'
logit <- function(p) {
  log(p / (1 - p))
}

#' Calculates the inverse logit of a number 
#'
#' @param x Numeric. positive and negative floats
#'
#' @return Numeric. inverse logit of x
#'
#' @examples
#' # inverse logit of 0.5
#' inv_logit(0.5)
#' 
#' # inverse logit of 10
#' inv_logit(10)
#'
inv_logit <- function(x) {
  exp(x) / (1 + exp(x))
}

#' Provides an interpretation of Bayes factor values
#' 
#'  Given a single value of a Bayes factor it writes a message given an
#'  interpretation of the strength of the evidence that the value provides for
#'  the model. Based on Jeffrey's scale 
#'  (https://en.wikipedia.org/wiki/Bayes_factor#Interpretation)
#'
#' @param k Numeric. a Bayes factor (a single value)
#' 
#' @examples
#' interpret_K(100)
#' interpret_K(10)
#' interpret_K(1)
#' interpret_K(0.1)
#' interpret_K(0.01)
#' 
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