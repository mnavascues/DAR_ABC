library(rcarbon)
library(zoo)
library(moments)
library(extraDistr)
library(weights)
library(Hmisc)

# simulate demographic values under a exponential model
get_exponential_model = function(lambda_min, lambda_max, time_range, expansion = FALSE){
  lambda = exp(runif(2,log(lambda_min),log(lambda_max)))
  if (expansion) lambda = sort(lambda)
  r = (log(lambda[2]) - log(lambda[1])) / (time_range[1] - time_range[2])
  t_values = seq_len(time_range[1]-time_range[2])
  lambda_t = lambda[1]*exp(r*t_values)
  return(list(lambda_t=lambda_t,
              lambda_0=lambda[1],
              lambda_f=lambda[2],
              rate=r))
}


# simulate demographic values under a logistic model
get_logistic_model = function(lambda_min, lambda_max, time_range, expansion = FALSE){
  lambda = exp(runif(2,log(lambda_min),log(lambda_max)))
  if (expansion) lambda = sort(lambda)
  K = max(lambda) + exp(runif(1,log(lambda_min),log(lambda_max)))
  while (K*lambda[1]/lambda[2]-lambda[1] <= 0 | K-lambda[1] <= 0){
    lambda = exp(runif(2,log(lambda_min),log(lambda_max)))
    if (expansion) lambda = sort(lambda)
    K = max(lambda) + exp(runif(1,log(lambda_min),log(lambda_max)))
  }
  r = (log(K-lambda[1]) - log(K*lambda[1]/lambda[2]-lambda[1])) / (time_range[1] - time_range[2])
  t_ = seq_len(time_range[1] - time_range[2])
  lambda_t = K*lambda[1] / (lambda[1] + (K - lambda[1])*exp(-r*t_) )
  return(list(lambda_t=lambda_t,
              lambda_0=lambda[1],
              lambda_f=lambda[2],
              lambda_K=K,
              rate=r))
}

get_dynamic_model = function(lambda_min, lambda_max,
                             b0_min, b0_max,
                             b1_min, b1_max,
                             b2_min, b2_max,
                             d, time_range){
  time_seq = seq(time_range[1], time_range[2] + 1,by = -1)
  lambda_t = array(NA, length(time_seq))

  while(lambda_t[length(lambda_t)]==0 | is.na(lambda_t[length(lambda_t)])){
    lambda_0 = exp(runif(1, log(lambda_min), log(lambda_max)))
    b0 = exp(runif(1, log(b0_min), log(b0_max)))
    b1 = -exp(runif(1, log(b1_min), log(b1_max)))
    b2 = -exp(runif(1, log(b2_min), log(b2_max)))
    ref_d = rnorm(1, mean = mean(d), sd = sd(d))
    
    D = d - ref_d
    for (i in seq_along(time_seq)){
      if (i == 1) {
        lambda_t[i] = lambda_0
      }else{
        lambda_t[i] = lambda_t[i-1] * exp(b0 + b1 * lambda_t[i-1]^2 + b2 * D[i]^2) 
      }
    }
  }

  return(list(lambda_t = lambda_t,
              lambda_0 = lambda_0,
              b0 = b0,
              b1 = b1,
              b2 = b2,
              d = ref_d))
}

# simulate demographic values under a piecewise demographic model
get_lambda_values_piecewise_model = function(num_of_periods,lambda_min,lambda_max){
  lambda_values = array(NA,num_of_periods)
  lambda_values[1] = exp(runif(1,log(lambda_min),log(lambda_max)))
  for (i in 2:(num_of_periods)){
    alpha = runif(1,log(0.1),log(10))
    lambda_values[i] = exp(max(min(log(lambda_values[i-1])+alpha,log(lambda_max)),log(lambda_min)))
  }
  return(lambda_values)
}
get_time_of_change = function(num_of_periods, time_range, intervals="dirichlet"){
  if (intervals=="dirichlet"){
    time_of_change = round(time_range[1] - cumsum(rdirichlet(1,rep(1,num_of_periods))) * (time_range[1] - time_range[2]))
    time_of_change = c(time_range[1], time_of_change)
    while (length(time_of_change)!=length(unique(time_of_change))){
      time_of_change = round(time_range[1] - cumsum(rdirichlet(1,rep(1,num_of_periods))) * (time_range[1] - time_range[2]))
      time_of_change = c(time_range[1], time_of_change)
    }
  }else if (intervals=="regular"){
    time_of_change = round(seq(time_range[1], time_range[2], - (time_range[1] - time_range[2]) / num_of_periods) )
  }else{
    stop("invalid intervals value. Values accepted: 'dirichlet' or 'regular'")
  }
  return(time_of_change)
}
get_growth_rates = function(num_of_periods,lambda_values,time_of_change){
  growth_rates = rep(NA,num_of_periods)
  for (i in seq_len(num_of_periods)){
    growth_rates[i] = (lambda_values[i+1]-lambda_values[i])/abs(time_of_change[i+1]-time_of_change[i])
  }
  return(growth_rates)
}
get_piecewise_constant_model = function(num_of_periods,lambda_min,lambda_max,time_range){
  lambda_values = get_lambda_values_piecewise_model(num_of_periods,lambda_min,lambda_max)
  time_of_change = get_time_of_change(num_of_periods,time_range)
  lambda_t = rep(lambda_values[1], time_of_change[1]-time_of_change[2])
  for (i in 2:num_of_periods){
    lambda_t = c(lambda_t, rep(lambda_values[1], time_of_change[1]-time_of_change[2]) )
  }
  length(lambda_t)
  return(lambda_t)
}
get_piecewise_linear_model = function(num_of_periods,lambda_min,lambda_max,time_range){
  lambda_values = get_lambda_values_piecewise_model(num_of_periods+1,lambda_min,lambda_max)
  time_of_change = get_time_of_change(num_of_periods,time_range)
  growth_rates = get_growth_rates(num_of_periods,lambda_values,time_of_change)
  lambda_t = numeric()
  for (i in 1:num_of_periods){
    lambda_t = c(lambda_t, lambda_values[i] + growth_rates[i]*seq_len(abs(time_of_change[i+1]-time_of_change[i])) )
  }
  return(list(lambda_t=lambda_t,
              lambda_values=lambda_values,
              time_of_change=time_of_change,
              growth_rates=growth_rates))
}
get_piecewise_exponential_model = function(num_of_periods,lambda_min,lambda_max,time_range, intervals="dirichlet", skyline=T){
  lambda_values = get_lambda_values_piecewise_model(num_of_periods+1,lambda_min,lambda_max)
  time_of_change = get_time_of_change(num_of_periods,time_range,intervals=intervals)
  growth_rates = get_growth_rates(num_of_periods,log(lambda_values),time_of_change)
  lambda_t = numeric()
  r_t = numeric()
  for (i in 1:num_of_periods){
    lambda_t = c(lambda_t,
                lambda_values[i] * exp(growth_rates[i]*seq_len(abs(time_of_change[i+1]-time_of_change[i])) ) )
    r_t = c(r_t, rep(growth_rates[i],abs(time_of_change[i+1]-time_of_change[i])))
  }
  if (intervals == "dirichlet" | skyline){
    skyline_years = c(1, seq(100, time_range[1] - time_range[2], by = 100))
    lambda_skyline = lambda_t[skyline_years]
    rate_skyline  = r_t[skyline_years]
  }else if (intervals == "regular" & !skyline){
    lambda_skyline = lambda_values
    rate_skyline = growth_rates
  }
  
  return(list(lambda_t       = lambda_t,
              lambda_skyline = lambda_skyline,
              rate_skyline  = rate_skyline))
}

get_piecewise_exponential_model_2_categories = function(num_of_periods,
                                                        lambda_min, lambda_max,
                                                        model = "independent", # "interdependent" "parallel"
                                                        alpha = NULL, beta = NULL, pi = NULL,
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
    
    lambda_values_A = get_lambda_values_piecewise_model(num_of_periods+1,lambda_min,lambda_max)
    growth_rates_A = get_growth_rates(num_of_periods,log(lambda_values_A),time_of_change)
    
    lambda_values_B = get_lambda_values_piecewise_model(num_of_periods+1,lambda_min,lambda_max)
    growth_rates_B = get_growth_rates(num_of_periods,log(lambda_values_B),time_of_change)
    
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

# simulate number of items and their age in the archaeological record giving a model with lambda_t
sim_dates = function(lambda_t,t_values){
  if (length(lambda_t)!=length(t_values)) {
    stop("lambda_t and t_values have to be of the same length")
  }
  n_t   = rpois(length(lambda_t), lambda_t)
  dates = rep(t_values, n_t)
  return(dates)
}
# simulate 14C measure
sim_14Cdates = function(dates,calCurves,errors){
  x = uncalibrate(dates, CRAerrors=sample(errors,size=length(dates),replace=T), calCurves=calCurves)
  C14Age = x$rCRA
  C14SD = x$rError
  return(list(C14Age=C14Age,C14SD=C14SD))
}
# simulation of SPD summary statistics
get_sumstats_SPD = function(dates, errors, calCurves, time_range, bins=NA, runm=100, window=100, window_reg=500){
  if (length(dates)>0){
    dates_calibrated = calibrate(x         = dates,
                                 errors    = errors,
                                 calCurves = calCurves)
    spd = spd(dates_calibrated, timeRange=time_range, bins=bins, runm=runm)
    reg_sumstats = rollapply(zoo(spd$grid),
                             width = window_reg, by = window,
                             FUN = function(Z){
                               t = lm(formula=PrDens~calBP, data = as.data.frame(Z), na.rm=T); 
                               return(-1*t$coef[2])
                             },
                             by.column=FALSE, align = "center")
    sumstats = rollapply(zoo(spd$grid$PrDens), width = window, by = window, FUN = mean, align = "center")
    mean_dates = time_range[1]-index(sumstats)
    sumstats = data.frame(t(sumstats))
    names(sumstats) = paste0("spd",mean_dates)
    mean_dates = time_range[1]-index(reg_sumstats)
    reg_sumstats =  data.frame(t(reg_sumstats))
    names(reg_sumstats) = paste0("slope",mean_dates)
  }else{
    t_values = seq(time_range[1],time_range[2]+1,by=-1)
    reg_sumstats = rollapply(zoo(data.frame(calBP=t_values,PrDens=rep(0,length(t_values)))),
                             width = window_reg, by = window,
                             FUN = function(Z){
                               t = lm(formula=PrDens~calBP, data = as.data.frame(Z), na.rm=T); 
                               return(-1*t$coef[2])
                             },
                             by.column=FALSE, align = "center")
    sumstats = rollapply(zoo(rep(0,length(t_values))), width = window, by = window, FUN = mean, align = "center")
    mean_dates = time_range[1]-index(sumstats)
    sumstats = data.frame(t(sumstats))
    names(sumstats) = paste0("spd",mean_dates)
    mean_dates = time_range[1]-index(reg_sumstats)
    reg_sumstats =  data.frame(t(reg_sumstats))
    names(reg_sumstats) = paste0("slope",mean_dates)
  }
  return(cbind(sumstats,reg_sumstats))
}
get_sumstats_uncalibrated = function(dates_uncalibrated, time_range, window = 100, w = NULL){
  x = dates_uncalibrated
  #x      = dates_uncalibrated$C14Age
  # errors = dates_uncalibrated$C14SD
  dates_on_range = intersect(which(x<time_range[1]), which(x>time_range[2]))
  x = x[dates_on_range]
  probs = seq(0, 1, 0.02)
  breaks = seq(time_range[1],time_range[2]-1,-window)
  dates = breaks[-1]+window/2
  num_of_bins = length(breaks)-1
  if (is.null(w)){
    n = length(x)
  } else{
    w = w[dates_on_range]
    n = sum(w)
  }
  if (n > 0){
    h   = rev(wtd.hist(x, breaks=breaks, plot=F, weight = w)$counts) ##
    dh  = diff(h) ##
    q   = wtd.quantile (x, probs = probs, weights=w)
    m   = wtd.mean(x, w)
    sd  = sqrt(wtd.var(x, w))
    # TODO -  Look for weighted version of:
    # gm  = exp(mean(log(x)))
    # gsd = exp(sqrt((length(x) - 1) / length(x)) * sd(log(x)))
    # hm  = 1 / mean(1 / x)
    # sk  = skewness(x)
    # k   = kurtosis(x)
    sumstats = c(n, h, dh, q, m, sd)
  }else{
    sumstats = rep(NA, 2 + length(probs) + num_of_bins*2 )
  }
  
  #sumstats = data.frame(t(sumstats))
  names(sumstats) = c("count",
                       paste0("hist",dates),
                       paste0("delta_h",dates[-1]),
                       paste0("quantile",probs),
                       "mean","sd")#,"g_mean","g_sd","h_mean","skew","kurt")
  return(as.data.frame(t(sumstats)))
}
get_sumstats_correlation = function(ssA,ssB){
  Hss = names(ssA)[grepl("hist",  names(ssA))]
  Dss = names(ssA)[grepl("delta",  names(ssA))]
  Qss = names(ssA)[grepl("quantile",  names(ssA))]
  HssB = names(ssB)[grepl("hist",  names(ssB))]
  DssB = names(ssB)[grepl("delta",  names(ssB))]
  QssB = names(ssB)[grepl("quantile",  names(ssB))]
  if (!all(Hss==HssB) | !all(Dss==DssB) | !all(Qss==QssB)){
    sumstats = rep(NA, 3*2*3)
  }else{
    HA = as.vector(t(ssA[,Hss]))
    HB = as.vector(t(ssB[,Hss]))
    DA = as.vector(t(ssA[,Dss]))
    DB = as.vector(t(ssB[,Dss]))
    QA = as.vector(t(ssA[,Qss]))
    QB = as.vector(t(ssB[,Qss]))
    
    PcorH = cor(HA, HB, method = "pearson")
    KcorH = cor(HA, HB, method = "kendall")
    ScorH = cor(HA, HB, method = "spearman")
    PcovH = cov(HA, HB, method = "pearson")
    KcovH = cov(HA, HB, method = "kendall")
    ScovH = cov(HA, HB, method = "spearman")
    
    PcorD = cor(DA, DB, method = "pearson")
    KcorD = cor(DA, DB, method = "kendall")
    ScorD = cor(DA, DB, method = "spearman")
    PcovD = cov(DA, DB, method = "pearson")
    KcovD = cov(DA, DB, method = "kendall")
    ScovD = cov(DA, DB, method = "spearman")
    
    PcorQ = cor(QA, QB, method = "pearson")
    KcorQ = cor(QA, QB, method = "kendall")
    ScorQ = cor(QA, QB, method = "spearman")
    PcovQ = cov(QA, QB, method = "pearson")
    KcovQ = cov(QA, QB, method = "kendall")
    ScovQ = cov(QA, QB, method = "spearman")

    ratioABtot  = ssA[,"count"]/ssB[,"count"]
    ratioABhist = (1+HA)/(1+HB)
    ratioABhist_names = character()
    for (i in seq_along(Hss)){
      ratioABhist_names = c(ratioABhist_names, paste0("ratioABhist",strsplit(Hss[i],"hist")[[1]][2]))
    }
        
    sumstats = c(PcorH, KcorH, ScorH, PcovH, KcovH, ScovH, 
                 PcorD, KcorD, ScorD, PcovD, KcovD, ScovD, 
                 PcorQ, KcorQ, ScorQ, PcovQ, KcovQ, ScovQ,
                 ratioABtot, ratioABhist)
  } 
  names(sumstats) = c("PcorH", "KcorH", "ScorH", "PcovH", "KcovH", "ScovH", 
                      "PcorD", "KcorD", "ScorD", "PcovD", "KcovD", "ScovD", 
                      "PcorQ", "KcorQ", "ScorQ", "PcovQ", "KcovQ", "ScovQ",
                      "ratioABtot", ratioABhist_names)
  return(as.data.frame(t(sumstats)))
}






# whole simulation
sim_all = function(lambda_t,t_values,SPD=T,calCurves="intcal20",errors=50,bins=NA,runm=100,window=100,window_reg=500,out_all=F){
  dates = sim_dates(lambda_t,t_values)
  time_range = c(t_values[1],t_values[length(t_values)])
  if (length(dates)>0){
    dates_uncalibrated = sim_14Cdates(dates,calCurves,errors=errors)
    sumstats = get_sumstats_uncalibrated(dates_uncalibrated$C14Age, time_range)
  }else{
    dates_uncalibrated = list(C14Age = numeric(), C14SD = numeric())
    sumstats = get_sumstats_uncalibrated(dates_uncalibrated$C14Age, time_range)
  }
  if (SPD){
    sumstats_spd = get_sumstats_SPD(dates_uncalibrated$C14Age,
                                    dates_uncalibrated$C14SD,
                                    calCurves, time_range, bins=bins, runm=runm, window=window, window_reg=window_reg)
    sumstats = cbind(sumstats, sumstats_spd)
  }
  if (out_all){
    return(list(dates              = dates,
                dates_uncalibrated = dates_uncalibrated,
                sumstats           = sumstats))
  }else{
    return(sumstats)
  }
}

# whole simulation with 2 categories
sim_2_categories = function(lambda_t_A, lambda_t_B, t_values,
                            calCurves = "intcal20", errors = 50, 
                            bins = NA,
                            runm = 100, window = 100, window_reg = 500, out_all = F){
  
  time_range = c(t_values[1],t_values[length(t_values)])
  
  for (i in 1:2){
    if (i == 1) lambda_t = lambda_t_A
    if (i == 2) lambda_t = lambda_t_B
    dates = sim_dates(lambda_t, t_values)
    if (length(dates)>0){
      dates_uncalibrated = sim_14Cdates(dates, calCurves, errors = errors)
    }else{
      dates_uncalibrated = list(C14Age = numeric(), C14SD = numeric())
    }
    sumstats = get_sumstats_uncalibrated(dates_uncalibrated$C14Age, time_range)
    if (i == 1) all_dates_uncalibrated = dates_uncalibrated
    if (i == 2) all_dates_uncalibrated = list(C14Age = c(dates_uncalibrated$C14Age,
                                                         all_dates_uncalibrated$C14Age),
                                              C14SD  = c(dates_uncalibrated$C14SD,
                                                         all_dates_uncalibrated$C14SD))
    if (i == 1) sumstats_14C_A = sumstats
    if (i == 2) sumstats_14C_B = sumstats
  }
  if (length(all_dates_uncalibrated$C14Age)==0){
    all_dates_uncalibrated = list(C14Age = numeric(), C14SD = numeric())
  }
  sumstats_14C_all = get_sumstats_uncalibrated(all_dates_uncalibrated$C14Age, time_range)
  sumstats_cor = get_sumstats_correlation(sumstats_14C_A, sumstats_14C_B)
  
  names(sumstats_14C_A) = paste0(names(sumstats_14C_A),"_A")
  names(sumstats_14C_B) = paste0(names(sumstats_14C_B),"_B")

  sumstats = cbind(sumstats_14C_A, sumstats_14C_B, sumstats_14C_all, sumstats_cor)
  
  return(sumstats)
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


