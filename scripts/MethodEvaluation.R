source("scripts/DARthABC.R")
require(rbenchmark)
require(abcrf)
require(ggplot2)
results_directory = "results/method_evaluation"
dir.create(results_directory)

# simulate target data and until X replicates have exactly n samples
# (where n is the expected number of samples under the model):
time_range = c(7000, 5000)
number_of_replicates = 300 # X
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
# p_t = transform_to_pdf(lambda_t)

simulated_target_file = paste0(results_directory, "/simulated_target.rda")
if (!file.exists(simulated_target_file)){
  counter = 0
  dates = sim_dates_lambda(lambda_t, time_range)
  dates_CRA = sim_CRA(dates, errors = errors)
  if (length(dates) == n){
    counter = counter + 1
    print(paste(counter,"replicates with",n,"samples"))
  }
  simulated_target_sumstats = get_sumstats(dates_CRA, time_range_CRA, window = c(10,50,100,500) )

  while (counter <= number_of_replicates){
    dates = sim_dates_lambda(lambda_t, time_range)
    dates_CRA = sim_CRA(dates, errors = errors)
    if (length(dates) == n){
      counter = counter + 1
      print(paste(counter,"replicates with",n,"samples"))
    }
    simulated_target_sumstats = rbind(simulated_target_sumstats,
                                      get_sumstats(dates_CRA, time_range_CRA, window = c(10,50,100,500)))
  }
  save(simulated_target_sumstats, file = simulated_target_file)
}else{
  #load(file = simulated_target_file)
}

# benchmark calculation of summary statistics based on SPD and histogram of CRA
bench_results_file = paste0(results_directory, "/bench_results.rda")
if (!file.exists(bench_results_file)){
  bench_results = rbenchmark::benchmark(
    "spd" = {
      dates = sample(seq(time_range[1],time_range[2],-1), n)
      dates_CRA = sim_CRA(dates, errors = errors)
      get_sumstats_spd(dates_CRA, time_range)
    },
    "hist" = {
      dates = sample(seq(time_range[1],time_range[2],-1), n)
      dates_CRA = sim_CRA(dates, errors = errors)
      get_sumstats_hist(dates_CRA, time_range)
    },
    "hist2" = {
      dates = sample(seq(time_range[1],time_range[2],-1), n)
      dates_CRA = sim_CRA(dates, errors = errors)
      get_sumstats(dates_CRA, time_range, window = c(10,50,100,500))
    },
    replications = 100)
  save(bench_results, file = bench_results_file)
}else{
  #load(file = bench_results_file)
}
#(bench_results)

# generate reference table for exponential model (probability distribution).
# reference table with two sets of summary stats (spd & hist)
rate_min = -0.01
rate_max = 0.01
reftable_file = paste0(results_directory, "/reftable_exponential_prob_model.rda")
if (!file.exists(reftable_file) ){
  num_of_sims = 20000
  ncores = 15
  require(doSNOW)
  require(doParallel)
  require(doRNG)
  # setup parallel computing
  cl <- makeCluster(ncores, type = "FORK")  
  registerDoParallel(cl)  
  registerDoRNG(seed = 33333)
  reftable <- foreach(sim = seq_len(num_of_sims), .combine = rbind) %dopar% {
    rm(params, lambda_t, p_t, dates, dates_CRA, lambda_target_sumstats,
       sumstats_spd, sumstats_hist, sumstats, p_target_sumstats, ss)
    gc()
    
    rate = runif(1, rate_min, rate_max)
    lambda_t = get_exponential_lambda_t(1, rate, time_range)
    p_t = transform_to_pdf(lambda_t)
    p_0 = p_t[1]
    p_f = tail(p_t, 1)

    dates = sim_dates_from_pdf(n, p_t, time_range)
    dates_CRA = sim_CRA(dates, errors = errors)
    sumstats_spd = get_sumstats_spd(dates_CRA, time_range)
    sumstats_hist = get_sumstats_hist(dates_CRA, time_range_CRA)
    sumstats_hist_2 = get_sumstats(dates_CRA, time_range_CRA, window = c(10,50,100,500))
    sumstats = cbind(sumstats_spd, sumstats_hist, sumstats_hist_2)
    params = cbind(rate, p_0, p_f)
    cbind(params, sumstats)
  }
  stopCluster(cl) 
  save(reftable, file = reftable_file)
}else{
  #load(file = reftable_file)
}

RF_model_p_hist_file = paste0(results_directory, "/RF_model_p_hist.rda")
if (!file.exists(RF_model_p_hist_file) ){
  load(file = reftable_file)
  sumstats = reftable[paste0("hist",seq(time_range_CRA[1], time_range_CRA[2], -1))]
  param    = reftable$rate
  RF_p_hist = regAbcrf(param~., data.frame(param,sumstats),
                                   ntree = 5000, paral = TRUE)
  save(RF_p_hist, file = RF_model_p_hist_file)
}else{
  # load(file = RF_model_p_hist_file)
}

RF_model_p_hist2_file = paste0(results_directory, "/RF_model_p_hist2.rda")
if (!file.exists(RF_model_p_hist2_file) ){
  load(file = reftable_file)
  sumstats = reftable[3994:4516]
  param    = reftable$rate
  RF_p_hist2 = regAbcrf(param~., data.frame(param,sumstats),
                       ntree = 5000, paral = TRUE)
  save(RF_p_hist2, file = RF_model_p_hist2_file)
}else{
  # load(file = RF_model_p_hist2_file)
}


RF_model_p_spd_file = paste0(results_directory, "/RF_model_p_spd.rda")
if (!file.exists(RF_model_p_spd_file) ){
  load(file = RF_model_p_hist_file)
  sumstats = reftable[paste0("spd",seq(time_range[1], time_range[2], -1))]
  param    = reftable$rate
  RF_p_spd = regAbcrf(param~., data.frame(param,sumstats),
                       ntree = 5000, paral = TRUE)
  save(RF_p_spd, file = RF_model_p_spd_file)
}else{
  #load(file = RF_model_p_spd_file)
}


########################################################################

# generate reference table for exponential model.
# reference table with two different models probability based and lambda based
lambda_min = 0.005
lambda_max = 5
rate_min = -0.005
rate_max = 0.005
reftable_file = paste0(results_directory, "/reftable_exponential_2_exponential_models.rda")
if (!file.exists(reftable_file) ){
  num_of_sims = 100000
  ncores = 20
  require(doSNOW)
  require(doParallel)
  require(doRNG)
  # setup parallel computing
  cl <- makeCluster(ncores, type = "FORK")  
  registerDoParallel(cl)  
  registerDoRNG(seed = 33333)
  reftable <- foreach(sim = seq_len(num_of_sims), .combine = rbind) %dopar% {
    rm(params, lambda_t, p_t, dates, dates_CRA, lambda_target_sumstats,
       sumstats_spd, sumstats_hist, sumstats, p_target_sumstats, ss)
    gc()
    
    #params = sample_exponential_parameters_from_priors(lambda_min, lambda_max, time_range)
    
    sum_lambda_t = 5001
    while (sum_lambda_t>5000){
      rate = runif(1, rate_min, rate_max)
      lambda_0 = 10^runif(1, log10(lambda_min), log10(lambda_max))
      lambda_t = get_exponential_lambda_t(lambda_0, rate, time_range)
      sum_lambda_t = sum(lambda_t)
    }
    lambda_f = tail(lambda_t, 1)
    p_t = transform_to_pdf(lambda_t)
    p_0 = p_t[1]
    p_f = tail(p_t, 1)
    params = cbind(lambda_0, lambda_f, rate, p_0, p_f)
    
    dates = sim_dates_lambda(lambda_t, time_range)
    dates_CRA = sim_CRA(dates, errors = errors)
    lambda_sumstats = get_sumstats(dates_CRA, time_range_CRA, window = c(10,50,100,500))
    
    dates = sim_dates_from_pdf(n, p_t, time_range)
    dates_CRA = sim_CRA(dates, errors = errors)
    p_sumstats = get_sumstats(dates_CRA, time_range_CRA, window = c(10,50,100,500))
    names(p_sumstats) = paste0("p_", names(p_sumstats))
    sumstats = cbind(lambda_sumstats, p_sumstats)
    
    cbind(params, sumstats)
  }
  stopCluster(cl) 
  save(reftable, file = reftable_file)
}else{
  #load(file = reftable_file)
}

num_of_param = 5

RF_model_p_file = paste0(results_directory, "/RF_model_p.rda")
if (!file.exists(RF_model_p_file) ){
  load(file = reftable_file)
  num_of_sumstats = (ncol(reftable) - num_of_param) / 2
  lambda_sumstats_names = names(reftable)[(1:num_of_sumstats) + num_of_param]
  p_sumstats_names = names(reftable)[(1:num_of_sumstats) + num_of_param + num_of_sumstats]
  sumstats = reftable[p_sumstats_names]
  param    = reftable$rate

  RF_p = regAbcrf(param~., data.frame(param,sumstats),
                      ntree = 5000, paral = TRUE)
  save(RF_p, file = RF_model_p_file)
}else{
  # load(file = RF_model_p_file)
}
results_RF_model_p_file = paste0(results_directory, "/results_RF_model_p.rda")
if (!file.exists(results_RF_model_p_file) ){
  load(file = reftable_file)
  num_of_sumstats = (ncol(reftable) - num_of_param) / 2
  lambda_sumstats_names = names(reftable)[(1:num_of_sumstats) + num_of_param]
  p_sumstats_names = names(reftable)[(1:num_of_sumstats) + num_of_param + num_of_sumstats]
  sumstats = reftable[p_sumstats_names]
  param    = reftable$rate
  load(file = RF_model_p_file)
  load(file = simulated_target_file)
  target_data = simulated_target_sumstats[which(simulated_target_sumstats$n==n),]
  names(target_data) = paste0("p_",names(target_data))
  results_p = predict(RF_p,
                      obs=target_data,
                      training=data.frame(param,sumstats),
                      quantiles = seq(0,1,0.005),
                      paral=T)
  save(results_p, file=results_RF_model_p_file)
}else{
  # load(file=results_RF_model_p_file)
}

RF_model_p_0_file = paste0(results_directory, "/RF_model_p_0.rda")
if (!file.exists(RF_model_p_0_file) ){
  load(file = reftable_file)
  num_of_sumstats = (ncol(reftable) - num_of_param) / 2
  lambda_sumstats_names = names(reftable)[(1:num_of_sumstats) + num_of_param]
  p_sumstats_names = names(reftable)[(1:num_of_sumstats) + num_of_param + num_of_sumstats]
  sumstats = reftable[p_sumstats_names]
  param    = log10(reftable$p_0)
  RF_p_0 = regAbcrf(param~., data.frame(param,sumstats),
                    ntree = 5000, paral = TRUE)
  save(RF_p_0, file = RF_model_p_0_file)
}else{
  # load(file = RF_model_p_0_file)
}
results_RF_model_p_0_file = paste0(results_directory, "/results_RF_model_p_0.rda")
if (!file.exists(results_RF_model_p_0_file) ){
  load(file = reftable_file)
  num_of_sumstats = (ncol(reftable) - num_of_param) / 2
  lambda_sumstats_names = names(reftable)[(1:num_of_sumstats) + num_of_param]
  p_sumstats_names = names(reftable)[(1:num_of_sumstats) + num_of_param + num_of_sumstats]
  sumstats = reftable[p_sumstats_names]
  param    = log10(reftable$p_0)
  load(file = RF_model_p_0_file)
  load(file = simulated_target_file)
  target_data = simulated_target_sumstats[which(simulated_target_sumstats$n == n),]
  names(target_data) = paste0("p_", names(target_data))
  results_p_0 = predict(RF_p_0,
                        obs = target_data,
                        training = data.frame(param,sumstats),
                        paral = T)
  save(results_p_0, file = results_RF_model_p_0_file)
}else{
  # load(file = results_RF_model_p_0_file)
}

RF_model_p_f_file = paste0(results_directory, "/RF_model_p_f.rda")
if (!file.exists(RF_model_p_f_file) ){
  load(file = reftable_file)
  num_of_sumstats = (ncol(reftable) - num_of_param) / 2
  lambda_sumstats_names = names(reftable)[(1:num_of_sumstats) + num_of_param]
  p_sumstats_names = names(reftable)[(1:num_of_sumstats) + num_of_param + num_of_sumstats]
  sumstats = reftable[p_sumstats_names]
  param    = reftable$p_f
  RF_p_f = regAbcrf(param~., data.frame(param,sumstats),
                    ntree = 5000, paral = TRUE)
  save(RF_p_f, file = RF_model_p_f_file)
}else{
  # load(file = RF_model_p_f_file)
}
results_RF_model_p_f_file = paste0(results_directory, "/results_RF_model_p_f.rda")
if (!file.exists(results_RF_model_p_f_file) ){
  load(file = reftable_file)
  num_of_sumstats = (ncol(reftable) - num_of_param) / 2
  lambda_sumstats_names = names(reftable)[(1:num_of_sumstats) + num_of_param]
  p_sumstats_names = names(reftable)[(1:num_of_sumstats) + num_of_param + num_of_sumstats]
  sumstats = reftable[p_sumstats_names]
  param    = reftable$p_f
  load(file = RF_model_p_f_file)
  load(file = simulated_target_file)
  target_data = simulated_target_sumstats[which(simulated_target_sumstats$n == n),]
  names(target_data) = paste0("p_", names(target_data))
  results_p_f = predict(RF_p_f,
                        obs = target_data,
                        training = data.frame(param,sumstats),
                        paral = T)
  save(results_p_f, file = results_RF_model_p_f_file)
}else{
  # load(file = results_RF_model_p_f_file)
}













RF_model_lambda_file = paste0(results_directory, "/RF_model_lambda.rda")
if (!file.exists(RF_model_lambda_file) ){
  load(file = reftable_file)
  num_of_sumstats = (ncol(reftable) - num_of_param) / 2
  lambda_sumstats_names = names(reftable)[(1:num_of_sumstats) + num_of_param]
  p_sumstats_names = names(reftable)[(1:num_of_sumstats) + num_of_param + num_of_sumstats]
  sumstats = reftable[lambda_sumstats_names]
  param    = reftable$rate
  RF_lambda = regAbcrf(param~., data.frame(param,sumstats),
                  ntree = 5000, paral = TRUE)
  save(RF_lambda, file = RF_model_lambda_file)
}else{
  # load(file = RF_model_lambda_file)
}
results_RF_model_lambda_file = paste0(results_directory, "/results_RF_model_lambda.rda")
if (!file.exists(results_RF_model_lambda_file) ){
  load(file = reftable_file)
  num_of_sumstats = (ncol(reftable) - num_of_param) / 2
  lambda_sumstats_names = names(reftable)[(1:num_of_sumstats) + num_of_param]
  p_sumstats_names = names(reftable)[(1:num_of_sumstats) + num_of_param + num_of_sumstats]
  sumstats = reftable[lambda_sumstats_names]
  param    = reftable$rate
  load(file = RF_model_lambda_file)
  load(file = simulated_target_file)
  target_data = simulated_target_sumstats[1:1000,]
  results_lambda = predict(RF_lambda,
                      obs=target_data,
                      training=data.frame(param,sumstats),
                      quantiles = seq(0,1,0.005),
                      paral=T)
  save(results_lambda, file=results_RF_model_lambda_file)
}else{
  #load(file=results_RF_model_lambda_file)
}



RF_model_lambda_0_file = paste0(results_directory, "/RF_model_lambda_0.rda")
if (!file.exists(RF_model_lambda_0_file) ){
  load(file = reftable_file)
  num_of_sumstats = (ncol(reftable) - num_of_param) / 2
  lambda_sumstats_names = names(reftable)[(1:num_of_sumstats) + num_of_param]
  sumstats = reftable[lambda_sumstats_names]
  param    = log10(reftable$lambda_0)
  RF_lambda_0 = regAbcrf(param~., data.frame(param,sumstats),
                       ntree = 5000, paral = TRUE)
  save(RF_lambda_0, file = RF_model_lambda_0_file)
}else{
  # load(file = RF_model_lambda_0_file)
}
results_RF_model_lambda_0_file = paste0(results_directory, "/results_RF_model_lambda_0.rda")
if (!file.exists(results_RF_model_lambda_0_file) ){
  load(file = reftable_file)
  num_of_sumstats = (ncol(reftable) - num_of_param) / 2
  lambda_sumstats_names = names(reftable)[(1:num_of_sumstats) + num_of_param]
  sumstats = reftable[lambda_sumstats_names]
  param    = log10(reftable$lambda_0)
  load(file = RF_model_lambda_0_file)
  load(file = simulated_target_file)
  target_data = simulated_target_sumstats[which(simulated_target_sumstats$n==n),]
  results_lambda_0 = predict(RF_lambda_0,
                           obs=target_data,
                           training=data.frame(param,sumstats),
                           paral=T)
  save(results_lambda_0, file=results_RF_model_lambda_0_file)
}else{
  #load(file=results_RF_model_lambda_0_file)
}



RF_model_lambda_f_file = paste0(results_directory, "/RF_model_lambda_f.rda")
if (!file.exists(RF_model_lambda_f_file) ){
  load(file = reftable_file)
  num_of_sumstats = (ncol(reftable) - num_of_param) / 2
  lambda_sumstats_names = names(reftable)[(1:num_of_sumstats) + num_of_param]
  sumstats = reftable[lambda_sumstats_names]
  param    = reftable$lambda_f
  RF_lambda_f = regAbcrf(param~., data.frame(param,sumstats),
                         ntree = 5000, paral = TRUE)
  save(RF_lambda_f, file = RF_model_lambda_f_file)
}else{
  # load(file = RF_model_lambda_f_file)
}
results_RF_model_lambda_f_file = paste0(results_directory, "/results_RF_model_lambda_f.rda")
if (!file.exists(results_RF_model_lambda_f_file) ){
  load(file = reftable_file)
  num_of_sumstats = (ncol(reftable) - num_of_param) / 2
  lambda_sumstats_names = names(reftable)[(1:num_of_sumstats) + num_of_param]
  sumstats = reftable[lambda_sumstats_names]
  param    = reftable$lambda_f
  load(file = RF_model_lambda_f_file)
  load(file = simulated_target_file)
  target_data = simulated_target_sumstats[which(simulated_target_sumstats$n==n),]
  results_lambda_f = predict(RF_lambda_f,
                             obs=target_data,
                             training=data.frame(param,sumstats),
                             paral=T)
  save(results_lambda_f, file=results_RF_model_lambda_f_file)
}else{
  #load(file=results_RF_model_lambda_f_file)
}

