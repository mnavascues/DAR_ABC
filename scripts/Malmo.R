source("scripts/DARthABC.R")
results_directory = "results/Malmo"
dir.create(results_directory)
library(doSNOW)
library(abcrf)
library(abc)

dates = read.csv("data/C14_projekt_Innovationsprocesser_master_230314.csv")
head(dates)

summary(dates$BP)
time_range_CRA = c(max(dates$BP),min(dates$BP))
time_range_calBP = c(6000,3200)

caldates_filename = paste0(results_directory, "/caldates.rda")
if (!file.exists(caldates_filename)){
  ncores = 30
  cl = makeCluster(ncores, type="SOCK")
  registerDoSNOW(cl)
  caldates = rcarbon::calibrate(x = dates$BP,
                                errors = dates$X.,
                                ids = dates$Id.nr,
                                calCurves = 'intcal20',
                                #timeRange = time_range_CRA,
                                normalised = FALSE,
                                ncores = ncores,
                                calMatrix = TRUE)
  stopCluster(cl)
  save(caldates, file = caldates_filename)
}else{
  load(file = caldates_filename)
}

spd_from_feature = rcarbon::spd(x = caldates[dates$Dating.from.feature],
                                timeRange = time_range_calBP, 
                                datenormalised = FALSE,
                                runm = 100)
spd_old_wood = rcarbon::spd(x = caldates[dates$Old.wood.effect.present & dates$Dating.from.feature],
                            timeRange = time_range_calBP, 
                            datenormalised = FALSE,
                            runm = 100)
spd_short_term = rcarbon::spd(x = caldates[dates$Short.term.human.impact & dates$Dating.from.feature],
                              timeRange = time_range_calBP, 
                              datenormalised = FALSE,
                              runm = 100)

plot(spd_from_feature$grid$calBP,
     spd_from_feature$grid$PrDens,
     xlim = time_range_calBP,
     ylim = c(0,0.65),
     #log = "y",
     type ="l",
     xlab = "Years cal BP",
     ylab = "SPD",
     lwd = 2)
lines(spd_old_wood$grid$calBP,
      spd_old_wood$grid$PrDens,
      col = "blue", lwd = 2)
lines(spd_short_term$grid$calBP,
      spd_short_term$grid$PrDens,
      col = "red", lwd = 2)

# generate reference table for piecewise exponential model.
lambda_min = 0.001
lambda_max = 2
num_of_periods = 10
num_of_sims = 100000
malmo_reftable = paste0(results_directory, "/malmo_reftable.rda")
if (!file.exists(malmo_reftable) ){
  ncores = 30
  require(doSNOW)
  require(doParallel)
  require(doRNG)
  # setup parallel computing
  cl <- makeCluster(ncores, type = "FORK")  
  registerDoParallel(cl)  
  registerDoRNG(seed = 1234567)
  reftable <- foreach(sim = seq_len(num_of_sims), .combine = rbind) %dopar% {
    rm(params, lambda_t, sim_dates, sim_dates_CRA, sumstats)
    gc()
    
    params = sample_exponential_piecewise_parameters_from_priors(lambda_min, lambda_max,
                                                                 time_range = time_range_calBP,
                                                                 num_of_periods)
    lambda_t = get_piecewise_exponential_lambda_t(as.numeric(params[seq_len(num_of_periods+1)]),
                                                  as.integer(params[1+num_of_periods+seq_len(num_of_periods+1)]),
                                                  as.numeric(params[2+2*num_of_periods+seq_len(num_of_periods)]))
    sim_dates = sim_dates_lambda(lambda_t, time_range_calBP)
    sim_dates_CRA = sim_CRA(sim_dates, errors = dates$X.[dates$Dating.from.feature])
    sumstats = get_sumstats(sim_dates_CRA, time_range_CRA,  window = c(25, 50, 100))

    cbind(params, sumstats)
  }
  stopCluster(cl) 
  save(reftable, file = malmo_reftable)
}else{
  load(file = malmo_reftable)
}

skyline_step = 25
skyline_points = 1 + (time_range_calBP[1]-time_range_calBP[2])/skyline_step
skyline_years = seq(time_range_calBP[1],time_range_calBP[2],-skyline_step)
malmo_reftable_skyline = paste0(results_directory, "/malmo_reftable_skyline.rda")
if (!file.exists(malmo_reftable_skyline) ){
  reftable_skyline_plot = data.frame(matrix(NA, nrow=nrow(reftable), ncol=skyline_points*2))
  names(reftable_skyline_plot) =  c(paste0("lambda",skyline_years),
                                    paste0("rate",skyline_years))
  
  for (sim in seq_len(num_of_sims)){
    print(paste("sim",sim))
    lambdas         = as.numeric(reftable[sim,seq_len(num_of_periods+1)])
    times_of_change = as.integer(reftable[sim,1+num_of_periods+seq_len(num_of_periods+1)])
    rates           = as.numeric(reftable[sim,2+2*num_of_periods+seq_len(num_of_periods)])
    
    lambda_skyline = get_piecewise_exponential_lambda_t(lambdas,
                                                        times_of_change,
                                                        rates)[seq(1,skyline_points*skyline_step,skyline_step)]
    
    
    rates_skyline = c(rep(rates,abs(diff(times_of_change))), rates[length(rates)])[seq(1,skyline_points*skyline_step,skyline_step)]
    
    reftable_skyline_plot[sim,] =  c(lambda_skyline,rates_skyline)
  }
  
  save(reftable_skyline_plot, file = malmo_reftable_skyline)
}else{
  load(file = malmo_reftable_skyline)
}


dates_CRA = list(CRA = dates$BP[dates$Dating.from.feature])
dates_CRA_NOWE = list(CRA = dates$BP[dates$Dating.from.feature & !dates$Old.wood.effect.present])
dates_CRA_SHORT = list(CRA = dates$BP[dates$Dating.from.feature & dates$Short.term.human.impact])
target_sumstats_all = get_sumstats(dates_CRA, time_range_CRA, window = c(25, 50, 100))
target_sumstats_nowe = get_sumstats(dates_CRA_NOWE, time_range_CRA, window = c(25, 50, 100))
target_sumstats_short = get_sumstats(dates_CRA_SHORT, time_range_CRA, window = c(25, 50, 100))

results_skyline_file = paste0(results_directory, "/skyline_results.rda")
if (!file.exists(results_skyline_file) ){
  
  skyline_all = list(median  = rep(NA, skyline_points),
                        q_0.025 = rep(NA, skyline_points),
                        q_0.975 = rep(NA, skyline_points),
                        error   = rep(NA, skyline_points))
  skyline_nowe = list(median  = rep(NA, skyline_points),
                        q_0.025 = rep(NA, skyline_points),
                        q_0.975 = rep(NA, skyline_points),
                        error   = rep(NA, skyline_points))
  skyline_short = list(median  = rep(NA, skyline_points),
                        q_0.025 = rep(NA, skyline_points),
                        q_0.975 = rep(NA, skyline_points),
                        error   = rep(NA, skyline_points))
  
  sumstats = reftable[, names(target_sumstats_all)]
  
  for (i in seq_len(skyline_points)){
    param = log10(reftable_skyline_plot[,i])
    RF_lambda = regAbcrf(param~., data.frame(param,sumstats),
                         ntree = 1000, paral = TRUE)
    posterior_lambda_all = predict(RF_lambda, target_sumstats_all,
                                   training = data.frame(param,sumstats),
                                   paral = TRUE, rf.weights = F)
    skyline_all$median[i] = posterior_lambda_all$med
    skyline_all$q_0.025[i] = posterior_lambda_all$quantiles[1]
    skyline_all$q_0.975[i] = posterior_lambda_all$quantiles[2]
    skyline_all$error[i] = RF_lambda$model.rf$prediction.error
    
    posterior_lambda_nowe = predict(RF_lambda, target_sumstats_nowe,
                                   training = data.frame(param,sumstats),
                                   paral = TRUE, rf.weights = F)
    skyline_nowe$median[i] = posterior_lambda_nowe$med
    skyline_nowe$q_0.025[i] = posterior_lambda_nowe$quantiles[1]
    skyline_nowe$q_0.975[i] = posterior_lambda_nowe$quantiles[2]
    skyline_nowe$error[i] = RF_lambda$model.rf$prediction.error
    
    posterior_lambda_short = predict(RF_lambda, target_sumstats_short,
                                    training = data.frame(param,sumstats),
                                    paral = TRUE, rf.weights = F)
    skyline_short$median[i] = posterior_lambda_short$med
    skyline_short$q_0.025[i] = posterior_lambda_short$quantiles[1]
    skyline_short$q_0.975[i] = posterior_lambda_short$quantiles[2]
    skyline_short$error[i] = RF_lambda$model.rf$prediction.error
    
    
    rm(param, RF_lambda, posterior_lambda_short, posterior_lambda_nowe, posterior_lambda_all)
    gc()
    save(skyline_short, skyline_nowe, skyline_all, file = results_skyline_file)
  }
  save(skyline_short, skyline_nowe, skyline_all, file = results_skyline_file)
}else{
  load(file = results_skyline_file)
  if (any(is.na(skyline_all$median))){
    sumstats = reftable[, names(target_sumstats_all)]
    
    initial_i = which(is.na(skyline_all$median))[1]
    for (i in initial_i:skyline_points){
      param = log10(reftable_skyline_plot[,i])
      RF_lambda = regAbcrf(param~., data.frame(param,sumstats),
                           ntree = 1000, paral = TRUE)
      posterior_lambda_all = predict(RF_lambda, target_sumstats_all,
                                     training = data.frame(param,sumstats),
                                     paral = TRUE, rf.weights = F)
      skyline_all$median[i] = posterior_lambda_all$med
      skyline_all$q_0.025[i] = posterior_lambda_all$quantiles[1]
      skyline_all$q_0.975[i] = posterior_lambda_all$quantiles[2]
      skyline_all$error[i] = RF_lambda$model.rf$prediction.error
      
      posterior_lambda_nowe = predict(RF_lambda, target_sumstats_nowe,
                                      training = data.frame(param,sumstats),
                                      paral = TRUE, rf.weights = F)
      skyline_nowe$median[i] = posterior_lambda_nowe$med
      skyline_nowe$q_0.025[i] = posterior_lambda_nowe$quantiles[1]
      skyline_nowe$q_0.975[i] = posterior_lambda_nowe$quantiles[2]
      skyline_nowe$error[i] = RF_lambda$model.rf$prediction.error
      
      posterior_lambda_short = predict(RF_lambda, target_sumstats_short,
                                       training = data.frame(param,sumstats),
                                       paral = TRUE, rf.weights = F)
      skyline_short$median[i] = posterior_lambda_short$med
      skyline_short$q_0.025[i] = posterior_lambda_short$quantiles[1]
      skyline_short$q_0.975[i] = posterior_lambda_short$quantiles[2]
      skyline_short$error[i] = RF_lambda$model.rf$prediction.error
      
      
      rm(param, RF_lambda, posterior_lambda_short, posterior_lambda_nowe, posterior_lambda_all)
      gc()
      save(skyline_short, skyline_nowe, skyline_all, file = results_skyline_file)
    }
    save(skyline_short, skyline_nowe, skyline_all, file = results_skyline_file)
  }
}



plot(skyline_years,
     skyline_all$median,
     #10^skyline_lambda$median,
     ylab = expression(log[10]*lambda),
     #ylab = expression(lambda),
     xlab = "cal BP",
     xlim = time_range_calBP,
     #ylim = c(0,1.5),
     ylim = c(-2,0.2),
     type="l",
     lwd=2)
lines(skyline_years,
      skyline_all$q_0.025,
      #10^skyline_lambda$q_0.025,
      lty=2)
lines(skyline_years,
      skyline_all$q_0.975,
      #10^skyline_lambda$q_0.975,
      lty=2)

lines(skyline_years,
      skyline_short$median,
     col="red",
     lwd=2)
lines(skyline_years,
      skyline_short$q_0.025,
      col="red",
      lty=2)
lines(skyline_years,
      skyline_short$q_0.975,
      col="red",
      lty=2)

lines(skyline_years,
      skyline_nowe$median,
      col="blue",
      lwd=2)
lines(skyline_years,
      skyline_nowe$q_0.025,
      col="blue",
      lty=2)
lines(skyline_years,
      skyline_nowe$q_0.975,
      col="blue",
      lty=2)


lines(spd_from_feature$grid$calBP,
      log10(spd_from_feature$grid$PrDens),
      #(spd_from_feature$grid$PrDens),
      col="gray")





sum(dates$Dating.from.feature)
old_wood_sample_size = sum(dates$Dating.from.feature & dates$Old.wood.effect.present)
short_term_sample_size = sum(dates$Dating.from.feature & dates$Short.term.human.impact)
proportion_old_wood = old_wood_sample_size/sum(dates$Dating.from.feature)
proportion_short_term = short_term_sample_size/sum(dates$Dating.from.feature)


# generate reference table for piecewise exponential model, 2 samples.
lambda_min = 0.001
lambda_max = 2
num_of_periods = 10
num_of_sims = 50000
malmo_reftable_no_old_wood = paste0(results_directory, "/malmo_reftable_no_old_wood_effect.rda")
if (!file.exists(malmo_reftable_no_old_wood) ){
  ncores = 10
  require(doSNOW)
  require(doParallel)
  require(doRNG)
  # setup parallel computing
  cl <- makeCluster(ncores, type = "FORK")  
  registerDoParallel(cl)  
  registerDoRNG(seed = 1234567)
  reftable <- foreach(sim = seq_len(num_of_sims), .combine = rbind) %dopar% {
  #for (i in 1:1000){
    rm(params, lambda_t, sim_dates, sim_dates_CRA, sumstats,
       old_wood_samples, sim_dates_old_wood, sim_dates_short_term,
       sim_dates_CRA_old_wood, sim_dates_CRA_short_term,
       sumstats_old_wood, sumstats_short_term, cor_sumstats)
    gc()
    
    params = sample_exponential_piecewise_parameters_from_priors(lambda_min, lambda_max,
                                                                 time_range = time_range_calBP,
                                                                 num_of_periods)
    prop_old_wood = rbeta(1, 45, 55)
    old_wood_effect = 0
    mean_old_wood_effect = 0
    lambda_t = get_piecewise_exponential_lambda_t(as.numeric(params[seq_len(num_of_periods+1)]),
                                                  as.integer(params[1+num_of_periods+seq_len(num_of_periods+1)]),
                                                  as.numeric(params[2+2*num_of_periods+seq_len(num_of_periods)]))
    sim_dates = sim_dates_lambda(lambda_t, time_range_calBP)
    old_wood_samples = sample(c(T,F),length(sim_dates),replace=T,prob=c(prop_old_wood,1-prop_old_wood))
    sim_dates_old_wood = sim_dates[old_wood_samples]
    sim_dates_short_term = sim_dates[!old_wood_samples]
    
    sim_dates_CRA_old_wood = sim_CRA(sim_dates_old_wood, errors = dates$X.[dates$Dating.from.feature])
    sim_dates_CRA_short_term = sim_CRA(sim_dates_short_term, errors = dates$X.[dates$Dating.from.feature])
    sumstats_old_wood = get_sumstats(sim_dates_CRA_old_wood, time_range_CRA,  window = c(25, 50, 100))
    sumstats_short_term = get_sumstats(sim_dates_CRA_short_term, time_range_CRA,  window = c(25, 50, 100))
    
    cor_sumstats = get_sumstats_correlation(sumstats_short_term,sumstats_old_wood, time_range_CRA,  window = c(25, 50, 100))
    
    names(sumstats_old_wood) = paste0(names(sumstats_old_wood),"_OW")
    names(sumstats_short_term) = paste0(names(sumstats_short_term),"_ST")
    
    sumstats = cbind(sumstats_old_wood, sumstats_short_term, cor_sumstats)
    params = cbind(params, prop_old_wood, mean_old_wood_effect, old_wood_effect)
    
    cbind(params, sumstats)
  }
  stopCluster(cl) 
  save(reftable, file = malmo_reftable_no_old_wood)
}else{
  load(file = malmo_reftable_no_old_wood)
}




# generate reference table for piecewise exponential model, 2 samples with old wood effect.
lambda_min = 0.001
lambda_max = 2
num_of_periods = 10
num_of_sims = 50000
malmo_reftable_old_wood = paste0(results_directory, "/malmo_reftable_old_wood_effect.rda")
if (!file.exists(malmo_reftable_old_wood) ){
  ncores = 15
  require(doSNOW)
  require(doParallel)
  require(doRNG)
  # setup parallel computing
  cl <- makeCluster(ncores, type = "FORK")  
  registerDoParallel(cl)  
  registerDoRNG(seed = 1234567)
  reftable <- foreach(sim = seq_len(num_of_sims), .combine = rbind) %dopar% {
    #for (i in 1:1000){
    rm(params, lambda_t, sim_dates, sim_dates_CRA, sumstats,
       old_wood_samples, sim_dates_old_wood, sim_dates_short_term,
       sim_dates_CRA_old_wood, sim_dates_CRA_short_term,
       sumstats_old_wood, sumstats_short_term, cor_sumstats)
    gc()
    
    params = sample_exponential_piecewise_parameters_from_priors(lambda_min, lambda_max,
                                                                 time_range = time_range_calBP,
                                                                 num_of_periods)
    prop_old_wood = rbeta(1, 45, 55)
    p_nbinom_owe = r_nbinom_owe = -1
    while ( p_nbinom_owe<=0 | p_nbinom_owe>=1 | r_nbinom_owe<=0 ){
      mean_old_wood_effect = runif(1,0,2000)
      var_old_wood_effect = 10^runif(1,0,10)
      p_nbinom_owe = mean_old_wood_effect/var_old_wood_effect
      r_nbinom_owe = mean_old_wood_effect*p_nbinom_owe/(1-p_nbinom_owe)
    }
    lambda_t = get_piecewise_exponential_lambda_t(as.numeric(params[seq_len(num_of_periods+1)]),
                                                  as.integer(params[1+num_of_periods+seq_len(num_of_periods+1)]),
                                                  as.numeric(params[2+2*num_of_periods+seq_len(num_of_periods)]))
    sim_dates = sim_dates_lambda(lambda_t, time_range_calBP)
    old_wood_samples = sample(c(T,F),length(sim_dates),replace=T,prob=c(prop_old_wood,1-prop_old_wood))
    sim_dates_old_wood = sim_dates[old_wood_samples]
    lag = rnbinom(length(sim_dates_old_wood), size = r_nbinom_owe, p = p_nbinom_owe)
    old_wood_effect = mean(lag)
    sim_dates_old_wood = sim_dates_old_wood + lag
    sim_dates_short_term = sim_dates[!old_wood_samples]
    
    sim_dates_CRA_old_wood = sim_CRA(sim_dates_old_wood, errors = dates$X.[dates$Dating.from.feature])
    sim_dates_CRA_short_term = sim_CRA(sim_dates_short_term, errors = dates$X.[dates$Dating.from.feature])
    sumstats_old_wood = get_sumstats(sim_dates_CRA_old_wood, time_range_CRA,  window = c(25, 50, 100))
    sumstats_short_term = get_sumstats(sim_dates_CRA_short_term, time_range_CRA,  window = c(25, 50, 100))
    
    cor_sumstats = get_sumstats_correlation(sumstats_short_term,sumstats_old_wood, time_range_CRA,  window = c(25, 50, 100))
    
    names(sumstats_old_wood) = paste0(names(sumstats_old_wood),"_OW")
    names(sumstats_short_term) = paste0(names(sumstats_short_term),"_ST")
    
    sumstats = cbind(sumstats_old_wood, sumstats_short_term, cor_sumstats)
    params = cbind(params, prop_old_wood, mean_old_wood_effect, var_old_wood_effect, r_nbinom_owe, p_nbinom_owe, old_wood_effect)
    
    cbind(params, sumstats)
  }
  stopCluster(cl) 
  save(reftable, file = malmo_reftable_old_wood)
}else{
  load(file = malmo_reftable_old_wood)
}

##### HERE


# samples marked as OLD WOOD as samples with OWE
dates_CRA_old_wood = list(CRA = dates$BP[dates$Dating.from.feature & dates$Old.wood.effect.present])
dates_CRA_short_term = list(CRA = dates$BP[dates$Dating.from.feature & !dates$Old.wood.effect.present])
sumstats_old_wood = get_sumstats(dates_CRA_old_wood, time_range_CRA,  window = c(25, 50, 100))
sumstats_short_term = get_sumstats(dates_CRA_short_term, time_range_CRA,  window = c(25, 50, 100))
cor_sumstats = get_sumstats_correlation(sumstats_short_term,sumstats_old_wood, time_range_CRA,  window = c(25, 50, 100))
names(sumstats_old_wood) = paste0(names(sumstats_old_wood),"_OW")
names(sumstats_short_term) = paste0(names(sumstats_short_term),"_ST")
target_sumstats_1 = cbind(sumstats_old_wood, sumstats_short_term, cor_sumstats)

# samples NOT marked as SHORT TERM as samples with OWE
dates_CRA_old_wood = list(CRA = dates$BP[dates$Dating.from.feature & !dates$Short.term.human.impact])
dates_CRA_short_term = list(CRA = dates$BP[dates$Dating.from.feature & dates$Short.term.human.impact])
sumstats_old_wood = get_sumstats(dates_CRA_old_wood, time_range_CRA,  window = c(25, 50, 100))
sumstats_short_term = get_sumstats(dates_CRA_short_term, time_range_CRA,  window = c(25, 50, 100))
cor_sumstats = get_sumstats_correlation(sumstats_short_term,sumstats_old_wood, time_range_CRA,  window = c(25, 50, 100))
names(sumstats_old_wood) = paste0(names(sumstats_old_wood),"_OW")
names(sumstats_short_term) = paste0(names(sumstats_short_term),"_ST")
target_sumstats_2 = cbind(sumstats_old_wood, sumstats_short_term, cor_sumstats)

# samples NOT marked as SHORT TERM as samples with OWE minus samples marked as OLD WOOD
dates_CRA_old_wood = list(CRA = dates$BP[dates$Dating.from.feature & !dates$Short.term.human.impact & !dates$Old.wood.effect.present])
dates_CRA_short_term = list(CRA = dates$BP[dates$Dating.from.feature & dates$Short.term.human.impact])
sumstats_old_wood = get_sumstats(dates_CRA_old_wood, time_range_CRA,  window = c(25, 50, 100))
sumstats_short_term = get_sumstats(dates_CRA_short_term, time_range_CRA,  window = c(25, 50, 100))
cor_sumstats = get_sumstats_correlation(sumstats_short_term,sumstats_old_wood, time_range_CRA,  window = c(25, 50, 100))
names(sumstats_old_wood) = paste0(names(sumstats_old_wood),"_OW")
names(sumstats_short_term) = paste0(names(sumstats_short_term),"_ST")
target_sumstats_3 = cbind(sumstats_old_wood, sumstats_short_term, cor_sumstats)



#### HERE

load(file = malmo_reftable_old_wood)
reftable_old_wood = reftable
load(file = malmo_reftable_no_old_wood)
reftable_no_old_wood = reftable
rm(reftable); gc()

reftable_sumstats = rbind(reftable_old_wood[names(target_sumstats_1)],
                          reftable_no_old_wood[names(target_sumstats_1)])

model = as.factor(c(rep("OW", nrow(reftable_old_wood)),
                    rep("NOW",nrow(reftable_no_old_wood))))

reftable = na.omit(data.frame(model=model, reftable_sumstats))


malmo_model_choice_file = paste0(results_directory, "/malmo_model_choice.rda")
if (!file.exists(malmo_model_choice_file) ){
  model_choice_RF = abcrf(model~., data = reftable, ntree = 2000, paral = T)
  model_choice_res_1 = predict(model_choice_RF, target_sumstats_1, reftable, ntree = 2000, paral = T)
  model_choice_res_2 = predict(model_choice_RF, target_sumstats_2, reftable, ntree = 2000, paral = T)
  model_choice_res_3 = predict(model_choice_RF, target_sumstats_3, reftable, ntree = 2000, paral = T)
  save(model_choice_RF, model_choice_res_1, model_choice_res_2, model_choice_res_3, file=malmo_model_choice_file)
}else{
  load(file=malmo_model_choice_file)
}
(model_choice_res_1)
interpret_K(model_choice_res_1$post.prob/(1-model_choice_res_1$post.prob))
(model_choice_res_2)
interpret_K(model_choice_res_2$post.prob/(1-model_choice_res_2$post.prob))
(model_choice_res_3)
interpret_K(model_choice_res_3$post.prob/(1-model_choice_res_3$post.prob))

##### HERE

owe = reftable_old_wood$mean_old_wood_effect
reftable_sumstats = reftable_old_wood[names(target_sumstats_1)]
reftable = na.omit(data.frame(owe=owe, reftable_sumstats))

malmo_old_wood_efect_estimate_file = paste0(results_directory, "/malmo_old_wood_effect_estimate.rda")
if (!file.exists(malmo_old_wood_efect_estimate_file) ){
  owe_RF = regAbcrf(owe~., reftable,
                    ntree = 1000, paral = TRUE)
  posterior_owe = predict(owe_RF, target_sumstats_2,
                             training = reftable,
                             paral = TRUE, rf.weights = T)

  save(owe_RF, posterior_owe, file=malmo_old_wood_efect_estimate_file)
}else{
  load(file=malmo_old_wood_efect_estimate_file)
}
(posterior_owe)
hist(owe)
weights::wtd.hist(x=reftable$owe, weight = posterior_owe$weights)


##### HERE



var_owe = log10(reftable_old_wood$var_old_wood_effect)
reftable = na.omit(data.frame(var_owe=var_owe, reftable_sumstats))



malmo_var_old_wood_efect_estimate_file = paste0(results_directory, "/malmo_var_old_wood_effect_estimate.rda")
if (!file.exists(malmo_var_old_wood_efect_estimate_file) ){
  var_owe_RF = regAbcrf(var_owe~., reftable,
                    ntree = 1000, paral = TRUE)
  posterior_var_owe = predict(var_owe_RF, target_sumstats_2,
                          training = reftable,
                          paral = TRUE, rf.weights = T)
  
  save(var_owe_RF, posterior_var_owe, file=malmo_var_old_wood_efect_estimate_file)
}else{
  load(file=malmo_var_old_wood_efect_estimate_file)
}
(posterior_var_owe)
hist(var_owe)
weights::wtd.hist(x=reftable$var_owe, weight = posterior_var_owe$weights)



mean_old_wood_effect_hat = posterior_owe$med
var_old_wood_effect_hat = 10^posterior_var_owe$med
p_hat = mean_old_wood_effect_hat/var_old_wood_effect_hat
r_hat = mean_old_wood_effect_hat*p_hat/(1-p_hat)
plot(dnbinom(1:2000,r_hat,p_hat),type="l",#log="x",
     xlab="years", ylab="probability density", main="old wood effect")
abline(h=0,col="grey",lty=2)


#### HERE



skyline_step = 25
skyline_points = 1 + (time_range_calBP[1]-time_range_calBP[2])/skyline_step
skyline_years = seq(time_range_calBP[1],time_range_calBP[2],-skyline_step)
malmo_reftable_skyline_owe = paste0(results_directory, "/malmo_reftable_skyline_owe.rda")
if (!file.exists(malmo_reftable_skyline_owe) ){
  reftable_skyline_plot_owe = data.frame(matrix(NA, nrow=nrow(reftable_old_wood), ncol=skyline_points*2))
  names(reftable_skyline_plot_owe) =  c(paste0("lambda",skyline_years),
                                    paste0("rate",skyline_years))
  
  for (sim in seq_len(num_of_sims)){
    print(paste("sim",sim))
    lambdas         = as.numeric(reftable_old_wood[sim,seq_len(num_of_periods+1)])
    times_of_change = as.integer(reftable_old_wood[sim,1+num_of_periods+seq_len(num_of_periods+1)])
    rates           = as.numeric(reftable_old_wood[sim,2+2*num_of_periods+seq_len(num_of_periods)])
    
    lambda_skyline = get_piecewise_exponential_lambda_t(lambdas,
                                                        times_of_change,
                                                        rates)[seq(1,skyline_points*skyline_step,skyline_step)]
    
    
    rates_skyline = c(rep(rates,abs(diff(times_of_change))), rates[length(rates)])[seq(1,skyline_points*skyline_step,skyline_step)]
    
    reftable_skyline_plot_owe[sim,] =  c(lambda_skyline,rates_skyline)
  }
  
  save(reftable_skyline_plot_owe, file = malmo_reftable_skyline_owe)
}else{
  load(file = malmo_reftable_skyline_owe)
}


##### HERE








results_skyline_owe_file = paste0(results_directory, "/skyline_owe_results.rda")
if (!file.exists(results_skyline_owe_file) ){
  
  skyline_1 = list(median  = rep(NA, skyline_points),
                   q_0.025 = rep(NA, skyline_points),
                   q_0.975 = rep(NA, skyline_points),
                   error   = rep(NA, skyline_points))
  skyline_2 = list(median  = rep(NA, skyline_points),
                   q_0.025 = rep(NA, skyline_points),
                   q_0.975 = rep(NA, skyline_points),
                   error   = rep(NA, skyline_points))
  skyline_3 = list(median  = rep(NA, skyline_points),
                   q_0.025 = rep(NA, skyline_points),
                   q_0.975 = rep(NA, skyline_points),
                   error   = rep(NA, skyline_points))

  sumstats = reftable_old_wood[, names(target_sumstats_1)]
  
  for (i in seq_len(skyline_points)){
    param = log10(reftable_skyline_plot_owe[,i])
    RF_lambda = regAbcrf(param~., data.frame(param,sumstats)[complete.cases(data.frame(param,sumstats)),],
                         ntree = 1000, paral = TRUE)
    posterior_lambda_1 = predict(RF_lambda, target_sumstats_1,
                                 training = data.frame(param,sumstats)[complete.cases(data.frame(param,sumstats)),],
                                 paral = TRUE, rf.weights = F)
    skyline_1$median[i] = posterior_lambda_1$med
    skyline_1$q_0.025[i] = posterior_lambda_1$quantiles[1]
    skyline_1$q_0.975[i] = posterior_lambda_1$quantiles[2]
    skyline_1$error[i] = RF_lambda$model.rf$prediction.error
    
    posterior_lambda_2 = predict(RF_lambda, target_sumstats_2,
                                 training = data.frame(param,sumstats)[complete.cases(data.frame(param,sumstats)),],
                                 paral = TRUE, rf.weights = F)
    skyline_2$median[i] = posterior_lambda_2$med
    skyline_2$q_0.025[i] = posterior_lambda_2$quantiles[1]
    skyline_2$q_0.975[i] = posterior_lambda_2$quantiles[2]
    skyline_2$error[i] = RF_lambda$model.rf$prediction.error
    
    posterior_lambda_3 = predict(RF_lambda, target_sumstats_3,
                                 training = data.frame(param,sumstats)[complete.cases(data.frame(param,sumstats)),],
                                 paral = TRUE, rf.weights = F)
    skyline_3$median[i] = posterior_lambda_3$med
    skyline_3$q_0.025[i] = posterior_lambda_3$quantiles[1]
    skyline_3$q_0.975[i] = posterior_lambda_3$quantiles[2]
    skyline_3$error[i] = RF_lambda$model.rf$prediction.error
    
    

    
    rm(param, RF_lambda, posterior_lambda_low_error, posterior_lambda_1, posterior_lambda_2, posterior_lambda_3)
    gc()
    save(skyline_1, skyline_2, skyline_3, file = results_skyline_owe_file)
  }
  save(skyline_1, skyline_2, skyline_3, file = results_skyline_owe_file)
}else{
  load(file = results_skyline_owe_file)
  
  if (any(is.na(skyline_3$median))){
    sumstats = reftable_old_wood[, names(target_sumstats_1)]
    initial_i = which(is.na(skyline_1$median))[1]
    for (i in initial_i:skyline_points){
      print(i)
      param = log10(reftable_skyline_plot_owe[,i])
      RF_lambda = regAbcrf(param~., data.frame(param,sumstats)[complete.cases(data.frame(param,sumstats)),],
                           ntree = 1000, paral = TRUE)
      posterior_lambda_1 = predict(RF_lambda, target_sumstats_1,
                                   training = data.frame(param,sumstats)[complete.cases(data.frame(param,sumstats)),],
                                   paral = TRUE, rf.weights = F)
      skyline_1$median[i] = posterior_lambda_1$med
      skyline_1$q_0.025[i] = posterior_lambda_1$quantiles[1]
      skyline_1$q_0.975[i] = posterior_lambda_1$quantiles[2]
      skyline_1$error[i] = RF_lambda$model.rf$prediction.error
      
      posterior_lambda_2 = predict(RF_lambda, target_sumstats_2,
                                   training = data.frame(param,sumstats)[complete.cases(data.frame(param,sumstats)),],
                                   paral = TRUE, rf.weights = F)
      skyline_2$median[i] = posterior_lambda_2$med
      skyline_2$q_0.025[i] = posterior_lambda_2$quantiles[1]
      skyline_2$q_0.975[i] = posterior_lambda_2$quantiles[2]
      skyline_2$error[i] = RF_lambda$model.rf$prediction.error
      
      posterior_lambda_3 = predict(RF_lambda, target_sumstats_3,
                                   training = data.frame(param,sumstats)[complete.cases(data.frame(param,sumstats)),],
                                   paral = TRUE, rf.weights = F)
      skyline_3$median[i] = posterior_lambda_3$med
      skyline_3$q_0.025[i] = posterior_lambda_3$quantiles[1]
      skyline_3$q_0.975[i] = posterior_lambda_3$quantiles[2]
      skyline_3$error[i] = RF_lambda$model.rf$prediction.error
      
      
      
      
      rm(param, RF_lambda, posterior_lambda_low_error, posterior_lambda_1, posterior_lambda_2, posterior_lambda_3)
      gc()
      save(skyline_1, skyline_2, skyline_3, file = results_skyline_owe_file)
    }
    save(skyline_1, skyline_2, skyline_3, file = results_skyline_owe_file)
  }
}


#### HERE


load(file = results_skyline_file)
plot(skyline_years,
     skyline_all$median,
     #10^skyline_lambda$median,
     ylab = expression(log[10]*lambda),
     #ylab = expression(lambda),
     xlab = "cal BP",
     xlim = time_range_calBP,
     #ylim = c(0,1.5),
     ylim = c(-2,0.2),
     type="l",
     lwd=2)
lines(skyline_years,
      skyline_all$q_0.025,
      #10^skyline_lambda$q_0.025,
      lty=2)
lines(skyline_years,
      skyline_all$q_0.975,
      #10^skyline_lambda$q_0.975,
      lty=2)

load(file = results_skyline_owe_file)

lines(skyline_years,
      skyline_1$median,
      col="blue",
      lwd=2)
lines(skyline_years,
      skyline_1$q_0.025,
      col="blue",
      lty=2)
lines(skyline_years,
      skyline_1$q_0.975,
      col="blue",
      lty=2)

plot(skyline_years,
     skyline_all$median,
     #10^skyline_lambda$median,
     ylab = expression(log[10]*lambda),
     #ylab = expression(lambda),
     xlab = "cal BP",
     xlim = time_range_calBP,
     #ylim = c(0,1.5),
     ylim = c(-2,0.2),
     type="l",
     lwd=2)
lines(skyline_years,
      skyline_all$q_0.025,
      #10^skyline_lambda$q_0.025,
      lty=2)
lines(skyline_years,
      skyline_all$q_0.975,
      #10^skyline_lambda$q_0.975,
      lty=2)

load(file = results_skyline_owe_file)

lines(skyline_years,
      skyline_2$median,
      col="blue",
      lwd=2)
lines(skyline_years,
      skyline_2$q_0.025,
      col="blue",
      lty=2)
lines(skyline_years,
      skyline_2$q_0.975,
      col="blue",
      lty=2)

plot(skyline_years,
     skyline_all$median,
     #10^skyline_lambda$median,
     ylab = expression(log[10]*lambda),
     #ylab = expression(lambda),
     xlab = "cal BP",
     xlim = time_range_calBP,
     #ylim = c(0,1.5),
     ylim = c(-2,0.2),
     type="l",
     lwd=2)
lines(skyline_years,
      skyline_all$q_0.025,
      #10^skyline_lambda$q_0.025,
      lty=2)
lines(skyline_years,
      skyline_all$q_0.975,
      #10^skyline_lambda$q_0.975,
      lty=2)

load(file = results_skyline_owe_file)

lines(skyline_years,
      skyline_3$median,
      col="blue",
      lwd=2)
lines(skyline_years,
      skyline_3$q_0.025,
      col="blue",
      lty=2)
lines(skyline_years,
      skyline_3$q_0.975,
      col="blue",
      lty=2)





