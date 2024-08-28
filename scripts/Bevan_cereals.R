########################################################################
# Analysis of real data from Britain and Ireland: wheat and barley
########################################################################
source("scripts/DARthABC.R")
require(abcrf)
require(doSNOW)
require(doParallel)
require(doRNG)
results_all_directory = "results/Bevan_all"
results_cereals_directory = "results/Bevan_cereals"
dir.create(results_cereals_directory)

# color
PCI_blue = rgb(44, 110.4, 148.3, 255, maxColorValue = 255)
PCI_t_blue = rgb(44, 110.4, 148.3, 100, maxColorValue = 255)
require(graphics)
cbPalette1 <- c("#000000", # Black
                "#E69F00", # Orange
                "#56B4E9", # Sky Blue
                "#009E73", # Bluish Green
                "#F0E442", # Yellow
                "#0072B2", # Blue
                "#D55E00", # Vermillion
                "#CC79A7") # Reddish Purple


# set time range for analysis
time_range_BP = c(6000,500)
t_ = seq(time_range_BP[1], time_range_BP[2] + 1,by = -1)

# set number of simulations
num_of_sims = 30000

# set prior parameters
lambda_min = 0.001
lambda_max = 2


# filter data for wheat and barley entries and calibrate data subset
dates_file = paste0(results_cereals_directory, "/dates.rda")
caldates_file = paste0(results_cereals_directory, "/caldates.rda")
spd_file = paste0(results_cereals_directory, "/spd.rda")
if (!file.exists(dates_file) |
    !file.exists(caldates_file) |
    !file.exists(spd_file)){

  load(file = paste0(results_all_directory,"/dates.rda"))
  load(file = paste0(results_all_directory,"/caldates.rda"))
  
  entries2exclude = c("Hordeum/Triticum",
                      "Triticum/Hordeum vulgare ",
                      "Hordeum vulgare/Triticum spelta",
                      "Hordeum vulgare/Triticum",
                      "Avena/Triticum")
  
  # Wheat
  triticum_entries = grepl("Triticum", dates$Species)
  triticum_entries = triticum_entries & !dates$Species %in% entries2exclude
  triticum_dates = dates[triticum_entries,]
  triticum_caldates = caldates[triticum_entries]
  triticum_bins = binPrep(sites = triticum_dates$SiteID, ages = triticum_dates$CRA, h = 100)
  triticum_spd = spd(x = triticum_caldates,
                     bins = triticum_bins, 
                     timeRange = time_range_BP, 
                     datenormalised = FALSE,
                     runm = 100)
  # Barley
  hordeum_entries = grepl("Hordeum",dates$Species)
  hordeum_entries = hordeum_entries & !dates$Species %in% entries2exclude
  hordeum_dates = dates[hordeum_entries,]
  hordeum_caldates = caldates[hordeum_entries]
  hordeum_bins = binPrep(sites = hordeum_dates$SiteID, ages = hordeum_dates$CRA, h = 100)
  hordeum_spd = spd(x = hordeum_caldates,
                    bins = hordeum_bins, 
                    timeRange = time_range_BP, 
                    datenormalised = FALSE,
                    runm = 100)
  # pooled together
  cereals_entries = c(which(triticum_entries),which(hordeum_entries))
  cereals_dates = dates[cereals_entries,]
  cereals_caldates = caldates[cereals_entries]
  cereals_spd = spd(x = cereals_caldates,
                    bins = c(triticum_bins, hordeum_bins),
                    timeRange = time_range_BP, 
                    datenormalised = FALSE,
                    runm = 100)
  plot(cereals_spd)
  
  save(cereals_dates, triticum_dates, hordeum_dates, file = dates_file)
  save(cereals_caldates, triticum_caldates, hordeum_caldates, file = caldates_file)
  save(cereals_spd, triticum_bins, triticum_spd, hordeum_bins, hordeum_spd, file = spd_file)
}
  


spd_plot_file = paste0(results_cereals_directory, "/spd.pdf")
hist_plot_file = paste0(results_cereals_directory, "/hist.pdf")
if (!file.exists(spd_plot_file) |
    !file.exists(hist_plot_file) ){

  load(file = spd_file)
  
  pdf(file=spd_plot_file, width=10, height=5)
  par(mar=c(4.5, 4.5, 1, 1) + 0.1)
  plot(hordeum_spd$grid$calBP, hordeum_spd$grid$PrDens,
       xlab="Years cal BP", ylab="density",
       xlim = time_range_BP, ylim=c(0,0.4),
       type="l",lwd=2, lty=2)
  lines(triticum_spd$grid$calBP, triticum_spd$grid$PrDens, col=PCI_blue, lwd=2, lty=1)
  legend("topleft",
         legend=c(expression(italic(Triticum)),expression(italic(Hordeum))),
         col=c(PCI_blue,"black"),
         lwd=2,
         lty=c(1,2))
  dev.off()
  
  
  load(file = dates_file)
  
  pdf(file=hist_plot_file, width=10, height=5)
  par(mar=c(4.5, 4.5, 1, 1) + 0.1)
  hordeum_weights = apply(as.array(hordeum_bins), 1, function(x) 1/sum(hordeum_bins == x))
  breaks = seq(6000,200,-100)
  HH = wtd.hist(hordeum_dates$CRA,
                breaks=breaks, 
                plot=T, weight = hordeum_weights,
                main="",
                xlim=c(6000,400), 
                xlab="Years BP (uncalibrated)",col="gray")
  
  triticum_weights = apply(as.array(triticum_bins), 1, function(x) 1/sum(triticum_bins == x))
  #breaks = seq(time_range_BP[1],time_range_BP[2]-1,-100)
  TH = wtd.hist(triticum_dates$CRA, breaks=breaks, plot=T, weight = triticum_weights,
                main="",xlim=c(6000,400), xlab="Years BP (uncalibrated)",col=PCI_t_blue,add=T)
  box()
  dev.off()
}

sumstats_windows = c(10,50,100,500)

# Calculate Summary statstics
sumstats_file = paste0(results_cereals_directory, "/sumstats.rda")
if (!file.exists(sumstats_file)){
  load(file = dates_file)
  
  triticum_sumstats_c14 = get_sumstats(triticum_dates["CRA"],
                                       time_range_BP,
                                       window = sumstats_windows,
                                       w=triticum_weights)
    
  hordeum_sumstats_c14 = get_sumstats(hordeum_dates["CRA"],
                                      time_range_BP, 
                                      window = sumstats_windows,
                                      w=hordeum_weights)
  cereals_sumstats_c14 = get_sumstats(cereals_dates["CRA"], 
                                      time_range_BP, 
                                      window = sumstats_windows,
                                      w=c(triticum_weights, hordeum_weights))
  cereals_sumstats_cor = get_sumstats_correlation(triticum_sumstats_c14,
                                                  hordeum_sumstats_c14,
                                                  time_range_BP,
                                                  window = sumstats_windows)
  names(triticum_sumstats_c14) = paste0(names(triticum_sumstats_c14),"_A")
  names(hordeum_sumstats_c14) = paste0(names(hordeum_sumstats_c14),"_B")
  
  cereals_sumstats_2_categories = cbind(triticum_sumstats_c14, hordeum_sumstats_c14,
                                        cereals_sumstats_c14, cereals_sumstats_cor)
  
  save(cereals_sumstats_2_categories, file=sumstats_file)
}









num_of_periods = round((time_range_BP[1] - time_range_BP[2]) / 100 / 4)
skyline_years = get_time_of_change(num_of_periods, time_range_BP, 
                                   intervals = "regular")

# make reference tables for piecewise exponential model : independent scenario
reftable_independent_file = paste0(results_cereals_directory, "/independent_scenario_reftable.rda")
if (!file.exists(reftable_independent_file) ){
  load(file = dates_file)
  # setup parallel computing
  ncores = 5
  cl <- makeCluster(ncores, type="FORK")  
  registerDoParallel(cl)  
  registerDoRNG(seed = sample(1:1000000,1) )
  reftable <- foreach(sim=seq_len(num_of_sims), .combine=rbind) %dopar% {
    gc()
    
    params = sample_exponential_piecewise_parameters_2_categories(lambda_min, lambda_max,
                                                                  time_range_BP,
                                                                  num_of_periods,
                                                                  intervals = "regular",
                                                                  scenario = "independent",
                                                                  alpha = NULL, beta = NULL)    
    
    lambda_t_A = get_piecewise_exponential_lambda_t(params$lambda_values_A,
                                                    params$time_of_change,
                                                    params$growth_rates_A)
    lambda_t_B = get_piecewise_exponential_lambda_t(params$lambda_values_B,
                                                    params$time_of_change,
                                                    params$growth_rates_B)
    
    names(params$lambda_values_A) = paste0("lambda_A_",seq_len(num_of_periods+1))
    names(params$lambda_values_B) = paste0("lambda_B_",seq_len(num_of_periods+1))
    names(params$pi) = paste0("pi_",seq_len(num_of_periods+1))
    names(params$time_of_change) = c("t_0",paste0("t_",seq_len(num_of_periods-1)),"t_f")
    names(params$growth_rates_A) = paste0("r_A_",seq_len(num_of_periods))
    names(params$growth_rates_B) = paste0("r_B_",seq_len(num_of_periods))
    params = data.frame(t(c(params$lambda_values_A,
                            params$lambda_values_B,
                            params$pi,
                            params$time_of_change,
                            params$growth_rates_A,
                            params$growth_rates_B)))
    
    na_in_lambda_t = any(c(is.na(lambda_t_A), is.na(lambda_t_B)))
    if (na_in_lambda_t){
      ss        = rep(NA,ncol(cereals_sumstats_2_categories))
      names(ss) = names(cereals_sumstats_2_categories)
    }else{
      sim_dates_A = sim_dates_lambda(lambda_t_A, time_range_BP)
      sim_dates_CRA_A = sim_CRA(sim_dates_A, errors = triticum_dates$Error)
      
      sim_dates_B = sim_dates_lambda(lambda_t_B, time_range_BP)
      sim_dates_CRA_B = sim_CRA(sim_dates_B, errors = hordeum_dates$Error)
      
      
      A_sumstats_c14 = get_sumstats(sim_dates_CRA_A["CRA"],
                                           time_range_BP,
                                           window = sumstats_windows)
      B_sumstats_c14 = get_sumstats(sim_dates_CRA_B["CRA"],
                                          time_range_BP, 
                                          window = sumstats_windows)
      AplusB_sumstats_c14 = get_sumstats(rbind(sim_dates_CRA_A,sim_dates_CRA_B)["CRA"], 
                                          time_range_BP, 
                                          window = sumstats_windows)
      sumstats_cor = get_sumstats_correlation(A_sumstats_c14,
                                              B_sumstats_c14,
                                              time_range_BP,
                                              window = sumstats_windows)
      names(A_sumstats_c14) = paste0(names(A_sumstats_c14),"_A")
      names(B_sumstats_c14) = paste0(names(B_sumstats_c14),"_B")
      
      ss = cbind(A_sumstats_c14, B_sumstats_c14, AplusB_sumstats_c14, sumstats_cor)

    }
    cbind(params,ss)
    
  }
  save(reftable, file=reftable_independent_file)
  stopCluster(cl)
  rm(reftable); gc()
}
# load( file=reftable_independent_file)



reftable_interdependent_file = paste0(results_cereals_directory, "/interdependent_scenario_reftable.rda")
if (!file.exists(reftable_interdependent_file) ){
  num_of_sims = 100000
  
  load(file = dates_file)
  # setup parallel computing
  ncores = 5
  cl <- makeCluster(ncores, type="FORK")  
  registerDoParallel(cl)  
  registerDoRNG(seed = sample(1:1000000,1) )
  reftable <- foreach(sim=seq_len(num_of_sims), .combine=rbind) %dopar% {
    gc()
    
    params = sample_exponential_piecewise_parameters_2_categories(lambda_min, lambda_max,
                                                                  time_range_BP,
                                                                  num_of_periods,
                                                                  intervals = "regular",
                                                                  scenario = "interdependent",
                                                                  alpha = 1, beta = 1)    
    
    lambda_t_A = get_piecewise_exponential_lambda_t(params$lambda_values_A,
                                                    params$time_of_change,
                                                    params$growth_rates_A)
    lambda_t_B = get_piecewise_exponential_lambda_t(params$lambda_values_B,
                                                    params$time_of_change,
                                                    params$growth_rates_B)
    
    names(params$lambda_values_A) = paste0("lambda_A_",seq_len(num_of_periods+1))
    names(params$lambda_values_B) = paste0("lambda_B_",seq_len(num_of_periods+1))
    names(params$pi) = paste0("pi_",seq_len(num_of_periods+1))
    names(params$time_of_change) = c("t_0",paste0("t_",seq_len(num_of_periods-1)),"t_f")
    names(params$growth_rates_A) = paste0("r_A_",seq_len(num_of_periods))
    names(params$growth_rates_B) = paste0("r_B_",seq_len(num_of_periods))
    params = data.frame(t(c(params$lambda_values_A,
                            params$lambda_values_B,
                            params$pi,
                            params$time_of_change,
                            params$growth_rates_A,
                            params$growth_rates_B)))
    
    na_in_lambda_t = any(c(is.na(lambda_t_A), is.na(lambda_t_B)))
    if (na_in_lambda_t){
      ss        = rep(NA,ncol(cereals_sumstats_2_categories))
      names(ss) = names(cereals_sumstats_2_categories)
    }else{
      sim_dates_A = sim_dates_lambda(lambda_t_A, time_range_BP)
      sim_dates_CRA_A = sim_CRA(sim_dates_A, errors = triticum_dates$Error)
      
      sim_dates_B = sim_dates_lambda(lambda_t_B, time_range_BP)
      sim_dates_CRA_B = sim_CRA(sim_dates_B, errors = hordeum_dates$Error)
      
      
      A_sumstats_c14 = get_sumstats(sim_dates_CRA_A["CRA"],
                                    time_range_BP,
                                    window = sumstats_windows)
      B_sumstats_c14 = get_sumstats(sim_dates_CRA_B["CRA"],
                                    time_range_BP, 
                                    window = sumstats_windows)
      AplusB_sumstats_c14 = get_sumstats(rbind(sim_dates_CRA_A,sim_dates_CRA_B)["CRA"], 
                                         time_range_BP, 
                                         window = sumstats_windows)
      sumstats_cor = get_sumstats_correlation(A_sumstats_c14,
                                              B_sumstats_c14,
                                              time_range_BP,
                                              window = sumstats_windows)
      names(A_sumstats_c14) = paste0(names(A_sumstats_c14),"_A")
      names(B_sumstats_c14) = paste0(names(B_sumstats_c14),"_B")
      
      ss = cbind(A_sumstats_c14, B_sumstats_c14, AplusB_sumstats_c14, sumstats_cor)
      
    }
    cbind(params,ss)
    
  }
  save(reftable, file=reftable_interdependent_file)
  stopCluster(cl) 
  rm(reftable); gc()
}
# load( file=reftable_interdependent_file)










reftable_parallel_file = paste0(results_cereals_directory, "/parallel_scenario_reftable.rda")
if (!file.exists(reftable_parallel_file) ){
  load(file = dates_file)
  # setup parallel computing
  ncores = 5
  cl <- makeCluster(ncores, type="FORK")  
  registerDoParallel(cl)  
  registerDoRNG(seed = sample(1:1000000,1) )
  reftable <- foreach(sim=seq_len(num_of_sims), .combine=rbind) %dopar% {
    gc()
    
    params = sample_exponential_piecewise_parameters_2_categories(lambda_min, lambda_max,
                                                                  time_range_BP,
                                                                  num_of_periods,
                                                                  intervals = "regular",
                                                                  scenario = "parallel",
                                                                  alpha = 1, beta = 1)    
    
    lambda_t_A = get_piecewise_exponential_lambda_t(params$lambda_values_A,
                                                    params$time_of_change,
                                                    params$growth_rates_A)
    lambda_t_B = get_piecewise_exponential_lambda_t(params$lambda_values_B,
                                                    params$time_of_change,
                                                    params$growth_rates_B)
    
    names(params$lambda_values_A) = paste0("lambda_A_",seq_len(num_of_periods+1))
    names(params$lambda_values_B) = paste0("lambda_B_",seq_len(num_of_periods+1))
    names(params$pi) = paste0("pi_",seq_len(num_of_periods+1))
    names(params$time_of_change) = c("t_0",paste0("t_",seq_len(num_of_periods-1)),"t_f")
    names(params$growth_rates_A) = paste0("r_A_",seq_len(num_of_periods))
    names(params$growth_rates_B) = paste0("r_B_",seq_len(num_of_periods))
    params = data.frame(t(c(params$lambda_values_A,
                            params$lambda_values_B,
                            params$pi,
                            params$time_of_change,
                            params$growth_rates_A,
                            params$growth_rates_B)))
    
    na_in_lambda_t = any(c(is.na(lambda_t_A), is.na(lambda_t_B)))
    if (na_in_lambda_t){
      ss        = rep(NA,ncol(cereals_sumstats_2_categories))
      names(ss) = names(cereals_sumstats_2_categories)
    }else{
      sim_dates_A = sim_dates_lambda(lambda_t_A, time_range_BP)
      sim_dates_CRA_A = sim_CRA(sim_dates_A, errors = triticum_dates$Error)
      
      sim_dates_B = sim_dates_lambda(lambda_t_B, time_range_BP)
      sim_dates_CRA_B = sim_CRA(sim_dates_B, errors = hordeum_dates$Error)
      
      
      A_sumstats_c14 = get_sumstats(sim_dates_CRA_A["CRA"],
                                    time_range_BP,
                                    window = sumstats_windows)
      B_sumstats_c14 = get_sumstats(sim_dates_CRA_B["CRA"],
                                    time_range_BP, 
                                    window = sumstats_windows)
      AplusB_sumstats_c14 = get_sumstats(rbind(sim_dates_CRA_A,sim_dates_CRA_B)["CRA"], 
                                         time_range_BP, 
                                         window = sumstats_windows)
      sumstats_cor = get_sumstats_correlation(A_sumstats_c14,
                                              B_sumstats_c14,
                                              time_range_BP,
                                              window = sumstats_windows)
      names(A_sumstats_c14) = paste0(names(A_sumstats_c14),"_A")
      names(B_sumstats_c14) = paste0(names(B_sumstats_c14),"_B")
      
      ss = cbind(A_sumstats_c14, B_sumstats_c14, AplusB_sumstats_c14, sumstats_cor)
      
    }
    cbind(params,ss)
    
  }
  save(reftable, file=reftable_parallel_file)
  stopCluster(cl) 
  rm(reftable); gc()
}
# load( file=reftable_parallel_file)



model_choice_file = paste0(results_cereals_directory, "/model_choice.rda")
# test 3 models
if ( !file.exists(model_choice_file) ){
  load( file=reftable_independent_file)
  reftable_independent = reftable
  load( file=reftable_interdependent_file)
  reftable_interdependent = reftable
  load( file=reftable_parallel_file)
  reftable_parallel = reftable
  rm(reftable);gc()
  
  load(file=sumstats_file)

  reftable_interdependent = reftable_interdependent[complete.cases(reftable_interdependent),]
  reftable_parallel = reftable_parallel[complete.cases(reftable_parallel),]
  
  number_of_sims = min(c(nrow(reftable_independent),
                         nrow(reftable_interdependent),
                         nrow(reftable_parallel)))
  
  
  model = as.factor(c(rep("independent",number_of_sims),
                      rep("interdependent",number_of_sims),
                      rep("parallel",number_of_sims)))
  sumstats = rbind(reftable_independent[seq_len(number_of_sims), names(cereals_sumstats_2_categories)],
                   reftable_interdependent[seq_len(number_of_sims), names(cereals_sumstats_2_categories)],
                   reftable_parallel[seq_len(number_of_sims), names(cereals_sumstats_2_categories)])
  
  RF_model_choice = abcrf(model~., data = data.frame(model,sumstats),
                          lda=F, ntree = 2000, paral = TRUE)
  
  posterior_model_choice =  predict(RF_model_choice, cereals_sumstats_2_categories,
                                    training = data.frame(model,sumstats),
                                    ntree = 2000, paral = TRUE) 
  save(RF_model_choice, posterior_model_choice, file=model_choice_file)
  
  RF_model_choice$model.rf$confusion.matrix
}
#load(file=model_choice_file)



model_checking_file = paste0(results_cereals_directory, "/model_checking.rda")
# PCA for model checking
if ( !file.exists(model_checking_file) ){
  load( file=reftable_independent_file)
  reftable_independent = reftable
  load( file=reftable_interdependent_file)
  reftable_interdependent = reftable
  load( file=reftable_parallel_file)
  reftable_parallel = reftable
  rm(reftable);gc()
  
  reftable_parallel = reftable_parallel[complete.cases(reftable_parallel),]
  
  number_of_sims = min(c(nrow(reftable_independent),
                         nrow(reftable_interdependent),
                         nrow(reftable_parallel)))
  
  
  model = as.factor(c(rep("independent",number_of_sims),
                      rep("interdependent",number_of_sims),
                      rep("parallel",number_of_sims)))
  sumstats = rbind(reftable_independent[seq_len(number_of_sims), names(cereals_sumstats_2_categories)],
                   reftable_interdependent[seq_len(number_of_sims), names(cereals_sumstats_2_categories)],
                   reftable_parallel[seq_len(number_of_sims), names(cereals_sumstats_2_categories)])
  
  PCA_stats  = princomp(sumstats)
  PCA_target = predict(PCA_stats, cereals_sumstats_2_categories)
  save(PCA_stats, PCA_target, file = model_checking_file)
  
  summary(PCA_stats)
  # proportion of variance explained by first 6 PC
  sum((PCA_stats$sdev[1:6])^2)/sum((PCA_stats$sdev)^2)
  
  #load(file=model_checking_file)
  
  points_per_model = 3000
  sims2plot = sort(c(sample(which(model==unique(model)[1]),points_per_model),
                     sample(which(model==unique(model)[2]),points_per_model),
                     sample(which(model==unique(model)[3]),points_per_model)))
  colors2plot = c(rep(cbPalette1[2],points_per_model),
                  rep(cbPalette1[3],points_per_model),
                  rep(cbPalette1[4],points_per_model))
  order2plot = sample(points_per_model*3)
  sims2plot = sims2plot[order2plot]
  colors2plot = colors2plot[order2plot]
  

  
  for (axis_pair in 1:3){
    
    if (axis_pair == 1){
      i = 1; j = 2
      xlim = c(-50000, 100000); ylim = c(-50000, 50000); id_plot = "a"; add_legend=TRUE
    }
    if (axis_pair == 2){
      i = 3; j = 4
      xlim = c(-50000,80000); ylim=c(-30000,30000); id_plot="b"; add_legend=FALSE
    }
    if (axis_pair == 3){
      i = 5; j = 6
      xlim=c(-20000,30000); ylim=c(-5000,10000); id_plot="c"; add_legend=FALSE
    }
    
    pdf_file_name = paste0(results_cereals_directory, "/PCA_PC",i,"_PC",j,".pdf")
    if ( !file.exists(pdf_file_name) ){
      pdf(file=pdf_file_name, width=4.5, height=4.5)
      par(mar=c(4.5, 4.5, 1, 1) + 0.1)
      plot(PCA_stats$scores[sims2plot,i],
           PCA_stats$scores[sims2plot,j],
           xlab = paste0("PC",i),
           ylab = paste0("PC",j),
           xlim=xlim,
           ylim=ylim,
           col=colors2plot, cex=0.5)
      points(PCA_target[,i],PCA_target[,j],pch="*",cex=3)
      text(x=xlim[1], y=ylim[2], label=id_plot, cex=2)
      if (add_legend) legend("bottomright",
                             c("Independent","Interdependent","Parallel"),
                             pch=1,
                             col=cbPalette1[c(2,3,4)],
                             bg = "white")
      dev.off()
    }
  }
}




parameter_inference_file = paste0(results_cereals_directory, "/interdependent_parameter_estimates.rda")
# Inference under interdependent scenario 
if ( !file.exists(parameter_inference_file) ){
  load(file = sumstats_file)
  load(file = reftable_interdependent_file)
  
  reftable = reftable[complete.cases(reftable),]
  
  sumstats = rbind(reftable[names(cereals_sumstats_2_categories)])

  # total
  lambda_error = rep(NA,length(skyline_years))
  lambda_hat = rep(NA,length(skyline_years))
  lambda_95low = rep(NA,length(skyline_years))
  lambda_95upp = rep(NA,length(skyline_years))
  for (i in seq_along(skyline_years)){
    results_file = paste0(results_cereals_directory, "/interdependent_posterior_log_lambda_",
                          skyline_years[i], ".rda")
    if ( !file.exists(results_file) ){
      param_name_A = paste0("lambda_A_",i)
      param_name_B = paste0("lambda_B_",i)
      param_index_A = which(names(reftable)==param_name_A)
      param_index_B = which(names(reftable)==param_name_B)
      
      param = log10(reftable[param_index_A]+reftable[param_index_B])
      names(param) = "param"
      RF_log_lambda = regAbcrf(param~., data.frame(param,sumstats),
                               ntree = 1000, paral = TRUE)
      
      posterior_log_lambda = predict(RF_log_lambda, cereals_sumstats_2_categories,
                                     training = data.frame(param,sumstats),
                                     paral = TRUE, rf.weights = TRUE) 
      save(RF_log_lambda, posterior_log_lambda,
           file = results_file)
    }else{load(file = results_file)}
    gc()
    lambda_error[i] = RF_log_lambda$model.rf$prediction.error
    lambda_hat[i] = 10^(posterior_log_lambda$med[1])
    lambda_95low[i] = 10^(posterior_log_lambda$quantiles[1])
    lambda_95upp[i] = 10^(posterior_log_lambda$quantiles[2])
    save(lambda_error, lambda_hat,
         lambda_95low, lambda_95upp,
         file = parameter_inference_file)
  }
  

  # pi
  pi_error = rep(NA,length(skyline_years))
  pi_hat = rep(NA,length(skyline_years))
  pi_95low = rep(NA,length(skyline_years))
  pi_95upp = rep(NA,length(skyline_years))
  for (i in seq_along(skyline_years)){
    results_file = paste0(results_cereals_directory,"/interdependent_posterior_pi_",
                          skyline_years[i], ".rda")
    if ( !file.exists(results_file) ){
      param_name = paste0("pi_",i)
      param_index = which(names(reftable)==param_name)
      
      p = reftable[param_index]
      param = log(p/(1-p))
      names(param) = "param"
      RF_logit_pi = regAbcrf(param~., data.frame(param,sumstats),
                             ntree = 1000, paral = TRUE)
      
      posterior_logit_pi = predict(RF_logit_pi, cereals_sumstats_2_categories,
                                   training = data.frame(param,sumstats),
                                   paral = TRUE, rf.weights = TRUE) 
      
      save(RF_logit_pi, posterior_logit_pi,
           file = results_file)
    }else{load(file = results_file)}
    gc()    
    pi_error[i] = RF_logit_pi$model.rf$prediction.error
    pi_hat[i] = exp(posterior_logit_pi$med[1])/(1+exp(posterior_logit_pi$med[1]))  
    pi_95low[i] = exp(posterior_logit_pi$quantiles[1])/(1+exp(posterior_logit_pi$quantiles[1]))  
    pi_95upp[i] = exp(posterior_logit_pi$quantiles[2])/(1+exp(posterior_logit_pi$quantiles[2]))   
  }
  
  save(lambda_error, lambda_hat,
       lambda_95low, lambda_95upp,
       pi_error, pi_hat,
       pi_95low, pi_95upp,
       file = parameter_inference_file)
  
  
  
}



plot_lambda_results_file = paste0(results_cereals_directory, "/interdependent_model_result_lambda.pdf")
if ( !file.exists(plot_lambda_results_file) ){
  load(file = spd_file)
  load(file = parameter_inference_file)
  
  pdf(file=plot_lambda_results_file, width=10, height=5)
  par(mar=c(4.5, 4.5, 1, 1) + 0.1)
  plot(hordeum_spd$grid$calBP, hordeum_spd$grid$PrDens, xlim = time_range_BP, ylim = c(0.0001, 1), log = "y",
       type="l", xlab="Years cal BP", ylab=expression(lambda), col="grey", lwd = 2)
  lines(triticum_spd$grid$calBP, triticum_spd$grid$PrDens, lty = 2, lwd = 2, col = "grey")
  lines(skyline_years, lambda_hat, col = PCI_blue, lwd = 2)
  lines(skyline_years, lambda_95low, lty = 2, lwd = 2, col = PCI_blue)
  lines(skyline_years, lambda_95upp, lty = 2, lwd = 2, col = PCI_blue)
  text(time_range_BP[1],1,"a",cex=2)
  dev.off()
  
}
  





# note that there was a change of notation from the initial code (pi)
# to the draft (q) regarding the proportion of Triticum samples
plot_q_results_file = paste0(results_cereals_directory, "/interdependent_model_result_q.pdf")
if ( !file.exists(plot_q_results_file) ){
  load(file = spd_file)
  load(file = parameter_inference_file)
  pdf(file=plot_q_results_file, width=10, height=5)
  par(mar=c(4.5, 4.5, 1, 1) + 0.1)
  plot(skyline_years, pi_hat, xlim = time_range_BP, ylim = c(0, 1), 
       type="l", xlab="Years cal BP", ylab=expression(italic(q)), col = PCI_blue, lwd = 2)
  lines(skyline_years, pi_95low, lty = 2, lwd = 2, col = PCI_blue)
  lines(skyline_years, pi_95upp, lty = 2, lwd = 2, col = PCI_blue)
  text(time_range_BP[1],1,"b",cex=2)
  dev.off()
  
  
}











