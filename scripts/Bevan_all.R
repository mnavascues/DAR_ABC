########################################################################
# Analysis of real data from Britain and Ireland
########################################################################
source("scripts/DARthABC.R")
require(abcrf)
require(doSNOW)
require(doParallel)
require(doRNG)
results_directory = "results/Bevan_all"
dir.create(results_directory)

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
  time_range_BP = c(10000, 500)
  t_ = seq(time_range_BP[1], time_range_BP[2] + 1,by = -1)

# Data processing 

# Load Bevan et al. 2017 data (obtained from DOI: 10.14324/000.ds.10025178)
data_file = paste0(results_directory, "/dates.rda")
if (!file.exists(data_file)){
  dates = read.csv("data/gbie14Csub/dates/archdates.csv",
                   header = TRUE,
                   stringsAsFactors = FALSE,
                   encoding = "UTF-8",
                   na.strings = c("NA", ""),
                   strip.white = TRUE)
  head(dates)
  save(dates, file = data_file)
}else{
  load(file=data_file)
}

# calculate bins and weights
weights_file = paste0(results_directory, "/weights.rda")
if (!file.exists(weights_file)){
  bins = binPrep(sites = dates$SiteID, ages = dates$CRA, h = 100)
  w = apply(as.array(bins), 1, function(x) 1/sum(bins == x))
  save(bins, w, file=weights_file)
}else{
  #load(file=weights_file)
}



# Make histogram plot
histogram_file = paste0(results_directory, "/hist.pdf")
if (!file.exists(histogram_file)){
  dates_on_range = intersect(which(dates$CRA<time_range_BP[1]), which(dates$CRA>time_range_BP[2]))
  window=100
  breaks = seq(time_range_BP[1],time_range_BP[2]-1,-window)
  pdf(file=histogram_file, width=10, height=5)
  par(mar=c(4.5, 4.5, 1, 1) + 0.1)
  wtd.hist(dates$CRA[dates_on_range], breaks=breaks, plot=T, weight = w[dates_on_range],
           main="",xlim=time_range_BP, xlab="Years BP (uncalibrated)",col=PCI_t_blue)
  box()
  dev.off()
}


# Calibrate dates
# (Difference with Bevan et al. 2017: using intcal20 instead of intcal13)
caldates_file = paste0(results_directory, "/caldates.rda")
if (!file.exists(caldates_file)){
  ncores = 30
  cl = makeCluster(ncores, type="SOCK")
  registerDoSNOW(cl)
  caldates = calibrate(x = dates$CRA,
                       errors = dates$Error,
                       ids = dates$DateID,
                       calCurves = 'intcal20',
                       timeRange = time_range_BP,
                       normalised = FALSE,
                       ncores = ncores,
                       calMatrix = TRUE)
  stopCluster(cl)
  save(caldates, file = caldates_file)
}else{
  #load(file=caldates_file)
}


# Calculate SPD
spd_file = paste0(results_directory, "/spd.rda")
if (!file.exists(spd_file)){
  allspd = spd(x = caldates, 
               bins = bins, 
               timeRange = time_range_BP, 
               datenormalised = FALSE,
               runm = 100)
  save(allspd, file = spd_file)
}else{
  #load(file=spd_file)
} 
# Make SPD plot
spd_plot_file = paste0(results_directory, "/spd.pdf")
if (!file.exists(spd_plot_file)){
  pdf(file=spd_plot_file, width=10, height=5)
  par(mar=c(4.5, 4.5, 1, 1) + 0.1)
  plot(allspd$grid$calBP, allspd$grid$PrDens, xlim = time_range_BP,
       type="l", xlab="Years cal BP", ylab="Sum of Probability Densities (SPD)", col=PCI_blue, lwd=2)
  dev.off()
}


# Calculate Summary statstics
sumstats_file = paste0(results_directory, "/sumstats.rda")
if (!file.exists(sumstats_file)){
  all_sumstats_c14 = get_sumstats(dates[c("CRA","Error")], time_range_BP, window = c(10,50,100,500), w=w)
  save(all_sumstats_c14, file=sumstats_file)
}

# set number of simulations
num_of_sims = 100000

# set prior parameters
lambda_min = 0.001
lambda_max = 12

set.seed(24)
ncores = 30

# make a reference table exponential model
reftable_exponential_file = paste0(results_directory, "/exponential_model_reftable.rda")
if (!file.exists(reftable_exponential_file)){    
  # setup parallel computing
  cl <- makeCluster(ncores, type="FORK")  
  registerDoParallel(cl)  
  registerDoRNG(seed = 73562846)
  reftable <- foreach(sim=seq_len(num_of_sims), .combine=rbind) %dopar% {
    gc()
    
    params = sample_exponential_parameters_from_priors(lambda_min, lambda_max, time_range_BP, force="expansion")
    
    lambda_t = get_exponential_lambda_t(params$lambda_0, params$rate, time_range_BP)
    
    sim_dates = sim_dates_lambda(lambda_t, time_range_BP)
    sim_dates_CRA = sim_CRA(sim_dates, errors = dates$Error)
    sumstats = get_sumstats(sim_dates_CRA, time_range_BP, window = c(10,50,100,500))
    
  
    cbind(params,sumstats)
  }
  stopCluster(cl) 
  save(reftable,file=reftable_exponential_file)
}



# make a reference table logistic model
reftable_logistic_file = paste0(results_directory, "/logistic_model_reftable.rda")
if (!file.exists(reftable_logistic_file)){    
  # setup parallel computing
  cl <- makeCluster(ncores, type="FORK")  
  registerDoParallel(cl)  
  registerDoRNG(seed = 73562846)
  reftable <- foreach(sim=seq_len(num_of_sims), .combine=rbind) %dopar% {
    gc()
    
    params = sample_logistic_parameters_from_priors(lambda_min, lambda_max, time_range_BP, force="expansion")
    
    lambda_t = get_logistic_lambda_t(params$lambda_0, params$K, params$rate, time_range_BP)
    
    sim_dates = sim_dates_lambda(lambda_t, time_range_BP)
    sim_dates_CRA = sim_CRA(sim_dates, errors = dates$Error)
    sumstats = get_sumstats(sim_dates_CRA, time_range_BP, window = c(10,50,100,500))
    
    
    cbind(params,sumstats)
  }
  stopCluster(cl) 
  save(reftable,file=reftable_logistic_file)
}




# make a reference table piecewise model
reftable_piecewise_file = paste0(results_directory, "/piecewise_model_reftable.rda")
num_of_periods = round((time_range_BP[1] - time_range_BP[2]) / 100 / 4)
if (!file.exists(reftable_piecewise_file)){    
  # setup parallel computing
  cl <- makeCluster(ncores, type="FORK")  
  registerDoParallel(cl)  
  registerDoRNG(seed = 73562846)
  reftable <- foreach(sim=seq_len(num_of_sims), .combine=rbind) %dopar% {
    gc()
    
    params = sample_exponential_piecewise_parameters_from_priors(lambda_min, lambda_max,
                                                                 time_range_BP,
                                                                 num_of_periods,
                                                                 intervals="regular")
      

    lambda_t = get_piecewise_exponential_lambda_t(as.numeric(params[seq_len(num_of_periods+1)]),
                                                  as.integer(params[1+num_of_periods+seq_len(num_of_periods+1)]),
                                                  as.numeric(params[2+2*num_of_periods+seq_len(num_of_periods)]))
    
    
    sim_dates = sim_dates_lambda(lambda_t, time_range_BP)
    sim_dates_CRA = sim_CRA(sim_dates, errors = dates$Error)
    sumstats = get_sumstats(sim_dates_CRA, time_range_BP, window = c(10,50,100,500))
    
    
    cbind(params,sumstats)
  }
  stopCluster(cl) 
  save(reftable,file=reftable_piecewise_file)
}


# MODEL CHOICE
model_choice_file = paste0(results_directory, "/model_choice.rda")
subset_size = 30000
if ( !file.exists(model_choice_file) ){
  
  load(file = sumstats_file)
  
  load(file = reftable_exponential_file)
  rows2keep = complete.cases(reftable)
  reftable_exponential = reftable[rows2keep,] #; head(reftable_exponential)
  
  load(file = reftable_logistic_file)
  rows2keep = complete.cases(reftable)
  reftable_logistic = reftable[rows2keep,] #; head(reftable_logistic)
  
  load(file = reftable_piecewise_file)
  rows2keep = complete.cases(reftable)
  reftable_piecewise = reftable[rows2keep,] #; head(reftable_piecewise)
  
  rm(reftable);gc()  
  
  
  model = as.factor(c(rep("E",subset_size),  
                      rep("L",subset_size),  
                      rep("P",subset_size))) 

  sumstats = rbind(reftable_exponential[seq_len(subset_size) , names(all_sumstats_c14)],
                   reftable_logistic[seq_len(subset_size) , names(all_sumstats_c14)],
                   reftable_piecewise[seq_len(subset_size) , names(all_sumstats_c14)])
  
  RF_model_choice = abcrf(model~., data = data.frame(model,sumstats),
                          lda=F, ntree = 2000, paral = TRUE)
  
  plot(RF_model_choice, data.frame(model,sumstats))
  
  posterior_model_choice =  predict(RF_model_choice, all_sumstats_c14,
                                    training = data.frame(model,sumstats),
                                    ntree = 2000, paral = TRUE)
  confusion_matrix = RF_model_choice$model.rf$confusion.matrix
  prediction_error = RF_model_choice$model.rf$prediction.error
  save(prediction_error,
       confusion_matrix,
       posterior_model_choice, file=model_choice_file)

  
  posterior_model_choice
  confusion_matrix
  prediction_error
  K = posterior_model_choice$post.prob * (1-subset_size/length(model)) / (1-posterior_model_choice$post.prob) / (subset_size / length(model))
  (K)
  interpret_K(K)

}else{
  load(file=model_choice_file)
}

# PCA for visual goodness of fit

PCA_file = paste0(results_directory, "/PCA.rda")
subset_size = 30000
model = as.factor(c(rep("E",subset_size),  # nrow(reftable_exponential)),
                    rep("L",subset_size),  # nrow(reftable_logistic)),
                    rep("P",subset_size))) # nrow(reftable_piecewise))))
if ( !file.exists(PCA_file) ){
  load(file = reftable_exponential_file)
  rows2keep = complete.cases(reftable)
  reftable_exponential = reftable[rows2keep,] #; head(reftable_exponential)
  
  load(file = reftable_logistic_file)
  rows2keep = complete.cases(reftable)
  reftable_logistic = reftable[rows2keep,] #; head(reftable_logistic)
  
  load(file = reftable_piecewise_file)
  rows2keep = complete.cases(reftable)
  reftable_piecewise = reftable[rows2keep,] #; head(reftable_piecewise)
  
  load(file = sumstats_file)
  
  
  sumstats = rbind(reftable_exponential[seq_len(subset_size) , names(all_sumstats_c14)],
                   reftable_logistic[seq_len(subset_size) , names(all_sumstats_c14)],
                   reftable_piecewise[seq_len(subset_size) , names(all_sumstats_c14)])
  PCA_stats  = princomp(sumstats)
  PCA_target = predict(PCA_stats, all_sumstats_c14)
  
  summary(PCA_stats)
  save(PCA_stats, PCA_target, file=PCA_file)
}else{
  #load(file=PCA_file)
  #sum((PCA_stats$sdev[1:6])^2)/sum((PCA_stats$sdev)^2)
}  

# plot PCA
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
    xlim = c(-10000, 60000); ylim = c(-10000, 10000); id_plot = "a"; add_legend=TRUE
  }
  if (axis_pair == 2){
    i = 3; j = 4
    xlim = c(-2000,4000); ylim=c(-2000,4000); id_plot="b"; add_legend=FALSE
  }
  if (axis_pair == 3){
    i = 5; j = 6
    xlim=c(-3000,2000); ylim=c(-2000,3000); id_plot="c"; add_legend=FALSE
  }
  
  pdf_file_name = paste0(results_directory, "/PCA_PC",i,"_PC",j,".pdf")
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
                           c("Exponential","Logistic","Piecewise exponential"),
                           pch=1,
                           col=cbPalette1[c(2,3,4)],
                           bg = "white")
    dev.off()
  }
}

exponential_model_inference_file = paste0(results_directory, "/paramter_estimates_exponential.rda")
if ( !file.exists(exponential_model_inference_file) ){
  load(file = reftable_exponential_file)
  reftable = reftable[complete.cases(reftable),] #; head(reftable_exponential)
  load(file = sumstats_file)
  
  sumstats = reftable[names(all_sumstats_c14)]
  
  param = log10(reftable$lambda_0)
  RF_log10lambda_0 = regAbcrf(param~., data.frame(param,sumstats),
                              ntree = 2000, paral = TRUE)
  posterior_log10lambda_0 = predict(RF_log10lambda_0, all_sumstats_c14,
                                    training = data.frame(param,sumstats),
                                    paral = TRUE, rf.weights = TRUE) 
  param = log10(reftable$lambda_f)
  RF_log10lambda_f = regAbcrf(param~., data.frame(param,sumstats),
                              ntree = 2000, paral = TRUE)
  posterior_log10lambda_f = predict(RF_log10lambda_f, all_sumstats_c14,
                                    training = data.frame(param,sumstats),
                                    paral = TRUE, rf.weights = TRUE) 
  save(RF_log10lambda_0, posterior_log10lambda_0,
       RF_log10lambda_f, posterior_log10lambda_f,
       file=exponential_model_inference_file)
}else{
  
  plot_posterior_lambda_0_file = paste0(results_directory, "/exponential_model_lambda_0_posterior.pdf")
  if ( !file.exists(plot_posterior_lambda_0_file) ){
    load(file = reftable_exponential_file)
    reftable = reftable[complete.cases(reftable),] #; head(reftable_exponential)

    load(file=exponential_model_inference_file)
    lambda_0_hat = 10^posterior_log10lambda_0$med[1]
    lambda_0_CI = 10^posterior_log10lambda_0$quantiles
    print(lambda_0_hat)
    print(lambda_0_CI)
    
    pdf(file = plot_posterior_lambda_0_file, width = 4, height = 4)
    par(mar = c(4.5, 4.5, 1, 1) + 0.1)
    breaks = seq(log10(lambda_min),log10(lambda_max)+0.05,0.05)
    hist(log10(reftable$lambda_0),
         breaks = breaks,
         main = "",#expression("Posterior probabilty number of periods ("*italic(m)*")"),
         xlab = expression(log[10](lambda[0])),
         ylim = c(0,3.5),
         col = adjustcolor("gray", alpha.f = 0.6), freq = F)
    text(x=log10(lambda_min), y=3.5, label="a", cex=2)
    wtd.hist(log10(reftable$lambda_0),
             breaks = breaks,
             col=PCI_t_blue,
             weight = posterior_log10lambda_0$weights,
             add=T, freq=F)
    box()
    dev.off()


    lambda_f_hat = 10^posterior_log10lambda_f$med[1]
    lambda_f_CI = 10^posterior_log10lambda_f$quantiles
    print(lambda_f_hat)
    print(lambda_f_CI)
    
    
    plot_posterior_lambda_f_file = paste0(results_directory, "/exponential_model_lambda_f_posterior.pdf")
    pdf(file = plot_posterior_lambda_f_file, width = 4, height = 4)
    par(mar = c(4.5, 4.5, 1, 1) + 0.1)
    breaks = seq(log10(lambda_min),log10(lambda_max)+0.05,0.05)
    hist(log10(reftable$lambda_f),
         breaks = breaks,
         main = "",#expression("Posterior probabilty number of periods ("*italic(m)*")"),
         xlab = expression(log[10](lambda[f])),
         ylim = c(0,12),
         col = adjustcolor("gray", alpha.f = 0.6), freq = F)
    text(x=log10(lambda_min), y=12, label="b", cex=2)
    wtd.hist(log10(reftable$lambda_f),
             breaks = breaks,
             col=PCI_t_blue,
             weight = posterior_log10lambda_f$weights,
             add=T, freq=F)
    box()
    dev.off()
    
        
  }  

}  






logistic_model_inference_file = paste0(results_directory, "/paramter_estimates_logistic.rda")
if ( !file.exists(logistic_model_inference_file) ){
  load(file = reftable_logistic_file)
  reftable = reftable[complete.cases(reftable),] #; head(reftable_exponential)
  load(file = sumstats_file)
  
  sumstats = reftable[names(all_sumstats_c14)]
  
  param = log10(reftable$lambda_0)
  RF_log10lambda_0 = regAbcrf(param~., data.frame(param,sumstats),
                              ntree = 2000, paral = TRUE)
  posterior_log10lambda_0 = predict(RF_log10lambda_0, all_sumstats_c14,
                                    training = data.frame(param,sumstats),
                                    paral = TRUE, rf.weights = TRUE) 
  param = log10(reftable$lambda_f)
  RF_log10lambda_f = regAbcrf(param~., data.frame(param,sumstats),
                              ntree = 2000, paral = TRUE)
  posterior_log10lambda_f = predict(RF_log10lambda_f, all_sumstats_c14,
                                    training = data.frame(param,sumstats),
                                    paral = TRUE, rf.weights = TRUE) 
  param = log10(reftable$K)
  RF_log10K = regAbcrf(param~., data.frame(param,sumstats),
                              ntree = 2000, paral = TRUE)
  posterior_log10K = predict(RF_log10K, all_sumstats_c14,
                                    training = data.frame(param,sumstats),
                                    paral = TRUE, rf.weights = TRUE) 
  save(RF_log10lambda_0, posterior_log10lambda_0,
       RF_log10lambda_f, posterior_log10lambda_f,
       RF_log10K, posterior_log10K,
       file=logistic_model_inference_file)
}else{
  
  plot_posterior_lambda_0_file = paste0(results_directory, "/logistic_model_lambda_0_posterior.pdf")
  if ( !file.exists(plot_posterior_lambda_0_file) ){
    load(file = reftable_logistic_file)
    reftable = reftable[complete.cases(reftable),] #; head(reftable_exponential)
    
    load(file=logistic_model_inference_file)
    lambda_0_hat = 10^posterior_log10lambda_0$med[1]
    lambda_0_CI = 10^posterior_log10lambda_0$quantiles
    print(lambda_0_hat)
    print(lambda_0_CI)
    
    pdf(file = plot_posterior_lambda_0_file, width = 4, height = 4)
    par(mar = c(4.5, 4.5, 1, 1) + 0.1)
    breaks = seq(log10(lambda_min),log10(lambda_max)+0.05,0.05)
    hist(log10(reftable$lambda_0),
         breaks = breaks,
         main = "",#expression("Posterior probabilty number of periods ("*italic(m)*")"),
         xlab = expression(log[10](lambda[0])),
         ylim = c(0,3.5),
         col = adjustcolor("gray", alpha.f = 0.6), freq = F)
    text(x=log10(lambda_min), y=3.5, label="a", cex=2)
    wtd.hist(log10(reftable$lambda_0),
             breaks = breaks,
             col=PCI_t_blue,
             weight = posterior_log10lambda_0$weights,
             add=T, freq=F)
    box()
    dev.off()
    
    
    lambda_f_hat = 10^posterior_log10lambda_f$med[1]
    lambda_f_CI = 10^posterior_log10lambda_f$quantiles
    print(lambda_f_hat)
    print(lambda_f_CI)
    
    
    plot_posterior_lambda_f_file = paste0(results_directory, "/logistic_model_lambda_f_posterior.pdf")
    pdf(file = plot_posterior_lambda_f_file, width = 4, height = 4)
    par(mar = c(4.5, 4.5, 1, 1) + 0.1)
    breaks = seq(log10(lambda_min),log10(lambda_max)+0.05,0.05)
    hist(log10(reftable$lambda_f),
         breaks = breaks,
         main = "",#expression("Posterior probabilty number of periods ("*italic(m)*")"),
         xlab = expression(log[10](lambda[f])),
         ylim = c(0,8),
         col = adjustcolor("gray", alpha.f = 0.6), freq = F)
    text(x=log10(lambda_min), y=8, label="b", cex=2)
    wtd.hist(log10(reftable$lambda_f),
             breaks = breaks,
             col=PCI_t_blue,
             weight = posterior_log10lambda_f$weights,
             add=T, freq=F)
    box()
    dev.off()
    
    plot_posterior_K_file = paste0(results_directory, "/logistic_model_K_posterior.pdf")
    pdf(file = plot_posterior_K_file, width = 4, height = 4)
    par(mar = c(4.5, 4.5, 1, 1) + 0.1)
    breaks = seq(log10(lambda_min),max(log10(reftable$K))+0.05,0.05)
    hist(log10(reftable$K),
         breaks = breaks,
         main = "",#expression("Posterior probabilty number of periods ("*italic(m)*")"),
         xlab = expression(log[10](italic(k))),
         ylim = c(0,5),
         col = adjustcolor("gray", alpha.f = 0.6), freq = F)
    text(x=log10(lambda_min), y=5, label="c", cex=2)
    wtd.hist(log10(reftable$K),
             breaks = breaks,
             col=PCI_t_blue,
             weight = posterior_log10K$weights,
             add=T, freq=F)
    box()
    dev.off()
    
  }  
  
}  


num_of_periods = round((time_range_BP[1] - time_range_BP[2]) / 100 / 4)
skyline_years = get_time_of_change(num_of_periods, time_range_BP, intervals="regular")


piecewise_model_inference_file = paste0(results_directory, "/parameter_estimates_piecewise.rda")
if ( !file.exists(piecewise_model_inference_file) ){
  
  load(file = reftable_piecewise_file)
  reftable = reftable[complete.cases(reftable),] #; head(reftable_exponential)
  load(file = sumstats_file)
  
  sumstats = reftable[names(all_sumstats_c14)]
  
  
  lambda_hat      = rep(NA,length(skyline_years))
  lambda_95CI_low = rep(NA,length(skyline_years))
  lambda_95CI_upp = rep(NA,length(skyline_years))
  lambda_error    = rep(NA,length(skyline_years))
  
  rate_error    = rep(NA,length(skyline_years)-1)
  rate_hat      = rep(NA,length(skyline_years)-1)
  rate_95CI_low = rep(NA,length(skyline_years)-1)
  rate_95CI_upp = rep(NA,length(skyline_years)-1)
  
  for (i in seq_len(length(skyline_years))){
    results_file = paste0(results_directory, "/piecewise_posterior_lambda_",i, ".rda")
    if ( !file.exists(results_file) ){
      param_name = paste0("lambda_",i)
      param_index = which(names(reftable)==param_name)
      param = log10(reftable[param_index])
      names(param) = "param"
      RF_log10lambda = regAbcrf(param~., data.frame(param,sumstats),
                                ntree = 2000, paral = TRUE)
      posterior_log10lambda = predict(RF_log10lambda, all_sumstats_c14,
                                      training = data.frame(param,sumstats),
                                      paral = TRUE, rf.weights = TRUE) 
      
      save(RF_log10lambda, posterior_log10lambda, file = results_file)
    }else{load(file = results_file)}
    lambda_error[i] = RF_log10lambda$model.rf$prediction.error
    rm(RF_log10lambda);gc()
    lambda_hat[i] = 10^(posterior_log10lambda$med[1])
    lambda_95CI_low[i] = 10^(posterior_log10lambda$quantiles[1])
    lambda_95CI_upp[i] = 10^(posterior_log10lambda$quantiles[2])
  }
  
  for (i in seq_len(length(skyline_years)-1)){
    results_file = paste0(results_directory,"/piecewise_posterior_r_",i, ".rda")
    if ( !file.exists(results_file) ){
      param_name = paste0("r_",i)
      param_index = which(names(reftable)==param_name)
      param = reftable[param_index]
      names(param) = "param"
      RF_rate = regAbcrf(param~., data.frame(param,sumstats),
                         ntree = 2000, paral = TRUE)
      posterior_rate = predict(RF_rate, all_sumstats_c14,
                               training = data.frame(param,sumstats),
                               paral = TRUE, rf.weights = TRUE) 
      
      save(RF_rate, posterior_rate, file = results_file)
    }else{load(file = results_file)}
    rate_error[i] = RF_rate$model.rf$prediction.error
    rm(RF_rate);gc()
    rate_hat[i] = (posterior_rate$med[1])
    rate_95CI_low[i] = (posterior_rate$quantiles[1])
    rate_95CI_upp[i] = (posterior_rate$quantiles[2])
  }
  
  
  
  save(lambda_hat,
       lambda_95CI_low, lambda_95CI_upp,
       lambda_error,
       rate_hat,
       rate_95CI_low, rate_95CI_upp,
       rate_error,
       file=piecewise_model_inference_file)
}else{
  #load(file=piecewise_model_inference_file)
}

piecewise_model_result_lambda_plot_file = paste0(results_directory, "/piecewise_model_result_lambda.pdf")
if (!file.exists(piecewise_model_result_lambda_plot_file)){
  load(file=spd_file)
  pdf(file=piecewise_model_result_lambda_plot_file, width=10, height=5)
  par(mar=c(4.5, 4.5, 1, 1) + 0.1)
  plot(allspd$grid$calBP, allspd$grid$PrDens, xlim = time_range_BP, ylim = c(0.04, 6), log = "y",
       type="l", xlab="Years cal BP", ylab=expression(lambda), col="grey", lwd = 2)
  lines(skyline_years, lambda_hat, col = PCI_blue, lwd = 2)
  lines(skyline_years, lambda_95CI_low, lty = 2, lwd = 2, col = PCI_blue)
  lines(skyline_years, lambda_95CI_upp, lty = 2, lwd = 2, col = PCI_blue)
  text(time_range_BP[1],6,"a",cex=2)
  dev.off()
}

piecewise_model_result_rate_plot_file = paste0(results_directory, "/piecewise_model_result_rate.pdf")
if (!file.exists(piecewise_model_result_rate_plot_file)){
  load(file=spd_file)
  step_wise_years = c(skyline_years[1],rep(skyline_years[2:(length(skyline_years)-1)],each=2),skyline_years[length(skyline_years)])
  pdf(file=piecewise_model_result_rate_plot_file, width=10, height=5)
  par(mar=c(4.5, 4.5, 1, 1) + 0.1)
  plot(step_wise_years, rep(rate_hat,each=2),
       xlab = "Years cal BP", ylab=expression(italic(r)),
       lwd = 2, type = "l", xlim = time_range_BP,
       ylim = c(-0.006, 0.006),
       #ylim=c(min(rate_95CI_low),max(rate_95CI_upp)),
       col=PCI_blue)
  lines(step_wise_years,rep(rate_95CI_low,each=2),lty=3, lwd=2,col=PCI_blue)
  lines(step_wise_years,rep(rate_95CI_upp,each=2),lty=3, lwd=2,col=PCI_blue)
  abline(h=0,col="gray")
  points( skyline_years[which(rate_95CI_low>0)]-196.5, rep(0,sum(rate_95CI_low>0)),
          pch="*", cex=2)
  points( skyline_years[which(rate_95CI_upp<0)]-196.5, rep(0,sum(rate_95CI_upp<0)),
          pch="*", cex=2)
  text(time_range_BP[1],0.006,"b",cex=2)
  dev.off()
}



pdf_file_name = paste0(results_directory, "/posterior_lambda_1_piecewise.pdf")
if (!file.exists(pdf_file_name)){
  load(file = reftable_piecewise_file)
  reftable = reftable[complete.cases(reftable),] #; head(reftable_exponential)
  load(file = sumstats_file)
  
  sumstats = reftable[names(all_sumstats_c14)]
  for (i in seq_len(length(skyline_years))){
    results_file = paste0(results_directory, "/piecewise_posterior_lambda_", i, ".rda")
    load(results_file)
    param_name = paste0("lambda_",i)
    param_index = which(names(reftable)==param_name)
    pdf_file_name = paste0(results_directory, "/posterior_lambda_", i, "_piecewise.pdf")
    pdf(file=pdf_file_name, width=4, height=4)
    par(mar=c(4.5, 4.5, 1, 1) + 0.1)
    breaks= seq(-3,1.1,0.05)
    
    hist(log10(t(reftable[param_index])),
         breaks = breaks,
         main = "",
         xlab = bquote(log[10]*(lambda[.(skyline_years[i])])),
         ylim = c(0,7),
         col = adjustcolor( "gray", alpha.f = 0.6), freq = F)
    wtd.hist(log10(t(reftable[param_index])),
             breaks = breaks,
             col=PCI_t_blue,
             weight = posterior_log10lambda$weights,
             add=T, freq=F)
    box()
    text(-3,7, letters[i],cex=2)
    dev.off()
  }
}









pdf_file_name = paste0(results_directory, "/posterior_rate_1_piecewise.pdf")
if (!file.exists(pdf_file_name)){
  load(file = reftable_piecewise_file)
  reftable = reftable[complete.cases(reftable),] #; head(reftable_exponential)
  load(file = sumstats_file)
  
  sumstats = reftable[names(all_sumstats_c14)]
  for (i in seq_len(length(skyline_years)-1)){
    results_file = paste0(results_directory, "/piecewise_posterior_r_", i, ".rda")
    load(results_file)
    param_name = paste0("r_",i)
    param_index = which(names(reftable)==param_name)
    pdf_file_name = paste0(results_directory, "/posterior_rate_", i, "_piecewise.pdf")
    pdf(file=pdf_file_name, width=4, height=4)
    par(mar=c(4.5, 4.5, 1, 1) + 0.1)
    breaks= seq(-0.006,0.006,0.0002)
    
    hist(t(reftable[param_index]),
         breaks = breaks,
         main = "",
         xlab = bquote(italic(r)["("*.(skyline_years[i])-.(skyline_years[i+1])*")"]),
         ylim = c(0,800),
         col = adjustcolor( "gray", alpha.f = 0.6), freq = F)
    wtd.hist(t(reftable[param_index]),
             breaks = breaks,
             col=PCI_t_blue,
             weight = posterior_rate$weights,
             add=T, freq=F)
    box()
    text(-0.006,800, letters[i],cex=2)
    dev.off()
  }
}

# Plot comparing models with SPD


SPD_plus_models_file_name = paste0(results_directory, "/SPD_plus_lambda_models.pdf")
if (!file.exists(SPD_plus_models_file_name)){
  load(file=spd_file)
  pdf(file=SPD_plus_models_file_name, width=10, height=5)
  
  par(mar=c(4.5, 4.5, 1, 1) + 0.1)
  plot(allspd$grid$calBP, allspd$grid$PrDens, xlim = time_range_BP, ylim = c(0.09,6), log = "y",
       type="l", xlab="Years cal BP", ylab=expression(lambda), col="grey", lwd=2)
  
  time_seq = seq(time_range_BP[1], time_range_BP[2] + 1,by = -1)

  load(file=paste0(results_directory,"/paramter_estimates_exponential.rda"))
  lambda_0_hat = 10^posterior_log10lambda_0$med[1]
  lambda_f_hat = 10^posterior_log10lambda_f$med[1]
  r_hat = (log(lambda_0_hat)-log(lambda_f_hat))/length(time_seq)
  lambda_t = lambda_0_hat*exp(-r_hat*seq_along(time_seq))
  lines(time_seq, lambda_t, col=cbPalette1[2], lwd=2, lty=2)
  
    
  load(file=paste0(results_directory,"/paramter_estimates_logistic.rda"))
  K_hat = 10^posterior_log10K$med[1]
  lambda_0_hat = 10^posterior_log10lambda_0$med[1]
  lambda_f_hat = 10^posterior_log10lambda_f$med[1]
  r_hat = -log( (K_hat-lambda_f_hat)*lambda_0_hat/(K_hat-lambda_0_hat)/lambda_f_hat )/length(time_seq)
  lambda_t = K_hat*lambda_0_hat / (lambda_0_hat + (K_hat - lambda_0_hat)*exp(-r_hat*seq_along(time_seq)) )
  lines(time_seq, lambda_t, col=cbPalette1[3], lwd=2, lty=3)
  
  
  load(file=paste0(results_directory,"/parameter_estimates_piecewise.rda"))
  num_of_periods = round((time_range_BP[1] - time_range_BP[2]) / 100 / 4)
  skyline_years = get_time_of_change(num_of_periods, time_range_BP, intervals="regular")
  lines(skyline_years, lambda_hat, col = cbPalette1[4], lwd = 2, lty=4)
  
  legend("bottomright",
         legend=c("SPD","Exponential","Logistic","Piecewise exponential"),
         col=c("grey",  cbPalette1[2:4]),
         lty=1:4,
         lwd=2)
  dev.off()
  
}
  
















