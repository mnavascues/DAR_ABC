source("scripts/DARthABC.R")
results_directory = "results/Coldsnap"
dir.create(results_directory)
library(abcrf)
library(abc)
require(graphics)
cbPalette1 <- c("#000000", # Black
                "#E69F00", # Orange
                "#56B4E9", # Sky Blue
                "#009E73", # Bluish Green
                "#F0E442", # Yellow
                "#0072B2", # Blue
                "#D55E00", # Vermillion
                "#CC79A7") # Reddish Purple
cbPalette1_transparent <- c(rgb(   0,   0,   0, 50, maxColorValue = 255), # Black
                            rgb( 230, 159,   0, 50, maxColorValue = 255), # Orange
                            rgb(  86, 158, 233, 50, maxColorValue = 255), # Sky Blue
                            rgb(   0, 158, 115, 50, maxColorValue = 255), # Bluish Green
                            rgb( 240, 228,  66, 50, maxColorValue = 255), # Yellow
                            rgb(   0, 114, 178, 50, maxColorValue = 255), # Blue
                            rgb( 213,  94,   0, 50, maxColorValue = 255), # Vermillion
                            rgb( 204, 121, 167, 50, maxColorValue = 255)) # Reddish Purple
                            
            


dates = read.csv("data/Dateringar_Sverige_Apel_12_3.csv")
# dates$Lab..no.[which(dates$Lab..no.=="missing")] = NA
head(dates)

summary(dates$C14.years.BP)
time_range_CRA = c(max(dates$C14.years.BP),min(dates$C14.years.BP))
time_range_calBP = c(10000, 5000)


caldates_filename = paste0(results_directory, "/caldates.rda")
if (!file.exists(caldates_filename)){
  ncores = 30
  cl = makeCluster(ncores, type="SOCK")
  registerDoSNOW(cl)
  caldates = rcarbon::calibrate(x = dates$C14.years.BP,
                                errors = dates$X.,
                                #ids = dates$Lab..no.,
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

summary(as.factor(dates$Province))

Eastern_Central_Sweden_samples = c(which(dates$Province=="Uppland"),
                                   which(dates$Province=="Västmanland"),
                                   which(dates$Province=="Södermanland"),
                                   which(dates$Province=="Närke"))

South_Eastern_Sweden_samples = c(which(dates$Province=="Blekinge"), #
                                 which(dates$Province=="Småland"),
                                 which(dates$Province=="Östergötland"))

Eastern_Islands_samples = c(which(dates$Province=="Öland"),
                            which(dates$Province=="Gotland"))

Western_Sweden_samples = c(which(dates$Province=="Halland"),
                           which(dates$Province=="Västergötland"),
                           which(dates$Province=="Bohuslän"), #
                           which(dates$Province=="Dalsland"))


Middle_Sweden_samples = c(which(dates$Province=="Värmland"),
                          which(dates$Province=="Dalarna"), #
                          which(dates$Province=="Gävleborg"), #
                          which(dates$Province=="Gästrikland")) #


Northern_Sweden_samples = c(which(dates$Province=="Hälsingland"),
                            which(dates$Province=="Härjedalen"),
                            which(dates$Province=="Medelpad"),
                            which(dates$Province=="Jämtland"),
                            which(dates$Province=="Ångermanland"), #
                            which(dates$Province=="Västerbotten"),
                            which(dates$Province=="Norrbotten"),
                            which(dates$Province=="Lappland"))

Southern_Sweden_samples = which(dates$Province=="Skåne")

set_of_regions = c("Eastern Central Sweden",
                   "South Eastern Sweden",
                   "Eastern Islands",
                   "Western Sweden",
                   "Middle Sweden",
                   "Northern Sweden",
                   "Southern Sweden")

# calculate bins and weights
weights_file = paste0(results_directory, "/weights.rda")
if (!file.exists(weights_file)){
  bins = binPrep(sites = dates$Site, ages = dates$C14.years.BP, h = 200)
  w = apply(as.array(bins), 1, function(x) 1/sum(bins == x))
  save(bins, w, file=weights_file)
}else{
  load(file=weights_file)
}




for (i in seq_along(set_of_regions)){
  spd_file = paste0(results_directory, "/", gsub(" ", "_",  set_of_regions[i]), "_spd.rda")
  if (!file.exists(spd_file)){
    samples = get(paste0(gsub(" ", "_",  set_of_regions[i]),"_samples"))
    spd = rcarbon::spd(x = caldates[samples],
                       bins = bins[samples], 
                       timeRange = time_range_calBP, 
                       datenormalised = FALSE,
                       runm = 100)
    assign(paste0(gsub(" ", "_",  set_of_regions[i]),"_spd"), spd)
    save(list=paste0(gsub(" ", "_",  set_of_regions[i]),"_spd"), file = spd_file)
  }
}

spd_plot_file = paste0(results_directory, "/spd.pdf")
if (!file.exists(spd_plot_file)){
  pdf(file=spd_plot_file)
  for (i in seq_along(set_of_regions)){
    spd = get(paste0(gsub(" ", "_",  set_of_regions[i]),"_spd"))
    if (i==1){
      plot(spd$grid$calBP,
           spd$grid$PrDens,
           xlim = time_range_calBP,
           ylim = c(0.000005,1),
           log="y",
           type ="l",
           xlab = "Years cal BP",
           ylab = "SPD",
           lwd = 2,
           col = cbPalette1[i])
    }else{
      lines(spd$grid$calBP,
            spd$grid$PrDens,
            col = cbPalette1[i], lwd = 2)   
    }
  }
  abline(v=8200, lty=2, col="grey")
  legend("bottomleft", legend=set_of_regions, lwd=2, col=cbPalette1[1:7])
  dev.off()
}  


# generate reference table for piecewise exponential model.
lambda_min = 0.00001
lambda_max = 2
num_of_periods = 50
num_of_sims = 100000
piecewise_exponential_reftable_file = paste0(results_directory, "/piecewise_exponential_reftable.rda")
if (!file.exists(piecewise_exponential_reftable_file) ){
  ncores = 15
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
    sim_dates_CRA = sim_CRA(sim_dates, errors = dates$X.)
    sumstats = get_sumstats(sim_dates_CRA, time_range_CRA,  window = c(25, 50, 100))

    cbind(params, sumstats)
  }
  stopCluster(cl) 
  save(reftable, file = piecewise_exponential_reftable_file)
}else{
  load(file = piecewise_exponential_reftable_file)
}


skyline_years = c(seq(time_range_calBP[1], 8700, -100),
                  seq(8680, 8020, -20),
                  seq(8000, time_range_calBP[2], -100))
skyline_points = length(skyline_years)

reftable_skyline_file = paste0(results_directory, "/skyline_reftable.rda")
if (!file.exists(reftable_skyline_file) ){
  #reftable_skyline_plot = data.frame(matrix(NA, nrow=nrow(reftable), ncol=skyline_points*2))
  #names(reftable_skyline_plot) =  c(paste0("lambda",skyline_years),
  #                                  paste0("rate",skyline_years))
  reftable_skyline_plot = data.frame(matrix(NA, nrow=nrow(reftable), ncol=skyline_points))
  names(reftable_skyline_plot) =  paste0("lambda",skyline_years)
  
  for (sim in seq_len(num_of_sims)){
    if (sim %% 100 == 0) print(paste("sim",sim))
    lambdas         = as.numeric(reftable[sim,seq_len(num_of_periods+1)])
    times_of_change = as.integer(reftable[sim,1+num_of_periods+seq_len(num_of_periods+1)])
    rates           = as.numeric(reftable[sim,2+2*num_of_periods+seq_len(num_of_periods)])
    
    lambda_skyline = get_piecewise_exponential_lambda_t(lambdas,
                                                        times_of_change,
                                                        rates)[match(skyline_years,seq(time_range_calBP[1],time_range_calBP[2],-1))]
    
    
    #rates_skyline = c(rep(rates,abs(diff(times_of_change))), rates[length(rates)])[seq(1,skyline_points*skyline_step,skyline_step)]
    
    reftable_skyline_plot[sim,] =  lambda_skyline
  }
  
  save(reftable_skyline_plot, file = reftable_skyline_file)
}else{
  #load(file = reftable_skyline_file)
}


# calculate observed summary statistics per region
for (i in seq_along(set_of_regions)){
  sumstats_file = paste0(results_directory, "/", gsub(" ", "_",  set_of_regions[i]), "_sumstats.rda")
  if (!file.exists(sumstats_file)){
    samples = get(paste0(gsub(" ", "_",  set_of_regions[i]),"_samples"))
    dates_CRA = list(CRA = dates$C14.years.BP[samples])
    target_sumstats = get_sumstats(dates_CRA, time_range_CRA, window = c(25, 50, 100), w=w[samples])
    assign(paste0(gsub(" ", "_",  set_of_regions[i]),"_sumstats"),
           target_sumstats)
    save(list = paste0(gsub(" ", "_",  set_of_regions[i]), "_sumstats"), file = sumstats_file)
  }
}


load(file=sumstats_file)
target_sumstats = get(paste0(gsub(" ", "_",  set_of_regions[i]),"_sumstats"))  
sumstats = reftable[, names(target_sumstats)]

i=1; skyline_file = paste0(results_directory, "/", gsub(" ", "_",  set_of_regions[i]), "_skyline.rda")
if (!file.exists(skyline_file)){
  for (j in seq_len(skyline_points)){
    print(paste0("Skyline point: ",skyline_years[j]))
    RF_file = paste0(results_directory, "/RF_skyline_",  skyline_years[j], ".rda")
    if (!file.exists(RF_file)){
      param = log10(reftable_skyline_plot[,j])
      RF_lambda = regAbcrf(param~., data.frame(param,sumstats),
                           ntree = 2000, paral = TRUE)
      save(RF_lambda, file=RF_file)
    }else{
      param = log10(reftable_skyline_plot[,j])
      load(file=RF_file)
    }
    for (i in seq_along(set_of_regions)){
      skyline_file = paste0(results_directory, "/", gsub(" ", "_",  set_of_regions[i]), "_skyline.rda")
      if (!file.exists(skyline_file)){
        print(paste0("     region: ",set_of_regions[i]))
        sumstats_file = paste0(results_directory, "/", gsub(" ", "_",  set_of_regions[i]), "_sumstats.rda")
        load(file = sumstats_file)
        target_sumstats = get(paste0(gsub(" ", "_",  set_of_regions[i]),"_sumstats"))
        posterior_lambda = predict(RF_lambda, target_sumstats,
                                   training = data.frame(param,sumstats),
                                   paral = TRUE, rf.weights = F)
        assign(paste0(gsub(" ", "_",  set_of_regions[i]),"_skyline_median_",j),
               posterior_lambda$med)
        assign(paste0(gsub(" ", "_",  set_of_regions[i]),"_skyline_q_0.025_",j),
               posterior_lambda$quantiles[1])
        assign(paste0(gsub(" ", "_",  set_of_regions[i]),"_skyline_q_0.975_",j),
               posterior_lambda$quantiles[2])
        assign(paste0(gsub(" ", "_",  set_of_regions[i]),"_skyline_error_",j),
               RF_lambda$model.rf$prediction.error)
        rm(posterior_lambda)
        gc()
      }  
    }
    rm(param, RF_lambda)
    gc()
  }
}

for (i in seq_along(set_of_regions)){
  skyline_file = paste0(results_directory, "/", gsub(" ", "_",  set_of_regions[i]), "_skyline.rda")
  if (!file.exists(skyline_file)){
    skyline = list(median  = rep(NA, skyline_points),
                   q_0.025 = rep(NA, skyline_points),
                   q_0.975 = rep(NA, skyline_points),
                   error   = rep(NA, skyline_points))
    for (j in seq_len(skyline_points)){
      skyline$median[j] = get(paste0(gsub(" ", "_",  set_of_regions[i]),"_skyline_median_",j))
      skyline$q_0.025[j] = get(paste0(gsub(" ", "_",  set_of_regions[i]),"_skyline_q_0.025_",j))
      skyline$q_0.975[j] = get(paste0(gsub(" ", "_",  set_of_regions[i]),"_skyline_q_0.975_",j))
      skyline$error[j] = get(paste0(gsub(" ", "_",  set_of_regions[i]),"_skyline_error_",j))
      
      rm(list = c(paste0(gsub(" ", "_",  set_of_regions[i]),"_skyline_median_",j),
                  paste0(gsub(" ", "_",  set_of_regions[i]),"_skyline_q_0.025_",j),
                  paste0(gsub(" ", "_",  set_of_regions[i]),"_skyline_q_0.975_",j),
                  paste0(gsub(" ", "_",  set_of_regions[i]),"_skyline_error_",j)))
    }
    
    assign(paste0(gsub(" ", "_",  set_of_regions[i]),"_skyline"),
           skyline)
    rm(skyline)
    save(list=paste0(gsub(" ", "_",  set_of_regions[i]),"_skyline"), file = skyline_file)
  }
}  


#skyline_years = seq(time_range_calBP[1], time_range_calBP[2], -100)  

skyline_plot_file = paste0(results_directory, "/skyline.pdf")
if (!file.exists(skyline_plot_file)){
  pdf(file = skyline_plot_file, width = 8.27, height = 11.69)

  par(mfrow = c(4, 2), mar = c(5, 5, 0.2, 0.2)) 
  
  for (i in seq_along(set_of_regions)){
    spd_file = paste0(results_directory, "/", gsub(" ", "_",  set_of_regions[i]), "_spd.rda")
    load(file = spd_file)
    spd = get(paste0(gsub(" ", "_",  set_of_regions[i]),"_spd"))
    if (i==1){
      plot(spd$grid$calBP,
           spd$grid$PrDens,
           xlim = time_range_calBP,
           ylim = c(0.000005,1.5),
           log="y",
           type ="l",
           xlab = "Years cal BP",
           ylab = "SPD",
           lwd = 2,
           col = cbPalette1[i])
    }else{
      lines(spd$grid$calBP,
            spd$grid$PrDens,
            col = cbPalette1[i], lwd = 2)   
    }
  }
  abline(v=8470, lty=2, col="grey")
  abline(v=8175, lty=2, col="grey")
  legend("bottomleft", legend=set_of_regions, lwd=2, col=cbPalette1[1:7])
  text(min(time_range_calBP), 1, "a", cex = 1.5)
  
  for (i in seq_along(set_of_regions)){
    skyline_file = paste0(results_directory, "/", gsub(" ", "_",  set_of_regions[i]), "_skyline.rda")
    load(file = skyline_file)
    skyline = get(paste0(gsub(" ", "_",  set_of_regions[i]),"_skyline"))
    plot(skyline_years,
         10^skyline$median,
         xlim = time_range_calBP,
         ylim = c(0.000005,1.5),
         log="y",
         type ="l",
         xlab = "Years cal BP",
         ylab = expression(lambda),
         lwd = 2,
         col = cbPalette1[i])
    polygon(c(skyline_years, rev(skyline_years)),
            c(10^skyline$q_0.025, rev(10^skyline$q_0.975)),
            col = cbPalette1_transparent[i], border = NA)
    text(min(time_range_calBP), 1, letters[i+1], cex = 1.5)
    abline(v=8470, lty=2, col="grey")
    abline(v=8175, lty=2, col="grey")
    text(mean(time_range_calBP), 1, set_of_regions[i], cex = 1.5)
  }  
  dev.off()
}  




