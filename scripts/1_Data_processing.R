library(rcarbon)
library(weights)
library(doSNOW)
library(doParallel)
source("../../Misc_R_tools/color.R")
source("scripts/sim14c.R")

# Load Bevan et al. 2017 data (obtained from DOI: 10.14324/000.ds.10025178)
if (!file.exists("results/Bevan_dates.rda")){
  dates = read.csv("data/gbie14Csub/dates/archdates.csv",
                 header = TRUE,
                 stringsAsFactors = FALSE,
                 encoding = "UTF-8",
                 na.strings = c("NA", ""),
                 strip.white = TRUE)
  head(dates)
  save(dates, file = "results/Bevan_dates.rda")
}

# Make histogram plot
if (!file.exists("results/Bevan_hist.pdf")){
  dates_on_range = intersect(which(dates$CRA<time_range_BP[1]), which(dates$CRA>time_range_BP[2]))
  window=100
  breaks = seq(time_range_BP[1],time_range_BP[2]-1,-window)
  bins = binPrep(sites = dates$SiteID, ages = dates$CRA, h = 100)
  w = apply(as.array(bins), 1, function(x) 1/sum(bins == x))
  pdf(file="results/Bevan_hist.pdf", width=10, height=5)
  par(mar=c(4.5, 4.5, 1, 1) + 0.1)
  wtd.hist(dates$CRA[dates_on_range], breaks=breaks, plot=T, weight = w[dates_on_range],
           main="",xlim=time_range_BP, xlab="Years BP (uncalibrated)",col=PCI_blue)
  box()
  dev.off()
  
  pdf(file="results/Bevan_hist_presentation.pdf", width=10, height=5)
  par(mar=c(4.5, 4.5, 1, 1) + 0.1)
  wtd.hist(dates$CRA[dates_on_range], breaks=breaks, plot=T, weight = w[dates_on_range],
           main="",xlim=time_range_BP, xlab="Years BP (uncalibrated)", col=INRAE_color)
  box()
  dev.off()
  
  pdf(file="results/calibration.pdf", width=5, height=5)
  random_date = sample(seq_along(dates$CRA), 1)
  calibrated_date_2_plot = calibrate(x=dates$CRA[random_date],
                                     errors=dates$Error[random_date],calCurves='intcal20')
  plot(calibrated_date_2_plot, HPD=TRUE,credMass=0.01)
  abline(h=dates$CRA[random_date], col="green", lty=2)
  dev.off()
}

# Calibrate dates
# (Difference with Bevan et al. 2017: using intcal20 instead of intcal13)
if (!file.exists("results/Bevan_caldates.rda")){
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
  save(caldates, file = "results/Bevan_caldates.rda")
}

# Make SPD plot
if (!file.exists("results/Bevan_spd.rda")){
  bins = binPrep(sites = dates$SiteID, ages = dates$CRA, h = 100)
  allspd = spd(x = caldates, 
               bins = bins, 
               timeRange = time_range_BP, 
               datenormalised = FALSE,
               runm = 100)
  save(bins, allspd, file = "results/Bevan_spd.rda")
  pdf(file="results/Bevan_spd.pdf", width=10, height=5)
  par(mar=c(4.5, 4.5, 1, 1) + 0.1)
  plot(allspd$grid$calBP, allspd$grid$PrDens, xlim = time_range_BP,
       type="l", xlab="Years cal BP", ylab="Sum of Probability Densities (SPD)", col=PCI_blue, lwd=2)
  dev.off()
  
  pdf(file="results/Bevan_spd_presentation.pdf", width=10, height=5)
  par(mar=c(4.5, 4.5, 1, 1) + 0.1)
  plot(allspd$grid$calBP, allspd$grid$PrDens, xlim = time_range_BP,
       type="l", xlab="Years cal BP", ylab="Sum of Probability Densities (SPD)", col=INRAE_color, lwd=2)
  dev.off()
}

# Calculate Summary statstics
if (!file.exists("results/Bevan_sumstats.rda")){
  bins = binPrep(sites = dates$SiteID, ages = dates$CRA, h = 100)
  w = apply(as.array(bins), 1, function(x) 1/sum(bins == x))
  all_sumstats_c14 = get_sumstats_uncalibrated(dates$CRA, time_range_BP, window = 100, w=w)
  all_sumstats_spd = get_sumstats_SPD(dates$CRA, dates$Error, calCurves = 'intcal20', bins = bins, time_range_BP)
  save(all_sumstats_c14, all_sumstats_spd, file="results/Bevan_sumstats.rda")
}


