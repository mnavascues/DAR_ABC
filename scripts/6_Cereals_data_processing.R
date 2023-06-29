library(rcarbon)
source("../../Misc_R_tools/color.R")
source("scripts/sim14c.R")

load(file = "results/dates.rda")
load(file = "results/caldates.rda")
load(file = "results/Cereals_time_range_BP.rda")

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


save(cereals_dates, triticum_dates, hordeum_dates, file = "results/Cereals_dates.rda")
save(cereals_caldates, triticum_caldates, hordeum_caldates, file = "results/Cereals_caldates.rda")
save(cereals_spd, triticum_bins, triticum_spd,
     hordeum_bins, hordeum_spd, file = "results/Cereals_spd.rda")



pdf(file="results/Cereals_spd.pdf", width=10, height=5)

par(mar=c(4.5, 4.5, 1, 1) + 0.1)
plot(hordeum_spd$grid$calBP, hordeum_spd$grid$PrDens,
     xlab="Years cal BP", ylab="density",
     xlim = time_range_BP, ylim=c(0,0.4),
     type="l",lwd=2)
lines(triticum_spd$grid$calBP, triticum_spd$grid$PrDens, col=PCI_blue,lwd=2)
legend(6000,0.4,
       legend=c(expression(italic(Triticum)),expression(italic(Hordeum))),
       col=c(PCI_blue,"black"),
       lwd=2)
dev.off()


pdf(file="results/Cereals_hist.pdf", width=10, height=5)
par(mar=c(4.5, 4.5, 1, 1) + 0.1)
hordeum_weights = apply(as.array(hordeum_bins), 1, function(x) 1/sum(hordeum_bins == x))
breaks = seq(6000,200,-100)
HH = wtd.hist(hordeum_dates$CRA,
              breaks=breaks, 
              plot=T, weight = hordeum_weights,
              main="",
              xlim=time_range_BP, 
              xlab="Years BP (uncalibrated)",col="gray")

triticum_weights = apply(as.array(triticum_bins), 1, function(x) 1/sum(triticum_bins == x))
#breaks = seq(time_range_BP[1],time_range_BP[2]-1,-100)
TH = wtd.hist(triticum_dates$CRA, breaks=breaks, plot=T, weight = triticum_weights,
         main="",xlim=time_range_BP, xlab="Years BP (uncalibrated)",col=PCI_t_blue,add=T)
box()
dev.off()


# Calculate Summary statstics
triticum_sumstats_c14 = get_sumstats_uncalibrated(triticum_dates$CRA, 
                                                  time_range_BP, 
                                                  window = 100, 
                                                  w=triticum_weights)
hordeum_sumstats_c14 = get_sumstats_uncalibrated(hordeum_dates$CRA, 
                                                 time_range_BP, 
                                                 window = 100, 
                                                 w=hordeum_weights)
cereals_sumstats_c14 = get_sumstats_uncalibrated(cereals_dates$CRA, 
                                                 time_range_BP, 
                                                 window = 100, 
                                                 w=c(triticum_weights,hordeum_weights))
cereals_sumstats_cor = get_sumstats_correlation(triticum_sumstats_c14,
                                                hordeum_sumstats_c14)

names(triticum_sumstats_c14) = paste0(names(triticum_sumstats_c14),"_A")
names(hordeum_sumstats_c14) = paste0(names(hordeum_sumstats_c14),"_B")

cereals_sumstats_2_categories = cbind(triticum_sumstats_c14, hordeum_sumstats_c14,
                                      cereals_sumstats_c14, cereals_sumstats_cor)

save(cereals_sumstats_2_categories, cereals_sumstats_c14, file="results/Cereals_sumstats.rda")

