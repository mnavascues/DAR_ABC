library(rcarbon)
library(dplyr, warn.conflicts = FALSE)

data = p3k14c::p3k14c_data[p3k14c::p3k14c_data$Continent=="Africa",]
# unique(data$Material)

p_millet = c(which(data$Material=="P. glaucum"),
             which(data$Material=="pearl millet"),
             which(data$Material=="Pennisestum americanum"),
             which(data$Material=="Pennisetu"),
             which(data$Material=="Pennisetum"),
             which(data$Material=="Pennisetum glaucum"),
             which(data$Material=="Pennisetum glaucum grain"),
             which(data$Material=="Pennisetum grain"),
             which(data$Material=="Pennisetum seeds"))
length(p_millet)

#sorghum = c(which(data$Material=="Sorghum"),
#            which(data$Material=="Sorghum seeds"))
#length(sorghum)

unique(data[p_millet,"Country"])

caldates = calibrate(x = pull(data,Age)[p_millet],
                     errors = pull(data,Error)[p_millet],
                     ids = pull(data,SiteName)[p_millet],
                     calCurves = 'intcal20',
                     timeRange = c(3500,1000),
                     normalised = FALSE,
                     calMatrix = TRUE)

bins = binPrep(sites = pull(data,SiteName)[p_millet], ages = pull(data,Age)[p_millet], h = 100)
allspd = spd(x = caldates, 
             bins = bins, 
             timeRange = c(3500,1000), 
             datenormalised = FALSE,
             runm = 100)
par(mar=c(4.5, 4.5, 1, 1) + 0.1)
plot(allspd$grid$calBP, allspd$grid$PrDens, xlim = c(3500,1000),
     type="l", xlab="Years cal BP", ylab="Sum of Probability Densities (SPD)", col="blue", lwd=2)



