source("scripts/sim14c.R")
source("../../Misc_R_tools/color.R")

load(file = "results/num_of_sims.rda")

# load target (i.e. observed) summary statistics
load(file = "results/sumstats.rda")

# lead reference tables for the 3 models
load(file = "results/logistic_model_reftable.rda")
rows2keep = complete.cases(reftable)
reftable_logistic = reftable[rows2keep,] #; head(reftable_logistic)
load(file = "results/piecewise_model_reftable.rda")
rows2keep = complete.cases(reftable)
reftable_piecewise = reftable[rows2keep,] #; head(reftable_piecewise)
load(file = "results/dynamic_model_reftable.rda")
rows2keep = complete.cases(reftable)
reftable_dynamic = reftable[rows2keep,] #; head(reftable_piecewise)
rm(reftable);gc()

sumstats = rbind(reftable_dynamic[names(all_sumstats_c14)[-c(2,3,97)]],
                 reftable_logistic[names(all_sumstats_c14)[-c(2,3,97)]],
                 reftable_piecewise[names(all_sumstats_c14)[-c(2,3,97)]])
model = c(rep.int("Dynamic", nrow(reftable_dynamic)),
          rep.int("Logistic", nrow(reftable_logistic)),
          rep.int("Piecewise exponential", nrow(reftable_piecewise)))

PCA_stats  = princomp(sumstats)
PCA_target = predict(PCA_stats, all_sumstats_c14)

summary(PCA_stats)

points_per_model = 2000
sims2plot = sort(c(sample(which(model==unique(model)[1]),points_per_model),
                   sample(which(model==unique(model)[2]),points_per_model),
                   sample(which(model==unique(model)[3]),points_per_model)))
colors2plot = c(rep(cbPalette1[2],points_per_model),
                rep(cbPalette1[3],points_per_model),
                rep(cbPalette1[4],points_per_model))
order2plot = sample(points_per_model*3)
sims2plot = sims2plot[order2plot]
colors2plot = colors2plot[order2plot]


i=1; j=2; xlim=c(-20000,100000); ylim=c(-30000,30000); id_plot="a"; add_legend=TRUE
i=3; j=4; xlim=c(-10000,10000); ylim=c(-10000,10000); id_plot="b"; add_legend=FALSE
i=5; j=6; xlim=c(-5000,5000); ylim=c(-5000,5000); id_plot="c"; add_legend=FALSE


# i=7; j=8; xlim=c(-5000,5000); ylim=c(-5000,5000); id_plot="d"; add_legend=FALSE
# i=9; j=10; xlim=c(-4000,4000); ylim=c(-4000,4000); id_plot="e"; add_legend=FALSE

#i=235; j=236; xlim=c(-6e-5,6e-5); ylim=c(-3e-4,3e-4); id_plot="d"; add_legend=FALSE

i=15; j=16

{
pdf_file_name = paste0("results/PCA_PC",i,"_PC",j,".pdf")
pdf(file=pdf_file_name, width=5, height=5)
par(mar=c(4.5, 4.5, 1, 1) + 0.1)
plot(PCA_stats$scores[sims2plot,i],
     PCA_stats$scores[sims2plot,j],
     xlab = paste0("PC",i),
     ylab = paste0("PC",j),
     #xlim=xlim,
     #ylim=ylim,
     col=colors2plot, cex=0.5)
points(PCA_target[,i],PCA_target[,j],pch="*",cex=3)
text(x=xlim[1], y=ylim[2], label=id_plot, cex=2)
if (add_legend) legend("bottomright", c("Dynamic","Logistic","Piecewise exponential"), pch=1, col=cbPalette1[2:4])
dev.off()
}
