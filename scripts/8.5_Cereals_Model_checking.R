source("scripts/sim14c.R")
source("../../Misc_R_tools/color.R")

load(file = "results/num_of_sims.rda")

# load target (i.e. observed) summary statistics
load(file = "results/Cereals_sumstats.rda")

# lead reference tables for the 3 models
load(file = "results/Cereals_independent_model_reftable.rda")
reftable = reftable[names(cereals_sumstats_2_categories)]
rows2keep = complete.cases(reftable)
reftable_independent = reftable[rows2keep,] #; head(reftable_independent)
load(file = "results/Cereals_interdependent_model_reftable.rda")
reftable = reftable[names(cereals_sumstats_2_categories)]
rows2keep = complete.cases(reftable)
reftable_interdependent = reftable[rows2keep,] #; head(reftable_interdependent)
load(file = "results/Cereals_parallel_model_reftable.rda")
reftable = reftable[names(cereals_sumstats_2_categories)]
rows2keep = complete.cases(reftable)
reftable_parallel = reftable[rows2keep,] #; head(reftable_interdependent)
rm(reftable);gc()

sumstats = rbind(reftable_parallel,
                 reftable_independent,
                 reftable_interdependent)
model = c(rep.int("Parallel", nrow(reftable_parallel)),
          rep.int("Independent", nrow(reftable_independent)),
          rep.int("Interdependent", nrow(reftable_interdependent)))

PCA_stats  = princomp(sumstats)
PCA_target = predict(PCA_stats, cereals_sumstats_2_categories)

summary(PCA_stats)
# proportion of variance explained by first 6 PC
sum((PCA_stats$sdev[1:6])^2)/sum((PCA_stats$sdev)^2)

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


i=1; j=2; xlim=c(-1000000,2000000); ylim=c(-25000,30000); id_plot="a"; add_legend=TRUE
i=3; j=4; xlim=c(-12000,12000); ylim=c(-6000,12000); id_plot="b"; add_legend=FALSE
i=5; j=6; xlim=c(-5000,5000); ylim=c(-5000,5000); id_plot="c"; add_legend=FALSE
i=7; j=8; xlim=c(-5000,5000); ylim=c(-5000,5000); id_plot="d"; add_legend=FALSE
i=9; j=10; xlim=c(-4000,4000); ylim=c(-4000,4000); id_plot="e"; add_legend=FALSE


pdf_file_name = paste0("results/Cereals_PCA_PC",i,"_PC",j,".pdf")
pdf(file=pdf_file_name, width=5, height=5)
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
if (add_legend) legend("bottomright", c("Parallel","Independent","Interdependent"), pch=1, col=cbPalette1[2:4])
dev.off()

i=561; j=562
i=562; j=563
i=100; j=101
i=200; j=201
i=300; j=301
i=400; j=401
i=500; j=501
par(mar=c(4.5, 4.5, 1, 1) + 0.1)
plot(PCA_stats$scores[sims2plot,i],
     PCA_stats$scores[sims2plot,j],
     xlab = paste0("PC",i),
     ylab = paste0("PC",j),
     col=colors2plot, cex=0.5)
points(PCA_target[,i],PCA_target[,j],pch="*",cex=4)
