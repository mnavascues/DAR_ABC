library(abcrf)
source("scripts/sim14c.R")

load(file = "results/Bevan_num_of_sims.rda")

# load target (i.e. observed) summary statistics
load(file = "results/Bevan_sumstats.rda")

# lead reference tables for the 2 models
load(file = "results/Bevan_logistic_model_reftable.rda")
reftable_logistic = reftable #; head(reftable_logistic)
load(file = "results/Bevan_piecewise_fixed_m_model_reftable.rda")
reftable_piecewise = reftable #; head(reftable_piecewise)
rm(reftable);gc()

# test logistic vs. piecewise
if ( !file.exists("results/Bevan_model_choice_4.rda") ){
  
  model = as.factor(c(rep("logistic",num_of_sims),
                      rep("piecewise",num_of_sims)))
  sumstats = rbind(reftable_logistic[names(all_sumstats_c14)],
                   reftable_piecewise[names(all_sumstats_c14)])
  
  RF_model_choice = abcrf(model~., data = data.frame(model,sumstats),
                          lda=F, ntree = 10000, paral = TRUE)
  
  posterior_model_choice =  predict(RF_model_choice, all_sumstats_c14,
                                    training = data.frame(model,sumstats),
                                    ntree = 10000, paral = TRUE) 
  save(RF_model_choice, posterior_model_choice, file="results/Bevan_model_choice_4.rda")
}
posterior_model_choice
K = posterior_model_choice$post.prob / (1-posterior_model_choice$post.prob)
(K)
interpret_K(K)

