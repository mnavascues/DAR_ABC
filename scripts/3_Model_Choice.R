library(abcrf)
source("scripts/sim14c.R")

load(file = "results/Bevan_num_of_sims.rda")

# load target (i.e. observed) summary statistics
load(file = "results/Bevan_sumstats.rda")

# lead reference tables for the 4 models
load(file = "results/Bevan_constant_model_reftable.rda")
reftable_constant = reftable #; head(reftable_constant)
load(file = "results/Bevan_exponential_model_reftable.rda")
reftable_exponential = reftable #; head(reftable_exponential)
load(file = "results/Bevan_logistic_model_reftable.rda")
reftable_logistic = reftable #; head(reftable_logistic)
load(file = "results/Bevan_piecewise_model_reftable.rda")
reftable_piecewise = reftable #; head(reftable_piecewise)
rm(reftable);gc()

# test constant vs. exponential
if ( !file.exists("results/Bevan_model_choice_1.rda") ){

  model = as.factor(c(rep("constant",num_of_sims),
                      rep("exponential",num_of_sims)))
  sumstats = rbind(reftable_constant[names(all_sumstats_c14)],
                   reftable_exponential[names(all_sumstats_c14)])
  
  RF_model_choice = abcrf(model~., data = data.frame(model,sumstats),
                          lda=F, ntree = 5000, paral = TRUE)
  
  posterior_model_choice =  predict(RF_model_choice, all_sumstats_c14,
                                    training = data.frame(model,sumstats),
                                    ntree = 5000, paral = TRUE) 
  save(RF_model_choice, posterior_model_choice, file="results/Bevan_model_choice_1.rda")
}
load(file="results/Bevan_model_choice_1.rda")
posterior_model_choice
K = posterior_model_choice$post.prob / (1-posterior_model_choice$post.prob)
(K)
interpret_K(K)

err.abcrf(RF_model_choice,data.frame(model,sumstats),paral=T)


# test exponential vs. logistic
if ( !file.exists("results/Bevan_model_choice_2.rda") ){
  
  model = as.factor(c(rep("exponential",num_of_sims),
                      rep("logistic",num_of_sims)))
  sumstats = rbind(reftable_exponential[names(all_sumstats_c14)],
                   reftable_logistic[names(all_sumstats_c14)])
  
  RF_model_choice = abcrf(model~., data = data.frame(model,sumstats),
                          lda=F, ntree = 5000, paral = TRUE)
  
  posterior_model_choice =  predict(RF_model_choice, all_sumstats_c14,
                                    training = data.frame(model,sumstats),
                                    ntree = 5000, paral = TRUE) 
  save(RF_model_choice, posterior_model_choice, file="results/Bevan_model_choice_2.rda")
}
load(file="results/Bevan_model_choice_2.rda")
posterior_model_choice
K = posterior_model_choice$post.prob / (1-posterior_model_choice$post.prob)
(K)
interpret_K(K)

# test logistic vs. piecewise
if ( !file.exists("results/Bevan_model_choice_3.rda") ){
  
  model = as.factor(c(rep("logistic",num_of_sims),
                      rep("piecewise",num_of_sims)))
  sumstats = rbind(reftable_logistic[names(all_sumstats_c14)],
                   reftable_piecewise[names(all_sumstats_c14)])
  
  RF_model_choice = abcrf(model~., data = data.frame(model,sumstats),
                          lda=F, ntree = 5000, paral = TRUE)
  
  posterior_model_choice =  predict(RF_model_choice, all_sumstats_c14,
                                    training = data.frame(model,sumstats),
                                    ntree = 5000, paral = TRUE) 
  save(RF_model_choice, posterior_model_choice, file="results/Bevan_model_choice_3.rda")
}
load(file="results/Bevan_model_choice_3.rda")
posterior_model_choice
K = posterior_model_choice$post.prob / (1-posterior_model_choice$post.prob)
(K)
interpret_K(K)

