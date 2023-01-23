library(abcrf)
source("scripts/sim14c.R")

load(file = "results/num_of_sims.rda")

# load target (i.e. observed) summary statistics
load(file = "results/Cereals_sumstats.rda")

# lead reference tables for the 3 models
load(file = "results/Cereals_independent_model_reftable.rda")
reftable = reftable[!is.na(reftable$count_A),]
reftable = reftable[reftable$count_A!=1,]
reftable = reftable[!is.na(reftable$count_B),]
reftable = reftable[reftable$count_B!=1,]
reftable_independent = reftable
load(file = "results/Cereals_correlated_model_reftable.rda")
reftable = reftable[!is.na(reftable$count_A),]
reftable = reftable[reftable$count_A!=1,]
reftable = reftable[!is.na(reftable$count_B),]
reftable = reftable[reftable$count_B!=1,]
reftable_correlated = reftable
load(file = "results/Cereals_proportional_model_reftable.rda")
reftable = reftable[!is.na(reftable$count_A),]
reftable = reftable[reftable$count_A!=1,]
reftable = reftable[!is.na(reftable$count_B),]
reftable = reftable[reftable$count_B!=1,]
reftable_proportional = reftable
rm(reftable);gc()




# test independent vs. proportional
if ( !file.exists("results/Cereals_model_choice_1.rda") ){

  model = as.factor(c(rep("independent",nrow(reftable_independent)),
                      rep("proportional",nrow(reftable_proportional))))
  sumstats = rbind(reftable_independent[names(cereals_sumstats_2_categories)],
                   reftable_proportional[names(cereals_sumstats_2_categories)])
  
  RF_model_choice = abcrf(model~., data = data.frame(model,sumstats),
                          lda=F, ntree = 1000, paral = TRUE)
  
  posterior_model_choice =  predict(RF_model_choice, cereals_sumstats_2_categories,
                                    training = data.frame(model,sumstats),
                                    ntree = 1000, paral = TRUE) 
  save(RF_model_choice, posterior_model_choice, file="results/Cereals_model_choice_1.rda")
}
#load(file="results/Cereals_model_choice_1.rda")
posterior_model_choice
K = posterior_model_choice$post.prob / (1-posterior_model_choice$post.prob)
(K)
interpret_K(K)

err.abcrf(RF_model_choice,data.frame(model,sumstats),paral=T)



# test correlated vs. proportional
if ( !file.exists("results/Cereals_model_choice_2.rda") ){
  
  model = as.factor(c(rep("correlated",nrow(reftable_correlated)),
                      rep("proportional",nrow(reftable_proportional))))
  sumstats = rbind(reftable_correlated[names(cereals_sumstats_2_categories)],
                   reftable_proportional[names(cereals_sumstats_2_categories)])
  
  RF_model_choice = abcrf(model~., data = data.frame(model,sumstats),
                          lda=F, ntree = 1000, paral = TRUE)
  
  posterior_model_choice =  predict(RF_model_choice, cereals_sumstats_2_categories,
                                    training = data.frame(model,sumstats),
                                    ntree = 1000, paral = TRUE) 
  save(RF_model_choice, posterior_model_choice, file="results/Cereals_model_choice_2.rda")
}
#load(file="results/Cereals_model_choice_2.rda")
posterior_model_choice
K = posterior_model_choice$post.prob / (1-posterior_model_choice$post.prob)
(K)
interpret_K(K)

err.abcrf(RF_model_choice,data.frame(model,sumstats),paral=T)






# test independent vs. correlated
if ( !file.exists("results/Cereals_model_choice_3.rda") ){
  
  model = as.factor(c(rep("independent",nrow(reftable_independent)),
                      rep("correlated",nrow(reftable_correlated))))
  sumstats = rbind(reftable_independent[names(cereals_sumstats_2_categories)],
                   reftable_correlated[names(cereals_sumstats_2_categories)])
  
  RF_model_choice = abcrf(model~., data = data.frame(model,sumstats),
                          lda=F, ntree = 1000, paral = TRUE)
  
  posterior_model_choice =  predict(RF_model_choice, cereals_sumstats_2_categories,
                                    training = data.frame(model,sumstats),
                                    ntree = 1000, paral = TRUE) 
  save(RF_model_choice, posterior_model_choice, file="results/Cereals_model_choice_3.rda")
}
#load(file="results/Cereals_model_choice_3.rda")
posterior_model_choice
K = posterior_model_choice$post.prob / (1-posterior_model_choice$post.prob)
(K)
interpret_K(K)

err.abcrf(RF_model_choice,data.frame(model,sumstats),paral=T)






# test 3 models
if ( !file.exists("results/Cereals_model_choice_4.rda") ){
  
  model = as.factor(c(rep("independent",nrow(reftable_independent)),
                      rep("correlated",nrow(reftable_correlated)),
                      rep("proportional",nrow(reftable_proportional))))
  sumstats = rbind(reftable_independent[names(cereals_sumstats_2_categories)],
                   reftable_correlated[names(cereals_sumstats_2_categories)],
                   reftable_proportional[names(cereals_sumstats_2_categories)])
  
  RF_model_choice = abcrf(model~., data = data.frame(model,sumstats),
                          lda=F, ntree = 1000, paral = TRUE)
  
  posterior_model_choice =  predict(RF_model_choice, cereals_sumstats_2_categories,
                                    training = data.frame(model,sumstats),
                                    ntree = 1000, paral = TRUE) 
  save(RF_model_choice, posterior_model_choice, file="results/Cereals_model_choice_4.rda")
}
#load(file="results/Cereals_model_choice_4.rda")
posterior_model_choice

err.abcrf(RF_model_choice,data.frame(model,sumstats),paral=T)



