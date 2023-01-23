library(abcrf)
source("scripts/sim14c.R")

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

# test all
if ( !file.exists("results/model_choice_0.rda") ){
  
  model = as.factor(c(rep("dynamic",nrow(reftable_dynamic)),
                      rep("logistic",nrow(reftable_logistic)),
                      rep("piecewise",nrow(reftable_piecewise))))
  # removing summary stats c(2,3,97) that "appear to be constant within groups" (for LDA)
  sumstats = rbind(reftable_dynamic[names(all_sumstats_c14)[-c(2,3,97)]],
                   reftable_logistic[names(all_sumstats_c14)[-c(2,3,97)]],
                   reftable_piecewise[names(all_sumstats_c14)[-c(2,3,97)]])
  
  RF_model_choice = abcrf(model~., data = data.frame(model,sumstats),
                          lda=T, ntree = 5000, paral = TRUE)
  
  plot(RF_model_choice, data.frame(model,sumstats))
  
  posterior_model_choice =  predict(RF_model_choice, all_sumstats_c14,
                                    training = data.frame(model,sumstats),
                                    ntree = 5000, paral = TRUE)
  confusion_matrix = RF_model_choice$model.rf$confusion.matrix
  prediction_error = RF_model_choice$model.rf$prediction.error
  save(prediction_error,
       confusion_matrix,
       posterior_model_choice, file="results/model_choice_0.rda")
}else{load(file="results/model_choice_0.rda")}
posterior_model_choice
confusion_matrix
prediction_error
K = posterior_model_choice$post.prob / (1-posterior_model_choice$post.prob)
(K)
interpret_K(K)

# test logistic vs. piecewise
if ( !file.exists("results/model_choice_1.rda") ){
  
  model = as.factor(c(rep("logistic", nrow(reftable_logistic)),
                      rep("piecewise", nrow(reftable_piecewise))))
  sumstats = rbind(reftable_logistic[names(all_sumstats_c14)[-c(2,3,97)]],
                   reftable_piecewise[names(all_sumstats_c14)[-c(2,3,97)]])
  
  RF_model_choice = abcrf(model~., data = data.frame(model,sumstats),
                          lda=T, ntree = 5000, paral = TRUE)
  
  posterior_model_choice =  predict(RF_model_choice, all_sumstats_c14,
                                    training = data.frame(model,sumstats),
                                    ntree = 5000, paral = TRUE) 
  confusion_matrix = RF_model_choice$model.rf$confusion.matrix
  prediction_error = RF_model_choice$model.rf$prediction.error
  save(prediction_error,
       confusion_matrix,
       posterior_model_choice, file="results/model_choice_1.rda")
}else{load(file="results/model_choice_1.rda")}
posterior_model_choice
K = posterior_model_choice$post.prob / (1-posterior_model_choice$post.prob)
(K)
interpret_K(K)


# test dynamic vs. piecewise
if ( !file.exists("results/model_choice_2.rda") ){
  
  model = as.factor(c(rep("dynamic", nrow(reftable_dynamic)),
                      rep("piecewise", nrow(reftable_piecewise))))
  sumstats = rbind(reftable_dynamic[names(all_sumstats_c14)],
                   reftable_piecewise[names(all_sumstats_c14)])
  
  RF_model_choice = abcrf(model~., data = data.frame(model,sumstats),
                          lda=F, ntree = 5000, paral = TRUE)
  
  posterior_model_choice =  predict(RF_model_choice, all_sumstats_c14,
                                    training = data.frame(model,sumstats),
                                    ntree = 5000, paral = TRUE) 
  confusion_matrix = RF_model_choice$model.rf$confusion.matrix
  prediction_error = RF_model_choice$model.rf$prediction.error
  save(prediction_error,
       confusion_matrix,
       posterior_model_choice, file="results/model_choice_2.rda")
}else{load(file="results/model_choice_2.rda")}
posterior_model_choice
K = posterior_model_choice$post.prob / (1-posterior_model_choice$post.prob)
(K)
interpret_K(K)



# test dynamic vs. logistic
if ( !file.exists("results/model_choice_3.rda") ){
  
  model = as.factor(c(rep("dynamic", nrow(reftable_dynamic)),
                      rep("logistic", nrow(reftable_logistic))))
  sumstats = rbind(reftable_dynamic[names(all_sumstats_c14)],
                   reftable_logistic[names(all_sumstats_c14)])
  
  RF_model_choice = abcrf(model~., data = data.frame(model,sumstats),
                          lda=F, ntree = 5000, paral = TRUE)
  
  posterior_model_choice =  predict(RF_model_choice, all_sumstats_c14,
                                    training = data.frame(model,sumstats),
                                    ntree = 5000, paral = TRUE) 
  confusion_matrix = RF_model_choice$model.rf$confusion.matrix
  prediction_error = RF_model_choice$model.rf$prediction.error
  save(prediction_error,
       confusion_matrix,
       posterior_model_choice, file="results/model_choice_3.rda")
}else{load(file="results/model_choice_3.rda")}
posterior_model_choice
K = posterior_model_choice$post.prob / (1-posterior_model_choice$post.prob)
(K)
interpret_K(K)

