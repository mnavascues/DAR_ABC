# set time range for analysis
time_range_BP = c(10000, 500)
t_ = seq(time_range_BP[1], time_range_BP[2] + 1,by = -1)
save(time_range_BP, t_, file = "results/time_range_BP.rda")
time_range_BP = c(6000,500)
t_ = seq(time_range_BP[1], time_range_BP[2] + 1,by = -1)
save(time_range_BP, t_, file = "results/Cereals_time_range_BP.rda")

# set number of simulations
num_of_sims = 100000
save(num_of_sims, file = "results/num_of_sims.rda")

# set prior parameters
lambda_min = 0.001
lambda_max = 12
save(lambda_min,lambda_max, file = "results/lambda_prior.rda")
lambda_min = 0.001
lambda_max = 2
save(lambda_min,lambda_max, file = "results/Cereals_lambda_prior.rda")


