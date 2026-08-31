##########
#### Multiple Outcomes Paper Simulations - Neyman Allocation Version
##########

### Author: Jasper Yang


#####
## Load packages
#####
library(dplyr)
library(tidyr)
library(optimall)
library(MASS)
library(survey)
library(parallel)
library(rprojroot)
library(optparse)
library(mice)

#####
## Source functions
#####
this_path <- normalizePath(".", mustWork = FALSE)
#proj_root <- rprojroot::find_root_file(criterion = ".projectile", path = this_path)

source(paste0(this_path, "/strategy_functions_priority_sensitivity.R"))
source(paste0(this_path, "/set_conditions.R"))
source(paste0(this_path, "/utils.R"))

simulations_df <- simulations_df1

#####
## Accept arguments
#####
parser <- OptionParser()
parser <- add_option(parser, "--scenario", type = "integer",
                     help = "scenario")
parser <- add_option(parser, "--pathname", type = "character",
                     help = "output path")
parser <- add_option(parser, "--n", type = "integer",
                     help = "output path")
parser <- add_option(parser, "--N", type = "integer",
                     help = "output path")
args <- parse_args(parser, convert_hyphens_to_underscores = TRUE)
print(args)

scenario <- args$scenario
pathname <- args$pathname
n <- args$n
N <- args$N

####
## Start simulations - Each scenario
####
print(simulations_df[,scenario])

# set first seed, for data - but set within functions when using cluster.
# set.seed(simulations_df["data_seed", scenario])

######
### Case 1: 01
######

# And iterate this 2500 times, storing the B-hat each time
cl <- makeCluster(detectCores() - 1)
clusterExport(cl, c("n", "simulations_df", "scenario", "Strategy_sensitivity", 
                    "inf_fun_logit","N"))

run_0 <- function(i) {
  set.seed(simulations_df["data_seed", scenario] + i)
  library(optimall)
  library(survey)
  library(dplyr)
  library(MASS)
  library(mice)
  run <- Strategy_sensitivity(n = n, simulations_df = simulations_df, scenario = scenario,
                   my_weights = c(0,1),
                   N = N)
  return(c(run[1], run[2], run[3], run[4], run[5], run[6], run[7],
           run[8], run[9], run[10], run[11], run[12],run[13],
           run[14], run[15], run[16], run[17], run[18]))
}

results <- parLapply(cl, 1:2500, run_0)
stopCluster(cl)

# Convert the results into a matrix and extract desired vectors
results_matrix <- do.call(rbind, results)
B_hat_Y1_01_GR <- results_matrix[, 1]
B_hat_Y2_01_GR <- results_matrix[, 2]
B_hat_Y1_01_IPW <- results_matrix[, 3]
B_hat_Y2_01_IPW <- results_matrix[, 4]
SE_hat_Y1_01_GR <- results_matrix[, 5]
SE_hat_Y2_01_GR <- results_matrix[, 6]
SE_hat_Y1_01_IPW <- results_matrix[, 7]
SE_hat_Y2_01_IPW <- results_matrix[, 8]
cover_Y1_01_IPW <- results_matrix[, 9]
cover_Y2_01_IPW <- results_matrix[, 10]
cover_Y1_01_GR <- results_matrix[, 11]
cover_Y2_01_GR <- results_matrix[, 12]
cor_Y1Y2_01 <- median(results_matrix[, 14])
prev_Y1_01 <- median(results_matrix[, 15])
prev_Y2_01 <- median(results_matrix[, 16])
B11_std_01 <- results_matrix[,17]
B12_std_01 <- results_matrix[,18]

# GR results
var_B_11_01_GR <- var(B_hat_Y1_01_GR)
var_B_12_01_GR <- var(B_hat_Y2_01_GR)
mean_B_11_01_GR <- mean(B_hat_Y1_01_GR)
mean_B_12_01_GR <- mean(B_hat_Y2_01_GR)
bias_B_11_01_GR <- mean(B_hat_Y1_01_GR - B11_std_01)
bias_B_12_01_GR <- mean(B_hat_Y2_01_GR - B12_std_01)
MSE_B_11_01_GR <- mean((B_hat_Y1_01_GR - 0.4)^2)
MSE_B_12_01_GR <- mean((B_hat_Y2_01_GR - 0.2)^2)
median_B_11_01_GR <- median(B_hat_Y1_01_GR)
median_B_12_01_GR <- median(B_hat_Y2_01_GR)
median_SE_hat_Y1_01_GR <- median(SE_hat_Y1_01_GR)
median_SE_hat_Y2_01_GR <- median(SE_hat_Y2_01_GR)
coverage_Y1_01_GR <- mean(cover_Y1_01_GR)
coverage_Y2_01_GR <- mean(cover_Y2_01_GR)

# IPW results
var_B_11_01_IPW <- var(B_hat_Y1_01_IPW)
var_B_12_01_IPW <- var(B_hat_Y2_01_IPW)
mean_B_11_01_IPW <- mean(B_hat_Y1_01_IPW)
mean_B_12_01_IPW <- mean(B_hat_Y2_01_IPW)
bias_B_11_01_IPW <- mean(B_hat_Y1_01_IPW - B11_std_01)
bias_B_12_01_IPW <- mean(B_hat_Y2_01_IPW - B12_std_01)
MSE_B_11_01_IPW <- mean((B_hat_Y1_01_IPW - 0.4)^2)
MSE_B_12_01_IPW <- mean((B_hat_Y2_01_IPW - 0.2)^2)
median_B_11_01_IPW <- median(B_hat_Y1_01_IPW)
median_B_12_01_IPW <- median(B_hat_Y2_01_IPW)
median_SE_hat_Y1_01_IPW <- median(SE_hat_Y1_01_IPW)
median_SE_hat_Y2_01_IPW <- median(SE_hat_Y2_01_IPW)
coverage_Y1_01_IPW <- mean(cover_Y1_01_IPW)
coverage_Y2_01_IPW <- mean(cover_Y2_01_IPW)

outpathname1 <- paste0(pathname,"/", names(simulations_df)[scenario], "01raw.csv")
write.csv(results_matrix, outpathname1, row.names = FALSE)

######
### 0.25 vs 0.75
######

# And iterate this 2500 times, storing the B-hat each time
cl <- makeCluster(detectCores() - 1)
clusterExport(cl, c("n", "simulations_df", "scenario", "Strategy_sensitivity", 
                    "inf_fun_logit","N"))

run_25 <- function(i) {
  set.seed(simulations_df["data_seed", scenario] + i)
  library(optimall)
  library(survey)
  library(dplyr)
  library(MASS)
  library(mice)
  run <- Strategy_sensitivity(n = n, simulations_df = simulations_df, scenario = scenario,
                              my_weights = c(0.25,0.75),
                              N = N)
  return(c(run[1], run[2], run[3], run[4], run[5], run[6], run[7],
           run[8], run[9], run[10], run[11], run[12],run[13],
           run[14], run[15], run[16], run[17], run[18]))
}

results <- parLapply(cl, 1:2500, run_25)
stopCluster(cl)

# Convert the results into a matrix and extract desired vectors
results_matrix <- do.call(rbind, results)
B_hat_Y1_25_GR <- results_matrix[, 1]
B_hat_Y2_25_GR <- results_matrix[, 2]
B_hat_Y1_25_IPW <- results_matrix[, 3]
B_hat_Y2_25_IPW <- results_matrix[, 4]
SE_hat_Y1_25_GR <- results_matrix[, 5]
SE_hat_Y2_25_GR <- results_matrix[, 6]
SE_hat_Y1_25_IPW <- results_matrix[, 7]
SE_hat_Y2_25_IPW <- results_matrix[, 8]
cover_Y1_25_IPW <- results_matrix[, 9]
cover_Y2_25_IPW <- results_matrix[, 10]
cover_Y1_25_GR <- results_matrix[, 11]
cover_Y2_25_GR <- results_matrix[, 12]
# Intentionally skip 13, which is n_oversampled
cor_Y1Y2_25 <- median(results_matrix[, 14])
prev_Y1_25 <- median(results_matrix[, 15])
prev_Y2_25 <- median(results_matrix[, 16])
B11_std_25 <- results_matrix[,17]
B12_std_25 <- results_matrix[,18]

# GR results
var_B_11_25_GR <- var(B_hat_Y1_25_GR)
var_B_12_25_GR <- var(B_hat_Y2_25_GR)
mean_B_11_25_GR <- mean(B_hat_Y1_25_GR)
mean_B_12_25_GR <- mean(B_hat_Y2_25_GR)
bias_B_11_25_GR <- mean(B_hat_Y1_25_GR - B11_std_25)
bias_B_12_25_GR <- mean(B_hat_Y2_25_GR - B12_std_25)
MSE_B_11_25_GR <- mean((B_hat_Y1_25_GR - 0.4)^2)
MSE_B_12_25_GR <- mean((B_hat_Y2_25_GR - 0.2)^2)
median_B_11_25_GR <- median(B_hat_Y1_25_GR)
median_B_12_25_GR <- median(B_hat_Y2_25_GR)
median_SE_hat_Y1_25_GR <- median(SE_hat_Y1_25_GR)
median_SE_hat_Y2_25_GR <- median(SE_hat_Y2_25_GR)
coverage_Y1_25_GR <- mean(cover_Y1_25_GR)
coverage_Y2_25_GR <- mean(cover_Y2_25_GR)

# IPW results
var_B_11_25_IPW <- var(B_hat_Y1_25_IPW)
var_B_12_25_IPW <- var(B_hat_Y2_25_IPW)
mean_B_11_25_IPW <- mean(B_hat_Y1_25_IPW)
mean_B_12_25_IPW <- mean(B_hat_Y2_25_IPW)
bias_B_11_25_IPW <- mean(B_hat_Y1_25_IPW - B11_std_25)
bias_B_12_25_IPW <- mean(B_hat_Y2_25_IPW - B12_std_25)
MSE_B_11_25_IPW <- mean((B_hat_Y1_25_IPW - 0.4)^2)
MSE_B_12_25_IPW <- mean((B_hat_Y2_25_IPW - 0.2)^2)
median_B_11_25_IPW <- median(B_hat_Y1_25_IPW)
median_B_12_25_IPW <- median(B_hat_Y2_25_IPW)
median_SE_hat_Y1_25_IPW <- median(SE_hat_Y1_25_IPW)
median_SE_hat_Y2_25_IPW <- median(SE_hat_Y2_25_IPW)
coverage_Y1_25_IPW <- mean(cover_Y1_25_IPW)
coverage_Y2_25_IPW <- mean(cover_Y2_25_IPW)

outpathname2 <- paste0(pathname,"/", names(simulations_df)[scenario], "25raw.csv")
write.csv(results_matrix, outpathname2, row.names = FALSE)

######
### 0.5 vs 0.5
######

# And iterate this 2500 times, storing the B-hat each time
cl <- makeCluster(detectCores() - 1)
clusterExport(cl, c("n", "simulations_df", "scenario", "Strategy_sensitivity", 
                    "inf_fun_logit","N"))

run_5 <- function(i) {
  set.seed(simulations_df["data_seed", scenario] + i)
  library(optimall)
  library(survey)
  library(dplyr)
  library(MASS)
  library(mice)
  run <- Strategy_sensitivity(n = n, simulations_df = simulations_df, scenario = scenario,
                              my_weights = c(0.5,0.5),
                              N = N)
  return(c(run[1], run[2], run[3], run[4], run[5], run[6], run[7],
           run[8], run[9], run[10], run[11], run[12],run[13],
           run[14], run[15], run[16], run[17], run[18]))
}

results <- parLapply(cl, 1:2500, run_5)
stopCluster(cl)

# Convert the results into a matrix and extract desired vectors
results_matrix <- do.call(rbind, results)
B_hat_Y1_5_GR <- results_matrix[, 1]
B_hat_Y2_5_GR <- results_matrix[, 2]
B_hat_Y1_5_IPW <- results_matrix[, 3]
B_hat_Y2_5_IPW <- results_matrix[, 4]
SE_hat_Y1_5_GR <- results_matrix[, 5]
SE_hat_Y2_5_GR <- results_matrix[, 6]
SE_hat_Y1_5_IPW <- results_matrix[, 7]
SE_hat_Y2_5_IPW <- results_matrix[, 8]
cover_Y1_5_IPW <- results_matrix[, 9]
cover_Y2_5_IPW <- results_matrix[, 10]
cover_Y1_5_GR <- results_matrix[, 11]
cover_Y2_5_GR <- results_matrix[, 12]
cor_Y1Y2_5 <- median(results_matrix[, 14])
prev_Y1_5 <- median(results_matrix[, 15])
prev_Y2_5 <- median(results_matrix[, 16])
B11_std_5 <- results_matrix[,17]
B12_std_5 <- results_matrix[,18]

# GR results
var_B_11_5_GR <- var(B_hat_Y1_5_GR)
var_B_12_5_GR <- var(B_hat_Y2_5_GR)
mean_B_11_5_GR <- mean(B_hat_Y1_5_GR)
mean_B_12_5_GR <- mean(B_hat_Y2_5_GR)
bias_B_11_5_GR <- mean(B_hat_Y1_5_GR - B11_std_5)
bias_B_12_5_GR <- mean(B_hat_Y2_5_GR - B12_std_5)
MSE_B_11_5_GR <- mean((B_hat_Y1_5_GR - 0.4)^2)
MSE_B_12_5_GR <- mean((B_hat_Y2_5_GR - 0.2)^2)
median_B_11_5_GR <- median(B_hat_Y1_5_GR)
median_B_12_5_GR <- median(B_hat_Y2_5_GR)
median_SE_hat_Y1_5_GR <- median(SE_hat_Y1_5_GR)
median_SE_hat_Y2_5_GR <- median(SE_hat_Y2_5_GR)
coverage_Y1_5_GR <- mean(cover_Y1_5_GR)
coverage_Y2_5_GR <- mean(cover_Y2_5_GR)

# IPW results
var_B_11_5_IPW <- var(B_hat_Y1_5_IPW)
var_B_12_5_IPW <- var(B_hat_Y2_5_IPW)
mean_B_11_5_IPW <- mean(B_hat_Y1_5_IPW)
mean_B_12_5_IPW <- mean(B_hat_Y2_5_IPW)
bias_B_11_5_IPW <- mean(B_hat_Y1_5_IPW - B11_std_5)
bias_B_12_5_IPW <- mean(B_hat_Y2_5_IPW - B12_std_5)
MSE_B_11_5_IPW <- mean((B_hat_Y1_5_IPW - 0.4)^2)
MSE_B_12_5_IPW <- mean((B_hat_Y2_5_IPW - 0.2)^2)
median_B_11_5_IPW <- median(B_hat_Y1_5_IPW)
median_B_12_5_IPW <- median(B_hat_Y2_5_IPW)
median_SE_hat_Y1_5_IPW <- median(SE_hat_Y1_5_IPW)
median_SE_hat_Y2_5_IPW <- median(SE_hat_Y2_5_IPW)
coverage_Y1_5_IPW <- mean(cover_Y1_5_IPW)
coverage_Y2_5_IPW <- mean(cover_Y2_5_IPW)

outpathname3 <- paste0(pathname,"/", names(simulations_df)[scenario], "5raw.csv")
write.csv(results_matrix, outpathname3, row.names = FALSE)

######
### 0.75 vs 0.25
######

# And iterate this 2500 times, storing the B-hat each time
cl <- makeCluster(detectCores() - 1)
clusterExport(cl, c("n", "simulations_df", "scenario", "Strategy_sensitivity", 
                    "inf_fun_logit","N"))

run_75 <- function(i) {
  set.seed(simulations_df["data_seed", scenario] + i)
  library(optimall)
  library(survey)
  library(dplyr)
  library(MASS)
  library(mice)
  run <- Strategy_sensitivity(n = n, simulations_df = simulations_df, scenario = scenario,
                              my_weights = c(0.75,0.25),
                              N = N)
  return(c(run[1], run[2], run[3], run[4], run[5], run[6], run[7],
           run[8], run[9], run[10], run[11], run[12],run[13],
           run[14], run[15], run[16], run[17], run[18]))
}

results <- parLapply(cl, 1:2500, run_75)
stopCluster(cl)

# Convert the results into a matrix and extract desired vectors
results_matrix <- do.call(rbind, results)
B_hat_Y1_75_GR <- results_matrix[, 1]
B_hat_Y2_75_GR <- results_matrix[, 2]
B_hat_Y1_75_IPW <- results_matrix[, 3]
B_hat_Y2_75_IPW <- results_matrix[, 4]
SE_hat_Y1_75_GR <- results_matrix[, 5]
SE_hat_Y2_75_GR <- results_matrix[, 6]
SE_hat_Y1_75_IPW <- results_matrix[, 7]
SE_hat_Y2_75_IPW <- results_matrix[, 8]
cover_Y1_75_IPW <- results_matrix[, 9]
cover_Y2_75_IPW <- results_matrix[, 10]
cover_Y1_75_GR <- results_matrix[, 11]
cover_Y2_75_GR <- results_matrix[, 12]
cor_Y1Y2_75 <- median(results_matrix[, 14])
prev_Y1_75 <- median(results_matrix[, 15])
prev_Y2_75 <- median(results_matrix[, 16])
B11_std_75 <- results_matrix[,17]
B12_std_75 <- results_matrix[,18]

# GR results
var_B_11_75_GR <- var(B_hat_Y1_75_GR)
var_B_12_75_GR <- var(B_hat_Y2_75_GR)
mean_B_11_75_GR <- mean(B_hat_Y1_75_GR)
mean_B_12_75_GR <- mean(B_hat_Y2_75_GR)
bias_B_11_75_GR <- mean(B_hat_Y1_75_GR - B11_std_75)
bias_B_12_75_GR <- mean(B_hat_Y2_75_GR - B12_std_75)
MSE_B_11_75_GR <- mean((B_hat_Y1_75_GR - 0.4)^2)
MSE_B_12_75_GR <- mean((B_hat_Y2_75_GR - 0.2)^2)
median_B_11_75_GR <- median(B_hat_Y1_75_GR)
median_B_12_75_GR <- median(B_hat_Y2_75_GR)
median_SE_hat_Y1_75_GR <- median(SE_hat_Y1_75_GR)
median_SE_hat_Y2_75_GR <- median(SE_hat_Y2_75_GR)
coverage_Y1_75_GR <- mean(cover_Y1_75_GR)
coverage_Y2_75_GR <- mean(cover_Y2_75_GR)

# IPW results
var_B_11_75_IPW <- var(B_hat_Y1_75_IPW)
var_B_12_75_IPW <- var(B_hat_Y2_75_IPW)
mean_B_11_75_IPW <- mean(B_hat_Y1_75_IPW)
mean_B_12_75_IPW <- mean(B_hat_Y2_75_IPW)
bias_B_11_75_IPW <- mean(B_hat_Y1_75_IPW - B11_std_75)
bias_B_12_75_IPW <- mean(B_hat_Y2_75_IPW - B12_std_75)
MSE_B_11_75_IPW <- mean((B_hat_Y1_75_IPW - 0.4)^2)
MSE_B_12_75_IPW <- mean((B_hat_Y2_75_IPW - 0.2)^2)
median_B_11_75_IPW <- median(B_hat_Y1_75_IPW)
median_B_12_75_IPW <- median(B_hat_Y2_75_IPW)
median_SE_hat_Y1_75_IPW <- median(SE_hat_Y1_75_IPW)
median_SE_hat_Y2_75_IPW <- median(SE_hat_Y2_75_IPW)
coverage_Y1_75_IPW <- mean(cover_Y1_75_IPW)
coverage_Y2_75_IPW <- mean(cover_Y2_75_IPW)

outpathname4 <- paste0(pathname,"/", names(simulations_df)[scenario], "75raw.csv")
write.csv(results_matrix, outpathname4, row.names = FALSE)

######
### 1 vs 0 
######

# And iterate this 2500 times, storing the B-hat each time
cl <- makeCluster(detectCores() - 1)
clusterExport(cl, c("n", "simulations_df", "scenario", "Strategy_sensitivity", 
                    "inf_fun_logit","N"))

run_100 <- function(i) {
  set.seed(simulations_df["data_seed", scenario] + i)
  library(optimall)
  library(survey)
  library(dplyr)
  library(MASS)
  library(mice)
  run <- Strategy_sensitivity(n = n, simulations_df = simulations_df, scenario = scenario,
                              my_weights = c(1,0),
                              N = N)
  return(c(run[1], run[2], run[3], run[4], run[5], run[6], run[7],
           run[8], run[9], run[10], run[11], run[12],run[13],
           run[14], run[15], run[16], run[17], run[18]))
}

results <- parLapply(cl, 1:2500, run_100)
stopCluster(cl)

# Convert the results into a matrix and extract desired vectors
results_matrix <- do.call(rbind, results)
B_hat_Y1_100_GR <- results_matrix[, 1]
B_hat_Y2_100_GR <- results_matrix[, 2]
B_hat_Y1_100_IPW <- results_matrix[, 3]
B_hat_Y2_100_IPW <- results_matrix[, 4]
SE_hat_Y1_100_GR <- results_matrix[, 5]
SE_hat_Y2_100_GR <- results_matrix[, 6]
SE_hat_Y1_100_IPW <- results_matrix[, 7]
SE_hat_Y2_100_IPW <- results_matrix[, 8]
cover_Y1_100_IPW <- results_matrix[, 9]
cover_Y2_100_IPW <- results_matrix[, 10]
cover_Y1_100_GR <- results_matrix[, 11]
cover_Y2_100_GR <- results_matrix[, 12]
# intentionally skip 13, which is number oversampled.
cor_Y1Y2_100 <- median(results_matrix[, 14])
prev_Y1_100 <- median(results_matrix[, 15])
prev_Y2_100 <- median(results_matrix[, 16])
B11_std_100 <- results_matrix[,17]
B12_std_100 <- results_matrix[,18]

# GR results
var_B_11_100_GR <- var(B_hat_Y1_100_GR)
var_B_12_100_GR <- var(B_hat_Y2_100_GR)
mean_B_11_100_GR <- mean(B_hat_Y1_100_GR)
mean_B_12_100_GR <- mean(B_hat_Y2_100_GR)
bias_B_11_100_GR <- mean(B_hat_Y1_100_GR - B11_std_100)
bias_B_12_100_GR <- mean(B_hat_Y2_100_GR - B12_std_100)
MSE_B_11_100_GR <- mean((B_hat_Y1_100_GR - 0.4)^2)
MSE_B_12_100_GR <- mean((B_hat_Y2_100_GR - 0.2)^2)
median_B_11_100_GR <- median(B_hat_Y1_100_GR)
median_B_12_100_GR <- median(B_hat_Y2_100_GR)
median_SE_hat_Y1_100_GR <- median(SE_hat_Y1_100_GR)
median_SE_hat_Y2_100_GR <- median(SE_hat_Y2_100_GR)
coverage_Y1_100_GR <- mean(cover_Y1_100_GR)
coverage_Y2_100_GR <- mean(cover_Y2_100_GR)

# IPW results
var_B_11_100_IPW <- var(B_hat_Y1_100_IPW)
var_B_12_100_IPW <- var(B_hat_Y2_100_IPW)
mean_B_11_100_IPW <- mean(B_hat_Y1_100_IPW)
mean_B_12_100_IPW <- mean(B_hat_Y2_100_IPW)
bias_B_11_100_IPW <- mean(B_hat_Y1_100_IPW - B11_std_100)
bias_B_12_100_IPW <- mean(B_hat_Y2_100_IPW - B12_std_100)
MSE_B_11_100_IPW <- mean((B_hat_Y1_100_IPW - 0.4)^2)
MSE_B_12_100_IPW <- mean((B_hat_Y2_100_IPW - 0.2)^2)
median_B_11_100_IPW <- median(B_hat_Y1_100_IPW)
median_B_12_100_IPW <- median(B_hat_Y2_100_IPW)
median_SE_hat_Y1_100_IPW <- median(SE_hat_Y1_100_IPW)
median_SE_hat_Y2_100_IPW <- median(SE_hat_Y2_100_IPW)
coverage_Y1_100_IPW <- mean(cover_Y1_100_IPW)
coverage_Y2_100_IPW <- mean(cover_Y2_100_IPW)

outpathname5 <- paste0(pathname,"/", names(simulations_df)[scenario], "100raw.csv")
write.csv(results_matrix, outpathname5, row.names = FALSE)


#####
## View results
#####

results <- data.frame(Scenario = c("01",".25",".5",".75","1"),
                      "true_B_11" = c(simulations_df["B11", scenario], 0,0,0,0),
                      "true_B_12" = c(simulations_df["B12", scenario], 0,0,0,0),
                      "true_B_11_std" = c(mean(B11_std_01),
                                          mean(B11_std_25),
                                          mean(B11_std_5),
                                          mean(B11_std_75),
                                          mean(B11_std_100)),
                      "true_B_12_std" = c(mean(B12_std_01),
                                          mean(B12_std_25),
                                          mean(B12_std_5),
                                          mean(B12_std_75),
                                          mean(B12_std_100)),
                      "cor_Y1_Y2" = c(cor_Y1Y2_01, cor_Y1Y2_25,
                                      cor_Y1Y2_5,cor_Y1Y2_75,
                                      cor_Y1Y2_100),
                      "prev_Y1" = c(prev_Y1_01, prev_Y1_25,
                                    prev_Y1_5, prev_Y1_75,
                                    prev_Y1_100),
                      "prev_Y2" = c(prev_Y2_01, prev_Y2_25,
                                    prev_Y2_5, prev_Y2_75,
                                    prev_Y2_100),
                      "mean(B_11)_GR" = c(mean_B_11_01_GR,
                                          mean_B_11_25_GR,
                                          mean_B_11_5_GR,
                                          mean_B_11_75_GR,
                                          mean_B_11_100_GR),
                      "mean(B_12)_GR" = c(mean_B_12_01_GR,
                                          mean_B_12_25_GR,
                                          mean_B_12_5_GR,
                                          mean_B_12_75_GR,
                                          mean_B_12_100_GR),
                      "median(B_11)_GR" = c(median_B_11_01_GR,
                                            median_B_11_25_GR,
                                            median_B_11_5_GR,
                                            median_B_11_75_GR,
                                            median_B_11_100_GR),
                      "median(B_12)_GR" = c(median_B_12_01_GR,
                                            median_B_12_25_GR,
                                            median_B_12_5_GR,
                                            median_B_12_75_GR,
                                            median_B_12_100_GR),
                      "var(B_11)_GR" = c(var_B_11_01_GR,
                                         var_B_11_25_GR,
                                         var_B_11_5_GR,
                                         var_B_11_75_GR,
                                         var_B_11_100_GR),
                      "var(B_12)_GR" = c(var_B_12_01_GR,
                                         var_B_12_25_GR,
                                         var_B_12_5_GR,
                                         var_B_12_75_GR,
                                         var_B_12_100_GR),
                      "mean(B_11)_IPW" = c(mean_B_11_01_IPW,
                                           mean_B_11_25_IPW,
                                           mean_B_11_5_IPW,
                                           mean_B_11_75_IPW,
                                           mean_B_11_100_IPW),
                      "mean(B_12)_IPW" = c(mean_B_12_01_IPW,
                                           mean_B_12_25_IPW,
                                           mean_B_12_5_IPW,
                                           mean_B_12_75_IPW,
                                           mean_B_12_100_IPW),
                      "median(B_11)_IPW" = c(median_B_11_01_IPW,
                                             median_B_11_25_IPW,
                                             median_B_11_5_IPW,
                                             median_B_11_75_IPW,
                                             median_B_11_100_IPW),
                      "median(B_12)_IPW" = c(median_B_12_01_IPW,
                                             median_B_12_25_IPW,
                                             median_B_12_5_IPW,
                                             median_B_12_75_IPW,
                                             median_B_12_100_IPW),
                      "var(B_11)_IPW" = c(var_B_11_01_IPW,
                                          var_B_11_25_IPW,
                                          var_B_11_5_IPW,
                                          var_B_11_75_IPW,
                                          var_B_11_100_IPW),
                      "var(B_12)_IPW" = c(var_B_12_01_IPW,
                                          var_B_12_25_IPW,
                                          var_B_12_5_IPW,
                                          var_B_12_75_IPW,
                                          var_B_12_100_IPW),
                      "coverage_Y1_IPW" = c(coverage_Y1_01_IPW,
                                            coverage_Y1_25_IPW,
                                            coverage_Y1_5_IPW,
                                            coverage_Y1_75_IPW,
                                            coverage_Y1_100_IPW),
                      "coverage_Y2_IPW" = c(coverage_Y2_01_IPW,
                                            coverage_Y2_25_IPW,
                                            coverage_Y2_5_IPW,
                                            coverage_Y2_75_IPW,
                                            coverage_Y2_100_IPW),
                      "coverage_Y1_GR" = c(coverage_Y1_01_GR,
                                           coverage_Y1_25_GR,
                                           coverage_Y1_5_GR,
                                           coverage_Y1_75_GR,
                                           coverage_Y1_100_GR),
                      "coverage_Y2_GR" = c(coverage_Y2_01_GR,
                                           coverage_Y2_25_GR,
                                           coverage_Y2_5_GR,
                                           coverage_Y2_75_GR,
                                           coverage_Y2_100_GR),
                      "Bias_Y1_GR" = c(bias_B_11_01_GR,
                                       bias_B_11_25_GR,
                                       bias_B_11_5_GR,
                                       bias_B_11_75_GR,
                                       bias_B_11_100_GR),
                      "Bias_Y2_GR" = c(bias_B_12_01_GR,
                                       bias_B_12_25_GR,
                                       bias_B_12_5_GR,
                                       bias_B_12_75_GR,
                                       bias_B_12_100_GR),
                      "Bias_Y1_IPW" = c(bias_B_11_01_IPW,
                                        bias_B_11_25_IPW,
                                        bias_B_11_5_IPW,
                                        bias_B_11_75_IPW,
                                        bias_B_11_100_IPW),
                      "Bias_Y2_IPW" = c(bias_B_12_01_IPW,
                                        bias_B_12_25_IPW,
                                        bias_B_12_5_IPW,
                                        bias_B_12_75_IPW,
                                        bias_B_12_100_IPW),
                      "median_ASE_Y1_IPW" = c(median_SE_hat_Y1_01_IPW,
                                              median_SE_hat_Y1_25_IPW,
                                              median_SE_hat_Y1_5_IPW,
                                              median_SE_hat_Y1_75_IPW,
                                              median_SE_hat_Y1_100_IPW),
                      "median_ASE_Y2_IPW" = c(median_SE_hat_Y2_01_IPW,
                                              median_SE_hat_Y2_25_IPW,
                                              median_SE_hat_Y2_5_IPW,
                                              median_SE_hat_Y2_75_IPW,
                                              median_SE_hat_Y2_100_IPW),
                      "median_ASE_Y1_GR" = c(median_SE_hat_Y1_01_GR,
                                             median_SE_hat_Y1_25_GR,
                                             median_SE_hat_Y1_5_GR,
                                             median_SE_hat_Y1_75_GR,
                                             median_SE_hat_Y1_100_GR),
                      "median_ASE_Y2_GR" = c(median_SE_hat_Y2_01_GR,
                                             median_SE_hat_Y2_25_GR,
                                             median_SE_hat_Y2_5_GR,
                                             median_SE_hat_Y2_75_GR,
                                             median_SE_hat_Y2_100_GR),
                      "MSE_Y1_IPW" = c(MSE_B_11_01_IPW,
                                       MSE_B_11_25_IPW,
                                       MSE_B_11_5_IPW,
                                       MSE_B_11_75_IPW,
                                       MSE_B_11_100_IPW),
                      "MSE_Y2_IPW" = c(MSE_B_12_01_IPW,
                                       MSE_B_12_25_IPW,
                                       MSE_B_12_5_IPW,
                                       MSE_B_12_75_IPW,
                                       MSE_B_12_100_IPW),
                      "MSE_Y1_GR" = c(MSE_B_11_01_GR,
                                       MSE_B_11_25_GR,
                                       MSE_B_11_5_GR,
                                       MSE_B_11_75_GR,
                                       MSE_B_11_100_GR),
                      "MSE_Y2_GR" = c(MSE_B_12_01_GR,
                                       MSE_B_12_25_GR,
                                       MSE_B_12_5_GR,
                                       MSE_B_12_75_GR,
                                       MSE_B_12_100_GR)
                      
)

print(results)


# save
outpathname <- paste0(pathname,"/", names(simulations_df)[scenario], "sensitivity", ".csv")
write.csv(results, outpathname, row.names = FALSE)

