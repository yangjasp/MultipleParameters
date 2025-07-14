##########
#### Multiple Outcomes Paper Simulations - Neyman Allocation Version
##########

### Author: Jasper Yang

### Note: These sims are same as dissertation, and strata are defined on X as in
### dissertation. The median is reported here.
### These simulations also have Strategy 3.5

###
###
### In this version, exact sampling probabilities are used to calculate the
### weights

###
###
### Update 7/19: Data generation is now inside the loops - every simulation
### regenerates new data.


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

source(paste0(this_path, "/strategy_functions_scenario2.R"))
source(paste0(this_path, "/set_conditions.R"))
source(paste0(this_path, "/utils.R"))

simulations_df <- simulations_df2

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
### Phase 2 Sampling - Strategy 1: Case-control
######

# And iterate this 1000 times, storing the B-hat each time
cl <- makeCluster(detectCores() - 1)
clusterExport(cl, c("n", "simulations_df", "scenario", "Strategy1",
                    "inf_fun_logit", "N"))

run_Strategy1 <- function(i) {
  set.seed(simulations_df["data_seed", scenario] + i)
  library(optimall)
  library(survey)
  library(dplyr)
  library(MASS)
  library(mice)
  run <- Strategy1(n = n, simulations_df = simulations_df, 
                   scenario = scenario, N = N)
  return(c(run[1], run[2], run[3], run[4], run[5], run[6], run[7],
           run[8], run[9], run[10], run[11], run[12],run[13],
           run[14], run[15], run[16], run[17]))
}

results <- parLapply(cl, 1:2500, run_Strategy1)
stopCluster(cl)

# Convert the results into a matrix and extract desired vectors
results_matrix <- do.call(rbind, results)
B_hat_X1_strat1_GR <- results_matrix[, 1]
B_hat_X2_strat1_GR <- results_matrix[, 2]
B_hat_X1_strat1_IPW <- results_matrix[, 3]
B_hat_X2_strat1_IPW <- results_matrix[, 4]
SE_hat_X1_strat1_GR <- results_matrix[, 5]
SE_hat_X2_strat1_GR <- results_matrix[, 6]
SE_hat_X1_strat1_IPW <- results_matrix[, 7]
SE_hat_X2_strat1_IPW <- results_matrix[, 8]
cover_X1_strat1_GR <- results_matrix[, 9]
cover_X2_strat1_GR <- results_matrix[, 10]
cover_X1_strat1_IPW <- results_matrix[, 11]
cover_X2_strat1_IPW <- results_matrix[, 12]
cor_X1X2_strat1 <- median(results_matrix[, 13])
prev_Y_strat1 <- median(results_matrix[, 14])
prev_Z_strat1 <- median(results_matrix[, 15])
BX1_std_strat1 <- results_matrix[,16]
BX2_std_strat1 <- results_matrix[,17]

# GR results
var_BX1_strat1_GR <- var(B_hat_X1_strat1_GR)
var_BX2_strat1_GR <- var(B_hat_X2_strat1_GR)
mean_BX1_strat1_GR <- mean(B_hat_X1_strat1_GR)
mean_BX2_strat1_GR <- mean(B_hat_X2_strat1_GR)
bias_BX1_strat1_GR <- mean(B_hat_X1_strat1_GR - BX1_std_strat1)
bias_BX2_strat1_GR <- mean(B_hat_X2_strat1_GR - BX2_std_strat1)
MSE_BX1_strat1_GR <- mean((B_hat_X1_strat1_GR - 0.3)^2)
MSE_BX2_strat1_GR <- mean((B_hat_X2_strat1_GR - 0.7)^2)
median_BX1_strat1_GR <- median(B_hat_X1_strat1_GR)
median_BX2_strat1_GR <- median(B_hat_X2_strat1_GR)
median_SE_hat_X1_strat1_GR <- median(SE_hat_X1_strat1_GR)
median_SE_hat_X2_strat1_GR <- median(SE_hat_X2_strat1_GR)
coverage_X1_strat1_GR <- mean(cover_X1_strat1_GR)
coverage_X2_strat1_GR <- mean(cover_X2_strat1_GR)

# IPW results
var_BX1_strat1_IPW <- var(B_hat_X1_strat1_IPW)
var_BX2_strat1_IPW <- var(B_hat_X2_strat1_IPW)
mean_BX1_strat1_IPW <- mean(B_hat_X1_strat1_IPW)
mean_BX2_strat1_IPW <- mean(B_hat_X2_strat1_IPW)
bias_BX1_strat1_IPW <- mean(B_hat_X1_strat1_IPW - BX1_std_strat1)
bias_BX2_strat1_IPW <- mean(B_hat_X2_strat1_IPW - BX2_std_strat1)
MSE_BX1_strat1_IPW <- mean((B_hat_X1_strat1_IPW - 0.3)^2)
MSE_BX2_strat1_IPW <- mean((B_hat_X2_strat1_IPW - 0.7)^2)
median_BX1_strat1_IPW <- median(B_hat_X1_strat1_IPW)
median_BX2_strat1_IPW <- median(B_hat_X2_strat1_IPW)
median_SE_hat_X1_strat1_IPW <- median(SE_hat_X1_strat1_IPW)
median_SE_hat_X2_strat1_IPW <- median(SE_hat_X2_strat1_IPW)
coverage_X1_strat1_IPW <- mean(cover_X1_strat1_IPW)
coverage_X2_strat1_IPW <- mean(cover_X2_strat1_IPW)

outpathname1 <- paste0(pathname,"/", names(simulations_df)[scenario], "Strat1raw.csv")
write.csv(results_matrix, outpathname1, row.names = FALSE)

######
### Phase 2 Sampling - Strategy 2: Independent, uniform strategy over 4 waves
######

## Here, each wave chooses n/4 samples where n/8 are chosen to be optimal for
## BX1 and n/8 to be optimal for BX2. We calculate Wave 1 
## allocation outside of the loop because allocation (but not sampling)
## is dependent on the fixed Phase 1 data.

# And iterate this 1000 times, storing the B-hat each time
cl <- makeCluster(detectCores() - 1)
clusterExport(cl, c("n", "simulations_df", "scenario", "Strategy2", 
                    "inf_fun_logit", "N"))

run_Strategy2 <- function(i) {
  set.seed(simulations_df["data_seed", scenario] + i)
  library(optimall)
  library(survey)
  library(dplyr)
  library(MASS)
  library(mice)
  run <- suppressWarnings(Strategy2(n = n, simulations_df = simulations_df,
                                    scenario = scenario, N = N))
  return(c(run[1], run[2], run[3], run[4], run[5], run[6], run[7],
           run[8], run[9], run[10], run[11], run[12], run[13],
           run[14], run[15], run[16], run[17], run[18], run[19]))
}

results <- parLapply(cl, 1:2500, run_Strategy2)
stopCluster(cl)

# Convert the results into a matrix and extract desired vectors
results_matrix <- do.call(rbind, results)
B_hat_X1_strat2_GR <- results_matrix[, 1]
B_hat_X2_strat2_GR <- results_matrix[, 2]
B_hat_X1_strat2_IPW <- results_matrix[, 3]
B_hat_X2_strat2_IPW <- results_matrix[, 4]
SE_hat_X1_strat2_GR <- results_matrix[, 5]
SE_hat_X2_strat2_GR <- results_matrix[, 6]
SE_hat_X1_strat2_IPW <- results_matrix[, 7]
SE_hat_X2_strat2_IPW <- results_matrix[, 8]
cover_X1_strat2_GR <- results_matrix[, 9]
cover_X2_strat2_GR <- results_matrix[, 10]
cover_X1_strat2_IPW <- results_matrix[, 11]
cover_X2_strat2_IPW <- results_matrix[, 12]
n_oversampled_X1_strat2 <- sum(results_matrix[,13])
n_oversampled_X2_strat2 <- sum(results_matrix[,14])
cor_X1X2_strat2 <- median(results_matrix[, 15])
prev_Y_strat2 <- median(results_matrix[, 16])
prev_Z_strat2 <- median(results_matrix[, 17])
BX1_std_strat2 <- results_matrix[,18]
BX2_std_strat2 <- results_matrix[,19]

# GR results
var_BX1_strat2_GR <- var(B_hat_X1_strat2_GR)
var_BX2_strat2_GR <- var(B_hat_X2_strat2_GR)
mean_BX1_strat2_GR <- mean(B_hat_X1_strat2_GR)
mean_BX2_strat2_GR <- mean(B_hat_X2_strat2_GR)
bias_BX1_strat2_GR <- mean(B_hat_X1_strat2_GR - BX1_std_strat2)
bias_BX2_strat2_GR <- mean(B_hat_X2_strat2_GR - BX2_std_strat2)
MSE_BX1_strat2_GR <- mean((B_hat_X1_strat2_GR - 0.3)^2)
MSE_BX2_strat2_GR <- mean((B_hat_X2_strat2_GR - 0.7)^2)
median_BX1_strat2_GR <- median(B_hat_X1_strat2_GR)
median_BX2_strat2_GR <- median(B_hat_X2_strat2_GR)
median_SE_hat_X1_strat2_GR <- median(SE_hat_X1_strat2_GR)
median_SE_hat_X2_strat2_GR <- median(SE_hat_X2_strat2_GR)
coverage_X1_strat2_GR <- mean(cover_X1_strat2_GR)
coverage_X2_strat2_GR <- mean(cover_X2_strat2_GR)

# IPW results
var_BX1_strat2_IPW <- var(B_hat_X1_strat2_IPW)
var_BX2_strat2_IPW <- var(B_hat_X2_strat2_IPW)
mean_BX1_strat2_IPW <- mean(B_hat_X1_strat2_IPW)
mean_BX2_strat2_IPW <- mean(B_hat_X2_strat2_IPW)
bias_BX1_strat2_IPW <- mean(B_hat_X1_strat2_IPW - BX1_std_strat2)
bias_BX2_strat2_IPW <- mean(B_hat_X2_strat2_IPW - BX2_std_strat2)
MSE_BX1_strat2_IPW <- mean((B_hat_X1_strat2_IPW - 0.3)^2)
MSE_BX2_strat2_IPW <- mean((B_hat_X2_strat2_IPW - 0.7)^2)
median_BX1_strat2_IPW <- median(B_hat_X1_strat2_IPW)
median_BX2_strat2_IPW <- median(B_hat_X2_strat2_IPW)
median_SE_hat_X1_strat2_IPW <- median(SE_hat_X1_strat2_IPW)
median_SE_hat_X2_strat2_IPW <- median(SE_hat_X2_strat2_IPW)
coverage_X1_strat2_IPW <- mean(cover_X1_strat2_IPW)
coverage_X2_strat2_IPW <- mean(cover_X2_strat2_IPW)

# Print raw results
outpathname2 <- paste0(pathname,"/", names(simulations_df)[scenario], "Strat2raw.csv")
write.csv(results_matrix, outpathname2, row.names = FALSE)

######
######
### Phase 2 Sampling - Strategy 3: Independent sequential strategy over 4 waves
######
######

## Here, 1st and 2nd waves (n/4 samples each) are optimized for BX1 and 
## 3rd and 4th last n/4 are optimized for BX2. We calculate Wave 1 
## allocation outside of the loop because allocation (but not sampling)
## is dependent on the fixed Phase 1 data.

# Run the function for 1,000 iterations,  storing the B-hat each time
cl <- makeCluster(detectCores() - 1)
clusterExport(cl, c("n", "simulations_df", "scenario", "Strategy3",
                    "inf_fun_logit", "N"))

run_Strategy3 <- function(i) {
  set.seed(simulations_df["data_seed", scenario] + i)
  library(optimall)
  library(survey)
  library(dplyr)
  library(MASS)
  library(mice)
  my_it <<- i
  run <- suppressWarnings(Strategy3(n = n, simulations_df = simulations_df,
                                    scenario = scenario, N = N))
  return(c(run[1], run[2], run[3], run[4], run[5], run[6], run[7],
           run[8], run[9], run[10], run[11], run[12], run[13],
           run[14], run[15], run[16], run[17], run[18]))
}

results <- parLapply(cl, 1:2500, run_Strategy3)
stopCluster(cl)


# Convert the results into a matrix and extract desired vectors
results_matrix <- do.call(rbind, results)
B_hat_X1_strat3_GR <- results_matrix[, 1]
B_hat_X2_strat3_GR <- results_matrix[, 2]
B_hat_X1_strat3_IPW <- results_matrix[, 3]
B_hat_X2_strat3_IPW <- results_matrix[, 4]
SE_hat_X1_strat3_GR <- results_matrix[, 5]
SE_hat_X2_strat3_GR <- results_matrix[, 6]
SE_hat_X1_strat3_IPW <- results_matrix[, 7]
SE_hat_X2_strat3_IPW <- results_matrix[, 8]
cover_X1_strat3_GR <- results_matrix[, 9]
cover_X2_strat3_GR <- results_matrix[, 10]
cover_X1_strat3_IPW <- results_matrix[, 11]
cover_X2_strat3_IPW <- results_matrix[, 12]
n_oversampled_X2_strat3 <- sum(results_matrix[,13])
cor_X1X2_strat3 <- median(results_matrix[, 14])
prev_Y_strat3 <- median(results_matrix[, 15])
prev_Z_strat3 <- median(results_matrix[, 16])
BX1_std_strat3 <- results_matrix[,17]
BX2_std_strat3 <- results_matrix[,18]

# GR results
var_BX1_strat3_GR <- var(B_hat_X1_strat3_GR)
var_BX2_strat3_GR <- var(B_hat_X2_strat3_GR)
mean_BX1_strat3_GR <- mean(B_hat_X1_strat3_GR)
mean_BX2_strat3_GR <- mean(B_hat_X2_strat3_GR)
bias_BX1_strat3_GR <- mean(B_hat_X1_strat3_GR - BX1_std_strat3)
bias_BX2_strat3_GR <- mean(B_hat_X2_strat3_GR - BX2_std_strat3)
MSE_BX1_strat3_GR <- mean((B_hat_X1_strat3_GR - 0.3)^2)
MSE_BX2_strat3_GR <- mean((B_hat_X2_strat3_GR - 0.7)^2)
median_BX1_strat3_GR <- median(B_hat_X1_strat3_GR)
median_BX2_strat3_GR <- median(B_hat_X2_strat3_GR)
median_SE_hat_X1_strat3_GR <- median(SE_hat_X1_strat3_GR)
median_SE_hat_X2_strat3_GR <- median(SE_hat_X2_strat3_GR)
coverage_X1_strat3_GR <- mean(cover_X1_strat3_GR)
coverage_X2_strat3_GR <- mean(cover_X2_strat3_GR)

# IPW results
var_BX1_strat3_IPW <- var(B_hat_X1_strat3_IPW)
var_BX2_strat3_IPW <- var(B_hat_X2_strat3_IPW)
mean_BX1_strat3_IPW <- mean(B_hat_X1_strat3_IPW)
mean_BX2_strat3_IPW <- mean(B_hat_X2_strat3_IPW)
bias_BX1_strat3_IPW <- mean(B_hat_X1_strat3_IPW - BX1_std_strat3)
bias_BX2_strat3_IPW <- mean(B_hat_X2_strat3_IPW - BX2_std_strat3)
MSE_BX1_strat3_IPW <- mean((B_hat_X1_strat3_IPW - 0.3)^2)
MSE_BX2_strat3_IPW <- mean((B_hat_X2_strat3_IPW - 0.7)^2)
median_BX1_strat3_IPW <- median(B_hat_X1_strat3_IPW)
median_BX2_strat3_IPW <- median(B_hat_X2_strat3_IPW)
median_SE_hat_X1_strat3_IPW <- median(SE_hat_X1_strat3_IPW)
median_SE_hat_X2_strat3_IPW <- median(SE_hat_X2_strat3_IPW)
coverage_X1_strat3_IPW <- mean(cover_X1_strat3_IPW)
coverage_X2_strat3_IPW <- mean(cover_X2_strat3_IPW)

#Save raw results
outpathname3 <- paste0(pathname,"/", names(simulations_df)[scenario], "Strat3raw.csv")
write.csv(results_matrix, outpathname3, row.names = FALSE)


######
######
### Phase 2 Sampling - Strategy 3.5: Independent sequential strategy over 4 waves
### Order flipped so rare outcome second
######
######

## Here, 1st and 2nd waves (n/4 samples each) are optimized for BX2 and 
## 3rd and 4th last n/4 are optimized for BX2. We calculate Wave 1 
## allocation outside of the loop because allocation (but not sampling)
## is dependent on the fixed Phase 1 data.

# Run the function for 1,000 iterations, storing the B-hat each time
cl <- makeCluster(detectCores() - 1)
clusterExport(cl, c("n", "simulations_df", "scenario", "Strategy3.5", 
                    "inf_fun_logit", "N"))

run_Strategy3.5 <- function(i) {
  set.seed(simulations_df["data_seed", scenario] + i)
  library(optimall)
  library(survey)
  library(dplyr)
  library(MASS)
  library(mice)
  my_it <<- i
  run <- suppressWarnings(Strategy3.5(n = n, simulations_df = simulations_df,
                                      scenario = scenario, N = N))
  return(c(run[1], run[2], run[3], run[4], run[5], run[6], run[7],
           run[8], run[9], run[10], run[11], run[12], run[13],
           run[14], run[15], run[16], run[17], run[18]))
}

results <- parLapply(cl, 1:2500, run_Strategy3.5)
stopCluster(cl)


# Convert the results into a matrix and extract desired vectors
results_matrix <- do.call(rbind, results)
B_hat_X1_strat3.5_GR <- results_matrix[, 1]
B_hat_X2_strat3.5_GR <- results_matrix[, 2]
B_hat_X1_strat3.5_IPW <- results_matrix[, 3]
B_hat_X2_strat3.5_IPW <- results_matrix[, 4]
SE_hat_X1_strat3.5_GR <- results_matrix[, 5]
SE_hat_X2_strat3.5_GR <- results_matrix[, 6]
SE_hat_X1_strat3.5_IPW <- results_matrix[, 7]
SE_hat_X2_strat3.5_IPW <- results_matrix[, 8]
cover_X1_strat3.5_GR <- results_matrix[, 9]
cover_X2_strat3.5_GR <- results_matrix[, 10]
cover_X1_strat3.5_IPW <- results_matrix[, 11]
cover_X2_strat3.5_IPW <- results_matrix[, 12]
n_oversampled_X1_strat3.5 <- sum(results_matrix[,13])
cor_X1X2_strat3.5 <- median(results_matrix[, 14])
prev_Y_strat3.5 <- median(results_matrix[, 15])
prev_Z_strat3.5 <- median(results_matrix[, 16])
BX1_std_strat3.5 <- results_matrix[,17]
BX2_std_strat3.5 <- results_matrix[,18]

# GR results
var_BX1_strat3.5_GR <- var(B_hat_X1_strat3.5_GR)
var_BX2_strat3.5_GR <- var(B_hat_X2_strat3.5_GR)
mean_BX1_strat3.5_GR <- mean(B_hat_X1_strat3.5_GR)
mean_BX2_strat3.5_GR <- mean(B_hat_X2_strat3.5_GR)
bias_BX1_strat3.5_GR <- mean(B_hat_X1_strat3.5_GR - BX1_std_strat3.5)
bias_BX2_strat3.5_GR <- mean(B_hat_X2_strat3.5_GR - BX2_std_strat3.5)
MSE_BX1_strat3.5_GR <- mean((B_hat_X1_strat3.5_GR - 0.3)^2)
MSE_BX2_strat3.5_GR <- mean((B_hat_X2_strat3.5_GR - 0.7)^2)
median_BX1_strat3.5_GR <- median(B_hat_X1_strat3.5_GR)
median_BX2_strat3.5_GR <- median(B_hat_X2_strat3.5_GR)
median_SE_hat_X1_strat3.5_GR <- median(SE_hat_X1_strat3.5_GR)
median_SE_hat_X2_strat3.5_GR <- median(SE_hat_X2_strat3.5_GR)
coverage_X1_strat3.5_GR <- mean(cover_X1_strat3.5_GR)
coverage_X2_strat3.5_GR <- mean(cover_X2_strat3.5_GR)

# IPW results
var_BX1_strat3.5_IPW <- var(B_hat_X1_strat3.5_IPW)
var_BX2_strat3.5_IPW <- var(B_hat_X2_strat3.5_IPW)
mean_BX1_strat3.5_IPW <- mean(B_hat_X1_strat3.5_IPW)
mean_BX2_strat3.5_IPW <- mean(B_hat_X2_strat3.5_IPW)
bias_BX1_strat3.5_IPW <- mean(B_hat_X1_strat3.5_IPW - BX1_std_strat3.5)
bias_BX2_strat3.5_IPW <- mean(B_hat_X2_strat3.5_IPW - BX2_std_strat3.5)
MSE_BX1_strat3.5_IPW <- mean((B_hat_X1_strat3.5_IPW - 0.3)^2)
MSE_BX2_strat3.5_IPW <- mean((B_hat_X2_strat3.5_IPW - 0.7)^2)
median_BX1_strat3.5_IPW <- median(B_hat_X1_strat3.5_IPW)
median_BX2_strat3.5_IPW <- median(B_hat_X2_strat3.5_IPW)
median_SE_hat_X1_strat3.5_IPW <- median(SE_hat_X1_strat3.5_IPW)
median_SE_hat_X2_strat3.5_IPW <- median(SE_hat_X2_strat3.5_IPW)
coverage_X1_strat3.5_IPW <- mean(cover_X1_strat3.5_IPW)
coverage_X2_strat3.5_IPW <- mean(cover_X2_strat3.5_IPW)

# Save raw results
outpathname4 <- paste0(pathname,"/", names(simulations_df)[scenario], "Strat4raw.csv")
write.csv(results_matrix, outpathname4, row.names = FALSE)

######
######
### Phase 2 Sampling - Strategy 4: Equally weighted A-optimality over 4 waves
######
######


# Run the function for 1,000 iterations,  storing the B-hat each time
cl <- makeCluster(detectCores() - 1)
clusterExport(cl, c("n", "simulations_df", "scenario", "Strategy4", 
                    "a_optimum_allocation",
                    "a_optimal_allocate_wave", "inf_fun_logit", "N"))

run_Strategy4 <- function(i) {
  set.seed(simulations_df["data_seed", scenario] + i)
  library(optimall)
  library(survey)
  library(dplyr)
  library(MASS)
  library(mice)
  my_it <<- i
  run <- suppressWarnings(Strategy4(n = n, simulations_df = simulations_df,
                                    scenario = scenario, N = N))
  return(c(run[1], run[2], run[3], run[4], run[5], run[6], run[7],
           run[8], run[9], run[10], run[11], run[12], run[13],
           run[14], run[15], run[16], run[17], run[18]))
}

results <- parLapply(cl, 1:2500, run_Strategy4)
stopCluster(cl)


# Convert the results into a matrix and extract desired vectors
results_matrix <- do.call(rbind, results)
B_hat_X1_strat4_GR <- results_matrix[, 1]
B_hat_X2_strat4_GR <- results_matrix[, 2]
B_hat_X1_strat4_IPW <- results_matrix[, 3]
B_hat_X2_strat4_IPW <- results_matrix[, 4]
SE_hat_X1_strat4_GR <- results_matrix[, 5]
SE_hat_X2_strat4_GR <- results_matrix[, 6]
SE_hat_X1_strat4_IPW <- results_matrix[, 7]
SE_hat_X2_strat4_IPW <- results_matrix[, 8]
cover_X1_strat4_GR <- results_matrix[, 9]
cover_X2_strat4_GR <- results_matrix[, 10]
cover_X1_strat4_IPW <- results_matrix[, 11]
cover_X2_strat4_IPW <- results_matrix[, 12]
n_oversampled_strat4 <- sum(results_matrix[,13])
cor_X1X2_strat4 <- median(results_matrix[, 14])
prev_Y_strat4 <- median(results_matrix[, 15])
prev_Z_strat4 <- median(results_matrix[, 16])
BX1_std_strat4 <- results_matrix[,17]
BX2_std_strat4 <- results_matrix[,18]

# GR results
var_BX1_strat4_GR <- var(B_hat_X1_strat4_GR)
var_BX2_strat4_GR <- var(B_hat_X2_strat4_GR)
mean_BX1_strat4_GR <- mean(B_hat_X1_strat4_GR)
mean_BX2_strat4_GR <- mean(B_hat_X2_strat4_GR)
bias_BX1_strat4_GR <- mean(B_hat_X1_strat4_GR - BX1_std_strat4)
bias_BX2_strat4_GR <- mean(B_hat_X2_strat4_GR - BX2_std_strat4)
MSE_BX1_strat4_GR <- mean((B_hat_X1_strat4_GR - 0.3)^2)
MSE_BX2_strat4_GR <- mean((B_hat_X2_strat4_GR - 0.7)^2)
median_BX1_strat4_GR <- median(B_hat_X1_strat4_GR)
median_BX2_strat4_GR <- median(B_hat_X2_strat4_GR)
median_SE_hat_X1_strat4_GR <- median(SE_hat_X1_strat4_GR)
median_SE_hat_X2_strat4_GR <- median(SE_hat_X2_strat4_GR)
coverage_X1_strat4_GR <- mean(cover_X1_strat4_GR)
coverage_X2_strat4_GR <- mean(cover_X2_strat4_GR)

# IPW results
var_BX1_strat4_IPW <- var(B_hat_X1_strat4_IPW)
var_BX2_strat4_IPW <- var(B_hat_X2_strat4_IPW)
mean_BX1_strat4_IPW <- mean(B_hat_X1_strat4_IPW)
mean_BX2_strat4_IPW <- mean(B_hat_X2_strat4_IPW)
bias_BX1_strat4_IPW <- mean(B_hat_X1_strat4_IPW - BX1_std_strat4)
bias_BX2_strat4_IPW <- mean(B_hat_X2_strat4_IPW - BX2_std_strat4)
MSE_BX1_strat4_IPW <- mean((B_hat_X1_strat4_IPW - 0.3)^2)
MSE_BX2_strat4_IPW <- mean((B_hat_X2_strat4_IPW - 0.7)^2)
median_BX1_strat4_IPW <- median(B_hat_X1_strat4_IPW)
median_BX2_strat4_IPW <- median(B_hat_X2_strat4_IPW)
median_SE_hat_X1_strat4_IPW <- median(SE_hat_X1_strat4_IPW)
median_SE_hat_X2_strat4_IPW <- median(SE_hat_X2_strat4_IPW)
coverage_X1_strat4_IPW <- mean(cover_X1_strat4_IPW)
coverage_X2_strat4_IPW <- mean(cover_X2_strat4_IPW)

# Save raw results
outpathname5 <- paste0(pathname,"/", names(simulations_df)[scenario], "Strat5raw.csv")
write.csv(results_matrix, outpathname5, row.names = FALSE)

######
######
### Phase 2 Sampling - Strategy 5: True A-optimality
######
######


# Run the function for 1,000 iterations,  storing the B-hat each time
cl <- makeCluster(detectCores() - 1)
clusterExport(cl, c("n", "simulations_df", "scenario", "Strategy5", 
                    "a_optimum_allocation",
                    "a_optimal_allocate_wave", "inf_fun_logit", "N"))

run_Strategy5 <- function(i) {
  set.seed(simulations_df["data_seed", scenario] + i)
  library(optimall)
  library(survey)
  library(dplyr)
  library(MASS)
  library(mice)
  my_it <<- i
  run <- suppressWarnings(Strategy5(n = n, simulations_df = simulations_df,
                                    scenario = scenario, N = N))
  return(c(run[1], run[2], run[3], run[4], run[5], run[6], run[7],
           run[8], run[9], run[10], run[11], run[12], run[13],
           run[14], run[15], run[16], run[17], run[18]))
}

results <- parLapply(cl, 1:2500, run_Strategy5)
stopCluster(cl)


# Convert the results into a matrix and extract desired vectors
results_matrix <- do.call(rbind, results)
B_hat_X1_strat5_GR <- results_matrix[, 1]
B_hat_X2_strat5_GR <- results_matrix[, 2]
B_hat_X1_strat5_IPW <- results_matrix[, 3]
B_hat_X2_strat5_IPW <- results_matrix[, 4]
SE_hat_X1_strat5_GR <- results_matrix[, 5]
SE_hat_X2_strat5_GR <- results_matrix[, 6]
SE_hat_X1_strat5_IPW <- results_matrix[, 7]
SE_hat_X2_strat5_IPW <- results_matrix[, 8]
cover_X1_strat5_GR <- results_matrix[, 9]
cover_X2_strat5_GR <- results_matrix[, 10]
cover_X1_strat5_IPW <- results_matrix[, 11]
cover_X2_strat5_IPW <- results_matrix[, 12]
cor_X1X2_strat5 <- median(results_matrix[, 13])
prev_Y_strat5 <- median(results_matrix[, 14])
prev_Z_strat5 <- median(results_matrix[, 15])
BX1_std_strat5 <- results_matrix[,16]
BX2_std_strat5 <- results_matrix[,17]

# GR results
var_BX1_strat5_GR <- var(B_hat_X1_strat5_GR)
var_BX2_strat5_GR <- var(B_hat_X2_strat5_GR)
mean_BX1_strat5_GR <- mean(B_hat_X1_strat5_GR)
mean_BX2_strat5_GR <- mean(B_hat_X2_strat5_GR)
bias_BX1_strat5_GR <- mean(B_hat_X1_strat5_GR - BX1_std_strat5)
bias_BX2_strat5_GR <- mean(B_hat_X2_strat5_GR - BX2_std_strat5)
MSE_BX1_strat5_GR <- mean((B_hat_X1_strat5_GR - 0.3)^2)
MSE_BX2_strat5_GR <- mean((B_hat_X2_strat5_GR - 0.7)^2)
median_BX1_strat5_GR <- median(B_hat_X1_strat5_GR)
median_BX2_strat5_GR <- median(B_hat_X2_strat5_GR)
median_SE_hat_X1_strat5_GR <- median(SE_hat_X1_strat5_GR)
median_SE_hat_X2_strat5_GR <- median(SE_hat_X2_strat5_GR)
coverage_X1_strat5_GR <- mean(cover_X1_strat5_GR)
coverage_X2_strat5_GR <- mean(cover_X2_strat5_GR)

# IPW results
var_BX1_strat5_IPW <- var(B_hat_X1_strat5_IPW)
var_BX2_strat5_IPW <- var(B_hat_X2_strat5_IPW)
mean_BX1_strat5_IPW <- mean(B_hat_X1_strat5_IPW)
mean_BX2_strat5_IPW <- mean(B_hat_X2_strat5_IPW)
bias_BX1_strat5_IPW <- mean(B_hat_X1_strat5_IPW - BX1_std_strat5)
bias_BX2_strat5_IPW <- mean(B_hat_X2_strat5_IPW - BX2_std_strat5)
MSE_BX1_strat5_IPW <- mean((B_hat_X1_strat5_IPW - 0.3)^2)
MSE_BX2_strat5_IPW <- mean((B_hat_X2_strat5_IPW - 0.7)^2)
median_BX1_strat5_IPW <- median(B_hat_X1_strat5_IPW)
median_BX2_strat5_IPW <- median(B_hat_X2_strat5_IPW)
median_SE_hat_X1_strat5_IPW <- median(SE_hat_X1_strat5_IPW)
median_SE_hat_X2_strat5_IPW <- median(SE_hat_X2_strat5_IPW)
coverage_X1_strat5_IPW <- mean(cover_X1_strat5_IPW)
coverage_X2_strat5_IPW <- mean(cover_X2_strat5_IPW)

# Save raw results
outpathname6 <- paste0(pathname,"/", names(simulations_df)[scenario], "Strat6raw.csv")
write.csv(results_matrix, outpathname6, row.names = FALSE)

#####
## View results
#####

results <- data.frame(Strategy = c(1,2,3,3.5,4,5),
                      "true_BX1" = c(simulations_df["B1", scenario], 0,0,0,0,0),
                      "true_BX2" = c(simulations_df["B2", scenario]*
                                       sqrt(simulations_df["prevX2", scenario]*
                                              (1-simulations_df["prevX2", 
                                                                scenario])), 
                                     0,0,0,0,0),
                      "true_BX1_std" = c(mean(BX1_std_strat1),
                                         mean(BX1_std_strat2),
                                         mean(BX1_std_strat3),
                                         mean(BX1_std_strat3.5),
                                         mean(BX1_std_strat4),
                                         mean(BX1_std_strat5)),
                      "true_BX2_std" = c(mean(BX2_std_strat1),
                                         mean(BX2_std_strat2),
                                         mean(BX2_std_strat3),
                                         mean(BX2_std_strat3.5),
                                         mean(BX2_std_strat4),
                                         mean(BX2_std_strat5)),
                      "cor_X1_X2" = c(cor_X1X2_strat1, cor_X1X2_strat2,
                                      cor_X1X2_strat3,cor_X1X2_strat3.5,
                                      cor_X1X2_strat4, cor_X1X2_strat5),
                      "prev_Y" = c(prev_Y_strat1, prev_Y_strat2,
                                    prev_Y_strat3, prev_Y_strat3.5,
                                    prev_Y_strat4, prev_Y_strat5),
                      "prev_Z" = c(prev_Z_strat1, prev_Z_strat2,
                                    prev_Z_strat3, prev_Z_strat3.5,
                                    prev_Z_strat4, prev_Z_strat5),
                      "mean(BX1)_GR" = c(mean_BX1_strat1_GR,
                                         mean_BX1_strat2_GR,
                                         mean_BX1_strat3_GR,
                                         mean_BX1_strat3.5_GR,
                                         mean_BX1_strat4_GR,
                                         mean_BX1_strat5_GR),
                      "mean(BX2)_GR" = c(mean_BX2_strat1_GR,
                                         mean_BX2_strat2_GR,
                                         mean_BX2_strat3_GR,
                                         mean_BX2_strat3.5_GR,
                                         mean_BX2_strat4_GR,
                                         mean_BX2_strat5_GR),
                      "median(BX1)_GR" = c(median_BX1_strat1_GR,
                                           median_BX1_strat2_GR,
                                           median_BX1_strat3_GR,
                                           median_BX1_strat3.5_GR,
                                           median_BX1_strat4_GR,
                                           median_BX1_strat5_GR),
                      "median(BX2)_GR" = c(median_BX2_strat1_GR,
                                           median_BX2_strat2_GR,
                                           median_BX2_strat3_GR,
                                           median_BX2_strat3.5_GR,
                                           median_BX2_strat4_GR,
                                           median_BX2_strat5_GR),
                      "var(BX1)_GR" = c(var_BX1_strat1_GR,
                                        var_BX1_strat2_GR,
                                        var_BX1_strat3_GR,
                                        var_BX1_strat3.5_GR,
                                        var_BX1_strat4_GR,
                                        var_BX1_strat5_GR),
                      "var(BX2)_GR" = c(var_BX2_strat1_GR,
                                        var_BX2_strat2_GR,
                                        var_BX2_strat3_GR,
                                        var_BX2_strat3.5_GR,
                                        var_BX2_strat4_GR,
                                        var_BX2_strat5_GR),
                      "mean(BX1)_IPW" = c(mean_BX1_strat1_IPW,
                                          mean_BX1_strat2_IPW,
                                          mean_BX1_strat3_IPW,
                                          mean_BX1_strat3.5_IPW,
                                          mean_BX1_strat4_IPW,
                                          mean_BX1_strat5_IPW),
                      "mean(BX2)_IPW" = c(mean_BX2_strat1_IPW,
                                          mean_BX2_strat2_IPW,
                                          mean_BX2_strat3_IPW,
                                          mean_BX2_strat3.5_IPW,
                                          mean_BX2_strat4_IPW,
                                          mean_BX2_strat5_IPW),
                      "median(BX1)_IPW" = c(median_BX1_strat1_IPW,
                                            median_BX1_strat2_IPW,
                                            median_BX1_strat3_IPW,
                                            median_BX1_strat3.5_IPW,
                                            median_BX1_strat4_IPW,
                                            median_BX1_strat5_IPW),
                      "median(BX2)_IPW" = c(median_BX2_strat1_IPW,
                                            median_BX2_strat2_IPW,
                                            median_BX2_strat3_IPW,
                                            median_BX2_strat3.5_IPW,
                                            median_BX2_strat4_IPW,
                                            median_BX2_strat5_IPW),
                      "var(BX1)_IPW" = c(var_BX1_strat1_IPW,
                                         var_BX1_strat2_IPW,
                                         var_BX1_strat3_IPW,
                                         var_BX1_strat3.5_IPW,
                                         var_BX1_strat4_IPW,
                                         var_BX1_strat5_IPW),
                      "var(BX2)_IPW" = c(var_BX2_strat1_IPW,
                                         var_BX2_strat2_IPW,
                                         var_BX2_strat3_IPW,
                                         var_BX2_strat3.5_IPW,
                                         var_BX2_strat4_IPW,
                                         var_BX2_strat5_IPW),
                      "n_oversampled_X1" = c(NA, n_oversampled_X1_strat2,
                                             NA, n_oversampled_X1_strat3.5, NA, NA),
                      "n_oversampled_X2" = c(NA, n_oversampled_X2_strat2,
                                             n_oversampled_X2_strat3, NA, NA, NA),
                      "n_oversampled_overall" = c(NA, NA, NA, NA,
                                                  n_oversampled_strat4, NA),
                      "coverage_X1_IPW" = c(coverage_X1_strat1_IPW,
                                            coverage_X1_strat2_IPW,
                                            coverage_X1_strat3_IPW,
                                            coverage_X1_strat3.5_IPW,
                                            coverage_X1_strat4_IPW,
                                            coverage_X1_strat5_IPW),
                      "coverage_X2_IPW" = c(coverage_X2_strat1_IPW,
                                            coverage_X2_strat2_IPW,
                                            coverage_X2_strat3_IPW,
                                            coverage_X2_strat3.5_IPW,
                                            coverage_X2_strat4_IPW,
                                            coverage_X2_strat5_IPW),
                      "coverage_X1_GR" = c(coverage_X1_strat1_GR,
                                           coverage_X1_strat2_GR,
                                           coverage_X1_strat3_GR,
                                           coverage_X1_strat3.5_GR,
                                           coverage_X1_strat4_GR,
                                           coverage_X1_strat5_GR),
                      "coverage_X2_GR" = c(coverage_X2_strat1_GR,
                                           coverage_X2_strat2_GR,
                                           coverage_X2_strat3_GR,
                                           coverage_X2_strat3.5_GR,
                                           coverage_X2_strat4_GR,
                                           coverage_X2_strat5_GR),
                      "Bias_X1_GR" = c(bias_BX1_strat1_GR,
                                       bias_BX1_strat2_GR,
                                       bias_BX1_strat3_GR,
                                       bias_BX1_strat3.5_GR,
                                       bias_BX1_strat4_GR,
                                       bias_BX1_strat5_GR),
                      "Bias_X2_GR" = c(bias_BX2_strat1_GR,
                                       bias_BX2_strat2_GR,
                                       bias_BX2_strat3_GR,
                                       bias_BX2_strat3.5_GR,
                                       bias_BX2_strat4_GR,
                                       bias_BX2_strat5_GR),
                      "Bias_X1_IPW" = c(bias_BX1_strat1_IPW,
                                        bias_BX1_strat2_IPW,
                                        bias_BX1_strat3_IPW,
                                        bias_BX1_strat3.5_IPW,
                                        bias_BX1_strat4_IPW,
                                        bias_BX1_strat5_IPW),
                      "Bias_X2_IPW" = c(bias_BX2_strat1_IPW,
                                        bias_BX2_strat2_IPW,
                                        bias_BX2_strat3_IPW,
                                        bias_BX2_strat3.5_IPW,
                                        bias_BX2_strat4_IPW,
                                        bias_BX2_strat5_IPW),
                      "median_ASE_X1_IPW" = c(median_SE_hat_X1_strat1_IPW,
                                              median_SE_hat_X1_strat2_IPW,
                                              median_SE_hat_X1_strat3_IPW,
                                              median_SE_hat_X1_strat3.5_IPW,
                                              median_SE_hat_X1_strat4_IPW,
                                              median_SE_hat_X1_strat5_IPW),
                      "median_ASE_X2_IPW" = c(median_SE_hat_X2_strat1_IPW,
                                              median_SE_hat_X2_strat2_IPW,
                                              median_SE_hat_X2_strat3_IPW,
                                              median_SE_hat_X2_strat3.5_IPW,
                                              median_SE_hat_X2_strat4_IPW,
                                              median_SE_hat_X2_strat5_IPW),
                      "median_ASE_X1_GR" = c(median_SE_hat_X1_strat1_GR,
                                             median_SE_hat_X1_strat2_GR,
                                             median_SE_hat_X1_strat3_GR,
                                             median_SE_hat_X1_strat3.5_GR,
                                             median_SE_hat_X1_strat4_GR,
                                             median_SE_hat_X1_strat4_GR),
                      "median_ASE_X2_GR" = c(median_SE_hat_X2_strat1_GR,
                                             median_SE_hat_X2_strat2_GR,
                                             median_SE_hat_X2_strat3_GR,
                                             median_SE_hat_X2_strat3.5_GR,
                                             median_SE_hat_X2_strat4_GR,
                                             median_SE_hat_X2_strat5_GR),
                      "MSE_X1_IPW" = c(MSE_BX1_strat1_IPW,
                                       MSE_BX1_strat2_IPW,
                                       MSE_BX1_strat3_IPW,
                                       MSE_BX1_strat3.5_IPW,
                                       MSE_BX1_strat4_IPW,
                                       MSE_BX1_strat5_IPW),
                      "MSE_X2_IPW" = c(MSE_BX2_strat1_IPW,
                                       MSE_BX2_strat2_IPW,
                                       MSE_BX2_strat3_IPW,
                                       MSE_BX2_strat3.5_IPW,
                                       MSE_BX2_strat4_IPW,
                                       MSE_BX2_strat5_IPW),
                      "MSE_X1_GR" = c(MSE_BX1_strat1_GR,
                                      MSE_BX1_strat2_GR,
                                      MSE_BX1_strat3_GR,
                                      MSE_BX1_strat3.5_GR,
                                      MSE_BX1_strat4_GR,
                                      MSE_BX1_strat5_GR),
                      "MSE_X2_GR" = c(MSE_BX2_strat1_GR,
                                      MSE_BX2_strat2_GR,
                                      MSE_BX2_strat3_GR,
                                      MSE_BX2_strat3.5_GR,
                                      MSE_BX2_strat4_GR,
                                      MSE_BX2_strat5_GR)
)

print(results)

# save
outpathname <- paste0(pathname,"/", names(simulations_df)[scenario], ".csv")
write.csv(results, outpathname, row.names = FALSE)

