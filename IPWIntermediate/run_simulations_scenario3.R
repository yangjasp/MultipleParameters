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

source(paste0(this_path, "/strategy_functions_scenario3.R"))
source(paste0(this_path, "/set_conditions.R"))
source(paste0(this_path, "/utils.R"))

simulations_df <- simulations_df3

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
  #### To run this manually (outside of .bat), set above args manually
  #### according to simulation scenario

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
                    "inf_fun_logit","N"))

run_Strategy1 <- function(i) {
  set.seed(simulations_df["data_seed", scenario] + i)
  library(optimall)
  library(survey)
  library(dplyr)
  library(MASS)
  library(mice)
  run <- Strategy1(n = n, simulations_df = simulations_df, scenario = scenario,
                   N = N)
  return(c(run[1], run[2], run[3], run[4], run[5], run[6], run[7],
           run[8], run[9], run[10], run[11], run[12],run[13],
           run[14], run[15], run[16], run[17], run[18], run[19], run[20],
           run[21], run[22], run[23], run[24], run[25],run[26],
           run[27], run[28], run[29], run[30], run[31], run[32]))
}

results <- parLapply(cl, 1:2500, run_Strategy1)
stopCluster(cl)

# Convert the results into a matrix and extract desired vectors
results_matrix <- do.call(rbind, results)
B_hat_Y1_X1_strat1_GR <- results_matrix[, 1]
B_hat_Y1_X2_strat1_GR <- results_matrix[, 2]
B_hat_Y2_X1_strat1_GR <- results_matrix[, 3]
B_hat_Y2_X2_strat1_GR <- results_matrix[, 4]
B_hat_Y1_X1_strat1_IPW <- results_matrix[, 5]
B_hat_Y1_X2_strat1_IPW <- results_matrix[, 6]
B_hat_Y2_X1_strat1_IPW <- results_matrix[, 7]
B_hat_Y2_X2_strat1_IPW <- results_matrix[, 8]
SE_hat_Y1_X1_strat1_GR <- results_matrix[, 9]
SE_hat_Y1_X2_strat1_GR <- results_matrix[, 10]
SE_hat_Y2_X1_strat1_GR <- results_matrix[, 11]
SE_hat_Y2_X2_strat1_GR <- results_matrix[, 12]
SE_hat_Y1_X1_strat1_IPW <- results_matrix[, 13]
SE_hat_Y1_X2_strat1_IPW <- results_matrix[, 14]
SE_hat_Y2_X1_strat1_IPW <- results_matrix[, 15]
SE_hat_Y2_X2_strat1_IPW <- results_matrix[, 16]
cover_Y1_X1_strat1_GR <- results_matrix[, 17]
cover_Y1_X2_strat1_GR <- results_matrix[, 18]
cover_Y2_X1_strat1_GR <- results_matrix[, 19]
cover_Y2_X2_strat1_GR <- results_matrix[, 20]
cover_Y1_X1_strat1_IPW <- results_matrix[, 21]
cover_Y1_X2_strat1_IPW <- results_matrix[, 22]
cover_Y2_X1_strat1_IPW <- results_matrix[, 23]
cover_Y2_X2_strat1_IPW <- results_matrix[, 24]
cor_Y1Y2_strat1 <- median(results_matrix[, 25])
prev_Y1_strat1 <- median(results_matrix[, 26])
prev_Y2_strat1 <- median(results_matrix[, 27])
B11_std_strat1 <- results_matrix[,28]
B21_std_strat1 <- results_matrix[,29]
B12_std_strat1 <- results_matrix[,30]
B22_std_strat1 <- results_matrix[,31]
cor_X1X2_strat1 <- results_matrix[,32]

# GR results
var_B_hat_Y1_X1_strat1_GR <- var(B_hat_Y1_X1_strat1_GR)
var_B_hat_Y1_X2_strat1_GR <- var(B_hat_Y1_X2_strat1_GR)
var_B_hat_Y2_X1_strat1_GR <- var(B_hat_Y2_X1_strat1_GR)
var_B_hat_Y2_X2_strat1_GR <- var(B_hat_Y2_X2_strat1_GR)
mean_B_hat_Y1_X1_strat1_GR <- mean(B_hat_Y1_X1_strat1_GR)
mean_B_hat_Y1_X2_strat1_GR <- mean(B_hat_Y1_X2_strat1_GR)
mean_B_hat_Y2_X1_strat1_GR <- mean(B_hat_Y2_X1_strat1_GR)
mean_B_hat_Y2_X2_strat1_GR <- mean(B_hat_Y2_X2_strat1_GR)
bias_B_11_strat1_GR <- mean(B_hat_Y1_X1_strat1_GR - B11_std_strat1)
bias_B_21_strat1_GR <- mean(B_hat_Y1_X2_strat1_GR - B21_std_strat1)
bias_B_12_strat1_GR <- mean(B_hat_Y2_X1_strat1_GR - B12_std_strat1)
bias_B_22_strat1_GR <- mean(B_hat_Y2_X2_strat1_GR - B22_std_strat1)
MSE_B_11_strat1_GR <- mean((B_hat_Y1_X1_strat1_GR - B11_std_strat1)^2)
MSE_B_21_strat1_GR <- mean((B_hat_Y1_X2_strat1_GR - B21_std_strat1)^2)
MSE_B_12_strat1_GR <- mean((B_hat_Y2_X1_strat1_GR - B12_std_strat1)^2)
MSE_B_22_strat1_GR <- mean((B_hat_Y2_X2_strat1_GR - B22_std_strat1)^2)
median_B_hat_Y1_X1_strat1_GR <- median(B_hat_Y1_X1_strat1_GR)
median_B_hat_Y1_X2_strat1_GR <- median(B_hat_Y1_X2_strat1_GR)
median_B_hat_Y2_X1_strat1_GR <- median(B_hat_Y2_X1_strat1_GR)
median_B_hat_Y2_X2_strat1_GR <- median(B_hat_Y2_X2_strat1_GR)
median_SE_hat_Y1_X1_strat1_GR <- median(SE_hat_Y1_X1_strat1_GR)
median_SE_hat_Y1_X2_strat1_GR <- median(SE_hat_Y1_X2_strat1_GR)
median_SE_hat_Y2_X1_strat1_GR <- median(SE_hat_Y2_X1_strat1_GR)
median_SE_hat_Y2_X2_strat1_GR <- median(SE_hat_Y2_X2_strat1_GR)
coverage_Y1_X1_strat1_GR <- mean(cover_Y1_X1_strat1_GR)
coverage_Y1_X2_strat1_GR <- mean(cover_Y1_X2_strat1_GR)
coverage_Y2_X1_strat1_GR <- mean(cover_Y2_X1_strat1_GR)
coverage_Y2_X2_strat1_GR <- mean(cover_Y2_X2_strat1_GR)

# IPW results
var_B_hat_Y1_X1_strat1_IPW <- var(B_hat_Y1_X1_strat1_IPW)
var_B_hat_Y1_X2_strat1_IPW <- var(B_hat_Y1_X2_strat1_IPW)
var_B_hat_Y2_X1_strat1_IPW <- var(B_hat_Y2_X1_strat1_IPW)
var_B_hat_Y2_X2_strat1_IPW <- var(B_hat_Y2_X2_strat1_IPW)
mean_B_hat_Y1_X1_strat1_IPW <- mean(B_hat_Y1_X1_strat1_IPW)
mean_B_hat_Y1_X2_strat1_IPW <- mean(B_hat_Y1_X2_strat1_IPW)
mean_B_hat_Y2_X1_strat1_IPW <- mean(B_hat_Y2_X1_strat1_IPW)
mean_B_hat_Y2_X2_strat1_IPW <- mean(B_hat_Y2_X2_strat1_IPW)
bias_B_11_strat1_IPW <- mean(B_hat_Y1_X1_strat1_IPW - B11_std_strat1)
bias_B_21_strat1_IPW <- mean(B_hat_Y1_X2_strat1_IPW - B21_std_strat1)
bias_B_12_strat1_IPW <- mean(B_hat_Y2_X1_strat1_IPW - B12_std_strat1)
bias_B_22_strat1_IPW <- mean(B_hat_Y2_X2_strat1_IPW - B22_std_strat1)
MSE_B_11_strat1_IPW <- mean((B_hat_Y1_X1_strat1_IPW - B11_std_strat1)^2)
MSE_B_21_strat1_IPW <- mean((B_hat_Y1_X2_strat1_IPW - B21_std_strat1)^2)
MSE_B_12_strat1_IPW <- mean((B_hat_Y2_X1_strat1_IPW - B12_std_strat1)^2)
MSE_B_22_strat1_IPW <- mean((B_hat_Y2_X2_strat1_IPW - B22_std_strat1)^2)
median_B_hat_Y1_X1_strat1_IPW <- median(B_hat_Y1_X1_strat1_IPW)
median_B_hat_Y1_X2_strat1_IPW <- median(B_hat_Y1_X2_strat1_IPW)
median_B_hat_Y2_X1_strat1_IPW <- median(B_hat_Y2_X1_strat1_IPW)
median_B_hat_Y2_X2_strat1_IPW <- median(B_hat_Y2_X2_strat1_IPW)
median_SE_hat_Y1_X1_strat1_IPW <- median(SE_hat_Y1_X1_strat1_IPW)
median_SE_hat_Y1_X2_strat1_IPW <- median(SE_hat_Y1_X2_strat1_IPW)
median_SE_hat_Y2_X1_strat1_IPW <- median(SE_hat_Y2_X1_strat1_IPW)
median_SE_hat_Y2_X2_strat1_IPW <- median(SE_hat_Y2_X2_strat1_IPW)
coverage_Y1_X1_strat1_IPW <- mean(cover_Y1_X1_strat1_IPW)
coverage_Y1_X2_strat1_IPW <- mean(cover_Y1_X2_strat1_IPW)
coverage_Y2_X1_strat1_IPW <- mean(cover_Y2_X1_strat1_IPW)
coverage_Y2_X2_strat1_IPW <- mean(cover_Y2_X2_strat1_IPW)

outpathname1 <- paste0(pathname,"/", names(simulations_df)[scenario], "Strat1raw.csv")
write.csv(results_matrix, outpathname1, row.names = FALSE)

######
### Phase 2 Sampling - Strategy 2: Independent, uniform strategy over 4 waves
######

## Here, each wave chooses n/4 samples where n/8 are chosen to be optimal for
## B_11 and n/8 to be optimal for B_12. We calculate Wave 1 
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
           run[14], run[15], run[16], run[17], run[18], run[19],
           run[20], run[21], run[22], run[23],
           run[24], run[25], run[26], run[27], run[28], run[29],
           run[30], run[31], run[32], run[33]))
}

results <- parLapply(cl, 1:2500, run_Strategy2)
stopCluster(cl)

# Convert the results into a matrix and extract desired vectors
results_matrix <- do.call(rbind, results)
B_hat_Y1_X1_strat2_GR <- results_matrix[, 1]
B_hat_Y1_X2_strat2_GR <- results_matrix[, 2]
B_hat_Y2_X1_strat2_GR <- results_matrix[, 3]
B_hat_Y2_X2_strat2_GR <- results_matrix[, 4]
B_hat_Y1_X1_strat2_IPW <- results_matrix[, 5]
B_hat_Y1_X2_strat2_IPW <- results_matrix[, 6]
B_hat_Y2_X1_strat2_IPW <- results_matrix[, 7]
B_hat_Y2_X2_strat2_IPW <- results_matrix[, 8]
SE_hat_Y1_X1_strat2_GR <- results_matrix[, 9]
SE_hat_Y1_X2_strat2_GR <- results_matrix[, 10]
SE_hat_Y2_X1_strat2_GR <- results_matrix[, 11]
SE_hat_Y2_X2_strat2_GR <- results_matrix[, 12]
SE_hat_Y1_X1_strat2_IPW <- results_matrix[, 13]
SE_hat_Y1_X2_strat2_IPW <- results_matrix[, 14]
SE_hat_Y2_X1_strat2_IPW <- results_matrix[, 15]
SE_hat_Y2_X2_strat2_IPW <- results_matrix[, 16]
cover_Y1_X1_strat2_GR <- results_matrix[, 17]
cover_Y1_X2_strat2_GR <- results_matrix[, 18]
cover_Y2_X1_strat2_GR <- results_matrix[, 19]
cover_Y2_X2_strat2_GR <- results_matrix[, 20]
cover_Y1_X1_strat2_IPW <- results_matrix[, 21]
cover_Y1_X2_strat2_IPW <- results_matrix[, 22]
cover_Y2_X1_strat2_IPW <- results_matrix[, 23]
cover_Y2_X2_strat2_IPW <- results_matrix[, 24]
cor_Y1Y2_strat2 <- median(results_matrix[, 27])
prev_Y1_strat2 <- median(results_matrix[, 28])
prev_Y2_strat2 <- median(results_matrix[, 29])
B11_std_strat2 <- results_matrix[,30]
B21_std_strat2 <- results_matrix[,31]
B12_std_strat2 <- results_matrix[,32]
B22_std_strat2 <- results_matrix[,33]

# GR results
var_B_hat_Y1_X1_strat2_GR <- var(B_hat_Y1_X1_strat2_GR)
var_B_hat_Y1_X2_strat2_GR <- var(B_hat_Y1_X2_strat2_GR)
var_B_hat_Y2_X1_strat2_GR <- var(B_hat_Y2_X1_strat2_GR)
var_B_hat_Y2_X2_strat2_GR <- var(B_hat_Y2_X2_strat2_GR)
mean_B_hat_Y1_X1_strat2_GR <- mean(B_hat_Y1_X1_strat2_GR)
mean_B_hat_Y1_X2_strat2_GR <- mean(B_hat_Y1_X2_strat2_GR)
mean_B_hat_Y2_X1_strat2_GR <- mean(B_hat_Y2_X1_strat2_GR)
mean_B_hat_Y2_X2_strat2_GR <- mean(B_hat_Y2_X2_strat2_GR)
bias_B_11_strat2_GR <- mean(B_hat_Y1_X1_strat2_GR - B11_std_strat2)
bias_B_21_strat2_GR <- mean(B_hat_Y1_X2_strat2_GR - B21_std_strat2)
bias_B_12_strat2_GR <- mean(B_hat_Y2_X1_strat2_GR - B12_std_strat2)
bias_B_22_strat2_GR <- mean(B_hat_Y2_X2_strat2_GR - B22_std_strat2)
MSE_B_11_strat2_GR <- mean((B_hat_Y1_X1_strat2_GR - B11_std_strat1)^2)
MSE_B_21_strat2_GR <- mean((B_hat_Y1_X2_strat2_GR - B21_std_strat1)^2)
MSE_B_12_strat2_GR <- mean((B_hat_Y2_X1_strat2_GR - B12_std_strat1)^2)
MSE_B_22_strat2_GR <- mean((B_hat_Y2_X2_strat2_GR - B22_std_strat1)^2)
median_B_hat_Y1_X1_strat2_GR <- median(B_hat_Y1_X1_strat2_GR)
median_B_hat_Y1_X2_strat2_GR <- median(B_hat_Y1_X2_strat2_GR)
median_B_hat_Y2_X1_strat2_GR <- median(B_hat_Y2_X1_strat2_GR)
median_B_hat_Y2_X2_strat2_GR <- median(B_hat_Y2_X2_strat2_GR)
median_SE_hat_Y1_X1_strat2_GR <- median(SE_hat_Y1_X1_strat2_GR)
median_SE_hat_Y1_X2_strat2_GR <- median(SE_hat_Y1_X2_strat2_GR)
median_SE_hat_Y2_X1_strat2_GR <- median(SE_hat_Y2_X1_strat2_GR)
median_SE_hat_Y2_X2_strat2_GR <- median(SE_hat_Y2_X2_strat2_GR)
coverage_Y1_X1_strat2_GR <- mean(cover_Y1_X1_strat2_GR)
coverage_Y1_X2_strat2_GR <- mean(cover_Y1_X2_strat2_GR)
coverage_Y2_X1_strat2_GR <- mean(cover_Y2_X1_strat2_GR)
coverage_Y2_X2_strat2_GR <- mean(cover_Y2_X2_strat2_GR)

# IPW results
var_B_hat_Y1_X1_strat2_IPW <- var(B_hat_Y1_X1_strat2_IPW)
var_B_hat_Y1_X2_strat2_IPW <- var(B_hat_Y1_X2_strat2_IPW)
var_B_hat_Y2_X1_strat2_IPW <- var(B_hat_Y2_X1_strat2_IPW)
var_B_hat_Y2_X2_strat2_IPW <- var(B_hat_Y2_X2_strat2_IPW)
mean_B_hat_Y1_X1_strat2_IPW <- mean(B_hat_Y1_X1_strat2_IPW)
mean_B_hat_Y1_X2_strat2_IPW <- mean(B_hat_Y1_X2_strat2_IPW)
mean_B_hat_Y2_X1_strat2_IPW <- mean(B_hat_Y2_X1_strat2_IPW)
mean_B_hat_Y2_X2_strat2_IPW <- mean(B_hat_Y2_X2_strat2_IPW)
bias_B_11_strat2_IPW <- mean(B_hat_Y1_X1_strat2_IPW - B11_std_strat2)
bias_B_21_strat2_IPW <- mean(B_hat_Y1_X2_strat2_IPW - B21_std_strat2)
bias_B_12_strat2_IPW <- mean(B_hat_Y2_X1_strat2_IPW - B12_std_strat2)
bias_B_22_strat2_IPW <- mean(B_hat_Y2_X2_strat2_IPW - B22_std_strat2)
MSE_B_11_strat2_IPW <- mean((B_hat_Y1_X1_strat2_IPW - B11_std_strat1)^2)
MSE_B_21_strat2_IPW <- mean((B_hat_Y1_X2_strat2_IPW - B21_std_strat1)^2)
MSE_B_12_strat2_IPW <- mean((B_hat_Y2_X1_strat2_IPW - B12_std_strat1)^2)
MSE_B_22_strat2_IPW <- mean((B_hat_Y2_X2_strat2_IPW - B22_std_strat1)^2)
median_B_hat_Y1_X1_strat2_IPW <- median(B_hat_Y1_X1_strat2_IPW)
median_B_hat_Y1_X2_strat2_IPW <- median(B_hat_Y1_X2_strat2_IPW)
median_B_hat_Y2_X1_strat2_IPW <- median(B_hat_Y2_X1_strat2_IPW)
median_B_hat_Y2_X2_strat2_IPW <- median(B_hat_Y2_X2_strat2_IPW)
median_SE_hat_Y1_X1_strat2_IPW <- median(SE_hat_Y1_X1_strat2_IPW)
median_SE_hat_Y1_X2_strat2_IPW <- median(SE_hat_Y1_X2_strat2_IPW)
median_SE_hat_Y2_X1_strat2_IPW <- median(SE_hat_Y2_X1_strat2_IPW)
median_SE_hat_Y2_X2_strat2_IPW <- median(SE_hat_Y2_X2_strat2_IPW)
coverage_Y1_X1_strat2_IPW <- mean(cover_Y1_X1_strat2_IPW)
coverage_Y1_X2_strat2_IPW <- mean(cover_Y1_X2_strat2_IPW)
coverage_Y2_X1_strat2_IPW <- mean(cover_Y2_X1_strat2_IPW)
coverage_Y2_X2_strat2_IPW <- mean(cover_Y2_X2_strat2_IPW)

# Print raw results
outpathname2 <- paste0(pathname,"/", names(simulations_df)[scenario], "Strat2raw.csv")
write.csv(results_matrix, outpathname2, row.names = FALSE)

######
######
### Phase 2 Sampling - Strategy 3: Independent sequential strategy over 4 waves
######
######

## Here, 1st and 2nd waves (n/4 samples each) are optimized for B_11 and 
## 3rd and 4th last n/4 are optimized for B_12. We calculate Wave 1 
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
           run[14], run[15], run[16], run[17], run[18],run[19],
           run[20], run[21], run[22], run[23],
           run[24], run[25], run[26], run[27], run[28], run[29],
           run[30], run[31], run[32]))
}

results <- parLapply(cl, 1:2500, run_Strategy3)
stopCluster(cl)


# Convert the results into a matrix and extract desired vectors
results_matrix <- do.call(rbind, results)
B_hat_Y1_X1_strat3_GR <- results_matrix[, 1]
B_hat_Y1_X2_strat3_GR <- results_matrix[, 2]
B_hat_Y2_X1_strat3_GR <- results_matrix[, 3]
B_hat_Y2_X2_strat3_GR <- results_matrix[, 4]
B_hat_Y1_X1_strat3_IPW <- results_matrix[, 5]
B_hat_Y1_X2_strat3_IPW <- results_matrix[, 6]
B_hat_Y2_X1_strat3_IPW <- results_matrix[, 7]
B_hat_Y2_X2_strat3_IPW <- results_matrix[, 8]
SE_hat_Y1_X1_strat3_GR <- results_matrix[, 9]
SE_hat_Y1_X2_strat3_GR <- results_matrix[, 10]
SE_hat_Y2_X1_strat3_GR <- results_matrix[, 11]
SE_hat_Y2_X2_strat3_GR <- results_matrix[, 12]
SE_hat_Y1_X1_strat3_IPW <- results_matrix[, 13]
SE_hat_Y1_X2_strat3_IPW <- results_matrix[, 14]
SE_hat_Y2_X1_strat3_IPW <- results_matrix[, 15]
SE_hat_Y2_X2_strat3_IPW <- results_matrix[, 16]
cover_Y1_X1_strat3_GR <- results_matrix[, 17]
cover_Y1_X2_strat3_GR <- results_matrix[, 18]
cover_Y2_X1_strat3_GR <- results_matrix[, 19]
cover_Y2_X2_strat3_GR <- results_matrix[, 20]
cover_Y1_X1_strat3_IPW <- results_matrix[, 21]
cover_Y1_X2_strat3_IPW <- results_matrix[, 22]
cover_Y2_X1_strat3_IPW <- results_matrix[, 23]
cover_Y2_X2_strat3_IPW <- results_matrix[, 24]
cor_Y1Y2_strat3 <- median(results_matrix[, 26])
prev_Y1_strat3 <- median(results_matrix[, 27])
prev_Y2_strat3 <- median(results_matrix[, 28])
B11_std_strat3 <- results_matrix[,29]
B21_std_strat3 <- results_matrix[,30]
B12_std_strat3 <- results_matrix[,31]
B22_std_strat3 <- results_matrix[,32]

# GR results
var_B_hat_Y1_X1_strat3_GR <- var(B_hat_Y1_X1_strat3_GR)
var_B_hat_Y1_X2_strat3_GR <- var(B_hat_Y1_X2_strat3_GR)
var_B_hat_Y2_X1_strat3_GR <- var(B_hat_Y2_X1_strat3_GR)
var_B_hat_Y2_X2_strat3_GR <- var(B_hat_Y2_X2_strat3_GR)
mean_B_hat_Y1_X1_strat3_GR <- mean(B_hat_Y1_X1_strat3_GR)
mean_B_hat_Y1_X2_strat3_GR <- mean(B_hat_Y1_X2_strat3_GR)
mean_B_hat_Y2_X1_strat3_GR <- mean(B_hat_Y2_X1_strat3_GR)
mean_B_hat_Y2_X2_strat3_GR <- mean(B_hat_Y2_X2_strat3_GR)
bias_B_11_strat3_GR <- mean(B_hat_Y1_X1_strat3_GR - B11_std_strat3)
bias_B_21_strat3_GR <- mean(B_hat_Y1_X2_strat3_GR - B21_std_strat3)
bias_B_12_strat3_GR <- mean(B_hat_Y2_X1_strat3_GR - B12_std_strat3)
bias_B_22_strat3_GR <- mean(B_hat_Y2_X2_strat3_GR - B22_std_strat3)
MSE_B_11_strat3_GR <- mean((B_hat_Y1_X1_strat3_GR - B11_std_strat1)^2)
MSE_B_21_strat3_GR <- mean((B_hat_Y1_X2_strat3_GR - B21_std_strat1)^2)
MSE_B_12_strat3_GR <- mean((B_hat_Y2_X1_strat3_GR - B12_std_strat1)^2)
MSE_B_22_strat3_GR <- mean((B_hat_Y2_X2_strat3_GR - B22_std_strat1)^2)
median_B_hat_Y1_X1_strat3_GR <- median(B_hat_Y1_X1_strat3_GR)
median_B_hat_Y1_X2_strat3_GR <- median(B_hat_Y1_X2_strat3_GR)
median_B_hat_Y2_X1_strat3_GR <- median(B_hat_Y2_X1_strat3_GR)
median_B_hat_Y2_X2_strat3_GR <- median(B_hat_Y2_X2_strat3_GR)
median_SE_hat_Y1_X1_strat3_GR <- median(SE_hat_Y1_X1_strat3_GR)
median_SE_hat_Y1_X2_strat3_GR <- median(SE_hat_Y1_X2_strat3_GR)
median_SE_hat_Y2_X1_strat3_GR <- median(SE_hat_Y2_X1_strat3_GR)
median_SE_hat_Y2_X2_strat3_GR <- median(SE_hat_Y2_X2_strat3_GR)
coverage_Y1_X1_strat3_GR <- mean(cover_Y1_X1_strat3_GR)
coverage_Y1_X2_strat3_GR <- mean(cover_Y1_X2_strat3_GR)
coverage_Y2_X1_strat3_GR <- mean(cover_Y2_X1_strat3_GR)
coverage_Y2_X2_strat3_GR <- mean(cover_Y2_X2_strat3_GR)

# IPW results
var_B_hat_Y1_X1_strat3_IPW <- var(B_hat_Y1_X1_strat3_IPW)
var_B_hat_Y1_X2_strat3_IPW <- var(B_hat_Y1_X2_strat3_IPW)
var_B_hat_Y2_X1_strat3_IPW <- var(B_hat_Y2_X1_strat3_IPW)
var_B_hat_Y2_X2_strat3_IPW <- var(B_hat_Y2_X2_strat3_IPW)
mean_B_hat_Y1_X1_strat3_IPW <- mean(B_hat_Y1_X1_strat3_IPW)
mean_B_hat_Y1_X2_strat3_IPW <- mean(B_hat_Y1_X2_strat3_IPW)
mean_B_hat_Y2_X1_strat3_IPW <- mean(B_hat_Y2_X1_strat3_IPW)
mean_B_hat_Y2_X2_strat3_IPW <- mean(B_hat_Y2_X2_strat3_IPW)
bias_B_11_strat3_IPW <- mean(B_hat_Y1_X1_strat3_IPW - B11_std_strat3)
bias_B_21_strat3_IPW <- mean(B_hat_Y1_X2_strat3_IPW - B21_std_strat3)
bias_B_12_strat3_IPW <- mean(B_hat_Y2_X1_strat3_IPW - B12_std_strat3)
bias_B_22_strat3_IPW <- mean(B_hat_Y2_X2_strat3_IPW - B22_std_strat3)
MSE_B_11_strat3_IPW <- mean((B_hat_Y1_X1_strat3_IPW - B11_std_strat1)^2)
MSE_B_21_strat3_IPW <- mean((B_hat_Y1_X2_strat3_IPW - B21_std_strat1)^2)
MSE_B_12_strat3_IPW <- mean((B_hat_Y2_X1_strat3_IPW - B12_std_strat1)^2)
MSE_B_22_strat3_IPW <- mean((B_hat_Y2_X2_strat3_IPW - B22_std_strat1)^2)
median_B_hat_Y1_X1_strat3_IPW <- median(B_hat_Y1_X1_strat3_IPW)
median_B_hat_Y1_X2_strat3_IPW <- median(B_hat_Y1_X2_strat3_IPW)
median_B_hat_Y2_X1_strat3_IPW <- median(B_hat_Y2_X1_strat3_IPW)
median_B_hat_Y2_X2_strat3_IPW <- median(B_hat_Y2_X2_strat3_IPW)
median_SE_hat_Y1_X1_strat3_IPW <- median(SE_hat_Y1_X1_strat3_IPW)
median_SE_hat_Y1_X2_strat3_IPW <- median(SE_hat_Y1_X2_strat3_IPW)
median_SE_hat_Y2_X1_strat3_IPW <- median(SE_hat_Y2_X1_strat3_IPW)
median_SE_hat_Y2_X2_strat3_IPW <- median(SE_hat_Y2_X2_strat3_IPW)
coverage_Y1_X1_strat3_IPW <- mean(cover_Y1_X1_strat3_IPW)
coverage_Y1_X2_strat3_IPW <- mean(cover_Y1_X2_strat3_IPW)
coverage_Y2_X1_strat3_IPW <- mean(cover_Y2_X1_strat3_IPW)
coverage_Y2_X2_strat3_IPW <- mean(cover_Y2_X2_strat3_IPW)

#Save raw results
outpathname3 <- paste0(pathname,"/", names(simulations_df)[scenario], "Strat3raw.csv")
write.csv(results_matrix, outpathname3, row.names = FALSE)

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
           run[14], run[15], run[16], run[17], run[18],run[19],
           run[20], run[21], run[22], run[23],
           run[24], run[25], run[26], run[27], run[28], run[29],
           run[30], run[31], run[32]))
}

results <- parLapply(cl, 1:2500, run_Strategy4)
stopCluster(cl)


# Convert the results into a matrix and extract desired vectors
results_matrix <- do.call(rbind, results)
B_hat_Y1_X1_strat4_GR <- results_matrix[, 1]
B_hat_Y1_X2_strat4_GR <- results_matrix[, 2]
B_hat_Y2_X1_strat4_GR <- results_matrix[, 3]
B_hat_Y2_X2_strat4_GR <- results_matrix[, 4]
B_hat_Y1_X1_strat4_IPW <- results_matrix[, 5]
B_hat_Y1_X2_strat4_IPW <- results_matrix[, 6]
B_hat_Y2_X1_strat4_IPW <- results_matrix[, 7]
B_hat_Y2_X2_strat4_IPW <- results_matrix[, 8]
SE_hat_Y1_X1_strat4_GR <- results_matrix[, 9]
SE_hat_Y1_X2_strat4_GR <- results_matrix[, 10]
SE_hat_Y2_X1_strat4_GR <- results_matrix[, 11]
SE_hat_Y2_X2_strat4_GR <- results_matrix[, 12]
SE_hat_Y1_X1_strat4_IPW <- results_matrix[, 13]
SE_hat_Y1_X2_strat4_IPW <- results_matrix[, 14]
SE_hat_Y2_X1_strat4_IPW <- results_matrix[, 15]
SE_hat_Y2_X2_strat4_IPW <- results_matrix[, 16]
cover_Y1_X1_strat4_GR <- results_matrix[, 17]
cover_Y1_X2_strat4_GR <- results_matrix[, 18]
cover_Y2_X1_strat4_GR <- results_matrix[, 19]
cover_Y2_X2_strat4_GR <- results_matrix[, 20]
cover_Y1_X1_strat4_IPW <- results_matrix[, 21]
cover_Y1_X2_strat4_IPW <- results_matrix[, 22]
cover_Y2_X1_strat4_IPW <- results_matrix[, 23]
cover_Y2_X2_strat4_IPW <- results_matrix[, 24]
cor_Y1Y2_strat4 <- median(results_matrix[, 26])
prev_Y1_strat4 <- median(results_matrix[, 27])
prev_Y2_strat4 <- median(results_matrix[, 28])
B11_std_strat4 <- results_matrix[,29]
B21_std_strat4 <- results_matrix[,30]
B12_std_strat4 <- results_matrix[,31]
B22_std_strat4 <- results_matrix[,32]

# GR results
var_B_hat_Y1_X1_strat4_GR <- var(B_hat_Y1_X1_strat4_GR)
var_B_hat_Y1_X2_strat4_GR <- var(B_hat_Y1_X2_strat4_GR)
var_B_hat_Y2_X1_strat4_GR <- var(B_hat_Y2_X1_strat4_GR)
var_B_hat_Y2_X2_strat4_GR <- var(B_hat_Y2_X2_strat4_GR)
mean_B_hat_Y1_X1_strat4_GR <- mean(B_hat_Y1_X1_strat4_GR)
mean_B_hat_Y1_X2_strat4_GR <- mean(B_hat_Y1_X2_strat4_GR)
mean_B_hat_Y2_X1_strat4_GR <- mean(B_hat_Y2_X1_strat4_GR)
mean_B_hat_Y2_X2_strat4_GR <- mean(B_hat_Y2_X2_strat4_GR)
bias_B_11_strat4_GR <- mean(B_hat_Y1_X1_strat4_GR - B11_std_strat4)
bias_B_21_strat4_GR <- mean(B_hat_Y1_X2_strat4_GR - B21_std_strat4)
bias_B_12_strat4_GR <- mean(B_hat_Y2_X1_strat4_GR - B12_std_strat4)
bias_B_22_strat4_GR <- mean(B_hat_Y2_X2_strat4_GR - B22_std_strat4)
MSE_B_11_strat4_GR <- mean((B_hat_Y1_X1_strat4_GR - B11_std_strat1)^2)
MSE_B_21_strat4_GR <- mean((B_hat_Y1_X2_strat4_GR - B21_std_strat1)^2)
MSE_B_12_strat4_GR <- mean((B_hat_Y2_X1_strat4_GR - B12_std_strat1)^2)
MSE_B_22_strat4_GR <- mean((B_hat_Y2_X2_strat4_GR - B22_std_strat1)^2)
median_B_hat_Y1_X1_strat4_GR <- median(B_hat_Y1_X1_strat4_GR)
median_B_hat_Y1_X2_strat4_GR <- median(B_hat_Y1_X2_strat4_GR)
median_B_hat_Y2_X1_strat4_GR <- median(B_hat_Y2_X1_strat4_GR)
median_B_hat_Y2_X2_strat4_GR <- median(B_hat_Y2_X2_strat4_GR)
median_SE_hat_Y1_X1_strat4_GR <- median(SE_hat_Y1_X1_strat4_GR)
median_SE_hat_Y1_X2_strat4_GR <- median(SE_hat_Y1_X2_strat4_GR)
median_SE_hat_Y2_X1_strat4_GR <- median(SE_hat_Y2_X1_strat4_GR)
median_SE_hat_Y2_X2_strat4_GR <- median(SE_hat_Y2_X2_strat4_GR)
coverage_Y1_X1_strat4_GR <- mean(cover_Y1_X1_strat4_GR)
coverage_Y1_X2_strat4_GR <- mean(cover_Y1_X2_strat4_GR)
coverage_Y2_X1_strat4_GR <- mean(cover_Y2_X1_strat4_GR)
coverage_Y2_X2_strat4_GR <- mean(cover_Y2_X2_strat4_GR)

# IPW results
var_B_hat_Y1_X1_strat4_IPW <- var(B_hat_Y1_X1_strat4_IPW)
var_B_hat_Y1_X2_strat4_IPW <- var(B_hat_Y1_X2_strat4_IPW)
var_B_hat_Y2_X1_strat4_IPW <- var(B_hat_Y2_X1_strat4_IPW)
var_B_hat_Y2_X2_strat4_IPW <- var(B_hat_Y2_X2_strat4_IPW)
mean_B_hat_Y1_X1_strat4_IPW <- mean(B_hat_Y1_X1_strat4_IPW)
mean_B_hat_Y1_X2_strat4_IPW <- mean(B_hat_Y1_X2_strat4_IPW)
mean_B_hat_Y2_X1_strat4_IPW <- mean(B_hat_Y2_X1_strat4_IPW)
mean_B_hat_Y2_X2_strat4_IPW <- mean(B_hat_Y2_X2_strat4_IPW)
bias_B_11_strat4_IPW <- mean(B_hat_Y1_X1_strat4_IPW - B11_std_strat4)
bias_B_21_strat4_IPW <- mean(B_hat_Y1_X2_strat4_IPW - B21_std_strat4)
bias_B_12_strat4_IPW <- mean(B_hat_Y2_X1_strat4_IPW - B12_std_strat4)
bias_B_22_strat4_IPW <- mean(B_hat_Y2_X2_strat4_IPW - B22_std_strat4)
MSE_B_11_strat4_IPW <- mean((B_hat_Y1_X1_strat4_IPW - B11_std_strat1)^2)
MSE_B_21_strat4_IPW <- mean((B_hat_Y1_X2_strat4_IPW - B21_std_strat1)^2)
MSE_B_12_strat4_IPW <- mean((B_hat_Y2_X1_strat4_IPW - B12_std_strat1)^2)
MSE_B_22_strat4_IPW <- mean((B_hat_Y2_X2_strat4_IPW - B22_std_strat1)^2)
median_B_hat_Y1_X1_strat4_IPW <- median(B_hat_Y1_X1_strat4_IPW)
median_B_hat_Y1_X2_strat4_IPW <- median(B_hat_Y1_X2_strat4_IPW)
median_B_hat_Y2_X1_strat4_IPW <- median(B_hat_Y2_X1_strat4_IPW)
median_B_hat_Y2_X2_strat4_IPW <- median(B_hat_Y2_X2_strat4_IPW)
median_SE_hat_Y1_X1_strat4_IPW <- median(SE_hat_Y1_X1_strat4_IPW)
median_SE_hat_Y1_X2_strat4_IPW <- median(SE_hat_Y1_X2_strat4_IPW)
median_SE_hat_Y2_X1_strat4_IPW <- median(SE_hat_Y2_X1_strat4_IPW)
median_SE_hat_Y2_X2_strat4_IPW <- median(SE_hat_Y2_X2_strat4_IPW)
coverage_Y1_X1_strat4_IPW <- mean(cover_Y1_X1_strat4_IPW)
coverage_Y1_X2_strat4_IPW <- mean(cover_Y1_X2_strat4_IPW)
coverage_Y2_X1_strat4_IPW <- mean(cover_Y2_X1_strat4_IPW)
coverage_Y2_X2_strat4_IPW <- mean(cover_Y2_X2_strat4_IPW)

# Save raw results
outpathname5 <- paste0(pathname,"/", names(simulations_df)[scenario], "Strat4raw.csv")
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
           run[14], run[15], run[16], run[17], run[18],run[19],
           run[20], run[21], run[22], run[23],
           run[24], run[25], run[26], run[27], run[28], run[29],
           run[30], run[31], run[32]))
}

results <- parLapply(cl, 1:2500, run_Strategy5)
stopCluster(cl)


# Convert the results into a matrix and extract desired vectors
results_matrix <- do.call(rbind, results)
B_hat_Y1_X1_strat5_GR <- results_matrix[, 1]
B_hat_Y1_X2_strat5_GR <- results_matrix[, 2]
B_hat_Y2_X1_strat5_GR <- results_matrix[, 3]
B_hat_Y2_X2_strat5_GR <- results_matrix[, 4]
B_hat_Y1_X1_strat5_IPW <- results_matrix[, 5]
B_hat_Y1_X2_strat5_IPW <- results_matrix[, 6]
B_hat_Y2_X1_strat5_IPW <- results_matrix[, 7]
B_hat_Y2_X2_strat5_IPW <- results_matrix[, 8]
SE_hat_Y1_X1_strat5_GR <- results_matrix[, 9]
SE_hat_Y1_X2_strat5_GR <- results_matrix[, 10]
SE_hat_Y2_X1_strat5_GR <- results_matrix[, 11]
SE_hat_Y2_X2_strat5_GR <- results_matrix[, 12]
SE_hat_Y1_X1_strat5_IPW <- results_matrix[, 13]
SE_hat_Y1_X2_strat5_IPW <- results_matrix[, 14]
SE_hat_Y2_X1_strat5_IPW <- results_matrix[, 15]
SE_hat_Y2_X2_strat5_IPW <- results_matrix[, 16]
cover_Y1_X1_strat5_GR <- results_matrix[, 17]
cover_Y1_X2_strat5_GR <- results_matrix[, 18]
cover_Y2_X1_strat5_GR <- results_matrix[, 19]
cover_Y2_X2_strat5_GR <- results_matrix[, 20]
cover_Y1_X1_strat5_IPW <- results_matrix[, 21]
cover_Y1_X2_strat5_IPW <- results_matrix[, 22]
cover_Y2_X1_strat5_IPW <- results_matrix[, 23]
cover_Y2_X2_strat5_IPW <- results_matrix[, 24]
cor_Y1Y2_strat5 <- median(results_matrix[, 25])
prev_Y1_strat5 <- median(results_matrix[, 26])
prev_Y2_strat5 <- median(results_matrix[, 27])
B11_std_strat5 <- results_matrix[,28]
B21_std_strat5 <- results_matrix[,29]
B12_std_strat5 <- results_matrix[,30]
B22_std_strat5 <- results_matrix[,31]

# GR results
var_B_hat_Y1_X1_strat5_GR <- var(B_hat_Y1_X1_strat5_GR)
var_B_hat_Y1_X2_strat5_GR <- var(B_hat_Y1_X2_strat5_GR)
var_B_hat_Y2_X1_strat5_GR <- var(B_hat_Y2_X1_strat5_GR)
var_B_hat_Y2_X2_strat5_GR <- var(B_hat_Y2_X2_strat5_GR)
mean_B_hat_Y1_X1_strat5_GR <- mean(B_hat_Y1_X1_strat5_GR)
mean_B_hat_Y1_X2_strat5_GR <- mean(B_hat_Y1_X2_strat5_GR)
mean_B_hat_Y2_X1_strat5_GR <- mean(B_hat_Y2_X1_strat5_GR)
mean_B_hat_Y2_X2_strat5_GR <- mean(B_hat_Y2_X2_strat5_GR)
bias_B_11_strat5_GR <- mean(B_hat_Y1_X1_strat5_GR - B11_std_strat5)
bias_B_21_strat5_GR <- mean(B_hat_Y1_X2_strat5_GR - B21_std_strat5)
bias_B_12_strat5_GR <- mean(B_hat_Y2_X1_strat5_GR - B12_std_strat5)
bias_B_22_strat5_GR <- mean(B_hat_Y2_X2_strat5_GR - B22_std_strat5)
MSE_B_11_strat5_GR <- mean((B_hat_Y1_X1_strat5_GR - B11_std_strat1)^2)
MSE_B_21_strat5_GR <- mean((B_hat_Y1_X2_strat5_GR - B21_std_strat1)^2)
MSE_B_12_strat5_GR <- mean((B_hat_Y2_X1_strat5_GR - B12_std_strat1)^2)
MSE_B_22_strat5_GR <- mean((B_hat_Y2_X2_strat5_GR - B22_std_strat1)^2)
median_B_hat_Y1_X1_strat5_GR <- median(B_hat_Y1_X1_strat5_GR)
median_B_hat_Y1_X2_strat5_GR <- median(B_hat_Y1_X2_strat5_GR)
median_B_hat_Y2_X1_strat5_GR <- median(B_hat_Y2_X1_strat5_GR)
median_B_hat_Y2_X2_strat5_GR <- median(B_hat_Y2_X2_strat5_GR)
median_SE_hat_Y1_X1_strat5_GR <- median(SE_hat_Y1_X1_strat5_GR)
median_SE_hat_Y1_X2_strat5_GR <- median(SE_hat_Y1_X2_strat5_GR)
median_SE_hat_Y2_X1_strat5_GR <- median(SE_hat_Y2_X1_strat5_GR)
median_SE_hat_Y2_X2_strat5_GR <- median(SE_hat_Y2_X2_strat5_GR)
coverage_Y1_X1_strat5_GR <- mean(cover_Y1_X1_strat5_GR)
coverage_Y1_X2_strat5_GR <- mean(cover_Y1_X2_strat5_GR)
coverage_Y2_X1_strat5_GR <- mean(cover_Y2_X1_strat5_GR)
coverage_Y2_X2_strat5_GR <- mean(cover_Y2_X2_strat5_GR)

# IPW results
var_B_hat_Y1_X1_strat5_IPW <- var(B_hat_Y1_X1_strat5_IPW)
var_B_hat_Y1_X2_strat5_IPW <- var(B_hat_Y1_X2_strat5_IPW)
var_B_hat_Y2_X1_strat5_IPW <- var(B_hat_Y2_X1_strat5_IPW)
var_B_hat_Y2_X2_strat5_IPW <- var(B_hat_Y2_X2_strat5_IPW)
mean_B_hat_Y1_X1_strat5_IPW <- mean(B_hat_Y1_X1_strat5_IPW)
mean_B_hat_Y1_X2_strat5_IPW <- mean(B_hat_Y1_X2_strat5_IPW)
mean_B_hat_Y2_X1_strat5_IPW <- mean(B_hat_Y2_X1_strat5_IPW)
mean_B_hat_Y2_X2_strat5_IPW <- mean(B_hat_Y2_X2_strat5_IPW)
bias_B_11_strat5_IPW <- mean(B_hat_Y1_X1_strat5_IPW - B11_std_strat5)
bias_B_21_strat5_IPW <- mean(B_hat_Y1_X2_strat5_IPW - B21_std_strat5)
bias_B_12_strat5_IPW <- mean(B_hat_Y2_X1_strat5_IPW - B12_std_strat5)
bias_B_22_strat5_IPW <- mean(B_hat_Y2_X2_strat5_IPW - B22_std_strat5)
MSE_B_11_strat5_IPW <- mean((B_hat_Y1_X1_strat5_IPW - B11_std_strat1)^2)
MSE_B_21_strat5_IPW <- mean((B_hat_Y1_X2_strat5_IPW - B21_std_strat1)^2)
MSE_B_12_strat5_IPW <- mean((B_hat_Y2_X1_strat5_IPW - B12_std_strat1)^2)
MSE_B_22_strat5_IPW <- mean((B_hat_Y2_X2_strat5_IPW - B22_std_strat1)^2)
median_B_hat_Y1_X1_strat5_IPW <- median(B_hat_Y1_X1_strat5_IPW)
median_B_hat_Y1_X2_strat5_IPW <- median(B_hat_Y1_X2_strat5_IPW)
median_B_hat_Y2_X1_strat5_IPW <- median(B_hat_Y2_X1_strat5_IPW)
median_B_hat_Y2_X2_strat5_IPW <- median(B_hat_Y2_X2_strat5_IPW)
median_SE_hat_Y1_X1_strat5_IPW <- median(SE_hat_Y1_X1_strat5_IPW)
median_SE_hat_Y1_X2_strat5_IPW <- median(SE_hat_Y1_X2_strat5_IPW)
median_SE_hat_Y2_X1_strat5_IPW <- median(SE_hat_Y2_X1_strat5_IPW)
median_SE_hat_Y2_X2_strat5_IPW <- median(SE_hat_Y2_X2_strat5_IPW)
coverage_Y1_X1_strat5_IPW <- mean(cover_Y1_X1_strat5_IPW)
coverage_Y1_X2_strat5_IPW <- mean(cover_Y1_X2_strat5_IPW)
coverage_Y2_X1_strat5_IPW <- mean(cover_Y2_X1_strat5_IPW)
coverage_Y2_X2_strat5_IPW <- mean(cover_Y2_X2_strat5_IPW)

# Save raw results
outpathname5 <- paste0(pathname,"/", names(simulations_df)[scenario], "Strat5raw.csv")
write.csv(results_matrix, outpathname5, row.names = FALSE)

#####
## View results
#####

results <- data.frame(Strategy = c(1,2,3,4,5),
                      "true_B_11" = c(simulations_df["B11", scenario], 0,0,0,0),
                      "true_B_21" = c(simulations_df["B21", scenario], 0,0,0,0),
                      "true_B_12" = c(simulations_df["B12", scenario], 0,0,0,0),
                      "true_B_22" = c(simulations_df["B22", scenario], 0,0,0,0),
                      "true_B_11_std" = c(mean(B11_std_strat1),
                                          mean(B11_std_strat2),
                                          mean(B11_std_strat3),
                                          mean(B11_std_strat4),
                                          mean(B11_std_strat5)),
                      "true_B_21_std" = c(mean(B21_std_strat1),
                                          mean(B21_std_strat2),
                                          mean(B21_std_strat3),
                                          mean(B21_std_strat4),
                                          mean(B21_std_strat5)),
                      "true_B_12_std" = c(mean(B12_std_strat1),
                                          mean(B12_std_strat2),
                                          mean(B12_std_strat3),
                                          mean(B12_std_strat4),
                                          mean(B12_std_strat5)),
                      "true_B_22_std" = c(mean(B22_std_strat1),
                                          mean(B22_std_strat2),
                                          mean(B22_std_strat3),
                                          mean(B22_std_strat4),
                                          mean(B22_std_strat5)),
                      "cor_Y1_Y2" = c(cor_Y1Y2_strat1, cor_Y1Y2_strat2,
                                      cor_Y1Y2_strat3,
                                      cor_Y1Y2_strat4,
                                      cor_Y1Y2_strat5),
                      "cor_X1_X2" = c(cor_X1X2_strat1[1], "", "", "", ""),
                      "prev_Y1" = c(prev_Y1_strat1, prev_Y1_strat2,
                                    prev_Y1_strat3, prev_Y1_strat4,
                                    prev_Y1_strat5),
                      "prev_Y2" = c(prev_Y2_strat1, prev_Y2_strat2,
                                    prev_Y2_strat3,  prev_Y2_strat4,
                                    prev_Y1_strat5),
                      "mean(B_11)_GR" = c(mean_B_hat_Y1_X1_strat1_GR,
                                          mean_B_hat_Y1_X1_strat2_GR,
                                          mean_B_hat_Y1_X1_strat3_GR,
                                          mean_B_hat_Y1_X1_strat4_GR,
                                          mean_B_hat_Y1_X1_strat5_GR),
                      "mean(B_21)_GR" = c(mean_B_hat_Y1_X2_strat1_GR,
                                          mean_B_hat_Y1_X2_strat2_GR,
                                          mean_B_hat_Y1_X2_strat3_GR,
                                          mean_B_hat_Y1_X2_strat4_GR,
                                          mean_B_hat_Y1_X2_strat5_GR),
                      "mean(B_12)_GR" = c(mean_B_hat_Y2_X1_strat1_GR,
                                          mean_B_hat_Y2_X1_strat2_GR,
                                          mean_B_hat_Y2_X1_strat3_GR,
                                          mean_B_hat_Y2_X1_strat4_GR,
                                          mean_B_hat_Y2_X1_strat5_GR),
                      "mean(B_22)_GR" = c(mean_B_hat_Y2_X2_strat1_GR,
                                          mean_B_hat_Y2_X2_strat2_GR,
                                          mean_B_hat_Y2_X2_strat3_GR,
                                          mean_B_hat_Y2_X2_strat4_GR,
                                          mean_B_hat_Y2_X2_strat5_GR),
                      "median(B_11)_GR" = c(median_B_hat_Y1_X1_strat1_GR,
                                          median_B_hat_Y1_X1_strat2_GR,
                                          median_B_hat_Y1_X1_strat3_GR,
                                          median_B_hat_Y1_X1_strat4_GR,
                                          median_B_hat_Y1_X1_strat5_GR),
                      "median(B_21)_GR" = c(median_B_hat_Y1_X2_strat1_GR,
                                          median_B_hat_Y1_X2_strat2_GR,
                                          median_B_hat_Y1_X2_strat3_GR,
                                          median_B_hat_Y1_X2_strat4_GR,
                                          median_B_hat_Y1_X2_strat5_GR),
                      "median(B_12)_GR" = c(median_B_hat_Y2_X1_strat1_GR,
                                          median_B_hat_Y2_X1_strat2_GR,
                                          median_B_hat_Y2_X1_strat3_GR,
                                          median_B_hat_Y2_X1_strat4_GR,
                                          median_B_hat_Y2_X1_strat5_GR),
                      "median(B_22)_GR" = c(median_B_hat_Y2_X2_strat1_GR,
                                          median_B_hat_Y2_X2_strat2_GR,
                                          median_B_hat_Y2_X2_strat3_GR,
                                          median_B_hat_Y2_X2_strat4_GR,
                                          median_B_hat_Y2_X2_strat5_GR),
                      "var(B_11)_GR" = c(var_B_hat_Y1_X1_strat1_GR,
                                         var_B_hat_Y1_X1_strat2_GR,
                                         var_B_hat_Y1_X1_strat3_GR,
                                         var_B_hat_Y1_X1_strat4_GR,
                                         var_B_hat_Y1_X1_strat5_GR),
                      "var(B_21)_GR" = c(var_B_hat_Y1_X2_strat1_GR,
                                         var_B_hat_Y1_X2_strat2_GR,
                                         var_B_hat_Y1_X2_strat3_GR,
                                         var_B_hat_Y1_X2_strat4_GR,
                                         var_B_hat_Y1_X2_strat5_GR),
                      "var(B_12)_GR" = c(var_B_hat_Y2_X1_strat1_GR,
                                         var_B_hat_Y2_X1_strat2_GR,
                                         var_B_hat_Y2_X1_strat3_GR,
                                         var_B_hat_Y2_X1_strat4_GR,
                                         var_B_hat_Y2_X1_strat5_GR),
                      "var(B_22)_GR" = c(var_B_hat_Y2_X2_strat1_GR,
                                         var_B_hat_Y2_X2_strat2_GR,
                                         var_B_hat_Y2_X2_strat3_GR,
                                         var_B_hat_Y2_X2_strat4_GR,
                                         var_B_hat_Y2_X2_strat5_GR),
                      "mean(B_11)_IPW" = c(mean_B_hat_Y1_X1_strat1_IPW,
                                          mean_B_hat_Y1_X1_strat2_IPW,
                                          mean_B_hat_Y1_X1_strat3_IPW,
                                          mean_B_hat_Y1_X1_strat4_IPW,
                                          mean_B_hat_Y1_X1_strat5_IPW),
                      "mean(B_21)_IPW" = c(mean_B_hat_Y1_X2_strat1_IPW,
                                          mean_B_hat_Y1_X2_strat2_IPW,
                                          mean_B_hat_Y1_X2_strat3_IPW,
                                          mean_B_hat_Y1_X2_strat4_IPW,
                                          mean_B_hat_Y1_X2_strat5_IPW),
                      "mean(B_12)_IPW" = c(mean_B_hat_Y2_X1_strat1_IPW,
                                          mean_B_hat_Y2_X1_strat2_IPW,
                                          mean_B_hat_Y2_X1_strat3_IPW,
                                          mean_B_hat_Y2_X1_strat4_IPW,
                                          mean_B_hat_Y2_X1_strat5_IPW),
                      "mean(B_22)_IPW" = c(mean_B_hat_Y2_X2_strat1_IPW,
                                          mean_B_hat_Y2_X2_strat2_IPW,
                                          mean_B_hat_Y2_X2_strat3_IPW,
                                          mean_B_hat_Y2_X2_strat4_IPW,
                                          mean_B_hat_Y2_X2_strat5_IPW),
                      "median(B_11)_IPW" = c(median_B_hat_Y1_X1_strat1_IPW,
                                            median_B_hat_Y1_X1_strat2_IPW,
                                            median_B_hat_Y1_X1_strat3_IPW,
                                            median_B_hat_Y1_X1_strat4_IPW,
                                            median_B_hat_Y1_X1_strat5_IPW),
                      "median(B_21)_IPW" = c(median_B_hat_Y1_X2_strat1_IPW,
                                            median_B_hat_Y1_X2_strat2_IPW,
                                            median_B_hat_Y1_X2_strat3_IPW,
                                            median_B_hat_Y1_X2_strat4_IPW,
                                            median_B_hat_Y1_X2_strat5_IPW),
                      "median(B_12)_IPW" = c(median_B_hat_Y2_X1_strat1_IPW,
                                            median_B_hat_Y2_X1_strat2_IPW,
                                            median_B_hat_Y2_X1_strat3_IPW,
                                            median_B_hat_Y2_X1_strat4_IPW,
                                            median_B_hat_Y2_X1_strat5_IPW),
                      "median(B_22)_IPW" = c(median_B_hat_Y2_X2_strat1_IPW,
                                            median_B_hat_Y2_X2_strat2_IPW,
                                            median_B_hat_Y2_X2_strat3_IPW,
                                            median_B_hat_Y2_X2_strat4_IPW,
                                            median_B_hat_Y2_X2_strat5_IPW),
                      "var(B_11)_IPW" = c(var_B_hat_Y1_X1_strat1_IPW,
                                         var_B_hat_Y1_X1_strat2_IPW,
                                         var_B_hat_Y1_X1_strat3_IPW,
                                         var_B_hat_Y1_X1_strat4_IPW,
                                         var_B_hat_Y1_X1_strat5_IPW),
                      "var(B_21)_IPW" = c(var_B_hat_Y1_X2_strat1_IPW,
                                         var_B_hat_Y1_X2_strat2_IPW,
                                         var_B_hat_Y1_X2_strat3_IPW,
                                         var_B_hat_Y1_X2_strat4_IPW,
                                         var_B_hat_Y1_X2_strat5_IPW),
                      "var(B_12)_IPW" = c(var_B_hat_Y2_X1_strat1_IPW,
                                         var_B_hat_Y2_X1_strat2_IPW,
                                         var_B_hat_Y2_X1_strat3_IPW,
                                         var_B_hat_Y2_X1_strat4_IPW,
                                         var_B_hat_Y2_X1_strat5_IPW),
                      "var(B_22)_IPW" = c(var_B_hat_Y2_X2_strat1_IPW,
                                         var_B_hat_Y2_X2_strat2_IPW,
                                         var_B_hat_Y2_X2_strat3_IPW,
                                         var_B_hat_Y2_X2_strat4_IPW,
                                         var_B_hat_Y2_X2_strat5_IPW),
                      "coverage_B11_GR" = c(coverage_Y1_X1_strat1_GR,
                                             coverage_Y1_X1_strat2_GR,
                                             coverage_Y1_X1_strat3_GR,
                                             coverage_Y1_X1_strat4_GR,
                                            coverage_Y1_X1_strat5_GR),
                      "coverage_B21_GR" = c(coverage_Y1_X2_strat1_GR,
                                             coverage_Y1_X2_strat2_GR,
                                             coverage_Y1_X2_strat3_GR,
                                             coverage_Y1_X2_strat4_GR,
                                            coverage_Y1_X2_strat5_GR),
                      "coverage_B12_GR" = c(coverage_Y2_X1_strat1_GR,
                                             coverage_Y2_X1_strat2_GR,
                                             coverage_Y2_X1_strat3_GR,
                                             coverage_Y2_X1_strat4_GR,
                                            coverage_Y2_X1_strat5_GR),
                      "coverage_B22_GR" = c(coverage_Y2_X2_strat1_GR,
                                             coverage_Y2_X2_strat2_GR,
                                             coverage_Y2_X2_strat3_GR,
                                             coverage_Y2_X2_strat4_GR,
                                            coverage_Y2_X2_strat5_GR),
                      "coverage_B11_IPW" = c(coverage_Y1_X1_strat1_IPW,
                                            coverage_Y1_X1_strat2_IPW,
                                            coverage_Y1_X1_strat3_IPW,
                                            coverage_Y1_X1_strat4_IPW,
                                            coverage_Y1_X1_strat5_IPW),
                      "coverage_B21_IPW" = c(coverage_Y1_X2_strat1_IPW,
                                            coverage_Y1_X2_strat2_IPW,
                                            coverage_Y1_X2_strat3_IPW,
                                            coverage_Y1_X2_strat4_IPW,
                                            coverage_Y1_X2_strat5_IPW),
                      "coverage_B12_IPW" = c(coverage_Y2_X1_strat1_IPW,
                                             coverage_Y2_X1_strat2_IPW,
                                             coverage_Y2_X1_strat3_IPW,
                                             coverage_Y2_X1_strat4_IPW,
                                             coverage_Y2_X1_strat5_IPW),
                      "coverage_B22_IPW" = c(coverage_Y2_X2_strat1_IPW,
                                             coverage_Y2_X2_strat2_IPW,
                                             coverage_Y2_X2_strat3_IPW,
                                             coverage_Y2_X2_strat4_IPW,
                                             coverage_Y2_X2_strat5_IPW),
                      "bias_B11_GR" = c(bias_B_11_strat1_GR,
                                            bias_B_11_strat2_GR,
                                            bias_B_11_strat3_GR,
                                            bias_B_11_strat4_GR,
                                        bias_B_11_strat5_GR),
                      "bias_B21_GR" = c(bias_B_21_strat1_GR,
                                            bias_B_21_strat2_GR,
                                            bias_B_21_strat3_GR,
                                            bias_B_21_strat4_GR,
                                        bias_B_21_strat5_GR),
                      "bias_B12_GR" = c(bias_B_12_strat1_GR,
                                            bias_B_12_strat2_GR,
                                            bias_B_12_strat3_GR,
                                            bias_B_12_strat4_GR,
                                        bias_B_12_strat5_GR),
                      "bias_B22_GR" = c(bias_B_22_strat1_GR,
                                            bias_B_22_strat2_GR,
                                            bias_B_22_strat3_GR,
                                            bias_B_22_strat4_GR,
                                        bias_B_22_strat5_GR),
                      "bias_B11_IPW" = c(bias_B_11_strat1_IPW,
                                             bias_B_11_strat2_IPW,
                                             bias_B_11_strat3_IPW,
                                             bias_B_11_strat4_IPW,
                                         bias_B_11_strat5_IPW),
                      "bias_B21_IPW" = c(bias_B_21_strat1_IPW,
                                             bias_B_21_strat2_IPW,
                                             bias_B_21_strat3_IPW,
                                             bias_B_21_strat4_IPW,
                                         bias_B_21_strat5_IPW),
                      "bias_B12_IPW" = c(bias_B_12_strat1_IPW,
                                             bias_B_12_strat2_IPW,
                                             bias_B_12_strat3_IPW,
                                             bias_B_12_strat4_IPW,
                                         bias_B_12_strat5_IPW),
                      "bias_B22_IPW" = c(bias_B_22_strat1_IPW,
                                             bias_B_22_strat2_IPW,
                                             bias_B_22_strat3_IPW,
                                             bias_B_22_strat4_IPW,
                                         bias_B_22_strat5_IPW),
                      "median_ASE_B11_GR" = c(median_SE_hat_Y1_X1_strat1_GR,
                                               median_SE_hat_Y1_X1_strat2_GR,
                                               median_SE_hat_Y1_X1_strat3_GR,
                                               median_SE_hat_Y1_X1_strat4_GR,
                                              median_SE_hat_Y1_X1_strat5_GR),
                      "median_ASE_B21_GR" = c(median_SE_hat_Y1_X2_strat1_GR,
                                               median_SE_hat_Y1_X2_strat2_GR,
                                               median_SE_hat_Y1_X2_strat3_GR,
                                               median_SE_hat_Y1_X2_strat4_GR,
                                              median_SE_hat_Y1_X2_strat5_GR),
                      "median_ASE_B12_GR" = c(median_SE_hat_Y2_X1_strat1_GR,
                                               median_SE_hat_Y2_X1_strat2_GR,
                                               median_SE_hat_Y2_X1_strat3_GR,
                                               median_SE_hat_Y2_X1_strat4_GR,
                                              median_SE_hat_Y2_X1_strat5_GR),
                      "median_ASE_B22_GR" = c(median_SE_hat_Y2_X2_strat1_GR,
                                               median_SE_hat_Y2_X2_strat2_GR,
                                               median_SE_hat_Y2_X2_strat3_GR,
                                               median_SE_hat_Y2_X2_strat4_GR,
                                              median_SE_hat_Y2_X2_strat5_GR),
                      "median_ASE_B11_IPW" = c(median_SE_hat_Y1_X1_strat1_IPW,
                                              median_SE_hat_Y1_X1_strat2_IPW,
                                              median_SE_hat_Y1_X1_strat3_IPW,
                                              median_SE_hat_Y1_X1_strat4_IPW,
                                              median_SE_hat_Y1_X1_strat5_IPW),
                      "median_ASE_B21_IPW" = c(median_SE_hat_Y1_X2_strat1_IPW,
                                               median_SE_hat_Y1_X2_strat2_IPW,
                                               median_SE_hat_Y1_X2_strat3_IPW,
                                               median_SE_hat_Y1_X2_strat4_IPW,
                                               median_SE_hat_Y1_X2_strat5_IPW),
                      "median_ASE_B12_IPW" = c(median_SE_hat_Y2_X1_strat1_IPW,
                                               median_SE_hat_Y2_X1_strat2_IPW,
                                               median_SE_hat_Y2_X1_strat3_IPW,
                                               median_SE_hat_Y2_X1_strat4_IPW,
                                               median_SE_hat_Y2_X1_strat5_IPW),
                      "median_ASE_B22_IPW" = c(median_SE_hat_Y2_X2_strat1_IPW,
                                               median_SE_hat_Y2_X2_strat2_IPW,
                                               median_SE_hat_Y2_X2_strat3_IPW,
                                               median_SE_hat_Y2_X2_strat4_IPW,
                                               median_SE_hat_Y2_X2_strat5_IPW),
                      "MSE_B11_IPW" = c(MSE_B_11_strat1_IPW,
                                       MSE_B_11_strat2_IPW,
                                       MSE_B_11_strat3_IPW,
                                       MSE_B_11_strat4_IPW,
                                       MSE_B_11_strat5_IPW),
                      "MSE_B21_IPW" = c(MSE_B_21_strat1_IPW,
                                       MSE_B_21_strat2_IPW,
                                       MSE_B_21_strat3_IPW,
                                       MSE_B_21_strat4_IPW,
                                       MSE_B_21_strat5_IPW),
                      "MSE_B12_IPW" = c(MSE_B_12_strat1_IPW,
                                        MSE_B_12_strat2_IPW,
                                        MSE_B_12_strat3_IPW,
                                        MSE_B_12_strat4_IPW,
                                        MSE_B_12_strat5_IPW),
                      "MSE_B22_IPW" = c(MSE_B_22_strat1_IPW,
                                        MSE_B_22_strat2_IPW,
                                        MSE_B_22_strat3_IPW,
                                        MSE_B_22_strat4_IPW,
                                        MSE_B_22_strat5_IPW),
                      "MSE_B11_GR" = c(MSE_B_11_strat1_GR,
                                        MSE_B_11_strat2_GR,
                                        MSE_B_11_strat3_GR,
                                        MSE_B_11_strat4_GR,
                                        MSE_B_11_strat5_GR),
                      "MSE_B21_GR" = c(MSE_B_21_strat1_GR,
                                        MSE_B_21_strat2_GR,
                                        MSE_B_21_strat3_GR,
                                        MSE_B_21_strat4_GR,
                                        MSE_B_21_strat5_GR),
                      "MSE_B12_GR" = c(MSE_B_12_strat1_GR,
                                        MSE_B_12_strat2_GR,
                                        MSE_B_12_strat3_GR,
                                        MSE_B_12_strat4_GR,
                                        MSE_B_12_strat5_GR),
                      "MSE_B22_GR" = c(MSE_B_22_strat1_GR,
                                        MSE_B_22_strat2_GR,
                                        MSE_B_22_strat3_GR,
                                        MSE_B_22_strat4_GR,
                                        MSE_B_22_strat5_GR)
)

print(results)

# save
outpathname <- paste0(pathname,"/", names(simulations_df)[scenario], ".csv")
write.csv(results, outpathname, row.names = FALSE)