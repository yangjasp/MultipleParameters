#####
### Functions for Each Strategy
#####

####
####  
#### Four coefficients of interest: B_Y1_X1, B_Y1_X2, B_Y2_X1, B_Y2_X2

####
## Strategy 1 - Case-control sampling. n/4 cases and n/4 controls for each outcome
####

# @param n sample size
# @param data Phase 1 dataset with cols for Y1, Y2, X, X2, Z and
# Y1_obs, Y2_obs, X1_obs, X2_obs, Z_obs

Strategy1 <- function(n, simulations_df,
                      scenario, N = 10000){
  
  # Set up betas for scenario
  B01 <- simulations_df["B01", scenario]
  B11 <- simulations_df["B11", scenario]
  B21 <- simulations_df["B21", scenario]
  B31 <- simulations_df["B31", scenario]
  B41 <- simulations_df["B41", scenario]
  B02 <- simulations_df["B02", scenario]
  B12 <- simulations_df["B12", scenario]
  B22 <- simulations_df["B22", scenario]
  B32 <- simulations_df["B32", scenario] 
  
  # Set up of Xs
  corX1X2 <- simulations_df["corX1X2", scenario]
  sensZ <- simulations_df["sensZ", scenario]
  specZ <- simulations_df["specZ", scenario]
  prevZ <- simulations_df["prevZ", scenario]
  
  # Measurement error
  sensY1 <- simulations_df["sensY1", scenario]
  specY1 <- simulations_df["specY1", scenario]
  sensY2 <- simulations_df["sensY2", scenario]
  specY2 <- simulations_df["specY2", scenario]
  error_varX1 <- simulations_df["error_varX1", scenario]
  error_varX2 <- simulations_df["error_varX2", scenario]
  
  
  #####
  ## Generate Phase 1 Data
  #####
  
  id <- 1:N
  
  # Generate correlated covariates
  sigma <- matrix(c(1, corX1X2, 0.1, corX1X2, 1,  0.25, 0.1, 0.25, 1), nrow = 3)
  covs <- mvrnorm(N, mu = c(0,0,0), Sigma = sigma)
  X1 <- covs[,1]
  X2 <- covs[,2]
  Z <- ifelse(covs[,3] >= quantile(covs[,3], 1 - prevZ), 1, 0)
  cor(cbind(X1, X2, Z)) 
  
  Y2_probs <- exp(B02 + B12*X1 + B22*X2 + B32*Z)/(1 + exp(B02 + B12*X1 + B22*X2 + B32*Z))
  Y2 <- rbinom(N, 1, Y2_probs) # Compute realized value for Y2
  Y1_probs <- exp(B01 + B11*X1 + B21*X2 + B31*Z + B41*Y2)/
    (1 + exp(B01 + B11*X1 + B21*X2 + B31*Z + B41*Y2))
  Y1 <- rbinom(N, 1, Y1_probs) # Compute realized value for Y1
  
  full_data <- data.frame(id, X1, X2, Z, Y1, Y2)
  
  ## Test true model
  true_model1 <- glm(Y1 ~ X1 + X2 + Z + Y2, family = "binomial", data = full_data)
  true_model2 <- glm(Y2 ~ X1 + X2 + Z, family = "binomial", data = full_data)
  
  # Find observed correlation
  cor_Y1_Y2 <- cor(Y1, Y2) 
  
  # Prevalence
  prev_Y1 <- table(Y1)["1"]/length(Y1)
  prev_Y2 <- table(Y2)["1"]/length(Y2)
  
  # Find true coefficients
  true_modelY1 <- glm(Y1 ~ X1 + X2 + Z + Y2, family = "binomial", data = full_data) 
  true_modelY2 <- glm(Y2 ~ X1 + X2 + Z, family = "binomial", data = full_data) 
  true_B_Y1_X1 <- coef(true_modelY1)["X1"]
  true_B_Y1_X2 <- coef(true_modelY1)["X2"]
  true_B_Y2_X1 <- coef(true_modelY2)["X1"]
  true_B_Y2_X2 <- coef(true_modelY2)["X2"]
  
  # compute true influence functions
  full_data$infl_Y1_X1_true <- inf_fun_logit(true_modelY1)[,"X1"]
  full_data$infl_Y1_X2_true <- inf_fun_logit(true_modelY1)[,"X2"]
  full_data$infl_Y2_X1_true <- inf_fun_logit(true_modelY2)[,"X1"]
  full_data$infl_Y2_X2_true <- inf_fun_logit(true_modelY2)[,"X2"]
  
  
  #####
  ## Generate error-prone Phase 1 data
  ######
  X1X2ZY1Y2_obs <- mvrnorm(N, c(0,0,0,0,0),
                           matrix(c(error_varX1, 0.05, 0.04, 0.15, 0.04,
                                    0.05, error_varX2, 0.02, 0, 0.4,
                                    0.04, 0.02, 1, 0.01, 0,
                                    0.15, 0, 0.01, 1 , 0.15,
                                    0.04, 0.4, 0, 0.15, 1), nrow = 5))
  full_data$X1_obs <- full_data$X1 + X1X2ZY1Y2_obs[,1]
  full_data$X2_obs <- full_data$X2 + X1X2ZY1Y2_obs[,2]
  
  threshold_positiveZ <- qnorm((sensZ + 1)/2)
  threshold_negativeZ <- qnorm((specZ + 1)/2)
  threshold_positiveY1 <- qnorm((sensY1 + 1)/2)
  threshold_negativeY1 <- qnorm((specY1 + 1)/2)
  threshold_positiveY2 <- qnorm((sensY2 + 1)/2)
  threshold_negativeY2 <- qnorm((specY2 + 1)/2)
  
  full_data$Z_obs <- with(full_data, ifelse((Z > 0 & abs(X1X2ZY1Y2_obs[, 3]) < threshold_positiveZ) | 
                                               (Z < 0 & abs(X1X2ZY1Y2_obs[, 3]) > threshold_negativeZ), 
                                             1, 0))
  
  full_data$Y1_obs <- with(full_data, ifelse((Y1 == 1 & abs(X1X2ZY1Y2_obs[, 4]) < threshold_positiveY1) | 
                                               (Y1 == 0 & abs(X1X2ZY1Y2_obs[, 4]) > threshold_negativeY1), 1, 0))
  
  full_data$Y2_obs <- with(full_data, ifelse((Y2 == 1 & abs(X1X2ZY1Y2_obs[, 5]) < threshold_positiveY2) | 
                                               (Y2 == 0 & abs(X1X2ZY1Y2_obs[, 5]) > threshold_negativeY2),1, 0))
  
  
  #### Compute naive influence functions
  fitY1_phase1 <-  glm(Y1_obs ~ X1_obs + X2_obs + Z_obs + Y2_obs, 
                       family = "binomial", data = full_data)
  fitY2_phase1 <-  glm(Y2_obs ~ X1_obs + X2_obs + Z_obs, 
                       family = "binomial", data = full_data)
  
  full_data$inflB11_phase1 <- inf_fun_logit(fitY1_phase1)[,"X1_obs"]
  full_data$inflB21_phase1 <- inf_fun_logit(fitY1_phase1)[,"X2_obs"]
  full_data$inflB12_phase1 <- inf_fun_logit(fitY2_phase1)[,"X1_obs"]
  full_data$inflB22_phase1 <- inf_fun_logit(fitY2_phase1)[,"X2_obs"]
  
  #### Divide into strata based on observed phase 1 data X
  full_data <- full_data |> 
    dplyr::mutate(Y1_strat = Y1_obs, Y2_strat = Y2_obs, strata = interaction(Y1_strat, Y2_strat))
  
  data <- full_data
  
  #####
  #####
  ##### Phase 2
  
  ### Collect samples via case-control sampling
  design_Y1 <- data.frame("strata" = c(1,0), "n_to_sample" = c(n/4,n/4))
  design_Y2 <- data.frame("strata" = c(1,0), "n_to_sample" = c(n/4,n/4))
  
  data <- data |>
    optimall::sample_strata(strata = "Y1_strat", id = "id",
                            design_data = design_Y1) |>
    dplyr::mutate(sample1_indicator = sample_indicator) |>
    optimall::sample_strata(strata = "Y2_strat", id = "id",
                            design_data = design_Y2,
                            already_sampled = "sample_indicator") |>
    dplyr::mutate(sample2_indicator = sample_indicator,
                  sample_indicator = sample2_indicator + sample1_indicator)
  
  # Sample
  data$X1 <- ifelse(data$sample_indicator == 1, full_data$X1 , NA)
  data$X2 <- ifelse(data$sample_indicator == 1, full_data$X2 , NA)
  data$Z <- ifelse(data$sample_indicator == 1, full_data$Z , NA)
  
  # Generalized Raking
  twophase_design <- twophase(id = list(~1, ~1), 
                              strata = list(NULL, ~strata), 
                              subset = ~as.logical(sample_indicator), data = data)
  
  # Weights
  weightY1 <- svyglm(Y1 ~ X1 + X2 + Z + Y2, family = quasibinomial, 
                     design = twophase_design)
  weightY2 <- svyglm(Y2 ~ X1 + X2 + Z, family = quasibinomial, 
                     design = twophase_design)
  
  #########
  #### Raking 
  #########
  # Creating a survey object and calibrating the weights
  mydesign <- survey::twophase(
    id = list(~1, ~1), subset = ~as.logical(sample_indicator),
    strata = list(NULL, ~strata), data = data,
  )
  infcalY1 <- survey::calibrate(mydesign, formula = ~inflB11_phase1 +
                                  inflB21_phase1 + inflB12_phase1 +
                                  inflB22_phase1 + strata, phase = 2, calfun = "raking")
  infcalY2 <- survey::calibrate(mydesign, formula =  ~inflB11_phase1 +
                                  inflB21_phase1 + inflB12_phase1 +
                                  inflB22_phase1 + strata, phase = 2, calfun = "raking")
  
  # Fitting the outcome model:
  fitY1 <- survey::svyglm(Y1 ~ X1 + X2 + Z + Y2, design = infcalY1, family = quasibinomial)
  fitY2 <- survey::svyglm(Y2 ~ X1 + X2 + Z, design = infcalY2, family = quasibinomial)
  
  # Check if truth is contained in interval
  confintIPWB11 <- confint(weightY1)["X1",]
  confintIPWB21 <- confint(weightY1)["X2",]
  confintIPWB12 <- confint(weightY2)["X1",]
  confintIPWB22 <- confint(weightY2)["X2",]
  confintGRB11 <- confint(fitY1)["X1",]
  confintGRB21 <- confint(fitY1)["X2",]
  confintGRB12 <- confint(fitY2)["X1",]
  confintGRB22 <- confint(fitY2)["X2",]
  coverIPWB11 <- ifelse(B11 >= confintIPWB11[1] & B11 <= confintIPWB11[2],1,0)
  coverIPWB21 <- ifelse(B21 >= confintIPWB21[1] & B21 <= confintIPWB21[2],1,0)
  coverIPWB12 <- ifelse(B12 >= confintIPWB12[1] & B12 <= confintIPWB12[2],1,0)
  coverIPWB22 <- ifelse(B22 >= confintIPWB22[1] & B22 <= confintIPWB22[2],1,0)
  coverGRB11 <- ifelse(B11 >= confintGRB11[1] & B11 <= confintGRB11[2],1,0)
  coverGRB21 <- ifelse(B21 >= confintGRB21[1] & B21 <= confintGRB21[2],1,0)
  coverGRB12 <- ifelse(B12 >= confintGRB12[1] & B12 <= confintGRB12[2],1,0)
  coverGRB22 <- ifelse(B22 >= confintGRB22[1] & B22 <= confintGRB22[2],1,0)
  
  
  output <- c(fitY1$coefficients["X1"], fitY1$coefficients["X2"],
              fitY2$coefficients["X1"], fitY2$coefficients["X2"],
              coef(weightY1)["X1"], coef(weightY1)["X2"], 
              coef(weightY2)["X1"], coef(weightY2)["X2"],  
              SE(fitY1)["X1"], SE(fitY1)["X2"],
              SE(fitY2)["X1"], SE(fitY2)["X2"], 
              SE(weightY1)["X1"], SE(weightY1)["X2"], 
              SE(weightY2)["X1"], SE(weightY2)["X2"], 
              coverIPWB11, coverIPWB21,
              coverIPWB12, coverIPWB22, 
              coverGRB11, coverGRB21,
              coverGRB12,coverGRB22,
              cor_Y1_Y2,
              prev_Y1, prev_Y2, B11, B21, B12, B22)
  return(output)
}


#####
## Strategy 2- Simulteaneous Multi-wave sampling
#####

# This is looped part, so
Strategy2 <- function(n, simulations_df, scenario, N = 10000){
  
  # Set up betas for scenario
  B01 <- simulations_df["B01", scenario]
  B11 <- simulations_df["B11", scenario]
  B21 <- simulations_df["B21", scenario]
  B31 <- simulations_df["B31", scenario]
  B41 <- simulations_df["B41", scenario]
  B02 <- simulations_df["B02", scenario]
  B12 <- simulations_df["B12", scenario]
  B22 <- simulations_df["B22", scenario]
  B32 <- simulations_df["B32", scenario] 
  
  # Set up of Xs
  corX1X2 <- simulations_df["corX1X2", scenario]
  sensZ <- simulations_df["sensZ", scenario]
  specZ <- simulations_df["specZ", scenario]
  prevZ <- simulations_df["prevZ", scenario]
  
  # Measurement error
  sensY1 <- simulations_df["sensY1", scenario]
  specY1 <- simulations_df["specY1", scenario]
  sensY2 <- simulations_df["sensY2", scenario]
  specY2 <- simulations_df["specY2", scenario]
  error_varX1 <- simulations_df["error_varX1", scenario]
  error_varX2 <- simulations_df["error_varX2", scenario]
  
  
  #####
  ## Generate Phase 1 Data
  #####
  
  id <- 1:N
  
  # Generate correlated covariates
  sigma <- matrix(c(1, corX1X2, 0.1, corX1X2, 1,  0.25, 0.1, 0.25, 1), nrow = 3)
  covs <- mvrnorm(N, mu = c(0,0,0), Sigma = sigma)
  X1 <- covs[,1]
  X2 <- covs[,2]
  Z <- ifelse(covs[,3] >= quantile(covs[,3], 1 - prevZ), 1, 0)
  cor(cbind(X1, X2, Z)) 
  
  Y2_probs <- exp(B02 + B12*X1 + B22*X2 + B32*Z)/(1 + exp(B02 + B12*X1 + B22*X2 + B32*Z))
  Y2 <- rbinom(N, 1, Y2_probs) # Compute realized value for Y2
  Y1_probs <- exp(B01 + B11*X1 + B21*X2 + B31*Z + B41*Y2)/
    (1 + exp(B01 + B11*X1 + B21*X2 + B31*Z + B41*Y2))
  Y1 <- rbinom(N, 1, Y1_probs) # Compute realized value for Y1
  
  full_data <- data.frame(id, X1, X2, Z, Y1, Y2)
  
  ## Test true model
  true_model1 <- glm(Y1 ~ X1 + X2 + Z + Y2, family = "binomial", data = full_data)
  true_model2 <- glm(Y2 ~ X1 + X2 + Z, family = "binomial", data = full_data)
  
  # Find observed correlation
  cor_Y1_Y2 <- cor(Y1, Y2) 
  
  # Prevalence
  prev_Y1 <- table(Y1)["1"]/length(Y1)
  prev_Y2 <- table(Y2)["1"]/length(Y2)
  
  # Find true coefficients
  true_modelY1 <- glm(Y1 ~ X1 + X2 + Z + Y2, family = "binomial", data = full_data) 
  true_modelY2 <- glm(Y2 ~ X1 + X2 + Z, family = "binomial", data = full_data) 
  true_B_Y1_X1 <- coef(true_modelY1)["X1"]
  true_B_Y1_X2 <- coef(true_modelY1)["X2"]
  true_B_Y2_X1 <- coef(true_modelY2)["X1"]
  true_B_Y2_X2 <- coef(true_modelY2)["X2"]
  
  # compute true influence functions
  full_data$infl_Y1_X1_true <- inf_fun_logit(true_modelY1)[,"X1"]
  full_data$infl_Y1_X2_true <- inf_fun_logit(true_modelY1)[,"X2"]
  full_data$infl_Y2_X1_true <- inf_fun_logit(true_modelY2)[,"X1"]
  full_data$infl_Y2_X2_true <- inf_fun_logit(true_modelY2)[,"X2"]
  
  
  #####
  ## Generate error-prone Phase 1 data
  ######
  X1X2ZY1Y2_obs <- mvrnorm(N, c(0,0,0,0,0),
                           matrix(c(error_varX1, 0.05, 0.04, 0.15, 0.04,
                                    0.05, error_varX2, 0.02, 0.10, 0.4,
                                    0.04, 0.02, 1, 0.01, 0,
                                    0.15, 0.10, 0.01, 1 , 0.15,
                                    0.04, 0.4, 0, 0.15, 1), nrow = 5))
  full_data$X1_obs <- full_data$X1 + X1X2ZY1Y2_obs[,1]
  full_data$X2_obs <- full_data$X2 + X1X2ZY1Y2_obs[,2]
  
  threshold_positiveZ <- qnorm((sensZ + 1)/2)
  threshold_negativeZ <- qnorm((specZ + 1)/2)
  threshold_positiveY1 <- qnorm((sensY1 + 1)/2)
  threshold_negativeY1 <- qnorm((specY1 + 1)/2)
  threshold_positiveY2 <- qnorm((sensY2 + 1)/2)
  threshold_negativeY2 <- qnorm((specY2 + 1)/2)
  
  full_data$Z_obs <- with(full_data, ifelse((Z > 0 & abs(X1X2ZY1Y2_obs[, 3]) < threshold_positiveZ) | 
                                              (Z < 0 & abs(X1X2ZY1Y2_obs[, 3]) > threshold_negativeZ), 
                                            1, 0))
  
  full_data$Y1_obs <- with(full_data, ifelse((Y1 == 1 & abs(X1X2ZY1Y2_obs[, 4]) < threshold_positiveY1) | 
                                               (Y1 == 0 & abs(X1X2ZY1Y2_obs[, 4]) > threshold_negativeY1), 1, 0))
  
  full_data$Y2_obs <- with(full_data, ifelse((Y2 == 1 & abs(X1X2ZY1Y2_obs[, 5]) < threshold_positiveY2) | 
                                               (Y2 == 0 & abs(X1X2ZY1Y2_obs[, 5]) > threshold_negativeY2),1, 0))
  
  
  #### Compute naive influence functions
  fitY1_phase1 <-  glm(Y1_obs ~ X1_obs + X2_obs + Z_obs + Y2_obs, 
                       family = "binomial", data = full_data)
  fitY2_phase1 <-  glm(Y2_obs ~ X1_obs + X2_obs + Z_obs, 
                       family = "binomial", data = full_data)
  
  full_data$inflB11_phase1 <- inf_fun_logit(fitY1_phase1)[,"X1_obs"]
  full_data$inflB21_phase1 <- inf_fun_logit(fitY1_phase1)[,"X2_obs"]
  full_data$inflB12_phase1 <- inf_fun_logit(fitY2_phase1)[,"X1_obs"]
  full_data$inflB22_phase1 <- inf_fun_logit(fitY2_phase1)[,"X2_obs"]
  
  #### Compute naive influence functions
  fitY1_phase1 <-  glm(Y1_obs ~ X1_obs + X2_obs + Z_obs + Y2_obs, 
                       family = "binomial", data = full_data)
  fitY2_phase1 <-  glm(Y2_obs ~ X1_obs + X2_obs + Z_obs, 
                       family = "binomial", data = full_data)
  
  full_data$inflB11_phase1 <- inf_fun_logit(fitY1_phase1)[,"X1_obs"]
  full_data$inflB12_phase1 <- inf_fun_logit(fitY1_phase1)[,"X2_obs"]
  full_data$inflB21_phase1 <- inf_fun_logit(fitY2_phase1)[,"X1_obs"]
  full_data$inflB22_phase1 <- inf_fun_logit(fitY2_phase1)[,"X2_obs"]
  
  #### Divide into strata based on observed phase 1 data X
  full_data <- full_data |>
    dplyr::mutate(Y1_strat = Y1_obs, Y2_strat = Y2_obs) |>
    split_strata(strata = c("Y1_strat", "Y2_strat"),
                 split_var = "X1_obs",
                 type = "local quantile",
                 split_at = c(0.5),
                 trunc = "X1") |>
    split_strata(strata = "new_strata",
                 split_var = "X2_obs",
                 type = "local quantile",
                 split_at = c(0.5),
                 trunc = "X2")
  
  names(full_data)[names(full_data) == "new_strata"] <- "strata"
  
  data <- full_data
  
  #####
  #####
  ##### Phase 2
  
  ## Here, each wave chooses n/4 samples where n/16 are chosen to be optimal for
  ## B11, n/16 are chosen to be optimal for B21, n/16 to be optimal for B12,
  ## and n/16 are chosen to be optimal for B22. 
  
  ####
  ## Wave 1
  ####
  
  # Step 1: Initialize phase1_data
  phase1_data <- full_data
  
  # Step 2: Determine optimum allocation for wave 1 with Wright algorithm
  allocation1 <- optimum_allocation(phase1_data,
                                    strata = "strata",
                                    y = "inflB11_phase1",
                                    nsample = n/16, method = "WrightII")
  allocation2 <- optimum_allocation(phase1_data,
                                    strata = "strata",
                                    y = "inflB21_phase1",
                                    nsample = n/16, method = "WrightII")
  allocation3 <- optimum_allocation(phase1_data,
                                    strata = "strata",
                                    y = "inflB12_phase1",
                                    nsample = n/16, method = "WrightII")
  allocation4 <- optimum_allocation(phase1_data,
                                    strata = "strata",
                                    y = "inflB22_phase1",
                                    nsample = n/16, method = "WrightII")
  
  # Step 5: Combine allocations for total allocation of n/4
  wave1_allocation <- dplyr::left_join(allocation1, allocation2, 
                                       by = "strata") |>
    dplyr::mutate(stratum_size = stratum_size.x + stratum_size.y) |>
    dplyr::left_join(allocation3,
                     by = "strata") |>
    dplyr::mutate(stratum_size = stratum_size.x.x + stratum_size.y.y) |>
    dplyr::left_join(allocation4,
                     by = "strata") |>
    dplyr::mutate(stratum_size = stratum_size.x.x.x + stratum_size.y.y.y) 
  
  
  # First 1/4th according to wave1_allocation
  phase2_wave1 <- sample_strata(data = phase1_data,
                                strata = "strata",
                                id = "id",
                                design_data = wave1_allocation,
                                n_allocated = "stratum_size")
  names(phase2_wave1)[names(phase2_wave1) == "sample_indicator"] <- "sampled_wave1"
  
  # Sample
  phase2_wave1$X1 <- ifelse(phase2_wave1$sampled_wave1 == 1, full_data$X1 , NA)
  phase2_wave1$X2 <- ifelse(phase2_wave1$sampled_wave1 == 1, full_data$X2 , NA)
  phase2_wave1$Z <- ifelse(phase2_wave1$sampled_wave1 == 1, full_data$Z , NA)
  
  #####
  ## Wave 2
  #####
  
  # Update estimates for influence functions using generalized raking
  twophase_design_Y1 <- twophase(id = list(~id, ~id), 
                                 strata = list(NULL, ~strata), 
                                 subset = ~as.logical(sampled_wave1), 
                                 data = phase2_wave1,
                                 method = "simple")
  
  twophase_design_Y2 <- twophase(id = list(~id, ~id), 
                                 strata = list(NULL, ~strata), 
                                 subset = ~as.logical(sampled_wave1), 
                                 data = phase2_wave1,
                                 method = "simple")
  
  # Calibrate
  calibrated_twophase_Y1_wave1 <- calibrate(twophase_design_Y1,
                                            ~ inflB11_phase1 +
                                              inflB21_phase1 +
                                              inflB12_phase1 +
                                              inflB22_phase1 +
                                              strata,
                                            phase = 2,
                                            calfun = "raking")
  calibrated_twophase_Y2_wave1 <- calibrate(twophase_design_Y2,
                                            ~ inflB11_phase1 +
                                              inflB21_phase1 +
                                              inflB12_phase1 +
                                              inflB22_phase1 +
                                              strata,
                                            phase = 2,
                                            calfun = "raking")
  
  # Run models
  fitY1_wave1 <- svyglm(Y1 ~ X1 + X2 + Z + Y2, family = quasibinomial,
                        design = calibrated_twophase_Y1_wave1)
  fitY2_wave1 <- svyglm(Y2 ~ X1 + X2 + Z, family = quasibinomial, 
                        design = calibrated_twophase_Y2_wave1)
  
  # Get IFs
  infl_Y1_wave1 <- inf_fun_logit(fitY1_wave1)
  infl_Y2_wave1 <- inf_fun_logit(fitY2_wave1)
  infl_wave1 <- data.frame(as.numeric(rownames(infl_Y1_wave1)), 
                           infl_Y1_wave1[,"X1"],
                           infl_Y1_wave1[,"X2"],
                           infl_Y2_wave1[,"X1"],
                           infl_Y2_wave1[,"X2"])
  names(infl_wave1) <- c("id", "inflB11", "inflB21", "inflB12", "inflB22")
  phase2_wave1 <- dplyr::left_join(phase2_wave1, infl_wave1, by = "id")
  
  # Regress Latest IFs (computed above) on Phase 1 IFs
  resid_model_B11_wave1 <- lm(phase2_wave1$inflB11 ~ phase2_wave1$inflB11_phase1,
                              na.action = na.exclude)
  resid_model_B21_wave1 <- lm(phase2_wave1$inflB21 ~ phase2_wave1$inflB21_phase1,
                              na.action = na.exclude)
  resid_model_B12_wave1 <- lm(phase2_wave1$inflB12 ~ phase2_wave1$inflB12_phase1,
                              na.action = na.exclude)
  resid_model_B22_wave1 <- lm(phase2_wave1$inflB22 ~ phase2_wave1$inflB22_phase1,
                              na.action = na.exclude)
  
  # Get residuals
  phase2_wave1$residB11_wave1 <- resid(resid_model_B11_wave1)
  phase2_wave1$residB21_wave1 <- resid(resid_model_B21_wave1)
  phase2_wave1$residB12_wave1 <- resid(resid_model_B12_wave1)
  phase2_wave1$residB22_wave1 <- resid(resid_model_B22_wave1)
  
  
  # Re-calculate allocations
  allocation1 <- allocate_wave(phase2_wave1,
                               strata = "strata",
                               y = "residB11_wave1", method = "iterative",
                               already_sampled = "sampled_wave1",
                               nsample = n/16, 
                               allocation_method = "Neyman")
  allocation2 <- allocate_wave(phase2_wave1,
                               strata = "strata",
                               y = "residB21_wave1", method = "iterative",
                               already_sampled = "sampled_wave1",
                               nsample = n/16, 
                               allocation_method = "Neyman")
  allocation3 <- allocate_wave(phase2_wave1,
                               strata = "strata",
                               y = "residB12_wave1",
                               already_sampled = "sampled_wave1",
                               nsample = n/16, method = "iterative",
                               allocation_method = "Neyman")
  allocation4 <- allocate_wave(phase2_wave1,
                               strata = "strata",
                               y = "residB22_wave1",
                               already_sampled = "sampled_wave1",
                               nsample = n/16, method = "iterative",
                               allocation_method = "Neyman")
  
  # and combine
  wave2_allocation <- dplyr::left_join(allocation1, allocation2, 
                                       by = "strata") |>
    dplyr::mutate(n_to_sample = n_to_sample.x + n_to_sample.y) |>
    dplyr::left_join(allocation3,
                     by = "strata") |>
    dplyr::mutate(n_to_sample = n_to_sample.x.x + n_to_sample.y.y) |>
    dplyr::left_join(allocation4,
                     by = "strata") |>
    dplyr::mutate(n_to_sample = n_to_sample.x.x.x + n_to_sample.y.y.y)  
  
  
  # sample and merge data
  phase2_wave2 <- sample_strata(data = phase2_wave1,
                                strata = "strata",
                                id = "id",
                                already_sampled = "sampled_wave1",
                                design_data = wave2_allocation,
                                n_allocated = "n_to_sample")
  names(phase2_wave2)[names(phase2_wave2) == "sample_indicator"] <- "sampled_wave2"
  
  # Sample
  phase2_wave2$X1 <- ifelse(phase2_wave2$sampled_wave1 == 1|
                              phase2_wave2$sampled_wave2 == 1, full_data$X1 , NA)
  phase2_wave2$X2 <- ifelse(phase2_wave2$sampled_wave1 == 1|
                              phase2_wave2$sampled_wave2 == 1, full_data$X2 , NA)
  phase2_wave2$Z <- ifelse(phase2_wave2$sampled_wave1 == 1|
                             phase2_wave2$sampled_wave2 == 1, full_data$Z , NA)
  
  
  # Also, remove inflB11 and inflB12 vars
  phase2_wave2 <- subset(phase2_wave2, select = -c(inflB11, inflB21,
                                                   inflB12, inflB22))
  
  #####
  ## Wave 3
  #####
  
  # Already sampled indicator
  phase2_wave2$already_sampled <- phase2_wave2$sampled_wave1 +
    phase2_wave2$sampled_wave2
  
  # Update estimates for influence functions using generalized raking
  twophase_design_Y1 <- twophase(id = list(~id, ~id), 
                                 strata = list(NULL, ~strata), 
                                 subset = ~as.logical(already_sampled), 
                                 data = phase2_wave2,
                                 method = "simple")
  
  twophase_design_Y2 <- twophase(id = list(~id, ~id), 
                                 strata = list(NULL, ~strata), 
                                 subset = ~as.logical(already_sampled), 
                                 data = phase2_wave2,
                                 method = "simple")
  
  # Calibrate
  calibrated_twophase_Y1_wave2 <- calibrate(twophase_design_Y1,
                                            ~ inflB11_phase1 +
                                              inflB21_phase1 +
                                              inflB12_phase1 +
                                              inflB22_phase1 +
                                              strata,
                                            phase = 2,
                                            calfun = "raking")
  calibrated_twophase_Y2_wave2 <- calibrate(twophase_design_Y2,
                                            ~ inflB11_phase1 +
                                              inflB21_phase1 +
                                              inflB12_phase1 +
                                              inflB22_phase1 +
                                              strata,
                                            phase = 2,
                                            calfun = "raking")
  
  # Run models
  fitY1_wave2 <- svyglm(Y1 ~ X1 + X2 + Z + Y2, family = quasibinomial,
                        design = calibrated_twophase_Y1_wave2)
  fitY2_wave2 <- svyglm(Y2 ~ X1 + X2 + Z, family = quasibinomial, 
                        design = calibrated_twophase_Y2_wave2)
  
  # Get IFs
  infl_Y1_wave2 <- inf_fun_logit(fitY1_wave2)
  infl_Y2_wave2 <- inf_fun_logit(fitY2_wave2)
  infl_wave2 <- data.frame(as.numeric(rownames(infl_Y1_wave2)), 
                           infl_Y1_wave2[,"X1"],
                           infl_Y1_wave2[,"X2"],
                           infl_Y2_wave2[,"X1"],
                           infl_Y2_wave2[,"X2"])
  names(infl_wave2) <- c("id", "inflB11", "inflB21", "inflB12", "inflB22")
  phase2_wave2 <- dplyr::left_join(phase2_wave2, infl_wave2, by = "id")
  
  # Regress Latest IFs (computed above) on Phase 1 IFs
  resid_model_B11_wave2 <- lm(phase2_wave2$inflB11 ~ phase2_wave2$inflB11_phase1,
                              na.action = na.exclude)
  resid_model_B21_wave2 <- lm(phase2_wave2$inflB21 ~ phase2_wave2$inflB21_phase1,
                              na.action = na.exclude)
  resid_model_B12_wave2 <- lm(phase2_wave2$inflB12 ~ phase2_wave2$inflB12_phase1,
                              na.action = na.exclude)
  resid_model_B22_wave2 <- lm(phase2_wave2$inflB22 ~ phase2_wave2$inflB22_phase1,
                              na.action = na.exclude)
  
  # Get residuals
  phase2_wave2$residB11_wave2 <- resid(resid_model_B11_wave2)
  phase2_wave2$residB21_wave2 <- resid(resid_model_B21_wave2)
  phase2_wave2$residB12_wave2 <- resid(resid_model_B12_wave2)
  phase2_wave2$residB22_wave2 <- resid(resid_model_B22_wave2)
  
  
  # Re-calculate allocations
  allocation1 <- allocate_wave(phase2_wave2,
                               strata = "strata",
                               y = "residB11_wave2",
                               method = "iterative",
                               nsample = n/16,
                               already_sampled = "already_sampled",
                               allocation_method = "Neyman")
  allocation2 <- allocate_wave(phase2_wave2,
                               strata = "strata",
                               y = "residB21_wave2",
                               method = "iterative",
                               nsample = n/16,
                               already_sampled = "already_sampled",
                               allocation_method = "Neyman")
  allocation3 <- allocate_wave(phase2_wave2,
                               strata = "strata",
                               y = "residB12_wave2",
                               method = "iterative",
                               nsample = n/16, 
                               already_sampled = "already_sampled",
                               allocation_method = "Neyman")
  allocation4 <- allocate_wave(phase2_wave2,
                               strata = "strata",
                               y = "residB22_wave2",
                               method = "iterative",
                               nsample = n/16, 
                               already_sampled = "already_sampled",
                               allocation_method = "Neyman")
  
  # and combine
  wave3_allocation <- dplyr::left_join(allocation1, allocation2, 
                                       by = "strata") |>
    dplyr::mutate(n_to_sample = n_to_sample.x + n_to_sample.y) |>
    dplyr::left_join(allocation3,
                     by = "strata") |>
    dplyr::mutate(n_to_sample = n_to_sample.x.x + n_to_sample.y.y) |>
    dplyr::left_join(allocation4,
                     by = "strata") |>
    dplyr::mutate(n_to_sample = n_to_sample.x.x.x + n_to_sample.y.y.y)   
  
  
  # sample and merge data
  phase2_wave3 <- sample_strata(data = phase2_wave2,
                                strata = "strata",
                                id = "id",
                                already_sampled = "already_sampled",
                                design_data = wave3_allocation,
                                n_allocated = "n_to_sample")
  names(phase2_wave3)[names(phase2_wave3) == "sample_indicator"] <- "sampled_wave3"
  
  # Sample
  phase2_wave3$X1 <- ifelse(phase2_wave3$sampled_wave1 == 1|
                              phase2_wave3$sampled_wave2 == 1|
                              phase2_wave3$sampled_wave3 == 1, full_data$X1 , NA)
  phase2_wave3$X2 <- ifelse(phase2_wave3$sampled_wave1 == 1|
                              phase2_wave3$sampled_wave2 == 1|
                              phase2_wave3$sampled_wave3 == 1, full_data$X2 , NA)
  phase2_wave3$Z <- ifelse(phase2_wave3$sampled_wave1 == 1|
                             phase2_wave3$sampled_wave2 == 1|
                             phase2_wave3$sampled_wave3 == 1, full_data$Z , NA)
  
  
  # Also, remove inflB11 and inflB12 vars
  phase2_wave3 <- subset(phase2_wave3, select = -c(inflB11, inflB21,
                                                   inflB12, inflB22))
  
  #####
  ## Wave 4
  #####
  
  # Already sampled indicator
  phase2_wave3$already_sampled <- phase2_wave3$sampled_wave1 +
    phase2_wave3$sampled_wave2 + phase2_wave3$sampled_wave3
  
  # Update estimates for influence functions using generalized raking
  twophase_design_Y1 <- twophase(id = list(~id, ~id), 
                                 strata = list(NULL, ~strata), 
                                 subset = ~as.logical(already_sampled), 
                                 data = phase2_wave3,
                                 method = "simple")
  
  twophase_design_Y2 <- twophase(id = list(~id, ~id), 
                                 strata = list(NULL, ~strata), 
                                 subset = ~as.logical(already_sampled), 
                                 data = phase2_wave3,
                                 method = "simple")
  
  # Calibrate
  calibrated_twophase_Y1_wave3 <- calibrate(twophase_design_Y1,
                                            ~ inflB11_phase1 +
                                              inflB21_phase1 +
                                              inflB12_phase1 +
                                              inflB22_phase1 +
                                              strata,
                                            phase = 2,
                                            calfun = "raking")
  calibrated_twophase_Y2_wave3 <- calibrate(twophase_design_Y2,
                                            ~ inflB11_phase1 +
                                              inflB21_phase1 +
                                              inflB12_phase1 +
                                              inflB22_phase1 +
                                              strata,
                                            phase = 2,
                                            calfun = "raking")
  
  # Run models
  fitY1_wave3 <- svyglm(Y1 ~ X1 + X2 + Z + Y2, family = quasibinomial,
                        design = calibrated_twophase_Y1_wave3)
  fitY2_wave3 <- svyglm(Y2 ~ X1 + X2 + Z, family = quasibinomial, 
                        design = calibrated_twophase_Y2_wave3)
  
  # Get IFs
  infl_Y1_wave3 <- inf_fun_logit(fitY1_wave3)
  infl_Y2_wave3 <- inf_fun_logit(fitY2_wave3)
  infl_wave3 <- data.frame(as.numeric(rownames(infl_Y1_wave3)), 
                           infl_Y1_wave3[,"X1"],
                           infl_Y1_wave3[,"X2"],
                           infl_Y2_wave3[,"X1"],
                           infl_Y2_wave3[,"X2"])
  names(infl_wave3) <- c("id", "inflB11", "inflB21", "inflB12", "inflB22")
  phase2_wave3 <- dplyr::left_join(phase2_wave3, infl_wave3, by = "id")
  
  # Regress Latest IFs (computed above) on Phase 1 IFs
  resid_model_B11_wave3 <- lm(phase2_wave3$inflB11 ~ phase2_wave3$inflB11_phase1,
                              na.action = na.exclude)
  resid_model_B21_wave3 <- lm(phase2_wave3$inflB21 ~ phase2_wave3$inflB21_phase1,
                              na.action = na.exclude)
  resid_model_B12_wave3 <- lm(phase2_wave3$inflB12 ~ phase2_wave3$inflB12_phase1,
                              na.action = na.exclude)
  resid_model_B22_wave3 <- lm(phase2_wave3$inflB22 ~ phase2_wave3$inflB22_phase1,
                              na.action = na.exclude)
  
  # Get residuals
  phase2_wave3$residB11_wave3 <- resid(resid_model_B11_wave3)
  phase2_wave3$residB21_wave3 <- resid(resid_model_B21_wave3)
  phase2_wave3$residB12_wave3 <- resid(resid_model_B12_wave3)
  phase2_wave3$residB22_wave3 <- resid(resid_model_B22_wave3)
  
  # Re-calculate allocations
  allocation1 <- allocate_wave(phase2_wave3,
                               strata = "strata",
                               y = "residB11_wave3",
                               method = "iterative",
                               nsample = n/16, 
                               already_sampled = "already_sampled",
                               detailed = TRUE,
                               allocation_method = "Neyman")
  allocation2 <- allocate_wave(phase2_wave3,
                               strata = "strata",
                               y = "residB21_wave3",
                               method = "iterative",
                               nsample = n/16, 
                               already_sampled = "already_sampled",
                               detailed = TRUE,
                               allocation_method = "Neyman")
  allocation3 <- allocate_wave(phase2_wave3,
                               strata = "strata",
                               y = "residB12_wave3",
                               method = "iterative",
                               nsample = n/16, 
                               already_sampled = "already_sampled",
                               detailed = TRUE,
                               allocation_method = "Neyman")
  allocation4 <- allocate_wave(phase2_wave3,
                               strata = "strata",
                               y = "residB22_wave3",
                               method = "iterative",
                               nsample = n/16, 
                               already_sampled = "already_sampled",
                               detailed = TRUE,
                               allocation_method = "Neyman")
  
  # Set indicators for oversampling
  oversampled_Y1 <- ifelse(all(allocation1$nsample_optimal == 
                                 allocation1$nsample_actual),
                           0, 1)
  oversampled_Y2 <- ifelse(all(allocation2$nsample_optimal == 
                                 allocation2$nsample_actual),
                           0, 1)
  
  # and combine
  wave4_allocation <- dplyr::left_join(allocation1, allocation2, 
                                       by = "strata") |>
    dplyr::mutate(n_to_sample = n_to_sample.x + n_to_sample.y) |>
    dplyr::left_join(allocation3,
                     by = "strata") |>
    dplyr::mutate(n_to_sample = n_to_sample.x.x + n_to_sample.y.y) |>
    dplyr::left_join(allocation4,
                     by = "strata") |>
    dplyr::mutate(n_to_sample = n_to_sample.x.x.x + n_to_sample.y.y.y)   
  
  # sample and merge data
  phase2_wave4 <- sample_strata(data = phase2_wave3,
                                strata = "strata",
                                id = "id",
                                already_sampled = "already_sampled",
                                design_data = wave4_allocation,
                                n_allocated = "n_to_sample")
  names(phase2_wave4)[names(phase2_wave4) == "sample_indicator"] <- "sampled_wave4"
  
  # Sample wave4
  phase2_wave4$X1 <- ifelse(phase2_wave4$sampled_wave1 == 1 |
                              phase2_wave4$sampled_wave2 == 1|
                              phase2_wave4$sampled_wave3 == 1|
                              phase2_wave4$sampled_wave4 == 1, full_data$X1 , NA)
  phase2_wave4$X2 <- ifelse(phase2_wave4$sampled_wave1 == 1|
                              phase2_wave4$sampled_wave2 == 1|
                              phase2_wave4$sampled_wave3 == 1|
                              phase2_wave4$sampled_wave4 == 1, full_data$X2 , NA)
  phase2_wave4$Z <- ifelse(phase2_wave4$sampled_wave1 == 1|
                             phase2_wave4$sampled_wave2 == 1|
                             phase2_wave4$sampled_wave3 == 1|
                             phase2_wave4$sampled_wave4 == 1, full_data$Z , NA)
  
  # Remove inflB11 and inflB12 vars
  phase2_wave4 <- subset(phase2_wave4, select = -c(inflB11, inflB21,
                                                   inflB12, inflB22))
  
  ####
  ## Sampling done: Now calculate the Beta estimates with raking 
  ## using the survey package
  ####
  
  phase2_wave4$already_sampled <- phase2_wave4$sampled_wave1 +
    phase2_wave4$sampled_wave2 + phase2_wave4$sampled_wave3 +
    phase2_wave4$sampled_wave4
  twophase_design <- twophase(id = list(~1, ~1), 
                              strata = list(NULL, ~strata), 
                              subset = ~as.logical(already_sampled), 
                              data = phase2_wave4)
  
  # Weights
  weightY1 <- svyglm(Y1 ~ X1 + X2 + Z + Y2, family = quasibinomial, 
                     design = twophase_design)
  weightY2 <- svyglm(Y2 ~ X1 + X2 + Z , family = quasibinomial, 
                     design = twophase_design)
  
  #########
  #### Raking
  #########
  data <- phase2_wave4
  
  # Creating a survey object and calibrating the weights
  mydesign <- survey::twophase(
    id = list(~1, ~1), subset = ~as.logical(already_sampled),
    strata = list(NULL, ~strata), data = data,
  )
  infcalY1 <- survey::calibrate(mydesign, formula =  ~inflB11_phase1 +
                                  inflB21_phase1 + inflB12_phase1 +
                                  inflB22_phase1 + strata, phase = 2, calfun = "raking")
  infcalY2 <- survey::calibrate(mydesign, formula =  ~inflB11_phase1 +
                                  inflB21_phase1 + inflB12_phase1 +
                                  inflB22_phase1 + strata, phase = 2, calfun = "raking")
  
  # Fitting the outcome model: conditional treatment effect of interest.
  fitY1 <- survey::svyglm(Y1 ~ X1 + X2 + Z + Y2, design = infcalY1, family = quasibinomial)
  fitY2 <- survey::svyglm(Y2 ~ X1 + X2 + Z, design = infcalY2, family = quasibinomial)
  
  # Check if truth is contained in interval
  confintIPWB11 <- confint(weightY1)["X1",]
  confintIPWB21 <- confint(weightY1)["X2",]
  confintIPWB12 <- confint(weightY2)["X1",]
  confintIPWB22 <- confint(weightY2)["X2",]
  confintGRB11 <- confint(fitY1)["X1",]
  confintGRB21 <- confint(fitY1)["X2",]
  confintGRB12 <- confint(fitY2)["X1",]
  confintGRB22 <- confint(fitY2)["X2",]
  coverIPWB11 <- ifelse(B11 >= confintIPWB11[1] & B11 <= confintIPWB11[2],1,0)
  coverIPWB21 <- ifelse(B21 >= confintIPWB21[1] & B21 <= confintIPWB21[2],1,0)
  coverIPWB12 <- ifelse(B12 >= confintIPWB12[1] & B12 <= confintIPWB12[2],1,0)
  coverIPWB22 <- ifelse(B22 >= confintIPWB22[1] & B22 <= confintIPWB22[2],1,0)
  coverGRB11 <- ifelse(B11 >= confintGRB11[1] & B11 <= confintGRB11[2],1,0)
  coverGRB21 <- ifelse(B21 >= confintGRB21[1] & B21 <= confintGRB21[2],1,0)
  coverGRB12 <- ifelse(B12 >= confintGRB12[1] & B12 <= confintGRB12[2],1,0)
  coverGRB22 <- ifelse(B22 >= confintGRB22[1] & B22 <= confintGRB22[2],1,0)
  
  output <- c(fitY1$coefficients["X1"], fitY1$coefficients["X2"],
              fitY2$coefficients["X1"], fitY2$coefficients["X2"],
              coef(weightY1)["X1"], coef(weightY1)["X2"], 
              coef(weightY2)["X1"], coef(weightY2)["X2"],  
              SE(fitY1)["X1"], SE(fitY1)["X2"],
              SE(fitY2)["X1"], SE(fitY2)["X2"], 
              SE(weightY1)["X1"], SE(weightY1)["X2"], 
              SE(weightY2)["X1"], SE(weightY2)["X2"], 
              coverIPWB11, coverIPWB21,
              coverIPWB12, coverIPWB22, 
              coverGRB11, coverGRB21,
              coverGRB12, coverGRB22,
              oversampled_Y1, oversampled_Y2, cor_Y1_Y2,
              prev_Y1, prev_Y2, B11, B21, B12, B22)
  return(output)
}

#####
## Strategy 3 - Independent sequential
#####

####
#### For this scenario I only consider one sequence (not all four orders). 
#### Arbitrarily, I let the sequence be B11, B21, B12, B22

# This is looped part, so Wave 1 allocation is required as input
Strategy3 <- function(n, simulations_df, scenario, N = 10000){
  
  # Set up betas for scenario
  B01 <- simulations_df["B01", scenario]
  B11 <- simulations_df["B11", scenario]
  B21 <- simulations_df["B21", scenario]
  B31 <- simulations_df["B31", scenario]
  B41 <- simulations_df["B41", scenario]
  B02 <- simulations_df["B02", scenario]
  B12 <- simulations_df["B12", scenario]
  B22 <- simulations_df["B22", scenario]
  B32 <- simulations_df["B32", scenario] 
  
  # Set up of Xs
  corX1X2 <- simulations_df["corX1X2", scenario]
  sensZ <- simulations_df["sensZ", scenario]
  specZ <- simulations_df["specZ", scenario]
  prevZ <- simulations_df["prevZ", scenario]
  
  # Measurement error
  sensY1 <- simulations_df["sensY1", scenario]
  specY1 <- simulations_df["specY1", scenario]
  sensY2 <- simulations_df["sensY2", scenario]
  specY2 <- simulations_df["specY2", scenario]
  error_varX1 <- simulations_df["error_varX1", scenario]
  error_varX2 <- simulations_df["error_varX2", scenario]
  
  
  #####
  ## Generate Phase 1 Data
  #####
  
  id <- 1:N
  
  # Generate correlated covariates
  sigma <- matrix(c(1, corX1X2, 0.1, corX1X2, 1,  0.25, 0.1, 0.25, 1), nrow = 3)
  covs <- mvrnorm(N, mu = c(0,0,0), Sigma = sigma)
  X1 <- covs[,1]
  X2 <- covs[,2]
  Z <- ifelse(covs[,3] >= quantile(covs[,3], 1 - prevZ), 1, 0)
  cor(cbind(X1, X2, Z)) 
  
  Y2_probs <- exp(B02 + B12*X1 + B22*X2 + B32*Z)/(1 + exp(B02 + B12*X1 + B22*X2 + B32*Z))
  Y2 <- rbinom(N, 1, Y2_probs) # Compute realized value for Y2
  Y1_probs <- exp(B01 + B11*X1 + B21*X2 + B31*Z + B41*Y2)/
    (1 + exp(B01 + B11*X1 + B21*X2 + B31*Z + B41*Y2))
  Y1 <- rbinom(N, 1, Y1_probs) # Compute realized value for Y1
  
  full_data <- data.frame(id, X1, X2, Z, Y1, Y2)
  
  ## Test true model
  true_model1 <- glm(Y1 ~ X1 + X2 + Z + Y2, family = "binomial", data = full_data)
  true_model2 <- glm(Y2 ~ X1 + X2 + Z, family = "binomial", data = full_data)
  
  # Find observed correlation
  cor_Y1_Y2 <- cor(Y1, Y2) 
  
  # Prevalence
  prev_Y1 <- table(Y1)["1"]/length(Y1)
  prev_Y2 <- table(Y2)["1"]/length(Y2)
  
  # Find true coefficients
  true_modelY1 <- glm(Y1 ~ X1 + X2 + Z + Y2, family = "binomial", data = full_data) 
  true_modelY2 <- glm(Y2 ~ X1 + X2 + Z, family = "binomial", data = full_data) 
  true_B_Y1_X1 <- coef(true_modelY1)["X1"]
  true_B_Y1_X2 <- coef(true_modelY1)["X2"]
  true_B_Y2_X1 <- coef(true_modelY2)["X1"]
  true_B_Y2_X2 <- coef(true_modelY2)["X2"]
  
  # compute true influence functions
  full_data$infl_Y1_X1_true <- inf_fun_logit(true_modelY1)[,"X1"]
  full_data$infl_Y1_X2_true <- inf_fun_logit(true_modelY1)[,"X2"]
  full_data$infl_Y2_X1_true <- inf_fun_logit(true_modelY2)[,"X1"]
  full_data$infl_Y2_X2_true <- inf_fun_logit(true_modelY2)[,"X2"]
  
  
  #####
  ## Generate error-prone Phase 1 data
  ######
  X1X2ZY1Y2_obs <- mvrnorm(N, c(0,0,0,0,0),
                           matrix(c(error_varX1, 0.05, 0.04, 0.15, 0.04,
                                    0.05, error_varX2, 0.02, 0.10, 0.4,
                                    0.04, 0.02, 1, 0.01, 0,
                                    0.15, 0.10, 0.01, 1 , 0.15,
                                    0.04, 0.4, 0, 0.15, 1), nrow = 5))
  full_data$X1_obs <- full_data$X1 + X1X2ZY1Y2_obs[,1]
  full_data$X2_obs <- full_data$X2 + X1X2ZY1Y2_obs[,2]
  
  threshold_positiveZ <- qnorm((sensZ + 1)/2)
  threshold_negativeZ <- qnorm((specZ + 1)/2)
  threshold_positiveY1 <- qnorm((sensY1 + 1)/2)
  threshold_negativeY1 <- qnorm((specY1 + 1)/2)
  threshold_positiveY2 <- qnorm((sensY2 + 1)/2)
  threshold_negativeY2 <- qnorm((specY2 + 1)/2)
  
  full_data$Z_obs <- with(full_data, ifelse((Z > 0 & abs(X1X2ZY1Y2_obs[, 3]) < threshold_positiveZ) | 
                                              (Z < 0 & abs(X1X2ZY1Y2_obs[, 3]) > threshold_negativeZ), 
                                            1, 0))
  
  full_data$Y1_obs <- with(full_data, ifelse((Y1 == 1 & abs(X1X2ZY1Y2_obs[, 4]) < threshold_positiveY1) | 
                                               (Y1 == 0 & abs(X1X2ZY1Y2_obs[, 4]) > threshold_negativeY1), 1, 0))
  
  full_data$Y2_obs <- with(full_data, ifelse((Y2 == 1 & abs(X1X2ZY1Y2_obs[, 5]) < threshold_positiveY2) | 
                                               (Y2 == 0 & abs(X1X2ZY1Y2_obs[, 5]) > threshold_negativeY2),1, 0))
  
  
  #### Compute naive influence functions
  fitY1_phase1 <-  glm(Y1_obs ~ X1_obs + X2_obs + Z_obs + Y2_obs, 
                       family = "binomial", data = full_data)
  fitY2_phase1 <-  glm(Y2_obs ~ X1_obs + X2_obs + Z_obs, 
                       family = "binomial", data = full_data)
  
  full_data$inflB11_phase1 <- inf_fun_logit(fitY1_phase1)[,"X1_obs"]
  full_data$inflB21_phase1 <- inf_fun_logit(fitY1_phase1)[,"X2_obs"]
  full_data$inflB12_phase1 <- inf_fun_logit(fitY2_phase1)[,"X1_obs"]
  full_data$inflB22_phase1 <- inf_fun_logit(fitY2_phase1)[,"X2_obs"]
  
  #### Compute naive influence functions
  fitY1_phase1 <-  glm(Y1_obs ~ X1_obs + X2_obs + Z_obs + Y2_obs, 
                       family = "binomial", data = full_data)
  fitY2_phase1 <-  glm(Y2_obs ~ X1_obs + X2_obs + Z_obs, 
                       family = "binomial", data = full_data)
  
  full_data$inflB11_phase1 <- inf_fun_logit(fitY1_phase1)[,"X1_obs"]
  full_data$inflB12_phase1 <- inf_fun_logit(fitY1_phase1)[,"X2_obs"]
  full_data$inflB21_phase1 <- inf_fun_logit(fitY2_phase1)[,"X1_obs"]
  full_data$inflB22_phase1 <- inf_fun_logit(fitY2_phase1)[,"X2_obs"]
  
  #### Divide into strata based on observed phase 1 data X
  full_data <- full_data |>
    dplyr::mutate(Y1_strat = Y1_obs, Y2_strat = Y2_obs) |>
    split_strata(strata = c("Y1_strat", "Y2_strat"),
                 split_var = "X1_obs",
                 type = "local quantile",
                 split_at = c(0.5),
                 trunc = "X1") |>
    split_strata(strata = "new_strata",
                 split_var = "X2_obs",
                 type = "local quantile",
                 split_at = c(0.5),
                 trunc = "X2")
  
  names(full_data)[names(full_data) == "new_strata"] <- "strata"
  
  data <- full_data
  
  #####
  #####
  ##### Phase 2
  
  ####
  ## Wave 1
  ####
  
  # Step 1: Initialize Phase 1 data
  phase1_data <- full_data
  
  # Step 2: Determine optimum allocation for wave 1 with Wright algorithm
  wave1_allocation <- optimum_allocation(phase1_data,
                                         strata = "strata",
                                         y = "inflB11_phase1",
                                         nsample = n/4, method = "Neyman")
  
  
  # Sample First n/4 according to wave1_allocation
  phase2_wave1 <- sample_strata(data = phase1_data,
                                strata = "strata",
                                id = "id",
                                design_data = wave1_allocation,
                                n_allocated = "stratum_size")
  names(phase2_wave1)[names(phase2_wave1) == "sample_indicator"] <- "sampled_wave1"
  
  #####
  ## Wave 2
  #####
  
  # Allocation optimal wrt B21
  
  # Update estimates for influence functions using generalized raking
  twophase_design_Y1 <- twophase(id = list(~id, ~id), 
                                 strata = list(NULL, ~strata), 
                                 subset = ~as.logical(sampled_wave1), 
                                 data = phase2_wave1,
                                 method = "simple")
  
  # Calibrate
  calibrated_twophase_Y1_wave1 <- calibrate(twophase_design_Y1,
                                            ~ inflB11_phase1 +
                                              inflB21_phase1 +
                                              inflB12_phase1 +
                                              inflB22_phase1 +
                                              strata,
                                            phase = 2,
                                            calfun = "raking")
  
  # Run models
  fitY1_wave1 <- svyglm(Y1 ~ X1 + X2 + Z + Y2, family = quasibinomial,
                        design = calibrated_twophase_Y1_wave1)
  
  # Get IFs
  infl_Y1_wave1 <- inf_fun_logit(fitY1_wave1)
  infl_wave1 <- data.frame(as.numeric(rownames(infl_Y1_wave1)), 
                           infl_Y1_wave1[,"X2"])
  names(infl_wave1) <- c("id", "inflB21")
  phase2_wave1 <- dplyr::left_join(phase2_wave1, infl_wave1, by = "id")
  
  # Regress Latest IFs (computed above) on Phase 1 IFs
  resid_model_B21_wave1 <- lm(phase2_wave1$inflB21 ~ phase2_wave1$inflB21_phase1,
                              na.action = na.exclude)
  
  # Get residuals
  phase2_wave1$residB21_wave1 <- resid(resid_model_B21_wave1)
  
  
  # Calculate allocation for Y1
  wave2_allocation <- allocate_wave(phase2_wave1,
                                    strata = "strata",
                                    y = "residB21_wave1",
                                    nsample = n/4, method = "iterative",
                                    already_sampled = "sampled_wave1",
                                    allocation_method = "Neyman")
  
  
  # sample and merge data
  phase2_wave2 <- sample_strata(data = phase2_wave1,
                                strata = "strata",
                                id = "id",
                                already_sampled = "sampled_wave1",
                                design_data = wave2_allocation,
                                n_allocated = "n_to_sample")
  names(phase2_wave2)[names(phase2_wave2) == "sample_indicator"] <- "sampled_wave2"
  
  # Sample
  phase2_wave2$X1 <- ifelse(phase2_wave2$sampled_wave1 == 1|
                              phase2_wave2$sampled_wave2 == 1, full_data$X1 , NA)
  phase2_wave2$X2 <- ifelse(phase2_wave2$sampled_wave1 == 1|
                              phase2_wave2$sampled_wave2 == 1, full_data$X2 , NA)
  phase2_wave2$Z <- ifelse(phase2_wave2$sampled_wave1 == 1|
                             phase2_wave2$sampled_wave2 == 1, full_data$Z , NA)
  
  phase2_wave2 <- subset(phase2_wave2, select = -c(inflB21))
  
  #####
  ## Wave 3
  #####
  
  # Allocation optimal wrt B12
  
  # Indicator for already sampled
  phase2_wave2$already_sampled <- phase2_wave2$sampled_wave1 +
    phase2_wave2$sampled_wave2
  
  # Update estimates for influence functions using generalized raking
  
  twophase_design_Y2 <- twophase(id = list(~id, ~id), 
                                 strata = list(NULL, ~strata), 
                                 subset = ~as.logical(already_sampled), 
                                 data = phase2_wave2,
                                 method = "simple")
  
  # Calibrate
  calibrated_twophase_Y2_wave2 <- calibrate(twophase_design_Y2,
                                            ~ inflB11_phase1 +
                                              inflB21_phase1 +
                                              inflB12_phase1 +
                                              inflB22_phase1 +
                                              strata,
                                            phase = 2,
                                            calfun = "raking")
  
  # Run models
  fitY2_wave2 <- svyglm(Y2 ~ X1 + X2 + Z, family = quasibinomial, 
                        design = calibrated_twophase_Y2_wave2)
  
  # Get IFs
  infl_Y2_wave2 <- inf_fun_logit(fitY2_wave2)
  infl_wave2 <- data.frame(as.numeric(rownames(infl_Y2_wave2)),
                           infl_Y2_wave2[,"X1"])
  names(infl_wave2) <- c("id", "inflB12")
  phase2_wave2 <- dplyr::left_join(phase2_wave2, infl_wave2, by = "id")
  
  # Regress Latest IFs (computed above) on Phase 1 IFs
  resid_model_B12_wave2 <- lm(phase2_wave2$inflB12 ~ phase2_wave2$inflB12_phase1,
                              na.action = na.exclude)
  
  # Get residuals
  phase2_wave2$residB12_wave2 <- resid(resid_model_B12_wave2)
  
  
  # Re-calculate allocation
  wave3_allocation <- allocate_wave(phase2_wave2,
                                    strata = "strata",
                                    y = "residB12_wave2",
                                    nsample = n/4, method = "iterative",
                                    already_sampled = "already_sampled",
                                    allocation_method = "Neyman")
  
  
  # sample and merge data
  phase2_wave3 <- sample_strata(data = phase2_wave2,
                                strata = "strata",
                                id = "id",
                                already_sampled = "already_sampled",
                                design_data = wave3_allocation,
                                n_allocated = "n_to_sample")
  names(phase2_wave3)[names(phase2_wave3) == "sample_indicator"] <- "sampled_wave3"
  
  # Sample
  phase2_wave3$X1 <- ifelse(phase2_wave3$sampled_wave1 == 1|
                              phase2_wave3$sampled_wave2 == 1|
                              phase2_wave3$sampled_wave3 == 1, full_data$X1 , NA)
  phase2_wave3$X2 <- ifelse(phase2_wave3$sampled_wave1 == 1|
                              phase2_wave3$sampled_wave2 == 1|
                              phase2_wave3$sampled_wave3 == 1, full_data$X2 , NA)
  phase2_wave3$Z <- ifelse(phase2_wave3$sampled_wave1 == 1|
                             phase2_wave3$sampled_wave2 == 1|
                             phase2_wave3$sampled_wave3 == 1, full_data$Z , NA)
  
  phase2_wave3 <- subset(phase2_wave3, select = -c(inflB12))
  
  
  #####
  ## Wave 4
  #####
  
  # Allocation optimal wrt B22
  
  # Indicator for already sampled
  phase2_wave3$already_sampled <- phase2_wave3$sampled_wave1 +
    phase2_wave3$sampled_wave2 + phase2_wave3$sampled_wave3
  
  # Update estimates for influence functions using generalized raking
  twophase_design_Y2 <- twophase(id = list(~id, ~id), 
                                 strata = list(NULL, ~strata), 
                                 subset = ~as.logical(already_sampled), 
                                 data = phase2_wave3,
                                 method = "simple")
  
  # Calibrate
  calibrated_twophase_Y2_wave3 <- calibrate(twophase_design_Y2,
                                            ~ inflB11_phase1 +
                                              inflB21_phase1 +
                                              inflB12_phase1 +
                                              inflB22_phase1 +
                                              strata,
                                            phase = 2,
                                            calfun = "raking")
  
  # Run models
  fitY2_wave3 <- svyglm(Y2 ~ X1 + X2 + Z, family = quasibinomial, 
                        design = calibrated_twophase_Y2_wave3)
  
  # Get IFs
  infl_Y2_wave3 <- inf_fun_logit(fitY2_wave3)
  infl_wave3 <- data.frame(as.numeric(rownames(infl_Y2_wave3)),
                           infl_Y2_wave3[,"X2"])
  names(infl_wave3) <- c("id", "inflB22")
  phase2_wave3 <- dplyr::left_join(phase2_wave3, infl_wave3, by = "id")
  
  # Regress Latest IFs (computed above) on Phase 1 IFs
  resid_model_B22_wave3 <- lm(phase2_wave3$inflB22 ~ phase2_wave3$inflB22_phase1,
                              na.action = na.exclude)
  
  # Get residuals
  phase2_wave3$residB22_wave3 <- resid(resid_model_B22_wave3)
  
  # Re-calculate allocation
  wave4_allocation <- allocate_wave(phase2_wave3,
                                    strata = "strata",
                                    y = "residB22_wave3",
                                    nsample = n/4, method = "iterative",
                                    already_sampled = "already_sampled",
                                    detailed = TRUE,
                                    allocation_method = "Neyman")
  
  # Check for oversampling
  oversampled_Y2 <- ifelse(all(wave4_allocation$nsample_optimal == 
                                 wave4_allocation$nsample_actual),
                           0, 1)
  # no reason to check for Y1 here, as optimality in later waves with respect to Y2
  
  # sample and merge data
  phase2_wave4 <- sample_strata(data = phase2_wave3,
                                strata = "strata",
                                id = "id",
                                already_sampled = "already_sampled",
                                design_data = wave4_allocation,
                                n_allocated = "n_to_sample")
  names(phase2_wave4)[names(phase2_wave4) == "sample_indicator"] <- "sampled_wave4"
  
  # Sample wave4
  phase2_wave4$X1 <- ifelse(phase2_wave4$sampled_wave1 == 1 |
                              phase2_wave4$sampled_wave2 == 1|
                              phase2_wave4$sampled_wave3 == 1|
                              phase2_wave4$sampled_wave4 == 1, full_data$X1 , NA)
  phase2_wave4$X2 <- ifelse(phase2_wave4$sampled_wave1 == 1|
                              phase2_wave4$sampled_wave2 == 1|
                              phase2_wave4$sampled_wave3 == 1|
                              phase2_wave4$sampled_wave4 == 1, full_data$X2 , NA)
  phase2_wave4$Z <- ifelse(phase2_wave4$sampled_wave1 == 1|
                             phase2_wave4$sampled_wave2 == 1|
                             phase2_wave4$sampled_wave3 == 1|
                             phase2_wave4$sampled_wave4 == 1, full_data$Z , NA)
  
  phase2_wave4 <- subset(phase2_wave4, select = -c(inflB22))
  
  ####
  ## Sampling done: Now calculate the Beta estimates with raking 
  ## using the survey package
  ####
  
  phase2_wave4$already_sampled <- phase2_wave4$sampled_wave1 +
    phase2_wave4$sampled_wave2 + phase2_wave4$sampled_wave3 +
    phase2_wave4$sampled_wave4
  twophase_design <- twophase(id = list(~1, ~1), 
                              strata = list(NULL, ~strata), 
                              subset = ~as.logical(already_sampled), 
                              data = phase2_wave4)
  
  # Weights
  weightY1 <- svyglm(Y1 ~ X1 + X2 + Z + Y2, family = quasibinomial, 
                     design = twophase_design)
  weightY2 <- svyglm(Y2 ~ X1 + X2 + Z , family = quasibinomial, 
                     design = twophase_design)
  
  
  #########
  #### Raking 
  #########
  data <- phase2_wave4
  
  # Creating a survey object and calibrating the weights
  mydesign <- survey::twophase(
    id = list(~1, ~1), subset = ~as.logical(already_sampled),
    strata = list(NULL, ~strata), data = data,
  )
  infcalY1 <- survey::calibrate(mydesign, formula =  ~inflB11_phase1 +
                                  inflB21_phase1 + inflB12_phase1 +
                                  inflB22_phase1 + strata, phase = 2, calfun = "raking")
  infcalY2 <- survey::calibrate(mydesign, formula =  ~inflB11_phase1 +
                                  inflB21_phase1 + inflB12_phase1 +
                                  inflB22_phase1 + strata, phase = 2, calfun = "raking")
  
  # Fitting the outcome model: conditional treatment effect of interest.
  fitY1 <- survey::svyglm(Y1 ~ X1 + X2 + Z + Y2, design = infcalY1, family = quasibinomial)
  fitY2 <- survey::svyglm(Y2 ~ X1 + X2 + Z, design = infcalY2, family = quasibinomial)
  
  
  # Check if truth is contained in interval
  confintIPWB11 <- confint(weightY1)["X1",]
  confintIPWB21 <- confint(weightY1)["X2",]
  confintIPWB12 <- confint(weightY2)["X1",]
  confintIPWB22 <- confint(weightY2)["X2",]
  confintGRB11 <- confint(fitY1)["X1",]
  confintGRB21 <- confint(fitY1)["X2",]
  confintGRB12 <- confint(fitY2)["X1",]
  confintGRB22 <- confint(fitY2)["X2",]
  coverIPWB11 <- ifelse(B11 >= confintIPWB11[1] & B11 <= confintIPWB11[2],1,0)
  coverIPWB21 <- ifelse(B21 >= confintIPWB21[1] & B21 <= confintIPWB21[2],1,0)
  coverIPWB12 <- ifelse(B12 >= confintIPWB12[1] & B12 <= confintIPWB12[2],1,0)
  coverIPWB22 <- ifelse(B22 >= confintIPWB22[1] & B22 <= confintIPWB22[2],1,0)
  coverGRB11 <- ifelse(B11 >= confintGRB11[1] & B11 <= confintGRB11[2],1,0)
  coverGRB21 <- ifelse(B21 >= confintGRB21[1] & B21 <= confintGRB21[2],1,0)
  coverGRB12 <- ifelse(B12 >= confintGRB12[1] & B12 <= confintGRB12[2],1,0)
  coverGRB22 <- ifelse(B22 >= confintGRB22[1] & B22 <= confintGRB22[2],1,0)
  
  output <- c(fitY1$coefficients["X1"], fitY1$coefficients["X2"],
              fitY2$coefficients["X1"], fitY2$coefficients["X2"],
              coef(weightY1)["X1"], coef(weightY1)["X2"], 
              coef(weightY2)["X1"], coef(weightY2)["X2"],  
              SE(fitY1)["X1"], SE(fitY1)["X2"],
              SE(fitY2)["X1"], SE(fitY2)["X2"], 
              SE(weightY1)["X1"], SE(weightY1)["X2"], 
              SE(weightY2)["X1"], SE(weightY2)["X2"], 
              coverIPWB11, coverIPWB21,
              coverIPWB12, coverIPWB22, 
              coverGRB11, coverGRB21,
              coverGRB12,coverGRB22,
              oversampled_Y2,
              cor_Y1_Y2,
              prev_Y1, prev_Y2, B11, B21, B12, B22)
  return(output)
}

######
### Strategy 4 - A-optimal sampling
#####

# Step 5: Create function to estimate beta-hat at each iteration
Strategy4 <- function(n, simulations_df, scenario, N = 10000){
  
  # Set up betas for scenario
  B01 <- simulations_df["B01", scenario]
  B11 <- simulations_df["B11", scenario]
  B21 <- simulations_df["B21", scenario]
  B31 <- simulations_df["B31", scenario]
  B41 <- simulations_df["B41", scenario]
  B02 <- simulations_df["B02", scenario]
  B12 <- simulations_df["B12", scenario]
  B22 <- simulations_df["B22", scenario]
  B32 <- simulations_df["B32", scenario] 
  
  # Set up of Xs
  corX1X2 <- simulations_df["corX1X2", scenario]
  sensZ <- simulations_df["sensZ", scenario]
  specZ <- simulations_df["specZ", scenario]
  prevZ <- simulations_df["prevZ", scenario]
  
  # Measurement error
  sensY1 <- simulations_df["sensY1", scenario]
  specY1 <- simulations_df["specY1", scenario]
  sensY2 <- simulations_df["sensY2", scenario]
  specY2 <- simulations_df["specY2", scenario]
  error_varX1 <- simulations_df["error_varX1", scenario]
  error_varX2 <- simulations_df["error_varX2", scenario]
  
  
  #####
  ## Generate Phase 1 Data
  #####
  
  id <- 1:N
  
  # Generate correlated covariates
  sigma <- matrix(c(1, corX1X2, 0.1, corX1X2, 1,  0.25, 0.1, 0.25, 1), nrow = 3)
  covs <- mvrnorm(N, mu = c(0,0,0), Sigma = sigma)
  X1 <- covs[,1]
  X2 <- covs[,2]
  Z <- ifelse(covs[,3] >= quantile(covs[,3], 1 - prevZ), 1, 0)
  cor(cbind(X1, X2, Z)) 
  
  Y2_probs <- exp(B02 + B12*X1 + B22*X2 + B32*Z)/(1 + exp(B02 + B12*X1 + B22*X2 + B32*Z))
  Y2 <- rbinom(N, 1, Y2_probs) # Compute realized value for Y2
  Y1_probs <- exp(B01 + B11*X1 + B21*X2 + B31*Z + B41*Y2)/
    (1 + exp(B01 + B11*X1 + B21*X2 + B31*Z + B41*Y2))
  Y1 <- rbinom(N, 1, Y1_probs) # Compute realized value for Y1
  
  full_data <- data.frame(id, X1, X2, Z, Y1, Y2)
  
  ## Test true model
  true_model1 <- glm(Y1 ~ X1 + X2 + Z + Y2, family = "binomial", data = full_data)
  true_model2 <- glm(Y2 ~ X1 + X2 + Z, family = "binomial", data = full_data)
  
  # Find observed correlation
  cor_Y1_Y2 <- cor(Y1, Y2) 
  
  # Prevalence
  prev_Y1 <- table(Y1)["1"]/length(Y1)
  prev_Y2 <- table(Y2)["1"]/length(Y2)
  
  # Find true coefficients
  true_modelY1 <- glm(Y1 ~ X1 + X2 + Z + Y2, family = "binomial", data = full_data) 
  true_modelY2 <- glm(Y2 ~ X1 + X2 + Z, family = "binomial", data = full_data) 
  true_B_Y1_X1 <- coef(true_modelY1)["X1"]
  true_B_Y1_X2 <- coef(true_modelY1)["X2"]
  true_B_Y2_X1 <- coef(true_modelY2)["X1"]
  true_B_Y2_X2 <- coef(true_modelY2)["X2"]
  
  # compute true influence functions
  full_data$infl_Y1_X1_true <- inf_fun_logit(true_modelY1)[,"X1"]
  full_data$infl_Y1_X2_true <- inf_fun_logit(true_modelY1)[,"X2"]
  full_data$infl_Y2_X1_true <- inf_fun_logit(true_modelY2)[,"X1"]
  full_data$infl_Y2_X2_true <- inf_fun_logit(true_modelY2)[,"X2"]
  
  
  #####
  ## Generate error-prone Phase 1 data
  ######
  X1X2ZY1Y2_obs <- mvrnorm(N, c(0,0,0,0,0),
                           matrix(c(error_varX1, 0.05, 0.04, 0.15, 0.04,
                                    0.05, error_varX2, 0.02, 0.10, 0.4,
                                    0.04, 0.02, 1, 0.01, 0,
                                    0.15, 0.10, 0.01, 1 , 0.15,
                                    0.04, 0.4, 0, 0.15, 1), nrow = 5))
  full_data$X1_obs <- full_data$X1 + X1X2ZY1Y2_obs[,1]
  full_data$X2_obs <- full_data$X2 + X1X2ZY1Y2_obs[,2]
  
  threshold_positiveZ <- qnorm((sensZ + 1)/2)
  threshold_negativeZ <- qnorm((specZ + 1)/2)
  threshold_positiveY1 <- qnorm((sensY1 + 1)/2)
  threshold_negativeY1 <- qnorm((specY1 + 1)/2)
  threshold_positiveY2 <- qnorm((sensY2 + 1)/2)
  threshold_negativeY2 <- qnorm((specY2 + 1)/2)
  
  full_data$Z_obs <- with(full_data, ifelse((Z > 0 & abs(X1X2ZY1Y2_obs[, 3]) < threshold_positiveZ) | 
                                              (Z < 0 & abs(X1X2ZY1Y2_obs[, 3]) > threshold_negativeZ), 
                                            1, 0))
  
  full_data$Y1_obs <- with(full_data, ifelse((Y1 == 1 & abs(X1X2ZY1Y2_obs[, 4]) < threshold_positiveY1) | 
                                               (Y1 == 0 & abs(X1X2ZY1Y2_obs[, 4]) > threshold_negativeY1), 1, 0))
  
  full_data$Y2_obs <- with(full_data, ifelse((Y2 == 1 & abs(X1X2ZY1Y2_obs[, 5]) < threshold_positiveY2) | 
                                               (Y2 == 0 & abs(X1X2ZY1Y2_obs[, 5]) > threshold_negativeY2),1, 0))
  
  
  #### Compute naive influence functions
  fitY1_phase1 <-  glm(Y1_obs ~ X1_obs + X2_obs + Z_obs + Y2_obs, 
                       family = "binomial", data = full_data)
  fitY2_phase1 <-  glm(Y2_obs ~ X1_obs + X2_obs + Z_obs, 
                       family = "binomial", data = full_data)
  
  full_data$inflB11_phase1 <- inf_fun_logit(fitY1_phase1)[,"X1_obs"]
  full_data$inflB21_phase1 <- inf_fun_logit(fitY1_phase1)[,"X2_obs"]
  full_data$inflB12_phase1 <- inf_fun_logit(fitY2_phase1)[,"X1_obs"]
  full_data$inflB22_phase1 <- inf_fun_logit(fitY2_phase1)[,"X2_obs"]
  
  #### Compute naive influence functions
  fitY1_phase1 <-  glm(Y1_obs ~ X1_obs + X2_obs + Z_obs + Y2_obs, 
                       family = "binomial", data = full_data)
  fitY2_phase1 <-  glm(Y2_obs ~ X1_obs + X2_obs + Z_obs, 
                       family = "binomial", data = full_data)
  
  full_data$inflB11_phase1 <- inf_fun_logit(fitY1_phase1)[,"X1_obs"]
  full_data$inflB12_phase1 <- inf_fun_logit(fitY1_phase1)[,"X2_obs"]
  full_data$inflB21_phase1 <- inf_fun_logit(fitY2_phase1)[,"X1_obs"]
  full_data$inflB22_phase1 <- inf_fun_logit(fitY2_phase1)[,"X2_obs"]
  
  #### Divide into strata based on observed phase 1 data X
  full_data <- full_data |>
    dplyr::mutate(Y1_strat = Y1_obs, Y2_strat = Y2_obs) |>
    split_strata(strata = c("Y1_strat", "Y2_strat"),
                 split_var = "X1_obs",
                 type = "local quantile",
                 split_at = c(0.5),
                 trunc = "X1") |>
    split_strata(strata = "new_strata",
                 split_var = "X2_obs",
                 type = "local quantile",
                 split_at = c(0.5),
                 trunc = "X2")
  
  names(full_data)[names(full_data) == "new_strata"] <- "strata"
  
  data <- full_data
  
  ####
  ## Wave 1
  ####
  
  # Step 1: Initialize phase1_data
  phase1_data <- full_data
  
  # Step 2: Determine optimum allocation for wave 1 with Wright algorithm
  wave1_allocation <- a_optimum_allocation(phase1_data,
                                           strata = "strata",
                                           nsample = n/4,
                                           vars = c("inflB11_phase1", "inflB12_phase1",
                                                    "inflB21_phase1", "inflB22_phase1"),
                                           weights = c(0.25, 0.25,
                                                       0.25, 0.25),
                                           method = "Neyman")
  
  # First n/4 according to wave1_allocation
  phase2_wave1 <- sample_strata(data = phase1_data,
                                strata = "strata",
                                id = "id",
                                design_data = wave1_allocation,
                                n_allocated = "stratum_size")
  
  names(phase2_wave1)[names(phase2_wave1) == "sample_indicator"] <- "sampled_wave1"
  
  # Sample wave 1
  phase2_wave1$X1 <- ifelse(phase2_wave1$sampled_wave1 == 1, full_data$X1 , NA)
  phase2_wave1$X2 <- ifelse(phase2_wave1$sampled_wave1 == 1, full_data$X2 , NA)
  phase2_wave1$Z <- ifelse(phase2_wave1$sampled_wave1 == 1, full_data$Z , NA)
  
  #####
  ## Wave 2
  #####
  
  # Update estimates for influence functions using generalized raking
  twophase_design_Y1 <- twophase(id = list(~id, ~id), 
                                 strata = list(NULL, ~strata), 
                                 subset = ~as.logical(sampled_wave1), 
                                 data = phase2_wave1,
                                 method = "simple")
  
  twophase_design_Y2 <- twophase(id = list(~id, ~id), 
                                 strata = list(NULL, ~strata), 
                                 subset = ~as.logical(sampled_wave1), 
                                 data = phase2_wave1,
                                 method = "simple")
  
  # Calibrate
  calibrated_twophase_Y1_wave1 <- calibrate(twophase_design_Y1,
                                            ~ inflB11_phase1 +
                                              inflB21_phase1 +
                                              inflB12_phase1 +
                                              inflB22_phase1 +
                                              strata,
                                            phase = 2,
                                            calfun = "raking")
  calibrated_twophase_Y2_wave1 <- calibrate(twophase_design_Y2,
                                            ~ inflB11_phase1 +
                                              inflB21_phase1 +
                                              inflB12_phase1 +
                                              inflB22_phase1 +
                                              strata,
                                            phase = 2,
                                            calfun = "raking")
  
  # Run models
  fitY1_wave1 <- svyglm(Y1 ~ X1 + X2 + Z + Y2, family = quasibinomial,
                        design = calibrated_twophase_Y1_wave1)
  fitY2_wave1 <- svyglm(Y2 ~ X1 + X2 + Z, family = quasibinomial, 
                        design = calibrated_twophase_Y2_wave1)
  
  # Get IFs
  infl_Y1_wave1 <- inf_fun_logit(fitY1_wave1)
  infl_Y2_wave1 <- inf_fun_logit(fitY2_wave1)
  infl_wave1 <- data.frame(as.numeric(rownames(infl_Y1_wave1)), 
                           infl_Y1_wave1[,"X1"],
                           infl_Y1_wave1[,"X2"],
                           infl_Y2_wave1[,"X1"],
                           infl_Y2_wave1[,"X2"])
  names(infl_wave1) <- c("id", "inflB11", "inflB21", "inflB12", "inflB22")
  phase2_wave1 <- dplyr::left_join(phase2_wave1, infl_wave1, by = "id")
  
  # Regress Latest IFs (computed above) on Phase 1 IFs
  resid_model_B11_wave1 <- lm(phase2_wave1$inflB11 ~ phase2_wave1$inflB11_phase1,
                              na.action = na.exclude)
  resid_model_B21_wave1 <- lm(phase2_wave1$inflB21 ~ phase2_wave1$inflB21_phase1,
                              na.action = na.exclude)
  resid_model_B12_wave1 <- lm(phase2_wave1$inflB12 ~ phase2_wave1$inflB12_phase1,
                              na.action = na.exclude)
  resid_model_B22_wave1 <- lm(phase2_wave1$inflB22 ~ phase2_wave1$inflB22_phase1,
                              na.action = na.exclude)
  
  # Get residuals
  phase2_wave1$residB11_wave1 <- resid(resid_model_B11_wave1)
  phase2_wave1$residB21_wave1 <- resid(resid_model_B21_wave1)
  phase2_wave1$residB12_wave1 <- resid(resid_model_B12_wave1)
  phase2_wave1$residB22_wave1 <- resid(resid_model_B22_wave1)
  
  # Calculate allocation for wave 2
  wave2_allocation <- a_optimal_allocate_wave(phase2_wave1,
                                              strata = "strata",
                                              vars = c("residB11_wave1", 
                                                       "residB21_wave1",
                                                       "residB12_wave1", 
                                                       "residB22_wave1"),
                                              weights = c(0.25, 0.25,
                                                          0.25, 0.25),
                                              already_sampled = "sampled_wave1",
                                              nsample = n/4,
                                              method = "simple",
                                              allocation_method = "Neyman",
                                              detailed = TRUE)
  
  
  # sample and merge data
  phase2_wave2 <- sample_strata(data = phase2_wave1,
                                strata = "strata",
                                id = "id",
                                already_sampled = "sampled_wave1",
                                design_data = wave2_allocation,
                                n_allocated = "n_to_sample")
  names(phase2_wave2)[names(phase2_wave2) == "sample_indicator"] <- "sampled_wave2"
  
  # Sample wave2
  phase2_wave2$X1 <- ifelse(phase2_wave2$sampled_wave1 == 1 |
                              phase2_wave2$sampled_wave2 == 1, full_data$X1 , NA)
  phase2_wave2$X2 <- ifelse(phase2_wave2$sampled_wave1 == 1|
                              phase2_wave2$sampled_wave2 == 1, full_data$X2 , NA)
  phase2_wave2$Z <- ifelse(phase2_wave2$sampled_wave1 == 1|
                             phase2_wave2$sampled_wave2 == 1, full_data$Z , NA)
  
  phase2_wave2 <- subset(phase2_wave2, select = -c(inflB11, inflB21,
                                                   inflB12, inflB22))
  
  #####
  ## Wave 3
  #####
  
  # Indicator for already sampled
  phase2_wave2$already_sampled <- phase2_wave2$sampled_wave1 +
    phase2_wave2$sampled_wave2 
  
  # Update estimates for influence functions using generalized raking
  twophase_design_Y1 <- twophase(id = list(~id, ~id), 
                                 strata = list(NULL, ~strata), 
                                 subset = ~as.logical(already_sampled), 
                                 data = phase2_wave2,
                                 method = "simple")
  
  twophase_design_Y2 <- twophase(id = list(~id, ~id), 
                                 strata = list(NULL, ~strata), 
                                 subset = ~as.logical(already_sampled), 
                                 data = phase2_wave2,
                                 method = "simple")
  
  # Calibrate
  calibrated_twophase_Y1_wave2 <- calibrate(twophase_design_Y1,
                                            ~ inflB11_phase1 +
                                              inflB21_phase1 +
                                              inflB12_phase1 +
                                              inflB22_phase1 +
                                              strata,
                                            phase = 2,
                                            calfun = "raking")
  calibrated_twophase_Y2_wave2 <- calibrate(twophase_design_Y2,
                                            ~ inflB11_phase1 +
                                              inflB21_phase1 +
                                              inflB12_phase1 +
                                              inflB22_phase1 +
                                              strata,
                                            phase = 2,
                                            calfun = "raking")
  
  # Run models
  fitY1_wave2 <- svyglm(Y1 ~ X1 + X2 + Z + Y2, family = quasibinomial,
                        design = calibrated_twophase_Y1_wave2)
  fitY2_wave2 <- svyglm(Y2 ~ X1 + X2 + Z, family = quasibinomial, 
                        design = calibrated_twophase_Y2_wave2)
  
  # Get IFs
  infl_Y1_wave2 <- inf_fun_logit(fitY1_wave2)
  infl_Y2_wave2 <- inf_fun_logit(fitY2_wave2)
  infl_wave2 <- data.frame(as.numeric(rownames(infl_Y1_wave2)), 
                           infl_Y1_wave2[,"X1"],
                           infl_Y1_wave2[,"X2"],
                           infl_Y2_wave2[,"X1"],
                           infl_Y2_wave2[,"X2"])
  names(infl_wave2) <- c("id", "inflB11", "inflB21", "inflB12", "inflB22")
  phase2_wave2 <- dplyr::left_join(phase2_wave2, infl_wave2, by = "id")
  
  # Regress Latest IFs (computed above) on Phase 1 IFs
  resid_model_B11_wave2 <- lm(phase2_wave2$inflB11 ~ phase2_wave2$inflB11_phase1,
                              na.action = na.exclude)
  resid_model_B21_wave2 <- lm(phase2_wave2$inflB21 ~ phase2_wave2$inflB21_phase1,
                              na.action = na.exclude)
  resid_model_B12_wave2 <- lm(phase2_wave2$inflB12 ~ phase2_wave2$inflB12_phase1,
                              na.action = na.exclude)
  resid_model_B22_wave2 <- lm(phase2_wave2$inflB22 ~ phase2_wave2$inflB22_phase1,
                              na.action = na.exclude)
  
  # Get residuals
  phase2_wave2$residB11_wave2 <- resid(resid_model_B11_wave2)
  phase2_wave2$residB21_wave2 <- resid(resid_model_B21_wave2)
  phase2_wave2$residB12_wave2 <- resid(resid_model_B12_wave2)
  phase2_wave2$residB22_wave2 <- resid(resid_model_B22_wave2)
  
  # Re-calculate allocations
  wave3_allocation <- a_optimal_allocate_wave(phase2_wave2,
                                              strata = "strata",
                                              vars = c("residB11_wave2", 
                                                       "residB21_wave2",
                                                       "residB12_wave2", 
                                                       "residB22_wave2"),
                                              weights = c(0.25, 0.25,
                                                          0.25, 0.25),
                                              already_sampled = "already_sampled",
                                              nsample = n/4,
                                              method = "simple",
                                              detailed = TRUE,
                                              allocation_method = "Neyman")
  
  
  # sample and merge data
  phase2_wave3 <- sample_strata(data = phase2_wave2,
                                strata = "strata",
                                id = "id",
                                already_sampled = "already_sampled",
                                design_data = wave3_allocation,
                                n_allocated = "n_to_sample")
  names(phase2_wave3)[names(phase2_wave3) == "sample_indicator"] <- "sampled_wave3"
  
  # Sample wave3
  phase2_wave3$X1 <- ifelse(phase2_wave3$sampled_wave1 == 1 |
                              phase2_wave3$sampled_wave2 == 1|
                              phase2_wave3$sampled_wave3 == 1, full_data$X1 , NA)
  phase2_wave3$X2 <- ifelse(phase2_wave3$sampled_wave1 == 1|
                              phase2_wave3$sampled_wave2 == 1|
                              phase2_wave3$sampled_wave3 == 1, full_data$X2 , NA)
  phase2_wave3$Z <- ifelse(phase2_wave3$sampled_wave1 == 1|
                             phase2_wave3$sampled_wave2 == 1|
                             phase2_wave3$sampled_wave3 == 1, full_data$Z , NA)
  
  phase2_wave3 <- subset(phase2_wave3, select = -c(inflB11, inflB21,
                                                   inflB12, inflB22))
  
  #####
  ## Wave 4
  #####
  
  # Indicator for already sampled
  phase2_wave3$already_sampled <- phase2_wave3$sampled_wave1 +
    phase2_wave3$sampled_wave2 + phase2_wave3$sampled_wave3
  
  # Update estimates for influence functions using generalized raking
  twophase_design_Y1 <- twophase(id = list(~id, ~id), 
                                 strata = list(NULL, ~strata), 
                                 subset = ~as.logical(already_sampled), 
                                 data = phase2_wave3,
                                 method = "simple")
  
  twophase_design_Y2 <- twophase(id = list(~id, ~id), 
                                 strata = list(NULL, ~strata), 
                                 subset = ~as.logical(already_sampled), 
                                 data = phase2_wave3,
                                 method = "simple")
  
  # Calibrate
  calibrated_twophase_Y1_wave3 <- calibrate(twophase_design_Y1,
                                            ~ inflB11_phase1 +
                                              inflB21_phase1 +
                                              inflB12_phase1 +
                                              inflB22_phase1 +
                                              strata,
                                            phase = 2,
                                            calfun = "raking")
  calibrated_twophase_Y2_wave3 <- calibrate(twophase_design_Y2,
                                            ~ inflB11_phase1 +
                                              inflB21_phase1 +
                                              inflB12_phase1 +
                                              inflB22_phase1 +
                                              strata,
                                            phase = 2,
                                            calfun = "raking")
  
  # Run models
  fitY1_wave3 <- svyglm(Y1 ~ X1 + X2 + Z + Y2, family = quasibinomial,
                        design = calibrated_twophase_Y1_wave3)
  fitY2_wave3 <- svyglm(Y2 ~ X1 + X2 + Z, family = quasibinomial, 
                        design = calibrated_twophase_Y2_wave3)
  
  # Get IFs
  infl_Y1_wave3 <- inf_fun_logit(fitY1_wave3)
  infl_Y2_wave3 <- inf_fun_logit(fitY2_wave3)
  infl_wave3 <- data.frame(as.numeric(rownames(infl_Y1_wave3)), 
                           infl_Y1_wave3[,"X1"],
                           infl_Y1_wave3[,"X2"],
                           infl_Y2_wave3[,"X1"],
                           infl_Y2_wave3[,"X2"])
  names(infl_wave3) <- c("id", "inflB11", "inflB21", "inflB12", "inflB22")
  phase2_wave3 <- dplyr::left_join(phase2_wave3, infl_wave3, by = "id")
  
  # Regress Latest IFs (computed above) on Phase 1 IFs
  resid_model_B11_wave3 <- lm(phase2_wave3$inflB11 ~ phase2_wave3$inflB11_phase1,
                              na.action = na.exclude)
  resid_model_B21_wave3 <- lm(phase2_wave3$inflB21 ~ phase2_wave3$inflB21_phase1,
                              na.action = na.exclude)
  resid_model_B12_wave3 <- lm(phase2_wave3$inflB12 ~ phase2_wave3$inflB12_phase1,
                              na.action = na.exclude)
  resid_model_B22_wave3 <- lm(phase2_wave3$inflB22 ~ phase2_wave3$inflB22_phase1,
                              na.action = na.exclude)
  
  # Get residuals
  phase2_wave3$residB11_wave3 <- resid(resid_model_B11_wave3)
  phase2_wave3$residB21_wave3 <- resid(resid_model_B21_wave3)
  phase2_wave3$residB12_wave3 <- resid(resid_model_B12_wave3)
  phase2_wave3$residB22_wave3 <- resid(resid_model_B22_wave3)
  
  # Re-calculate allocation
  wave4_allocation <- a_optimal_allocate_wave(phase2_wave3,
                                              strata = "strata",
                                              vars = c("residB11_wave3", 
                                                       "residB21_wave3",
                                                       "residB12_wave3", 
                                                       "residB22_wave3"),
                                              weights = c(0.25, 0.25,
                                                          0.25, 0.25),
                                              already_sampled = "already_sampled",
                                              nsample = n/4,
                                              method = "simple",
                                              allocation_method = "Neyman",
                                              detailed = TRUE)
  
  # Check for oversampling
  oversampled_A_optimal <- ifelse(all(wave4_allocation$nsample_optimal == 
                                        wave4_allocation$nsample_actual),
                                  0, 1)
  
  
  # sample and merge data
  phase2_wave4 <- sample_strata(data = phase2_wave3,
                                strata = "strata",
                                id = "id",
                                design_data = wave4_allocation,
                                already_sampled = "already_sampled",
                                n_allocated = "n_to_sample")
  names(phase2_wave4)[names(phase2_wave4) == "sample_indicator"] <- "sampled_wave4"
  
  # Sample wave4
  phase2_wave4$X1 <- ifelse(phase2_wave4$sampled_wave1 == 1 |
                              phase2_wave4$sampled_wave2 == 1|
                              phase2_wave4$sampled_wave3 == 1|
                              phase2_wave4$sampled_wave4 == 1, full_data$X1 , NA)
  phase2_wave4$X2 <- ifelse(phase2_wave4$sampled_wave1 == 1|
                              phase2_wave4$sampled_wave2 == 1|
                              phase2_wave4$sampled_wave3 == 1|
                              phase2_wave4$sampled_wave4 == 1, full_data$X2 , NA)
  phase2_wave4$Z <- ifelse(phase2_wave4$sampled_wave1 == 1|
                             phase2_wave4$sampled_wave2 == 1|
                             phase2_wave4$sampled_wave3 == 1|
                             phase2_wave4$sampled_wave4 == 1, full_data$Z , NA)
  
  phase2_wave4 <- subset(phase2_wave4, select = -c(inflB11,
                                                   inflB12))
  
  ####
  ## Sampling done: Now calculate the Beta estimates with raking 
  ## using the survey package
  ####
  
  phase2_wave4$already_sampled <- phase2_wave4$sampled_wave1 + 
    phase2_wave4$sampled_wave2 + phase2_wave4$sampled_wave3 + 
    phase2_wave4$sampled_wave4
  twophase_design <- twophase(id = list(~1, ~1), 
                              strata = list(NULL, ~strata), 
                              subset = ~as.logical(already_sampled), 
                              data = phase2_wave4)
  
  # Weights
  weightY1 <- svyglm(Y1 ~ X1 + X2 + Z + Y2, family = quasibinomial, 
                     design = twophase_design)
  weightY2 <- svyglm(Y2 ~ X1 + X2 + Z , family = quasibinomial, 
                     design = twophase_design)
  
  #########
  #### Raking 
  #########
  data <- phase2_wave4
  
  # Creating a survey object and calibrating the weights
  mydesign <- survey::twophase(
    id = list(~1, ~1), subset = ~as.logical(already_sampled),
    strata = list(NULL, ~strata), data = data,
  )
  infcalY1 <- survey::calibrate(mydesign, formula =  ~inflB11_phase1 +
                                  inflB21_phase1 + inflB12_phase1 +
                                  inflB22_phase1 + strata, phase = 2, calfun = "raking")
  infcalY2 <- survey::calibrate(mydesign, formula =  ~inflB11_phase1 +
                                  inflB21_phase1 + inflB12_phase1 +
                                  inflB22_phase1 + strata, phase = 2, calfun = "raking")
  
  # Fitting the outcome model: conditional treatment effect of interest.
  fitY1 <- survey::svyglm(Y1 ~ X1 + X2 + Z + Y2, design = infcalY1, family = quasibinomial)
  fitY2 <- survey::svyglm(Y2 ~ X1 + X2 + Z, design = infcalY2, family = quasibinomial)
  
  
  # Check if truth is contained in interval
  confintIPWB11 <- confint(weightY1)["X1",]
  confintIPWB21 <- confint(weightY1)["X2",]
  confintIPWB12 <- confint(weightY2)["X1",]
  confintIPWB22 <- confint(weightY2)["X2",]
  confintGRB11 <- confint(fitY1)["X1",]
  confintGRB21 <- confint(fitY1)["X2",]
  confintGRB12 <- confint(fitY2)["X1",]
  confintGRB22 <- confint(fitY2)["X2",]
  coverIPWB11 <- ifelse(B11 >= confintIPWB11[1] & B11 <= confintIPWB11[2],1,0)
  coverIPWB21 <- ifelse(B21 >= confintIPWB21[1] & B21 <= confintIPWB21[2],1,0)
  coverIPWB12 <- ifelse(B12 >= confintIPWB12[1] & B12 <= confintIPWB12[2],1,0)
  coverIPWB22 <- ifelse(B22 >= confintIPWB22[1] & B22 <= confintIPWB22[2],1,0)
  coverGRB11 <- ifelse(B11 >= confintGRB11[1] & B11 <= confintGRB11[2],1,0)
  coverGRB21 <- ifelse(B21 >= confintGRB21[1] & B21 <= confintGRB21[2],1,0)
  coverGRB12 <- ifelse(B12 >= confintGRB12[1] & B12 <= confintGRB12[2],1,0)
  coverGRB22 <- ifelse(B22 >= confintGRB22[1] & B22 <= confintGRB22[2],1,0)
  
  output <- c(fitY1$coefficients["X1"], fitY1$coefficients["X2"],
              fitY2$coefficients["X1"], fitY2$coefficients["X2"],
              coef(weightY1)["X1"], coef(weightY1)["X2"], 
              coef(weightY2)["X1"], coef(weightY2)["X2"],  
              SE(fitY1)["X1"], SE(fitY1)["X2"],
              SE(fitY2)["X1"], SE(fitY2)["X2"], 
              SE(weightY1)["X1"], SE(weightY1)["X2"], 
              SE(weightY2)["X1"], SE(weightY2)["X2"], 
              coverIPWB11, coverIPWB21,
              coverIPWB12, coverIPWB22, 
              coverGRB11, coverGRB21,
              coverGRB12,coverGRB22,
              oversampled_A_optimal, cor_Y1_Y2,
              prev_Y1, prev_Y2, B11, B21, B12, B22)
  return(output)
}

######
### Strategy 5 - A-optimal sampling with true influence functions
#####

Strategy5 <- function(n, simulations_df, scenario, N = 10000){
  
  # Set up betas for scenario
  B01 <- simulations_df["B01", scenario]
  B11 <- simulations_df["B11", scenario]
  B21 <- simulations_df["B21", scenario]
  B31 <- simulations_df["B31", scenario]
  B41 <- simulations_df["B41", scenario]
  B02 <- simulations_df["B02", scenario]
  B12 <- simulations_df["B12", scenario]
  B22 <- simulations_df["B22", scenario]
  B32 <- simulations_df["B32", scenario] 
  
  # Set up of Xs
  corX1X2 <- simulations_df["corX1X2", scenario]
  sensZ <- simulations_df["sensZ", scenario]
  specZ <- simulations_df["specZ", scenario]
  prevZ <- simulations_df["prevZ", scenario]
  
  # Measurement error
  sensY1 <- simulations_df["sensY1", scenario]
  specY1 <- simulations_df["specY1", scenario]
  sensY2 <- simulations_df["sensY2", scenario]
  specY2 <- simulations_df["specY2", scenario]
  error_varX1 <- simulations_df["error_varX1", scenario]
  error_varX2 <- simulations_df["error_varX2", scenario]
  
  
  #####
  ## Generate Phase 1 Data
  #####
  
  id <- 1:N
  
  # Generate correlated covariates
  sigma <- matrix(c(1, corX1X2, 0.1, corX1X2, 1,  0.25, 0.1, 0.25, 1), nrow = 3)
  covs <- mvrnorm(N, mu = c(0,0,0), Sigma = sigma)
  X1 <- covs[,1]
  X2 <- covs[,2]
  Z <- ifelse(covs[,3] >= quantile(covs[,3], 1 - prevZ), 1, 0)
  cor(cbind(X1, X2, Z)) 
  
  Y2_probs <- exp(B02 + B12*X1 + B22*X2 + B32*Z)/(1 + exp(B02 + B12*X1 + B22*X2 + B32*Z))
  Y2 <- rbinom(N, 1, Y2_probs) # Compute realized value for Y2
  Y1_probs <- exp(B01 + B11*X1 + B21*X2 + B31*Z + B41*Y2)/
    (1 + exp(B01 + B11*X1 + B21*X2 + B31*Z + B41*Y2))
  Y1 <- rbinom(N, 1, Y1_probs) # Compute realized value for Y1
  
  full_data <- data.frame(id, X1, X2, Z, Y1, Y2)
  
  ## Test true model
  true_model1 <- glm(Y1 ~ X1 + X2 + Z + Y2, family = "binomial", data = full_data)
  true_model2 <- glm(Y2 ~ X1 + X2 + Z, family = "binomial", data = full_data)
  
  # Find observed correlation
  cor_Y1_Y2 <- cor(Y1, Y2) 
  
  # Prevalence
  prev_Y1 <- table(Y1)["1"]/length(Y1)
  prev_Y2 <- table(Y2)["1"]/length(Y2)
  
  # Find true coefficients
  true_modelY1 <- glm(Y1 ~ X1 + X2 + Z + Y2, family = "binomial", data = full_data) 
  true_modelY2 <- glm(Y2 ~ X1 + X2 + Z, family = "binomial", data = full_data) 
  true_B_Y1_X1 <- coef(true_modelY1)["X1"]
  true_B_Y1_X2 <- coef(true_modelY1)["X2"]
  true_B_Y2_X1 <- coef(true_modelY2)["X1"]
  true_B_Y2_X2 <- coef(true_modelY2)["X2"]
  
  # compute true influence functions
  full_data$infl_Y1_X1_true <- inf_fun_logit(true_modelY1)[,"X1"]
  full_data$infl_Y1_X2_true <- inf_fun_logit(true_modelY1)[,"X2"]
  full_data$infl_Y2_X1_true <- inf_fun_logit(true_modelY2)[,"X1"]
  full_data$infl_Y2_X2_true <- inf_fun_logit(true_modelY2)[,"X2"]
  
  
  #####
  ## Generate error-prone Phase 1 data
  ######
  X1X2ZY1Y2_obs <- mvrnorm(N, c(0,0,0,0,0),
                           matrix(c(error_varX1, 0.05, 0.04, 0.15, 0.04,
                                    0.05, error_varX2, 0.02, 0.10, 0.4,
                                    0.04, 0.02, 1, 0.01, 0,
                                    0.15, 0.10, 0.01, 1 , 0.15,
                                    0.04, 0.4, 0, 0.15, 1), nrow = 5))
  full_data$X1_obs <- full_data$X1 + X1X2ZY1Y2_obs[,1]
  full_data$X2_obs <- full_data$X2 + X1X2ZY1Y2_obs[,2]
  
  threshold_positiveZ <- qnorm((sensZ + 1)/2)
  threshold_negativeZ <- qnorm((specZ + 1)/2)
  threshold_positiveY1 <- qnorm((sensY1 + 1)/2)
  threshold_negativeY1 <- qnorm((specY1 + 1)/2)
  threshold_positiveY2 <- qnorm((sensY2 + 1)/2)
  threshold_negativeY2 <- qnorm((specY2 + 1)/2)
  
  full_data$Z_obs <- with(full_data, ifelse((Z > 0 & abs(X1X2ZY1Y2_obs[, 3]) < threshold_positiveZ) | 
                                              (Z < 0 & abs(X1X2ZY1Y2_obs[, 3]) > threshold_negativeZ), 
                                            1, 0))
  
  full_data$Y1_obs <- with(full_data, ifelse((Y1 == 1 & abs(X1X2ZY1Y2_obs[, 4]) < threshold_positiveY1) | 
                                               (Y1 == 0 & abs(X1X2ZY1Y2_obs[, 4]) > threshold_negativeY1), 1, 0))
  
  full_data$Y2_obs <- with(full_data, ifelse((Y2 == 1 & abs(X1X2ZY1Y2_obs[, 5]) < threshold_positiveY2) | 
                                               (Y2 == 0 & abs(X1X2ZY1Y2_obs[, 5]) > threshold_negativeY2),1, 0))
  
  
  #### Compute naive influence functions
  fitY1_phase1 <-  glm(Y1_obs ~ X1_obs + X2_obs + Z_obs + Y2_obs, 
                       family = "binomial", data = full_data)
  fitY2_phase1 <-  glm(Y2_obs ~ X1_obs + X2_obs + Z_obs, 
                       family = "binomial", data = full_data)
  
  full_data$inflB11_phase1 <- inf_fun_logit(fitY1_phase1)[,"X1_obs"]
  full_data$inflB21_phase1 <- inf_fun_logit(fitY1_phase1)[,"X2_obs"]
  full_data$inflB12_phase1 <- inf_fun_logit(fitY2_phase1)[,"X1_obs"]
  full_data$inflB22_phase1 <- inf_fun_logit(fitY2_phase1)[,"X2_obs"]
  
  #### Compute naive influence functions
  fitY1_phase1 <-  glm(Y1_obs ~ X1_obs + X2_obs + Z_obs + Y2_obs, 
                       family = "binomial", data = full_data)
  fitY2_phase1 <-  glm(Y2_obs ~ X1_obs + X2_obs + Z_obs, 
                       family = "binomial", data = full_data)
  
  full_data$inflB11_phase1 <- inf_fun_logit(fitY1_phase1)[,"X1_obs"]
  full_data$inflB12_phase1 <- inf_fun_logit(fitY1_phase1)[,"X2_obs"]
  full_data$inflB21_phase1 <- inf_fun_logit(fitY2_phase1)[,"X1_obs"]
  full_data$inflB22_phase1 <- inf_fun_logit(fitY2_phase1)[,"X2_obs"]
  
  #### Divide into strata based on observed phase 1 data X
  full_data <- full_data |>
    dplyr::mutate(Y1_strat = Y1_obs, Y2_strat = Y2_obs) |>
    split_strata(strata = c("Y1_strat", "Y2_strat"),
                 split_var = "X1_obs",
                 type = "local quantile",
                 split_at = c(0.5),
                 trunc = "X1") |>
    split_strata(strata = "new_strata",
                 split_var = "X2_obs",
                 type = "local quantile",
                 split_at = c(0.5),
                 trunc = "X2")
  
  names(full_data)[names(full_data) == "new_strata"] <- "strata"
  
  data <- full_data
  
  ####
  ## Wave 1
  ####
  
  # Step 1: Initialize phase 2 data
  data <- full_data
  
  # Step 2:  Calculate true influence functions using Tong's function
  data$infl_B11_true <- inf_fun_logit(true_modelY1)[,"X1"]
  data$infl_B21_true <- inf_fun_logit(true_modelY1)[,"X2"]
  data$infl_B12_true <- inf_fun_logit(true_modelY2)[,"X1"]
  data$infl_B22_true <- inf_fun_logit(true_modelY2)[,"X2"]
  
  # Step 4: Determine optimum allocation for phase 2
  allocation <- a_optimum_allocation(data,
                                           strata = "strata",
                                           nsample = n,
                                           vars = c("infl_B11_true", "infl_B12_true",
                                                    "infl_B21_true", "infl_B22_true"),
                                           weights = c(0.25, 0.25,
                                                       0.25, 0.25),
                                           method = "Neyman")
  
  #Sample
  data <- sample_strata(data = data,
                                strata = "strata",
                                id = "id",
                                design_data = allocation,
                                n_allocated = "stratum_size")
  
  # Sample 
  data$X1 <- ifelse(data$sample_indicator == 1, full_data$X1 , NA)
  data$X2 <- ifelse(data$sample_indicator == 1, full_data$X2 , NA)
  data$Z <- ifelse(data$sample_indicator == 1, full_data$Z , NA)
  
  # Generalized Raking
  twophase_design <- twophase(id = list(~1, ~1), 
                              strata = list(NULL, ~strata), 
                              subset = ~as.logical(sample_indicator), data = data)
  
  # Weights
  weightY1 <- svyglm(Y1 ~ X1 + X2 + Z + Y2, family = quasibinomial, 
                     design = twophase_design)
  weightY2 <- svyglm(Y2 ~ X1 + X2 + Z, family = quasibinomial, 
                     design = twophase_design)
  
  #########
  #### Raking
  #########
  # Creating a survey object and calibrating the weights
  mydesign <- survey::twophase(
    id = list(~1, ~1), subset = ~as.logical(sample_indicator),
    strata = list(NULL, ~strata), data = data,
  )
  infcalY1 <- survey::calibrate(mydesign, formula =  ~inflB11_phase1 +
                                  inflB21_phase1 + inflB12_phase1 +
                                  inflB22_phase1 + strata, phase = 2, calfun = "raking")
  infcalY2 <- survey::calibrate(mydesign, formula = ~inflB11_phase1 +
                                  inflB21_phase1 + inflB12_phase1 +
                                  inflB22_phase1 + strata, phase = 2, calfun = "raking")
  
  # Fitting the outcome model: conditional treatment effect of interest.
  fitY1 <- survey::svyglm(Y1 ~ X1 + X2 + Z + Y2, design = infcalY1, family = quasibinomial)
  fitY2 <- survey::svyglm(Y2 ~ X1 + X2 + Z, design = infcalY2, family = quasibinomial)
  
  # Check if truth is contained in interval
  confintIPWB11 <- confint(weightY1)["X1",]
  confintIPWB21 <- confint(weightY1)["X2",]
  confintIPWB12 <- confint(weightY2)["X1",]
  confintIPWB22 <- confint(weightY2)["X2",]
  confintGRB11 <- confint(fitY1)["X1",]
  confintGRB21 <- confint(fitY1)["X2",]
  confintGRB12 <- confint(fitY2)["X1",]
  confintGRB22 <- confint(fitY2)["X2",]
  coverIPWB11 <- ifelse(B11 >= confintIPWB11[1] & B11 <= confintIPWB11[2],1,0)
  coverIPWB21 <- ifelse(B21 >= confintIPWB21[1] & B21 <= confintIPWB21[2],1,0)
  coverIPWB12 <- ifelse(B12 >= confintIPWB12[1] & B12 <= confintIPWB12[2],1,0)
  coverIPWB22 <- ifelse(B22 >= confintIPWB22[1] & B22 <= confintIPWB22[2],1,0)
  coverGRB11 <- ifelse(B11 >= confintGRB11[1] & B11 <= confintGRB11[2],1,0)
  coverGRB21 <- ifelse(B21 >= confintGRB21[1] & B21 <= confintGRB21[2],1,0)
  coverGRB12 <- ifelse(B12 >= confintGRB12[1] & B12 <= confintGRB12[2],1,0)
  coverGRB22 <- ifelse(B22 >= confintGRB22[1] & B22 <= confintGRB22[2],1,0)
  
  
  output <- c(fitY1$coefficients["X1"], fitY1$coefficients["X2"],
              fitY2$coefficients["X1"], fitY2$coefficients["X2"],
              coef(weightY1)["X1"], coef(weightY1)["X2"], 
              coef(weightY2)["X1"], coef(weightY2)["X2"],  
              SE(fitY1)["X1"], SE(fitY1)["X2"],
              SE(fitY2)["X1"], SE(fitY2)["X2"], 
              SE(weightY1)["X1"], SE(weightY1)["X2"], 
              SE(weightY2)["X1"], SE(weightY2)["X2"], 
              coverIPWB11, coverIPWB21,
              coverIPWB12, coverIPWB22, 
              coverGRB11, coverGRB21,
              coverGRB12,coverGRB22,
              cor_Y1_Y2,
              prev_Y1, prev_Y2, B11, B21, B12, B22)
  return(output)
}

