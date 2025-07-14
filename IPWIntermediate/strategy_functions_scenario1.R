#####
### Functions for Each Strategy
#####

####
### Now, errors are correlated for X, Z1, and standardization occurs only on
### observed covariates (so phase1, then between each wave)
###
### I also now rake on strata
####

####
## Strategy 1 - Case-control sampling. n/4 cases and n/4 controls for each outcome
####

# @param n sample size
# @param data Phase 1 dataset with cols for Y1, Y2, X, Z1, Z2 and
# Y1_obs, Y2_obs, X_obs, Z1_obs, Z2_obs

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
  
  # Measurement error
  sensY1 <- simulations_df["sensY1", scenario]
  specY1 <- simulations_df["specY1", scenario]
  sensY2 <- simulations_df["sensY2", scenario]
  specY2 <- simulations_df["specY2", scenario]
  error_varX <- simulations_df["error_varX", scenario]
  
  # These stay across sims for now
  sensZ1 <- simulations_df["sensZ1", scenario]
  specZ1 <- simulations_df["specZ1", scenario]
  error_varZ2 <- simulations_df["error_varZ2", scenario]
  
  
  #####
  ## Generate Phase 1 Data
  #####
  
  id <- 1:N
  
  # Generate correlated covariates
  sigma <- matrix(c(1, 0.15, 0.1, 0.15, 1,  0.25, 0.1, 0.25, 1), nrow = 3)
  covs <- mvrnorm(N, mu = c(0,0,0), Sigma = sigma)
  X <- covs[,1]
  Z1 <- ifelse(covs[,2] >= quantile(covs[,2], 0.8), 1, 0)
  Z2 <- covs[,3]
  cor(cbind(X, Z1, Z2)) 
  
  Y2_probs <- exp(B02 + B12*X + B22*Z1 + B32*Z2)/(1 + exp(B02 + B12*X + B22*Z1 + B32*Z2))
  Y2 <- rbinom(N, 1, Y2_probs) # Compute realized value for Y2
  Y1_probs <- exp(B01 + B11*X + B21*Z1 + B31*Z2 + B41*Y2)/
    (1 + exp(B01 + B11*X+ B21*Z1 + B31*Z2 + B41*Y2))
  Y1 <- rbinom(N, 1, Y1_probs) # Compute realized value for Y1
  
  full_data <- data.frame(id, X, Z1, Z2, Y1, Y2)
  
  ## Test true model
  true_model1 <- glm(Y1 ~ X + Z1 + Z2 + Y2, family = "binomial", data = full_data)
  true_model2 <- glm(Y2 ~ X + Z1 + Z2, family = "binomial", data = full_data)
  
  # Find observed correlation
  cor_Y1_Y2 <- cor(Y1, Y2) 
  
  # Prevalence
  prev_Y1 <- table(Y1)["1"]/length(Y1)
  prev_Y2 <- table(Y2)["1"]/length(Y2)
  
  # Find true coefficients
  true_modelY1 <- glm(Y1 ~ X + Z1 + Z2 + Y2, family = "binomial", data = full_data) 
  true_modelY2 <- glm(Y2 ~ X + Z1 + Z2, family = "binomial", data = full_data) 
  true_B_11 <- coef(true_modelY1)["X"]
  true_B_12 <- coef(true_modelY2)["X"]
  
  
  #####
  ## Generate error-prone Phase 1 data
  ######
  
  Z1_prob_obs1 <- ifelse(full_data$Z1 == 1, sensZ1, 1 - specZ1)
  full_data$Z1_obs <- rbinom(N, 1, Z1_prob_obs1)
  XZY1Y2_obs <- mvrnorm(N, c(0,0,0,0),
                        matrix(c(error_varX, 0.03, 0.02, 0.025,
                                 0.03, error_varZ2, 0.01, 0,
                                 0.02, 0.01, 1, 0.25,
                                 0.025, 0, 0.25, 1), nrow = 4))
  full_data$X_obs <- full_data$X + XZY1Y2_obs[,1]
  full_data$Z2_obs <- full_data$Z2 + XZY1Y2_obs[,2]
  
  threshold_positiveY1 <- qnorm((sensY1 + 1)/2)
  threshold_negativeY1 <- qnorm((specY1 + 1)/2)
  threshold_positiveY2 <- qnorm((sensY2 + 1)/2)
  threshold_negativeY2 <- qnorm((specY2 + 1)/2)
  full_data$Y1_obs <- with(full_data, ifelse((Y1 == 1 & abs(XZY1Y2_obs[, 3]) < threshold_positiveY1) | 
                                               (Y1 == 0 & abs(XZY1Y2_obs[, 3]) > threshold_negativeY1), 1, 0))
  
  full_data$Y2_obs <- with(full_data, ifelse((Y2 == 1 & abs(XZY1Y2_obs[, 4]) < threshold_positiveY2) | 
                                               (Y2 == 0 & abs(XZY1Y2_obs[, 4]) > threshold_negativeY2),1, 0))
  
  #### Compute naive influence functions
  fitY1_phase1 <-  glm(Y1_obs ~ X_obs + Z1_obs + Z2_obs + Y2_obs, 
                       family = "binomial", data = full_data)
  fitY2_phase1 <-  glm(Y2_obs ~ X_obs + Z1_obs + Z2_obs, 
                       family = "binomial", data = full_data)
  
  full_data$inflB11_phase1 <- inf_fun_logit(fitY1_phase1)[,"X_obs"]
  full_data$inflB12_phase1 <- inf_fun_logit(fitY2_phase1)[,"X_obs"]
  
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
  data$X <- ifelse(data$sample_indicator == 1, full_data$X , NA)
  data$Z1 <- ifelse(data$sample_indicator == 1, full_data$Z1 , NA)
  data$Z2 <- ifelse(data$sample_indicator == 1, full_data$Z2 , NA)
  
  # Generalized Raking
  twophase_design <- twophase(id = list(~1, ~1), 
                              strata = list(NULL, ~strata), 
                              subset = ~as.logical(sample_indicator), data = data)
  
  # Weights
  weightY1 <- svyglm(Y1 ~ X + Z1 + Z2 + Y2, family = quasibinomial, 
                     design = twophase_design)
  weightY2 <- svyglm(Y2 ~ X + Z1 + Z2, family = quasibinomial, 
                     design = twophase_design)
  
  #########
  #### Raking
  #########
  
  # Creating a survey object and calibrating the weights
  mydesign <- survey::twophase(
    id = list(~1, ~1), subset = ~as.logical(sample_indicator),
    strata = list(NULL, ~strata), data = data
  )
  infcalY1 <- survey::calibrate(mydesign, formula = ~ inflB11_phase1 +
                                  inflB12_phase1 + strata, phase = 2, calfun = "raking")
  infcalY2 <- survey::calibrate(mydesign, formula = ~ inflB11_phase1 +
                                  inflB12_phase1 + strata, phase = 2, calfun = "raking")
  
  # Fitting the outcome model: conditional treatment effect of interest.
  fitY1 <- survey::svyglm(Y1 ~ X + Z1 + Z2 + Y2, design = infcalY1, family = quasibinomial)
  fitY2 <- survey::svyglm(Y2 ~ X + Z1 + Z2, design = infcalY2, family = quasibinomial)
  
  # Check if truth is contained in interval
  confintIPWB11 <- confint(weightY1)["X",]
  confintIPWB12 <- confint(weightY2)["X",]
  confintGRB11 <- confint(fitY1)["X",]
  confintGRB12 <- confint(fitY2)["X",]
  coverIPWB11 <- ifelse(B11 >= confintIPWB11[1] & B11 <= confintIPWB11[2],1,0)
  coverIPWB12 <- ifelse(B12 >= confintIPWB12[1] & B12 <= confintIPWB12[2],1,0)
  coverGRB11 <- ifelse(B11 >= confintGRB11[1] & B11 <= confintGRB11[2],1,0)
  coverGRB12 <- ifelse(B12 >= confintGRB12[1] & B12 <= confintGRB12[2],1,0)
  
  
  output <- c(fitY1$coefficients["X"], fitY2$coefficients["X"],
              coef(weightY1)["X"], coef(weightY2)["X"], SE(fitY1)["X"],
              SE(fitY2)["X"], SE(weightY1)["X"], SE(weightY2)["X"], 
              coverIPWB11, coverIPWB12, coverGRB11, coverGRB12, cor_Y1_Y2,
              prev_Y1, prev_Y2, true_B_11, true_B_12)
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
  
  # Measurement error
  sensY1 <- simulations_df["sensY1", scenario]
  specY1 <- simulations_df["specY1", scenario]
  sensY2 <- simulations_df["sensY2", scenario]
  specY2 <- simulations_df["specY2", scenario]
  error_varX <- simulations_df["error_varX", scenario]
  
  # These stay across sims for now
  sensZ1 <- simulations_df["sensZ1", scenario]
  specZ1 <- simulations_df["specZ1", scenario]
  error_varZ2 <- simulations_df["error_varZ2", scenario]
  
  
  #####
  ## Generate Phase 1 Data
  #####
  
  id <- 1:N
  
  # Generate correlated covariates
  sigma <- matrix(c(1, 0.15, 0.1, 0.15, 1,  0.25, 0.1, 0.25, 1), nrow = 3)
  covs <- mvrnorm(N, mu = c(0,0,0), Sigma = sigma)
  X <- covs[,1]
  Z1 <- ifelse(covs[,2] >= quantile(covs[,2], 0.8), 1, 0)
  Z2 <- covs[,3]
  cor(cbind(X, Z1, Z2)) 
  
  Y2_probs <- exp(B02 + B12*X + B22*Z1 + B32*Z2)/(1 + exp(B02 + B12*X + B22*Z1 + B32*Z2))
  Y2 <- rbinom(N, 1, Y2_probs) # Compute realized value for Y2
  Y1_probs <- exp(B01 + B11*X + B21*Z1 + B31*Z2 + B41*Y2)/
    (1 + exp(B01 + B11*X+ B21*Z1 + B31*Z2 + B41*Y2))
  Y1 <- rbinom(N, 1, Y1_probs) # Compute realized value for Y1
  
  full_data <- data.frame(id, X, Z1, Z2, Y1, Y2) 
  
  # Find observed correlation
  cor_Y1_Y2 <- cor(Y1, Y2) 
  
  # Prevalence
  prev_Y1 <- table(Y1)["1"]/length(Y1)
  prev_Y2 <- table(Y2)["1"]/length(Y2)
  
  # Find true coefficients
  true_modelY1 <- glm(Y1 ~ X + Z1 + Z2 + Y2, family = "binomial", data = full_data) 
  true_modelY2 <- glm(Y2 ~ X + Z1 + Z2, family = "binomial", data = full_data) 
  true_B_11 <- coef(true_modelY1)["X"]
  true_B_12 <- coef(true_modelY2)["X"]
  
  
  #####
  ## Generate error-prone Phase 1 data
  ######
  
  Z1_prob_obs1 <- ifelse(full_data$Z1 == 1, sensZ1, 1 - specZ1)
  full_data$Z1_obs <- rbinom(N, 1, Z1_prob_obs1)
  XZY1Y2_obs <- mvrnorm(N, c(0,0,0,0),
                        matrix(c(error_varX, 0.03, 0.02, 0.025,
                                 0.03, error_varZ2, 0.01, 0,
                                 0.02, 0.01, 1, 0.25,
                                 0.025, 0, 0.25, 1), nrow = 4))
  full_data$X_obs <- full_data$X + XZY1Y2_obs[,1]
  full_data$Z2_obs <- full_data$Z2 + XZY1Y2_obs[,2]
  
  threshold_positiveY1 <- qnorm((sensY1 + 1)/2)
  threshold_negativeY1 <- qnorm((specY1 + 1)/2)
  threshold_positiveY2 <- qnorm((sensY2 + 1)/2)
  threshold_negativeY2 <- qnorm((specY2 + 1)/2)
  full_data$Y1_obs <- with(full_data, ifelse((Y1 == 1 & abs(XZY1Y2_obs[, 3]) < threshold_positiveY1) | 
                                               (Y1 == 0 & abs(XZY1Y2_obs[, 3]) > threshold_negativeY1), 1, 0))
  
  full_data$Y2_obs <- with(full_data, ifelse((Y2 == 1 & abs(XZY1Y2_obs[, 4]) < threshold_positiveY2) | 
                                               (Y2 == 0 & abs(XZY1Y2_obs[, 4]) > threshold_negativeY2),1, 0))
  
  ####
  #### Divide into strata based on observed phase 1 data
  #full_data <- full_data |>
  #  dplyr::mutate(Y1_strat = Y1_obs, Y2_strat = Y2_obs, X_strat = X_obs) |>
  #  split_strata(strata = c("Y1_strat", "Y2_strat" ),
  #               split_var = "X_strat",
  #               type = "global quantile",
  #               split_at = c(0.25, 0.75),
  #               trunc = "X")
  
  #names(full_data)[names(full_data) == "new_strata"] <- "strata"
  
  
  
  #### Stratify on naive influence functions
  fitY1_phase1 <-  glm(Y1_obs ~ X_obs + Z1_obs + Z2_obs + Y2_obs, 
                       family = "binomial", data = full_data)
  fitY2_phase1 <-  glm(Y2_obs ~ X_obs + Z1_obs + Z2_obs, 
                       family = "binomial", data = full_data)
  
  full_data$inflB11_phase1 <- inf_fun_logit(fitY1_phase1)[,"X_obs"]
  full_data$inflB12_phase1 <- inf_fun_logit(fitY2_phase1)[,"X_obs"]
  
  #### Divide into strata based on observed phase 1 data X
  full_data <- full_data |>
    dplyr::mutate(Y1_strat = Y1_obs, Y2_strat = Y2_obs, X_strat = X_obs,
                  sum_abs_inf = abs(inflB11_phase1) + abs(inflB12_phase1)) |>
    split_strata(strata = c("Y1_strat", "Y2_strat"),
                 split_var = "X_obs",
                 type = "local quantile",
                 split_at = c(0.25, 0.75),
                 trunc = "X")
  
  names(full_data)[names(full_data) == "new_strata"] <- "strata"
  
  data <- full_data
  
  #####
  #####
  ##### Phase 2
  
  ## Here, each wave chooses n/4 samples where n/8 are chosen to be optimal for
  ## B_11 and n/8 to be optimal for B_12. 
  
  ####
  ## Wave 1
  ####
  
  # Step 1: Initialize phase1_data
  phase1_data <- full_data
  
  # Step 2: Logistic regression on Phase 1 data
  fitY1_phase1 <-  glm(Y1_obs ~ X_obs + Z1_obs + Z2_obs + Y2_obs, 
                       family = "binomial", data = phase1_data)
  fitY2_phase1 <-  glm(Y2_obs ~ X_obs + Z1_obs + Z2_obs, 
                       family = "binomial", data = phase1_data)
  
  # Step 3: Calculate influence functions using Tong's function
  phase1_data$inflB11 <- inf_fun_logit(fitY1_phase1)[,"X_obs"]
  phase1_data$inflB12 <- inf_fun_logit(fitY2_phase1)[,"X_obs"]
  
  # Step 4: Determine optimum allocation for wave 1 with Wright algorithm
  allocation1 <- optimum_allocation(phase1_data,
                                    strata = "strata",
                                    y = "inflB11",
                                    nsample = n/8, method = "Neyman")
  
  allocation2 <- optimum_allocation(phase1_data,
                                    strata = "strata",
                                    y = "inflB12",
                                    nsample = n/8, method = "Neyman")
  
  # Step 5: Combine allocations for total allocation of n/4
  wave1_allocation <- dplyr::left_join(allocation1, allocation2, by = "strata") |>
    dplyr::mutate(stratum_size = stratum_size.x + stratum_size.y)
  
  # Also, remove inflB11 and inflB12 vars
  phase1_data <- subset(phase1_data, select = -c(inflB11,
                                                 inflB12))
  
  # First 1/4th according to wave1_allocation
  phase2_wave1 <- sample_strata(data = phase1_data,
                                strata = "strata",
                                id = "id",
                                design_data = wave1_allocation,
                                n_allocated = "stratum_size")
  names(phase2_wave1)[names(phase2_wave1) == "sample_indicator"] <- "sampled_wave1"
  
  # Sample
  phase2_wave1$X <- ifelse(phase2_wave1$sampled_wave1 == 1, full_data$X , NA)
  phase2_wave1$Z1 <- ifelse(phase2_wave1$sampled_wave1 == 1, full_data$Z1 , NA)
  phase2_wave1$Z2 <- ifelse(phase2_wave1$sampled_wave1 == 1, full_data$Z2 , NA)
  
  #####
  ## Wave 2
  #####
  
  # Calculate influence functions using IPW
  twophase_design_wave1 <- twophase(id = list(~1, ~1), 
                                    strata = list(NULL, ~strata), 
                                    subset = ~as.logical(sampled_wave1), 
                                    data = phase2_wave1)
  fitY1_wave1 <- svyglm(Y1 ~ X + Z1 + Z2 + Y2, family = quasibinomial, 
                        design = twophase_design_wave1)
  fitY2_wave1 <- svyglm(Y2 ~ X + Z1 + Z2, family = quasibinomial, 
                        design = twophase_design_wave1)
  infl_wave1_Y1 <- inf_fun_logit(fitY1_wave1)
  infl_wave1_Y2 <- inf_fun_logit(fitY2_wave1)
  infl_wave1 <- data.frame(as.numeric(rownames(infl_wave1_Y1)), 
                           infl_wave1_Y1[,"X"],
                           infl_wave1_Y2[,"X"])
  names(infl_wave1) <- c("id", "inflB11", "inflB12")
  phase2_wave1 <- dplyr::left_join(phase2_wave1, infl_wave1, by = "id")
  
  
  # Re-calculate allocations
  allocation1 <- allocate_wave(phase2_wave1,
                               strata = "strata",
                               y = "inflB11", method = "iterative",
                               already_sampled = "sampled_wave1",
                               nsample = n/8, 
                               allocation_method = "Neyman")
  
  allocation2 <- allocate_wave(phase2_wave1,
                               strata = "strata",
                               y = "inflB12",
                               already_sampled = "sampled_wave1",
                               nsample = n/8, method = "iterative",
                               allocation_method = "Neyman")
  
  # and combine
  wave2_allocation <- dplyr::left_join(allocation1, allocation2, 
                                       by = "strata") |>
    dplyr::mutate(stratum_size = n_to_sample.x + n_to_sample.y)
  
  
  # sample and merge data
  phase2_wave2 <- sample_strata(data = phase2_wave1,
                                strata = "strata",
                                id = "id",
                                already_sampled = "sampled_wave1",
                                design_data = wave2_allocation,
                                n_allocated = "stratum_size")
  names(phase2_wave2)[names(phase2_wave2) == "sample_indicator"] <- "sampled_wave2"
  
  # Sample
  phase2_wave2$X <- ifelse(phase2_wave2$sampled_wave1 == 1|
                             phase2_wave2$sampled_wave2 == 1, full_data$X , NA)
  phase2_wave2$Z1 <- ifelse(phase2_wave2$sampled_wave1 == 1|
                              phase2_wave2$sampled_wave2 == 1, full_data$Z1 , NA)
  phase2_wave2$Z2 <- ifelse(phase2_wave2$sampled_wave1 == 1|
                              phase2_wave2$sampled_wave2 == 1, full_data$Z2 , NA)
  
  
  # Also, remove inflB11 and inflB12 vars
  phase2_wave2 <- subset(phase2_wave2, select = -c(inflB11,
                                                   inflB12))
  
  #####
  ## Wave 3
  #####
  
  # Already sampled indicator
  phase2_wave2$already_sampled <- phase2_wave2$sampled_wave1 +
    phase2_wave2$sampled_wave2
  
  # Calculate influence functions using IPW
  twophase_design_wave2 <- twophase(id = list(~1, ~1), 
                                    strata = list(NULL, ~strata), 
                                    subset = ~as.logical(already_sampled), 
                                    data = phase2_wave2)
  fitY1_wave2 <- svyglm(Y1 ~ X + Z1 + Z2 + Y2, family = quasibinomial, 
                        design = twophase_design_wave2)
  fitY2_wave2 <- svyglm(Y2 ~ X + Z1 + Z2, family = quasibinomial, 
                        design = twophase_design_wave2)
  infl_wave2_Y1 <- inf_fun_logit(fitY1_wave2)
  infl_wave2_Y2 <- inf_fun_logit(fitY2_wave2)
  infl_wave2 <- data.frame(as.numeric(rownames(infl_wave2_Y1)), 
                           infl_wave2_Y1[,"X"],
                           infl_wave2_Y2[,"X"])
  names(infl_wave2) <- c("id", "inflB11", "inflB12")
  phase2_wave2 <- dplyr::left_join(phase2_wave2, infl_wave2, by = "id")
  
  
  # Re-calculate allocations
  allocation1 <- allocate_wave(phase2_wave2,
                               strata = "strata",
                               y = "inflB11",
                               method = "iterative",
                               nsample = n/8,
                               already_sampled = "already_sampled",
                               allocation_method = "Neyman")
  
  allocation2 <- allocate_wave(phase2_wave2,
                               strata = "strata",
                               y = "inflB12",
                               method = "iterative",
                               nsample = n/8, 
                               already_sampled = "already_sampled",
                               allocation_method = "Neyman")
  
  # and combine
  wave3_allocation <- dplyr::left_join(allocation1, allocation2, 
                                       by = "strata") |>
    dplyr::mutate(stratum_size = n_to_sample.x + n_to_sample.y)
  
  
  # sample and merge data
  phase2_wave3 <- sample_strata(data = phase2_wave2,
                                strata = "strata",
                                id = "id",
                                already_sampled = "already_sampled",
                                design_data = wave3_allocation,
                                n_allocated = "stratum_size")
  names(phase2_wave3)[names(phase2_wave3) == "sample_indicator"] <- "sampled_wave3"
  
  # Sample
  phase2_wave3$X <- ifelse(phase2_wave3$sampled_wave1 == 1|
                             phase2_wave3$sampled_wave2 == 1|
                             phase2_wave3$sampled_wave3 == 1, full_data$X , NA)
  phase2_wave3$Z1 <- ifelse(phase2_wave3$sampled_wave1 == 1|
                              phase2_wave3$sampled_wave2 == 1|
                              phase2_wave3$sampled_wave3 == 1, full_data$Z1 , NA)
  phase2_wave3$Z2 <- ifelse(phase2_wave3$sampled_wave1 == 1|
                              phase2_wave3$sampled_wave2 == 1|
                              phase2_wave3$sampled_wave3 == 1, full_data$Z2 , NA)
  
  
  # Also, remove inflB11 and inflB12 vars
  phase2_wave3 <- subset(phase2_wave3, select = -c(inflB11,
                                                   inflB12))
  
  #####
  ## Wave 4
  #####
  
  # Already sampled indicator
  phase2_wave3$already_sampled <- phase2_wave3$sampled_wave1 +
    phase2_wave3$sampled_wave2 + phase2_wave3$sampled_wave3
  
  # Calculate influence functions using IPW
  twophase_design_wave3 <- twophase(id = list(~1, ~1), 
                                    strata = list(NULL, ~strata), 
                                    subset = ~as.logical(already_sampled), 
                                    data = phase2_wave3)
  fitY1_wave3 <- svyglm(Y1 ~ X + Z1 + Z2 + Y2, family = quasibinomial, 
                        design = twophase_design_wave3)
  fitY2_wave3 <- svyglm(Y2 ~ X + Z1 + Z2, family = quasibinomial, 
                        design = twophase_design_wave3)
  infl_wave3_Y1 <- inf_fun_logit(fitY1_wave3)
  infl_wave3_Y2 <- inf_fun_logit(fitY2_wave3)
  infl_wave3 <- data.frame(as.numeric(rownames(infl_wave3_Y1)), 
                           infl_wave3_Y1[,"X"],
                           infl_wave3_Y2[,"X"])
  names(infl_wave3) <- c("id", "inflB11", "inflB12")
  phase2_wave3 <- dplyr::left_join(phase2_wave3, infl_wave3, by = "id")
  
  
  # Re-calculate allocations
  allocation1 <- allocate_wave(phase2_wave3,
                               strata = "strata",
                               y = "inflB11",
                               method = "iterative",
                               nsample = n/8, 
                               already_sampled = "already_sampled",
                               detailed = TRUE,
                               allocation_method = "Neyman")
  
  allocation2 <- allocate_wave(phase2_wave3,
                               strata = "strata",
                               y = "inflB12",
                               method = "iterative",
                               nsample = n/8, 
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
    dplyr::mutate(stratum_size = n_to_sample.x + n_to_sample.y)
  
  # sample and merge data
  phase2_wave4 <- sample_strata(data = phase2_wave3,
                                strata = "strata",
                                id = "id",
                                already_sampled = "already_sampled",
                                design_data = wave4_allocation,
                                n_allocated = "stratum_size")
  names(phase2_wave4)[names(phase2_wave4) == "sample_indicator"] <- "sampled_wave4"
  
  # Sample wave4
  phase2_wave4$X <- ifelse(phase2_wave4$sampled_wave1 == 1 |
                             phase2_wave4$sampled_wave2 == 1|
                             phase2_wave4$sampled_wave3 == 1|
                             phase2_wave4$sampled_wave4 == 1, full_data$X , NA)
  phase2_wave4$Z1 <- ifelse(phase2_wave4$sampled_wave1 == 1|
                              phase2_wave4$sampled_wave2 == 1|
                              phase2_wave4$sampled_wave3 == 1|
                              phase2_wave4$sampled_wave4 == 1, full_data$Z1 , NA)
  phase2_wave4$Z2 <- ifelse(phase2_wave4$sampled_wave1 == 1|
                              phase2_wave4$sampled_wave2 == 1|
                              phase2_wave4$sampled_wave3 == 1|
                              phase2_wave4$sampled_wave4 == 1, full_data$Z2 , NA)
  
  # Remove inflB11 and inflB12 vars
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
  weightY1 <- svyglm(Y1 ~ X + Z1 + Z2 + Y2, family = quasibinomial, 
                     design = twophase_design)
  weightY2 <- svyglm(Y2 ~ X + Z1 + Z2 , family = quasibinomial, 
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
  infcalY1 <- survey::calibrate(mydesign, formula = ~  inflB11_phase1 +
                                  inflB12_phase1 + strata, phase = 2, calfun = "raking")
  infcalY2 <- survey::calibrate(mydesign, formula = ~ inflB11_phase1 +
                                  inflB12_phase1 + strata, phase = 2, calfun = "raking")
  
  # Fitting the outcome model: conditional treatment effect of interest.
  fitY1 <- survey::svyglm(Y1 ~ X + Z1 + Z2 + Y2, design = infcalY1, family = quasibinomial)
  fitY2 <- survey::svyglm(Y2 ~ X + Z1 + Z2, design = infcalY2, family = quasibinomial)
  
  # Check if truth is contained in interval
  confintIPWB11 <- confint(weightY1)["X",]
  confintIPWB12 <- confint(weightY2)["X",]
  confintGRB11 <- confint(fitY1)["X",]
  confintGRB12 <- confint(fitY2)["X",]
  coverIPWB11 <- ifelse(B11 >= confintIPWB11[1] & B11 <= confintIPWB11[2],1,0)
  coverIPWB12 <- ifelse(B12 >= confintIPWB12[1] & B12 <= confintIPWB12[2],1,0)
  coverGRB11 <- ifelse(B11 >= confintGRB11[1] & B11 <= confintGRB11[2],1,0)
  coverGRB12 <- ifelse(B12 >= confintGRB12[1] & B12 <= confintGRB12[2],1,0)
  
  output <- c(fitY1$coefficients["X"], fitY2$coefficients["X"],
              coef(weightY1)["X"], coef(weightY2)["X"], SE(fitY1)["X"],
              SE(fitY2)["X"], SE(weightY1)["X"], SE(weightY2)["X"], 
              coverIPWB11, coverIPWB12, coverGRB11, coverGRB12,
              oversampled_Y1, oversampled_Y2, cor_Y1_Y2,
              prev_Y1, prev_Y2, true_B_11, true_B_12)
  return(output)
}

#####
## Strategy 3 - Independent sequential
#####

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
  
  # Measurement error
  sensY1 <- simulations_df["sensY1", scenario]
  specY1 <- simulations_df["specY1", scenario]
  sensY2 <- simulations_df["sensY2", scenario]
  specY2 <- simulations_df["specY2", scenario]
  error_varX <- simulations_df["error_varX", scenario]
  
  # These stay across sims for now
  sensZ1 <- simulations_df["sensZ1", scenario]
  specZ1 <- simulations_df["specZ1", scenario]
  error_varZ2 <- simulations_df["error_varZ2", scenario]
  
  
  #####
  ## Generate Phase 1 Data
  #####
  
  id <- 1:N
  
  # Generate correlated covariates
  sigma <- matrix(c(1, 0.15, 0.1, 0.15, 1,  0.25, 0.1, 0.25, 1), nrow = 3)
  covs <- mvrnorm(N, mu = c(0,0,0), Sigma = sigma)
  X <- covs[,1]
  Z1 <- ifelse(covs[,2] >= quantile(covs[,2], 0.8), 1, 0)
  Z2 <- covs[,3]
  cor(cbind(X, Z1, Z2)) 
  
  Y2_probs <- exp(B02 + B12*X + B22*Z1 + B32*Z2)/(1 + exp(B02 + B12*X + B22*Z1 + B32*Z2))
  Y2 <- rbinom(N, 1, Y2_probs) # Compute realized value for Y2
  Y1_probs <- exp(B01 + B11*X + B21*Z1 + B31*Z2 + B41*Y2)/
    (1 + exp(B01 + B11*X+ B21*Z1 + B31*Z2 + B41*Y2))
  Y1 <- rbinom(N, 1, Y1_probs) # Compute realized value for Y1
  
  full_data <- data.frame(id, X, Z1, Z2, Y1, Y2) 
  
  # Find observed correlation
  cor_Y1_Y2 <- cor(Y1, Y2) 
  
  # Prevalence
  prev_Y1 <- table(Y1)["1"]/length(Y1)
  prev_Y2 <- table(Y2)["1"]/length(Y2)
  
  # Find true coefficients
  true_modelY1 <- glm(Y1 ~ X + Z1 + Z2 + Y2, family = "binomial", data = full_data) 
  true_modelY2 <- glm(Y2 ~ X + Z1 + Z2, family = "binomial", data = full_data) 
  true_B_11 <- coef(true_modelY1)["X"]
  true_B_12 <- coef(true_modelY2)["X"]
  
  
  #####
  ## Generate error-prone Phase 1 data
  ######
  
  Z1_prob_obs1 <- ifelse(full_data$Z1 == 1, sensZ1, 1 - specZ1)
  full_data$Z1_obs <- rbinom(N, 1, Z1_prob_obs1)
  XZY1Y2_obs <- mvrnorm(N, c(0,0,0,0),
                        matrix(c(error_varX, 0.03, 0.02, 0.025,
                                 0.03, error_varZ2, 0.01, 0,
                                 0.02, 0.01, 1, 0.25,
                                 0.025, 0, 0.25, 1), nrow = 4))
  full_data$X_obs <- full_data$X + XZY1Y2_obs[,1]
  full_data$Z2_obs <- full_data$Z2 + XZY1Y2_obs[,2]
  
  threshold_positiveY1 <- qnorm((sensY1 + 1)/2)
  threshold_negativeY1 <- qnorm((specY1 + 1)/2)
  threshold_positiveY2 <- qnorm((sensY2 + 1)/2)
  threshold_negativeY2 <- qnorm((specY2 + 1)/2)
  full_data$Y1_obs <- with(full_data, ifelse((Y1 == 1 & abs(XZY1Y2_obs[, 3]) < threshold_positiveY1) | 
                                               (Y1 == 0 & abs(XZY1Y2_obs[, 3]) > threshold_negativeY1), 1, 0))
  
  full_data$Y2_obs <- with(full_data, ifelse((Y2 == 1 & abs(XZY1Y2_obs[, 4]) < threshold_positiveY2) | 
                                               (Y2 == 0 & abs(XZY1Y2_obs[, 4]) > threshold_negativeY2),1, 0))
  
  ####
  #### Divide into strata based on observed phase 1 data
  #full_data <- full_data |>
  #  dplyr::mutate(Y1_strat = Y1_obs, Y2_strat = Y2_obs, X_strat = X_obs) |>
  #  split_strata(strata = c("Y1_strat", "Y2_strat" ),
  #               split_var = "X_strat",
  #               type = "global quantile",
  #               split_at = c(0.25, 0.75),
  #               trunc = "X")
  
  #names(full_data)[names(full_data) == "new_strata"] <- "strata"
  
  
  
  #### Stratify on naive influence functions
  fitY1_phase1 <-  glm(Y1_obs ~ X_obs + Z1_obs + Z2_obs + Y2_obs, 
                       family = "binomial", data = full_data)
  fitY2_phase1 <-  glm(Y2_obs ~ X_obs + Z1_obs + Z2_obs, 
                       family = "binomial", data = full_data)
  
  full_data$inflB11_phase1 <- inf_fun_logit(fitY1_phase1)[,"X_obs"]
  full_data$inflB12_phase1 <- inf_fun_logit(fitY2_phase1)[,"X_obs"]
  
  #### Divide into strata based on observed phase 1 data X
  full_data <- full_data |>
    dplyr::mutate(Y1_strat = Y1_obs, Y2_strat = Y2_obs, X_strat = X_obs,
                  sum_abs_inf = abs(inflB11_phase1) + abs(inflB12_phase1)) |>
    split_strata(strata = c("Y1_strat", "Y2_strat"),
                 split_var = "X_obs",
                 type = "local quantile",
                 split_at = c(0.25, 0.75),
                 trunc = "X")
  
  names(full_data)[names(full_data) == "new_strata"] <- "strata"
  
  data <- full_data
  
  #####
  #####
  ##### Phase 2
  
  ## Here, each wave chooses n/4 samples where n/8 are chosen to be optimal for
  ## B_11 and n/8 to be optimal for B_12.
  
  ####
  ## Wave 1
  ####
  
  # Step 1: Initialize Phase 1 data
  phase1_data <- full_data
  
  # Step 2: Logistic regression on Phase 1 data (Y1 only needed)
  fitY1_phase1 <-  glm(Y1_obs ~ X_obs + Z1_obs + Z2_obs + Y2_obs, 
                       family = "binomial", data = phase1_data)
  
  # Step 3: Calculate influence functions using Tong's function
  phase1_data$inflB11 <- inf_fun_logit(fitY1_phase1)[,"X_obs"]
  
  # Step 4: Determine optimum allocation for wave 1 with Wright algorithm
  wave1_allocation <- optimum_allocation(phase1_data,
                                         strata = "strata",
                                         y = "inflB11",
                                         nsample = n/4, method = "Neyman")
  
  # Step5: Now can remove inflB11 column, 
  phase1_data <- subset(phase1_data, select = -c(inflB11))
  
  
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
  
  # Calculate influence functions using IPW
  twophase_design_wave1 <- twophase(id = list(~1, ~1), 
                                    strata = list(NULL, ~strata), 
                                    subset = ~as.logical(sampled_wave1), 
                                    data = phase2_wave1)
  fitY1_wave1 <- svyglm(Y1 ~ X + Z1 + Z2 + Y2, family = quasibinomial, 
                        design = twophase_design_wave1)
  infl_wave1_Y1 <- inf_fun_logit(fitY1_wave1)
  infl_wave1 <- data.frame(as.numeric(rownames(infl_wave1_Y1)), 
                           infl_wave1_Y1[,"X"])
  names(infl_wave1) <- c("id", "inflB11")
  phase2_wave1 <- dplyr::left_join(phase2_wave1, infl_wave1, by = "id")
  
  
  # Calculate allocation for Y1
  wave2_allocation <- allocate_wave(phase2_wave1,
                                    strata = "strata",
                                    y = "inflB11",
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
  phase2_wave2$X <- ifelse(phase2_wave2$sampled_wave1 == 1|
                             phase2_wave2$sampled_wave2 == 1, full_data$X , NA)
  phase2_wave2$Z1 <- ifelse(phase2_wave2$sampled_wave1 == 1|
                              phase2_wave2$sampled_wave2 == 1, full_data$Z1 , NA)
  phase2_wave2$Z2 <- ifelse(phase2_wave2$sampled_wave1 == 1|
                              phase2_wave2$sampled_wave2 == 1, full_data$Z2 , NA)
  
  phase2_wave2 <- subset(phase2_wave2, select = -c(inflB11))
  
  #####
  ## Wave 3
  #####
  
  # Indicator for already sampled
  phase2_wave2$already_sampled <- phase2_wave2$sampled_wave1 +
    phase2_wave2$sampled_wave2
  
  # Calculate influence functions using IPW
  twophase_design_wave2 <- twophase(id = list(~1, ~1), 
                                    strata = list(NULL, ~strata), 
                                    subset = ~as.logical(already_sampled), 
                                    data = phase2_wave2)
  fitY2_wave2 <- svyglm(Y2 ~ X + Z1 + Z2, family = quasibinomial, 
                        design = twophase_design_wave2)
  infl_wave2_Y2 <- inf_fun_logit(fitY2_wave2)
  infl_wave2 <- data.frame(as.numeric(rownames(infl_wave2_Y2)), 
                           infl_wave2_Y2[,"X"])
  names(infl_wave2) <- c("id", "inflB12")
  phase2_wave2 <- dplyr::left_join(phase2_wave2, infl_wave2, by = "id")
  
  
  # Re-calculate allocation
  wave3_allocation <- allocate_wave(phase2_wave2,
                                    strata = "strata",
                                    y = "inflB12",
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
  phase2_wave3$X <- ifelse(phase2_wave3$sampled_wave1 == 1|
                             phase2_wave3$sampled_wave2 == 1|
                             phase2_wave3$sampled_wave3 == 1, full_data$X , NA)
  phase2_wave3$Z1 <- ifelse(phase2_wave3$sampled_wave1 == 1|
                              phase2_wave3$sampled_wave2 == 1|
                              phase2_wave3$sampled_wave3 == 1, full_data$Z1 , NA)
  phase2_wave3$Z2 <- ifelse(phase2_wave3$sampled_wave1 == 1|
                              phase2_wave3$sampled_wave2 == 1|
                              phase2_wave3$sampled_wave3 == 1, full_data$Z2 , NA)
  
  phase2_wave3 <- subset(phase2_wave3, select = -c(inflB12))
  
  
  #####
  ## Wave 4
  #####
  
  # Indicator for already sampled
  phase2_wave3$already_sampled <- phase2_wave3$sampled_wave1 +
    phase2_wave3$sampled_wave2 + phase2_wave3$sampled_wave3
  
  # Calculate influence functions using IPW
  twophase_design_wave3 <- twophase(id = list(~1, ~1), 
                                    strata = list(NULL, ~strata), 
                                    subset = ~as.logical(already_sampled), 
                                    data = phase2_wave3)
  fitY2_wave3 <- svyglm(Y2 ~ X + Z1 + Z2, family = quasibinomial, 
                        design = twophase_design_wave3)
  infl_wave3_Y2 <- inf_fun_logit(fitY2_wave3)
  infl_wave3 <- data.frame(as.numeric(rownames(infl_wave3_Y2)),
                           infl_wave3_Y2[,"X"])
  names(infl_wave3) <- c("id", "inflB12")
  phase2_wave3 <- dplyr::left_join(phase2_wave3, infl_wave3, by = "id")
  
  
  # Re-calculate allocation
  wave4_allocation <- allocate_wave(phase2_wave3,
                                    strata = "strata",
                                    y = "inflB12",
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
  
  # Sample
  phase2_wave4$X <- ifelse(phase2_wave4$sampled_wave1 == 1 |
                             phase2_wave4$sampled_wave2 == 1|
                             phase2_wave4$sampled_wave3 == 1|
                             phase2_wave4$sampled_wave4 == 1, full_data$X , NA)
  phase2_wave4$Z1 <- ifelse(phase2_wave4$sampled_wave1 == 1|
                              phase2_wave4$sampled_wave2 == 1|
                              phase2_wave4$sampled_wave3 == 1|
                              phase2_wave4$sampled_wave4 == 1, full_data$Z1 , NA)
  phase2_wave4$Z2 <- ifelse(phase2_wave4$sampled_wave1 == 1|
                              phase2_wave4$sampled_wave2 == 1|
                              phase2_wave4$sampled_wave3 == 1|
                              phase2_wave4$sampled_wave4 == 1, full_data$Z2 , NA)
  
  phase2_wave4 <- subset(phase2_wave4, select = -c(inflB12))
  
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
  weightY1 <- svyglm(Y1 ~ X + Z1 + Z2 + Y2, family = quasibinomial, 
                     design = twophase_design)
  weightY2 <- svyglm(Y2 ~ X + Z1 + Z2 , family = quasibinomial, 
                     design = twophase_design)
  
  
  #########
  #### Raking 
  #########
  data <- phase2_wave4
  
  mydesign <- survey::twophase(
    id = list(~1, ~1), subset = ~as.logical(already_sampled),
    strata = list(NULL, ~strata), data = data,
  )
  infcalY1 <- survey::calibrate(mydesign, formula = ~ inflB11_phase1 +
                                  inflB12_phase1 + strata, phase = 2, calfun = "raking")
  infcalY2 <- survey::calibrate(mydesign, formula = ~ inflB11_phase1 +
                                  inflB12_phase1 + strata, phase = 2, calfun = "raking")
  
  # Fitting the outcome model: conditional treatment effect of interest.
  fitY1 <- survey::svyglm(Y1 ~ X + Z1 + Z2 + Y2, design = infcalY1, family = quasibinomial)
  fitY2 <- survey::svyglm(Y2 ~ X + Z1 + Z2, design = infcalY2, family = quasibinomial)
  
  
  # Check if truth is contained in interval
  confintIPWB11 <- confint(weightY1)["X",]
  confintIPWB12 <- confint(weightY2)["X",]
  confintGRB11 <- confint(fitY1)["X",]
  confintGRB12 <- confint(fitY2)["X",]
  coverIPWB11 <- ifelse(B11 >= confintIPWB11[1] & B11 <= confintIPWB11[2],1,0)
  coverIPWB12 <- ifelse(B12 >= confintIPWB12[1] & B12 <= confintIPWB12[2],1,0)
  coverGRB11 <- ifelse(B11 >= confintGRB11[1] & B11 <= confintGRB11[2],1,0)
  coverGRB12 <- ifelse(B12 >= confintGRB12[1] & B12 <= confintGRB12[2],1,0)
  
  output <- c(fitY1$coefficients["X"], fitY2$coefficients["X"],
              coef(weightY1)["X"], coef(weightY2)["X"], SE(fitY1)["X"],
              SE(fitY2)["X"], SE(weightY1)["X"], SE(weightY2)["X"], 
              coverIPWB11, coverIPWB12, coverGRB11, coverGRB12, oversampled_Y2,
              cor_Y1_Y2,
              prev_Y1, prev_Y2, true_B_11, true_B_12)
  return(output)
}


#####
## Strategy 3.5 - Independent sequential with rare outcome first
#####

# This is looped part, so Wave 1 allocation is required as input
Strategy3.5 <- function(n, simulations_df, scenario, N = 10000){
  
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
  
  # Measurement error
  sensY1 <- simulations_df["sensY1", scenario]
  specY1 <- simulations_df["specY1", scenario]
  sensY2 <- simulations_df["sensY2", scenario]
  specY2 <- simulations_df["specY2", scenario]
  error_varX <- simulations_df["error_varX", scenario]
  
  # These stay across sims for now
  sensZ1 <- simulations_df["sensZ1", scenario]
  specZ1 <- simulations_df["specZ1", scenario]
  error_varZ2 <- simulations_df["error_varZ2", scenario]
  
  
  #####
  ## Generate Phase 1 Data
  #####
  
  id <- 1:N
  
  # Generate correlated covariates
  sigma <- matrix(c(1, 0.15, 0.1, 0.15, 1,  0.25, 0.1, 0.25, 1), nrow = 3)
  covs <- mvrnorm(N, mu = c(0,0,0), Sigma = sigma)
  X <- covs[,1]
  Z1 <- ifelse(covs[,2] >= quantile(covs[,2], 0.8), 1, 0)
  Z2 <- covs[,3]
  cor(cbind(X, Z1, Z2)) 
  
  Y2_probs <- exp(B02 + B12*X + B22*Z1 + B32*Z2)/(1 + exp(B02 + B12*X + B22*Z1 + B32*Z2))
  Y2 <- rbinom(N, 1, Y2_probs) # Compute realized value for Y2
  Y1_probs <- exp(B01 + B11*X + B21*Z1 + B31*Z2 + B41*Y2)/
    (1 + exp(B01 + B11*X+ B21*Z1 + B31*Z2 + B41*Y2))
  Y1 <- rbinom(N, 1, Y1_probs) # Compute realized value for Y1
  
  full_data <- data.frame(id, X, Z1, Z2, Y1, Y2) 
  
  # Find observed correlation
  cor_Y1_Y2 <- cor(Y1, Y2) 
  
  # Prevalence
  prev_Y1 <- table(Y1)["1"]/length(Y1)
  prev_Y2 <- table(Y2)["1"]/length(Y2)
  
  # Find true coefficients
  true_modelY1 <- glm(Y1 ~ X + Z1 + Z2 + Y2, family = "binomial", data = full_data) 
  true_modelY2 <- glm(Y2 ~ X + Z1 + Z2, family = "binomial", data = full_data) 
  true_B_11 <- coef(true_modelY1)["X"]
  true_B_12 <- coef(true_modelY2)["X"]
  
  
  #####
  ## Generate error-prone Phase 1 data
  ######
  
  Z1_prob_obs1 <- ifelse(full_data$Z1 == 1, sensZ1, 1 - specZ1)
  full_data$Z1_obs <- rbinom(N, 1, Z1_prob_obs1)
  XZY1Y2_obs <- mvrnorm(N, c(0,0,0,0),
                        matrix(c(error_varX, 0.03, 0.02, 0.025,
                                 0.03, error_varZ2, 0.01, 0,
                                 0.02, 0.01, 1, 0.25,
                                 0.025, 0, 0.25, 1), nrow = 4))
  full_data$X_obs <- full_data$X + XZY1Y2_obs[,1]
  full_data$Z2_obs <- full_data$Z2 + XZY1Y2_obs[,2]
  
  threshold_positiveY1 <- qnorm((sensY1 + 1)/2)
  threshold_negativeY1 <- qnorm((specY1 + 1)/2)
  threshold_positiveY2 <- qnorm((sensY2 + 1)/2)
  threshold_negativeY2 <- qnorm((specY2 + 1)/2)
  full_data$Y1_obs <- with(full_data, ifelse((Y1 == 1 & abs(XZY1Y2_obs[, 3]) < threshold_positiveY1) | 
                                               (Y1 == 0 & abs(XZY1Y2_obs[, 3]) > threshold_negativeY1), 1, 0))
  
  full_data$Y2_obs <- with(full_data, ifelse((Y2 == 1 & abs(XZY1Y2_obs[, 4]) < threshold_positiveY2) | 
                                               (Y2 == 0 & abs(XZY1Y2_obs[, 4]) > threshold_negativeY2),1, 0))
  
  
  ####
  #### Divide into strata based on observed phase 1 data
  #full_data <- full_data |>
  #  dplyr::mutate(Y1_strat = Y1_obs, Y2_strat = Y2_obs, X_strat = X_obs) |>
  #  split_strata(strata = c("Y1_strat", "Y2_strat" ),
  #               split_var = "X_strat",
  #               type = "global quantile",
  #               split_at = c(0.25, 0.75),
  #               trunc = "X")
  
  #names(full_data)[names(full_data) == "new_strata"] <- "strata"
  
  
  
  #### Stratify on naive influence functions
  fitY1_phase1 <-  glm(Y1_obs ~ X_obs + Z1_obs + Z2_obs + Y2_obs, 
                       family = "binomial", data = full_data)
  fitY2_phase1 <-  glm(Y2_obs ~ X_obs + Z1_obs + Z2_obs, 
                       family = "binomial", data = full_data)
  
  full_data$inflB11_phase1 <- inf_fun_logit(fitY1_phase1)[,"X_obs"]
  full_data$inflB12_phase1 <- inf_fun_logit(fitY2_phase1)[,"X_obs"]
  
  #### Divide into strata based on observed phase 1 data X
  full_data <- full_data |>
    dplyr::mutate(Y1_strat = Y1_obs, Y2_strat = Y2_obs, X_strat = X_obs,
                  sum_abs_inf = abs(inflB11_phase1) + abs(inflB12_phase1)) |>
    split_strata(strata = c("Y1_strat", "Y2_strat"),
                 split_var = "X_obs",
                 type = "local quantile",
                 split_at = c(0.25, 0.75),
                 trunc = "X")
  
  names(full_data)[names(full_data) == "new_strata"] <- "strata"
  
  data <- full_data
  
  ####
  ## Wave 1
  ####
  
  ## Here, 1st and 2nd waves (n/4 samples each) are optimized for B_12 and 
  ## 3rd and 4th last n/4 are optimized for B_12.
  
  # Step 1: Initialize Phase 1 data
  phase1_data <- full_data
  
  # Step 2: Logistic regression on Phase 1 data
  fitY2_phase1 <-  glm(Y2_obs ~ X_obs + Z1_obs + Z2_obs, 
                       family = "binomial", data = phase1_data)
  
  # Step 3: Calculate influence functions using Tong's function
  phase1_data$inflB12 <- inf_fun_logit(fitY1_phase1)[,"X_obs"]
  
  # Step 4: Determine optimum allocation for wave 1 with Wright algorithm
  wave1_allocation <- optimum_allocation(phase1_data,
                                         strata = "strata",
                                         y = "inflB12",
                                         nsample = n/4, method = "Neyman")
  
  # Step 6a: Now can remove inflB11 column, 
  phase1_data <- subset(phase1_data, select = -c(inflB12))
  
  
  # First n/4 according to wave1_allocation
  phase2_wave1 <- sample_strata(data = phase1_data,
                                strata = "strata",
                                id = "id",
                                design_data = wave1_allocation,
                                n_allocated = "stratum_size")
  names(phase2_wave1)[names(phase2_wave1) == "sample_indicator"] <- "sampled_wave1"
  
  #####
  ## Wave 2
  #####
  
  # Calculate influence functions using IPW
  twophase_design_wave1 <- twophase(id = list(~1, ~1), 
                                    strata = list(NULL, ~strata), 
                                    subset = ~as.logical(sampled_wave1), 
                                    data = phase2_wave1)
  fitY1_wave1 <- svyglm(Y2 ~ X + Z1 + Z2, family = quasibinomial, 
                        design = twophase_design_wave1)
  infl_wave1_Y2 <- inf_fun_logit(fitY1_wave1)
  infl_wave1 <- data.frame(as.numeric(rownames(infl_wave1_Y2)), 
                           infl_wave1_Y2[,"X"])
  names(infl_wave1) <- c("id", "inflB12")
  phase2_wave1 <- dplyr::left_join(phase2_wave1, infl_wave1, by = "id")
  
  
  # Calculate allocation for Y1
  wave2_allocation <- allocate_wave(phase2_wave1,
                                    strata = "strata",
                                    y = "inflB12",
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
  phase2_wave2$X <- ifelse(phase2_wave2$sampled_wave1 == 1|
                             phase2_wave2$sampled_wave2 == 1, full_data$X , NA)
  phase2_wave2$Z1 <- ifelse(phase2_wave2$sampled_wave1 == 1|
                              phase2_wave2$sampled_wave2 == 1, full_data$Z1 , NA)
  phase2_wave2$Z2 <- ifelse(phase2_wave2$sampled_wave1 == 1|
                              phase2_wave2$sampled_wave2 == 1, full_data$Z2 , NA)
  
  phase2_wave2 <- subset(phase2_wave2, select = -c(inflB12))
  
  #####
  ## Wave 3
  #####
  
  # Indicator for already sampled
  phase2_wave2$already_sampled <- phase2_wave2$sampled_wave1 +
    phase2_wave2$sampled_wave2
  
  # Calculate influence functions using IPW
  twophase_design_wave2 <- twophase(id = list(~1, ~1), 
                                    strata = list(NULL, ~strata), 
                                    subset = ~as.logical(already_sampled), 
                                    data = phase2_wave2)
  fitY2_wave2 <- svyglm(Y1 ~ X + Z1 + Z2 + Y2, family = quasibinomial, 
                        design = twophase_design_wave2)
  infl_wave2_Y1 <- inf_fun_logit(fitY2_wave2)
  infl_wave2 <- data.frame(as.numeric(rownames(infl_wave2_Y1)), 
                           infl_wave2_Y1[,"X"])
  names(infl_wave2) <- c("id", "inflB11")
  phase2_wave2 <- dplyr::left_join(phase2_wave2, infl_wave2, by = "id")
  
  
  # Re-calculate allocation
  wave3_allocation <- allocate_wave(phase2_wave2,
                                    strata = "strata",
                                    y = "inflB11",
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
  phase2_wave3$X <- ifelse(phase2_wave3$sampled_wave1 == 1|
                             phase2_wave3$sampled_wave2 == 1|
                             phase2_wave3$sampled_wave3 == 1, full_data$X , NA)
  phase2_wave3$Z1 <- ifelse(phase2_wave3$sampled_wave1 == 1|
                              phase2_wave3$sampled_wave2 == 1|
                              phase2_wave3$sampled_wave3 == 1, full_data$Z1 , NA)
  phase2_wave3$Z2 <- ifelse(phase2_wave3$sampled_wave1 == 1|
                              phase2_wave3$sampled_wave2 == 1|
                              phase2_wave3$sampled_wave3 == 1, full_data$Z2 , NA)
  
  phase2_wave3 <- subset(phase2_wave3, select = -c(inflB11))
  
  
  #####
  ## Wave 4
  #####
  
  # Indicator for already sampled
  phase2_wave3$already_sampled <- phase2_wave3$sampled_wave1 +
    phase2_wave3$sampled_wave2 + phase2_wave3$sampled_wave3
  
  # Calculate influence functions using IPW
  twophase_design_wave3 <- twophase(id = list(~1, ~1), 
                                    strata = list(NULL, ~strata), 
                                    subset = ~as.logical(already_sampled), 
                                    data = phase2_wave3)
  fitY1_wave3 <- svyglm(Y1 ~ X + Z1 + Z2 + Y2, family = quasibinomial, 
                        design = twophase_design_wave3)
  infl_wave3_Y1 <- inf_fun_logit(fitY1_wave3)
  infl_wave3 <- data.frame(as.numeric(rownames(infl_wave3_Y1)),
                           infl_wave3_Y1[,"X"])
  names(infl_wave3) <- c("id", "inflB11")
  phase2_wave3 <- dplyr::left_join(phase2_wave3, infl_wave3, by = "id")
  
  
  # Re-calculate allocation
  wave4_allocation <- allocate_wave(phase2_wave3,
                                    strata = "strata",
                                    y = "inflB11",
                                    nsample = n/4, method = "iterative",
                                    already_sampled = "already_sampled",
                                    detailed = TRUE,
                                    allocation_method = "Neyman")
  
  # Check for oversampling
  oversampled_Y1 <- ifelse(all(wave4_allocation$nsample_optimal == 
                                 wave4_allocation$nsample_actual),
                           0, 1)
  # no reason to check for Y2 here, as optimality in later waves with respect to Y1
  
  # sample and merge data
  phase2_wave4 <- sample_strata(data = phase2_wave3,
                                strata = "strata",
                                id = "id",
                                already_sampled = "already_sampled",
                                design_data = wave4_allocation,
                                n_allocated = "n_to_sample")
  names(phase2_wave4)[names(phase2_wave4) == "sample_indicator"] <- "sampled_wave4"
  
  # Sample wave4
  phase2_wave4$X <- ifelse(phase2_wave4$sampled_wave1 == 1 |
                             phase2_wave4$sampled_wave2 == 1|
                             phase2_wave4$sampled_wave3 == 1|
                             phase2_wave4$sampled_wave4 == 1, full_data$X , NA)
  phase2_wave4$Z1 <- ifelse(phase2_wave4$sampled_wave1 == 1|
                              phase2_wave4$sampled_wave2 == 1|
                              phase2_wave4$sampled_wave3 == 1|
                              phase2_wave4$sampled_wave4 == 1, full_data$Z1 , NA)
  phase2_wave4$Z2 <- ifelse(phase2_wave4$sampled_wave1 == 1|
                              phase2_wave4$sampled_wave2 == 1|
                              phase2_wave4$sampled_wave3 == 1|
                              phase2_wave4$sampled_wave4 == 1, full_data$Z2 , NA)
  
  phase2_wave4 <- subset(phase2_wave4, select = -c(inflB11))
  
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
  weightY1 <- svyglm(Y1 ~ X + Z1 + Z2 + Y2, family = quasibinomial, 
                     design = twophase_design)
  weightY2 <- svyglm(Y2 ~ X + Z1 + Z2 , family = quasibinomial, 
                     design = twophase_design)
  
  #########
  #### Raking 
  #########
  data <- phase2_wave4
  
  mydesign <- survey::twophase(
    id = list(~1, ~1), subset = ~as.logical(already_sampled),
    strata = list(NULL, ~strata), data = data,
  )
  infcalY1 <- survey::calibrate(mydesign, formula = ~ inflB11_phase1 +
                                  inflB12_phase1 + strata, phase = 2, calfun = "raking")
  infcalY2 <- survey::calibrate(mydesign, formula = ~ inflB11_phase1 +
                                  inflB12_phase1 + strata, phase = 2, calfun = "raking")
  
  # Fitting the outcome model: conditional treatment effect of interest.
  fitY1 <- survey::svyglm(Y1 ~ X + Z1 + Z2 + Y2, design = infcalY1, family = quasibinomial)
  fitY2 <- survey::svyglm(Y2 ~ X + Z1 + Z2, design = infcalY2, family = quasibinomial)
  
  # Check if truth is contained in interval
  confintIPWB11 <- confint(weightY1)["X",]
  confintIPWB12 <- confint(weightY2)["X",]
  confintGRB11 <- confint(fitY1)["X",]
  confintGRB12 <- confint(fitY2)["X",]
  coverIPWB11 <- ifelse(B11 >= confintIPWB11[1] & B11 <= confintIPWB11[2],1,0)
  coverIPWB12 <- ifelse(B12 >= confintIPWB12[1] & B12 <= confintIPWB12[2],1,0)
  coverGRB11 <- ifelse(B11 >= confintGRB11[1] & B11 <= confintGRB11[2],1,0)
  coverGRB12 <- ifelse(B12 >= confintGRB12[1] & B12 <= confintGRB12[2],1,0)
  
  output <- c(fitY1$coefficients["X"], fitY2$coefficients["X"],
              coef(weightY1)["X"], coef(weightY2)["X"], SE(fitY1)["X"],
              SE(fitY2)["X"], SE(weightY1)["X"], SE(weightY2)["X"], 
              coverIPWB11, coverIPWB12, coverGRB11, coverGRB12,
              oversampled_Y1, cor_Y1_Y2,
              prev_Y1, prev_Y2, true_B_11, true_B_12)
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
  
  # Measurement error
  sensY1 <- simulations_df["sensY1", scenario]
  specY1 <- simulations_df["specY1", scenario]
  sensY2 <- simulations_df["sensY2", scenario]
  specY2 <- simulations_df["specY2", scenario]
  error_varX <- simulations_df["error_varX", scenario]
  
  # These stay across sims for now
  sensZ1 <- simulations_df["sensZ1", scenario]
  specZ1 <- simulations_df["specZ1", scenario]
  error_varZ2 <- simulations_df["error_varZ2", scenario]
  
  
  #####
  ## Generate Phase 1 Data
  #####
  
  id <- 1:N
  
  # Generate correlated covariates
  sigma <- matrix(c(1, 0.15, 0.1, 0.15, 1,  0.25, 0.1, 0.25, 1), nrow = 3)
  covs <- mvrnorm(N, mu = c(0,0,0), Sigma = sigma)
  X <- covs[,1]
  Z1 <- ifelse(covs[,2] >= quantile(covs[,2], 0.8), 1, 0)
  Z2 <- covs[,3]
  cor(cbind(X, Z1, Z2)) 
  
  Y2_probs <- exp(B02 + B12*X + B22*Z1 + B32*Z2)/(1 + exp(B02 + B12*X + B22*Z1 + B32*Z2))
  Y2 <- rbinom(N, 1, Y2_probs) # Compute realized value for Y2
  Y1_probs <- exp(B01 + B11*X + B21*Z1 + B31*Z2 + B41*Y2)/
    (1 + exp(B01 + B11*X+ B21*Z1 + B31*Z2 + B41*Y2))
  Y1 <- rbinom(N, 1, Y1_probs) # Compute realized value for Y1
  
  full_data <- data.frame(id, X, Z1, Z2, Y1, Y2) 
  
  # Find observed correlation
  cor_Y1_Y2 <- cor(Y1, Y2) 
  
  # Prevalence
  prev_Y1 <- table(Y1)["1"]/length(Y1)
  prev_Y2 <- table(Y2)["1"]/length(Y2)
  
  # Find true coefficients
  true_modelY1 <- glm(Y1 ~ X + Z1 + Z2 + Y2, family = "binomial", data = full_data) 
  true_modelY2 <- glm(Y2 ~ X + Z1 + Z2, family = "binomial", data = full_data) 
  true_B_11 <- coef(true_modelY1)["X"]
  true_B_12 <- coef(true_modelY2)["X"]
  
  
  #####
  ## Generate error-prone Phase 1 data
  ######
  
  Z1_prob_obs1 <- ifelse(full_data$Z1 == 1, sensZ1, 1 - specZ1)
  full_data$Z1_obs <- rbinom(N, 1, Z1_prob_obs1)
  XZY1Y2_obs <- mvrnorm(N, c(0,0,0,0),
                        matrix(c(error_varX, 0.03, 0.02, 0.025,
                                 0.03, error_varZ2, 0.01, 0,
                                 0.02, 0.01, 1, 0.25,
                                 0.025, 0, 0.25, 1), nrow = 4))
  full_data$X_obs <- full_data$X + XZY1Y2_obs[,1]
  full_data$Z2_obs <- full_data$Z2 + XZY1Y2_obs[,2]
  
  threshold_positiveY1 <- qnorm((sensY1 + 1)/2)
  threshold_negativeY1 <- qnorm((specY1 + 1)/2)
  threshold_positiveY2 <- qnorm((sensY2 + 1)/2)
  threshold_negativeY2 <- qnorm((specY2 + 1)/2)
  full_data$Y1_obs <- with(full_data, ifelse((Y1 == 1 & abs(XZY1Y2_obs[, 3]) < threshold_positiveY1) | 
                                               (Y1 == 0 & abs(XZY1Y2_obs[, 3]) > threshold_negativeY1), 1, 0))
  
  full_data$Y2_obs <- with(full_data, ifelse((Y2 == 1 & abs(XZY1Y2_obs[, 4]) < threshold_positiveY2) | 
                                               (Y2 == 0 & abs(XZY1Y2_obs[, 4]) > threshold_negativeY2),1, 0))
  
  ####
  #### Divide into strata based on observed phase 1 data
  #full_data <- full_data |>
  #  dplyr::mutate(Y1_strat = Y1_obs, Y2_strat = Y2_obs, X_strat = X_obs) |>
  #  split_strata(strata = c("Y1_strat", "Y2_strat" ),
  #               split_var = "X_strat",
  #               type = "global quantile",
  #               split_at = c(0.25, 0.75),
  #               trunc = "X")
  
  #names(full_data)[names(full_data) == "new_strata"] <- "strata"
  
  
  
  #### Stratify on naive influence functions
  fitY1_phase1 <-  glm(Y1_obs ~ X_obs + Z1_obs + Z2_obs + Y2_obs, 
                       family = "binomial", data = full_data)
  fitY2_phase1 <-  glm(Y2_obs ~ X_obs + Z1_obs + Z2_obs, 
                       family = "binomial", data = full_data)
  
  full_data$inflB11_phase1 <- inf_fun_logit(fitY1_phase1)[,"X_obs"]
  full_data$inflB12_phase1 <- inf_fun_logit(fitY2_phase1)[,"X_obs"]
  
  #### Divide into strata based on observed phase 1 data X
  full_data <- full_data |>
    dplyr::mutate(Y1_strat = Y1_obs, Y2_strat = Y2_obs, X_strat = X_obs,
                  sum_abs_inf = abs(inflB11_phase1) + abs(inflB12_phase1)) |>
    split_strata(strata = c("Y1_strat", "Y2_strat"),
                 split_var = "X_obs",
                 type = "local quantile",
                 split_at = c(0.25, 0.75),
                 trunc = "X")
  
  names(full_data)[names(full_data) == "new_strata"] <- "strata"
  
  data <- full_data
  
  ####
  ## Wave 1
  ####
  
  # Step 1: Initialize phase1_data
  phase1_data <- full_data
  
  # Step 2: Logistic regression on Phase 1 data
  fitY1_phase1 <-  glm(Y1_obs ~ X_obs + Z1_obs + Z2_obs + Y2_obs, 
                       family = "binomial", data = phase1_data)
  fitY2_phase1 <-  glm(Y2_obs ~ X_obs + Z1_obs + Z2_obs, 
                       family = "binomial", data = phase1_data)
  
  # Step 3: Calculate influence functions using Tong's function
  phase1_data$inflB11 <- inf_fun_logit(fitY1_phase1)[,"X_obs"]
  phase1_data$inflB12 <- inf_fun_logit(fitY2_phase1)[,"X_obs"]
  
  # Step 4: Determine optimum allocation for wave 1 with Wright algorithm
  wave1_allocation <- a_optimum_allocation(phase1_data,
                                           strata = "strata",
                                           nsample = n/4,
                                           vars = c("inflB11", "inflB12"),
                                           weights = c(0.5, 0.5),
                                           method = "Neyman")
  
  # Also, remove inflB11 and inflB12 vars
  phase1_data <- subset(phase1_data, select = -c(inflB11,
                                                 inflB12))
  
  # First n/4 according to wave1_allocation
  phase2_wave1 <- sample_strata(data = phase1_data,
                                strata = "strata",
                                id = "id",
                                design_data = wave1_allocation,
                                n_allocated = "stratum_size")
  
  names(phase2_wave1)[names(phase2_wave1) == "sample_indicator"] <- "sampled_wave1"
  
  # Sample wave 1
  phase2_wave1$X <- ifelse(phase2_wave1$sampled_wave1 == 1, full_data$X , NA)
  phase2_wave1$Z1 <- ifelse(phase2_wave1$sampled_wave1 == 1, full_data$Z1 , NA)
  phase2_wave1$Z2 <- ifelse(phase2_wave1$sampled_wave1 == 1, full_data$Z2 , NA)
  
  #####
  ## Wave 2
  #####
  
  # Calculate influence functions using IPW
  twophase_design_wave1 <- twophase(id = list(~1, ~1), 
                                    strata = list(NULL, ~strata), 
                                    subset = ~as.logical(sampled_wave1), 
                                    data = phase2_wave1)
  fitY1_wave1 <- svyglm(Y1 ~ X + Z1 + Z2 + Y2, family = quasibinomial, 
                        design = twophase_design_wave1)
  fitY2_wave1 <- svyglm(Y2 ~ X + Z1 + Z2, family = quasibinomial, 
                        design = twophase_design_wave1)
  infl_wave1_Y1 <- inf_fun_logit(fitY1_wave1)
  infl_wave1_Y2 <- inf_fun_logit(fitY2_wave1)
  infl_wave1 <- data.frame(as.numeric(rownames(infl_wave1_Y1)), 
                           infl_wave1_Y1[,"X"],
                           infl_wave1_Y2[,"X"])
  names(infl_wave1) <- c("id", "inflB11", "inflB12")
  phase2_wave1 <- dplyr::left_join(phase2_wave1, infl_wave1, by = "id")
  
  
  # Calculate allocation for wave 2
  wave2_allocation <- a_optimal_allocate_wave(phase2_wave1,
                                              strata = "strata",
                                              vars = c("inflB11", "inflB12"),
                                              weights = c(0.5, 0.5),
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
  phase2_wave2$X <- ifelse(phase2_wave2$sampled_wave1 == 1 |
                             phase2_wave2$sampled_wave2 == 1, full_data$X , NA)
  phase2_wave2$Z1 <- ifelse(phase2_wave2$sampled_wave1 == 1|
                              phase2_wave2$sampled_wave2 == 1, full_data$Z1 , NA)
  phase2_wave2$Z2 <- ifelse(phase2_wave2$sampled_wave1 == 1|
                              phase2_wave2$sampled_wave2 == 1, full_data$Z2 , NA)
  
  phase2_wave2 <- subset(phase2_wave2, select = -c(inflB11,
                                                   inflB12))
  
  #####
  ## Wave 3
  #####
  
  # Calculate influence functions using IPW
  phase2_wave2$already_sampled <- phase2_wave2$sampled_wave1 +
    phase2_wave2$sampled_wave2
  twophase_design_wave2 <- twophase(id = list(~1, ~1), 
                                    strata = list(NULL, ~strata), 
                                    subset = ~as.logical(already_sampled), 
                                    data = phase2_wave2)
  fitY1_wave2 <- svyglm(Y1 ~ X + Z1 + Z2 + Y2, family = quasibinomial, 
                        design = twophase_design_wave2)
  fitY2_wave2 <- svyglm(Y2 ~ X + Z1 + Z2, family = quasibinomial, 
                        design = twophase_design_wave2)
  infl_wave2_Y1 <- inf_fun_logit(fitY1_wave2)
  infl_wave2_Y2 <- inf_fun_logit(fitY2_wave2)
  infl_wave2 <- data.frame(as.numeric(rownames(infl_wave2_Y1)), 
                           infl_wave2_Y1[,"X"],
                           infl_wave2_Y2[,"X"])
  names(infl_wave2) <- c("id", "inflB11", "inflB12")
  phase2_wave2 <- dplyr::left_join(phase2_wave2, infl_wave2, by = "id")
  
  
  # Re-calculate allocations
  wave3_allocation <- a_optimal_allocate_wave(phase2_wave2,
                                              strata = "strata",
                                              vars = c("inflB11", "inflB12"),
                                              weights = c(0.5, 0.5),
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
  phase2_wave3$X <- ifelse(phase2_wave3$sampled_wave1 == 1 |
                             phase2_wave3$sampled_wave2 == 1|
                             phase2_wave3$sampled_wave3 == 1, full_data$X , NA)
  phase2_wave3$Z1 <- ifelse(phase2_wave3$sampled_wave1 == 1|
                              phase2_wave3$sampled_wave2 == 1|
                              phase2_wave3$sampled_wave3 == 1, full_data$Z1 , NA)
  phase2_wave3$Z2 <- ifelse(phase2_wave3$sampled_wave1 == 1|
                              phase2_wave3$sampled_wave2 == 1|
                              phase2_wave3$sampled_wave3 == 1, full_data$Z2 , NA)
  
  phase2_wave3 <- subset(phase2_wave3, select = -c(inflB11, inflB12))
  
  #####
  ## Wave 4
  #####
  
  # Indicator for already sampled
  phase2_wave3$already_sampled <- phase2_wave3$sampled_wave1 +
    phase2_wave3$sampled_wave2 + phase2_wave3$sampled_wave3
  
  # Calculate influence functions using IPW
  twophase_design_wave3 <- twophase(id = list(~1, ~1), 
                                    strata = list(NULL, ~strata), 
                                    subset = ~as.logical(already_sampled), 
                                    data = phase2_wave3)
  fitY1_wave3 <- svyglm(Y1 ~ X + Z1 + Z2 + Y2, family = quasibinomial, 
                        design = twophase_design_wave3)
  fitY2_wave3 <- svyglm(Y2 ~ X + Z1 + Z2, family = quasibinomial, 
                        design = twophase_design_wave3)
  infl_wave3_Y1 <- inf_fun_logit(fitY1_wave3)
  infl_wave3_Y2 <- inf_fun_logit(fitY2_wave3)
  infl_wave3 <- data.frame(as.numeric(rownames(infl_wave3_Y1)), 
                           infl_wave3_Y1[,"X"],
                           infl_wave3_Y2[,"X"])
  names(infl_wave3) <- c("id", "inflB11", "inflB12")
  phase2_wave3 <- dplyr::left_join(phase2_wave3, infl_wave3, by = "id")
  
  
  # Re-calculate allocation
  wave4_allocation <- a_optimal_allocate_wave(phase2_wave3,
                                              strata = "strata",
                                              vars = c("inflB11", "inflB12"),
                                              weights = c(0.5, 0.5),
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
  phase2_wave4$X <- ifelse(phase2_wave4$sampled_wave1 == 1 |
                             phase2_wave4$sampled_wave2 == 1|
                             phase2_wave4$sampled_wave3 == 1|
                             phase2_wave4$sampled_wave4 == 1, full_data$X , NA)
  phase2_wave4$Z1 <- ifelse(phase2_wave4$sampled_wave1 == 1|
                              phase2_wave4$sampled_wave2 == 1|
                              phase2_wave4$sampled_wave3 == 1|
                              phase2_wave4$sampled_wave4 == 1, full_data$Z1 , NA)
  phase2_wave4$Z2 <- ifelse(phase2_wave4$sampled_wave1 == 1|
                              phase2_wave4$sampled_wave2 == 1|
                              phase2_wave4$sampled_wave3 == 1|
                              phase2_wave4$sampled_wave4 == 1, full_data$Z2 , NA)
  
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
  weightY1 <- svyglm(Y1 ~ X + Z1 + Z2 + Y2, family = quasibinomial, 
                     design = twophase_design)
  weightY2 <- svyglm(Y2 ~ X + Z1 + Z2 , family = quasibinomial, 
                     design = twophase_design)
  
  #########
  #### Raking 
  #########
  data <- phase2_wave4
  
  mydesign <- survey::twophase(
    id = list(~1, ~1), subset = ~as.logical(already_sampled),
    strata = list(NULL, ~strata), data = data,
  )
  infcalY1 <- survey::calibrate(mydesign, formula = ~ inflB11_phase1 +
                                  inflB12_phase1 + strata, phase = 2, calfun = "raking")
  infcalY2 <- survey::calibrate(mydesign, formula = ~ inflB11_phase1 +
                                  inflB12_phase1 + strata, phase = 2, calfun = "raking")
  
  # Fitting the outcome model: conditional treatment effect of interest.
  fitY1 <- survey::svyglm(Y1 ~ X + Z1 + Z2 + Y2, design = infcalY1, family = quasibinomial)
  fitY2 <- survey::svyglm(Y2 ~ X + Z1 + Z2, design = infcalY2, family = quasibinomial)
  
  
  # Check if truth is contained in interval
  confintIPWB11 <- confint(weightY1)["X",]
  confintIPWB12 <- confint(weightY2)["X",]
  confintGRB11 <- confint(fitY1)["X",]
  confintGRB12 <- confint(fitY2)["X",]
  coverIPWB11 <- ifelse(B11 >= confintIPWB11[1] & B11 <= confintIPWB11[2],1,0)
  coverIPWB12 <- ifelse(B12 >= confintIPWB12[1] & B12 <= confintIPWB12[2],1,0)
  coverGRB11 <- ifelse(B11 >= confintGRB11[1] & B11 <= confintGRB11[2],1,0)
  coverGRB12 <- ifelse(B12 >= confintGRB12[1] & B12 <= confintGRB12[2],1,0)
  
  output <- c(fitY1$coefficients["X"], fitY2$coefficients["X"],
              coef(weightY1)["X"], coef(weightY2)["X"], SE(fitY1)["X"],
              SE(fitY2)["X"], SE(weightY1)["X"], SE(weightY2)["X"], 
              coverIPWB11, coverIPWB12, coverGRB11, coverGRB12,
              oversampled_A_optimal, cor_Y1_Y2,
              prev_Y1, prev_Y2, true_B_11, true_B_12)
  return(output)
}

######
### Strategy 6 - True A-optimal design
#####

# Step 5: Create function to estimate beta-hat at each iteration
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
  
  # Measurement error
  sensY1 <- simulations_df["sensY1", scenario]
  specY1 <- simulations_df["specY1", scenario]
  sensY2 <- simulations_df["sensY2", scenario]
  specY2 <- simulations_df["specY2", scenario]
  error_varX <- simulations_df["error_varX", scenario]
  
  # These stay across sims for now
  sensZ1 <- simulations_df["sensZ1", scenario]
  specZ1 <- simulations_df["specZ1", scenario]
  error_varZ2 <- simulations_df["error_varZ2", scenario]
  
  
  #####
  ## Generate Phase 1 Data
  #####
  
  id <- 1:N
  
  # Generate correlated covariates
  sigma <- matrix(c(1, 0.15, 0.1, 0.15, 1,  0.25, 0.1, 0.25, 1), nrow = 3)
  covs <- mvrnorm(N, mu = c(0,0,0), Sigma = sigma)
  X <- covs[,1]
  Z1 <- ifelse(covs[,2] >= quantile(covs[,2], 0.8), 1, 0)
  Z2 <- covs[,3]
  cor(cbind(X, Z1, Z2)) 
  
  Y2_probs <- exp(B02 + B12*X + B22*Z1 + B32*Z2)/(1 + exp(B02 + B12*X + B22*Z1 + B32*Z2))
  Y2 <- rbinom(N, 1, Y2_probs) # Compute realized value for Y2
  Y1_probs <- exp(B01 + B11*X + B21*Z1 + B31*Z2 + B41*Y2)/
    (1 + exp(B01 + B11*X+ B21*Z1 + B31*Z2 + B41*Y2))
  Y1 <- rbinom(N, 1, Y1_probs) # Compute realized value for Y1
  
  full_data <- data.frame(id, X, Z1, Z2, Y1, Y2) 
  
  # Find observed correlation
  cor_Y1_Y2 <- cor(Y1, Y2) 
  
  # Prevalence
  prev_Y1 <- table(Y1)["1"]/length(Y1)
  prev_Y2 <- table(Y2)["1"]/length(Y2)
  
  # Find true coefficients
  true_modelY1 <- glm(Y1 ~ X + Z1 + Z2 + Y2, family = "binomial", data = full_data) 
  true_modelY2 <- glm(Y2 ~ X + Z1 + Z2, family = "binomial", data = full_data) 
  true_B_11 <- coef(true_modelY1)["X"]
  true_B_12 <- coef(true_modelY2)["X"]
  
  
  #####
  ## Generate error-prone Phase 1 data
  ######
  
  Z1_prob_obs1 <- ifelse(full_data$Z1 == 1, sensZ1, 1 - specZ1)
  full_data$Z1_obs <- rbinom(N, 1, Z1_prob_obs1)
  XZY1Y2_obs <- mvrnorm(N, c(0,0,0,0),
                        matrix(c(error_varX, 0.03, 0.02, 0.025,
                                 0.03, error_varZ2, 0.01, 0,
                                 0.02, 0.01, 1, 0.25,
                                 0.025, 0, 0.25, 1), nrow = 4))
  full_data$X_obs <- full_data$X + XZY1Y2_obs[,1]
  full_data$Z2_obs <- full_data$Z2 + XZY1Y2_obs[,2]
  
  threshold_positiveY1 <- qnorm((sensY1 + 1)/2)
  threshold_negativeY1 <- qnorm((specY1 + 1)/2)
  threshold_positiveY2 <- qnorm((sensY2 + 1)/2)
  threshold_negativeY2 <- qnorm((specY2 + 1)/2)
  full_data$Y1_obs <- with(full_data, ifelse((Y1 == 1 & abs(XZY1Y2_obs[, 3]) < threshold_positiveY1) | 
                                               (Y1 == 0 & abs(XZY1Y2_obs[, 3]) > threshold_negativeY1), 1, 0))
  
  full_data$Y2_obs <- with(full_data, ifelse((Y2 == 1 & abs(XZY1Y2_obs[, 4]) < threshold_positiveY2) | 
                                               (Y2 == 0 & abs(XZY1Y2_obs[, 4]) > threshold_negativeY2),1, 0))
  
  ####
  #### Divide into strata based on observed phase 1 data
  #full_data <- full_data |>
  #  dplyr::mutate(Y1_strat = Y1_obs, Y2_strat = Y2_obs, X_strat = X_obs) |>
  #  split_strata(strata = c("Y1_strat", "Y2_strat" ),
  #               split_var = "X_strat",
  #               type = "global quantile",
  #               split_at = c(0.25, 0.75),
  #               trunc = "X")
  
  #names(full_data)[names(full_data) == "new_strata"] <- "strata"
  
  
  
  #### Stratify on naive influence functions
  fitY1_phase1 <-  glm(Y1_obs ~ X_obs + Z1_obs + Z2_obs + Y2_obs, 
                       family = "binomial", data = full_data)
  fitY2_phase1 <-  glm(Y2_obs ~ X_obs + Z1_obs + Z2_obs, 
                       family = "binomial", data = full_data)
  
  full_data$inflB11_phase1 <- inf_fun_logit(fitY1_phase1)[,"X_obs"]
  full_data$inflB12_phase1 <- inf_fun_logit(fitY2_phase1)[,"X_obs"]
  
  #### Divide into strata based on observed phase 1 data X
  full_data <- full_data |>
    dplyr::mutate(Y1_strat = Y1_obs, Y2_strat = Y2_obs, X_strat = X_obs,
                  sum_abs_inf = abs(inflB11_phase1) + abs(inflB12_phase1)) |>
    split_strata(strata = c("Y1_strat", "Y2_strat"),
                 split_var = "X_obs",
                 type = "local quantile",
                 split_at = c(0.25, 0.75),
                 trunc = "X")
  
  names(full_data)[names(full_data) == "new_strata"] <- "strata"
  
  ####
  ## Sample
  ####
  
  #Step 1: Compute true influence functions
  full_data$inflB11_true <- inf_fun_logit(true_modelY1)[,"X"]
  full_data$inflB12_true <- inf_fun_logit(true_modelY2)[,"X"]
  data <- full_data
  
  # Step 3: Determine optimum allocation for phase 2 with Wright algorithm
  phase2_allocation <- a_optimum_allocation(data,
                                           strata = "strata",
                                           nsample = n,
                                           vars = c("inflB11_true",
                                                    "inflB12_true"),
                                           weights = c(0.5, 0.5),
                                           method = "Neyman")
  
  # Sample
  data <- sample_strata(data = data,
                                strata = "strata",
                                id = "id",
                                design_data = phase2_allocation,
                                n_allocated = "stratum_size")
  
  # Generalized Raking
  twophase_design <- twophase(id = list(~1, ~1), 
                              strata = list(NULL, ~strata), 
                              subset = ~as.logical(sample_indicator), data = data)
  
  # Weights
  weightY1 <- svyglm(Y1 ~ X + Z1 + Z2 + Y2, family = quasibinomial, 
                     design = twophase_design)
  weightY2 <- svyglm(Y2 ~ X + Z1 + Z2, family = quasibinomial, 
                     design = twophase_design)
  
  #########
  #### Raking 
  #########
  
  mydesign <- survey::twophase(
    id = list(~1, ~1), subset = ~as.logical(sample_indicator),
    strata = list(NULL, ~strata), data = data,
  )
  infcalY1 <- survey::calibrate(mydesign, formula = ~ inflB11_phase1 +
                                  inflB12_phase1 + strata, phase = 2, calfun = "raking")
  infcalY2 <- survey::calibrate(mydesign, formula = ~ inflB11_phase1 +
                                  inflB12_phase1 + strata, phase = 2, calfun = "raking")
  
  # Fitting the outcome model: conditional treatment effect of interest.
  fitY1 <- survey::svyglm(Y1 ~ X + Z1 + Z2 + Y2, design = infcalY1, family = quasibinomial)
  fitY2 <- survey::svyglm(Y2 ~ X + Z1 + Z2, design = infcalY2, family = quasibinomial)
  
  # Check if truth is contained in interval
  confintIPWB11 <- confint(weightY1)["X",]
  confintIPWB12 <- confint(weightY2)["X",]
  confintGRB11 <- confint(fitY1)["X",]
  confintGRB12 <- confint(fitY2)["X",]
  coverIPWB11 <- ifelse(B11 >= confintIPWB11[1] & B11 <= confintIPWB11[2],1,0)
  coverIPWB12 <- ifelse(B12 >= confintIPWB12[1] & B12 <= confintIPWB12[2],1,0)
  coverGRB11 <- ifelse(B11 >= confintGRB11[1] & B11 <= confintGRB11[2],1,0)
  coverGRB12 <- ifelse(B12 >= confintGRB12[1] & B12 <= confintGRB12[2],1,0)
  
  
  output <- c(fitY1$coefficients["X"], fitY2$coefficients["X"],
              coef(weightY1)["X"], coef(weightY2)["X"], SE(fitY1)["X"],
              SE(fitY2)["X"], SE(weightY1)["X"], SE(weightY2)["X"], 
              coverIPWB11, coverIPWB12, coverGRB11, coverGRB12, cor_Y1_Y2,
              prev_Y1, prev_Y2, true_B_11, true_B_12)
  return(output)
}
