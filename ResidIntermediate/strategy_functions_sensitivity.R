#####
### Sensitivity analysis of A-optimality with different priority weights
#####



######
###  A-optimal sampling
#####

# Step 5: Create function to estimate beta-hat at each iteration
Strategy_sensitivity <- function(n, simulations_df, scenario, my_weights, N = 10000){
  
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
  
  
  
  #### Compute naive influence functions
  fitY1_phase1 <-  glm(Y1_obs ~ X_obs + Z1_obs + Z2_obs + Y2_obs, 
                       family = "binomial", data = full_data)
  fitY2_phase1 <-  glm(Y2_obs ~ X_obs + Z1_obs + Z2_obs, 
                       family = "binomial", data = full_data)
  
  full_data$inflB11_phase1 <- inf_fun_logit(fitY1_phase1)[,"X_obs"]
  full_data$inflB12_phase1 <- inf_fun_logit(fitY2_phase1)[,"X_obs"]
  full_data$inflBZ1Y1_phase1 <- inf_fun_logit(fitY1_phase1)[,"Z1_obs"]
  full_data$inflBZ2Y1_phase1 <- inf_fun_logit(fitY1_phase1)[,"Z2_obs"]
  full_data$inflBZ1Y2_phase1 <- inf_fun_logit(fitY2_phase1)[,"Z1_obs"]
  full_data$inflBZ2Y2_phase1 <- inf_fun_logit(fitY2_phase1)[,"Z2_obs"]
  
  #### Divide into strata based on observed phase 1 data X
  full_data <- full_data |>
    dplyr::mutate(Y1_strat = Y1_obs, Y2_strat = Y2_obs, X_strat = X_obs) |>
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
  
  # Step 2: Determine optimum allocation for wave 1 with Wright algorithm
  wave1_allocation <- optimum_allocation(phase1_data,
                                           strata = "strata",
                                           nsample = n/4,
                                           y = c("inflB11_phase1", "inflB12_phase1"),
                                           weights = my_weights)
  
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
                                              inflB12_phase1 +
                                              inflBZ1Y1_phase1 +
                                              inflBZ2Y1_phase1 +
                                              inflBZ1Y2_phase1 +
                                              inflBZ2Y2_phase1 +
                                              strata,
                                            phase = 2,
                                            calfun = "raking")
  calibrated_twophase_Y2_wave1 <- calibrate(twophase_design_Y2,
                                            ~ inflB11_phase1 +
                                              inflB12_phase1 +
                                              inflBZ1Y1_phase1 +
                                              inflBZ2Y1_phase1 +
                                              inflBZ1Y2_phase1 +
                                              inflBZ2Y2_phase1 +
                                              strata,
                                            phase = 2,
                                            calfun = "raking")
  
  # Run models
  fitY1_wave1 <- svyglm(Y1 ~ X + Z1 + Z2 + Y2, family = quasibinomial,
                        design = calibrated_twophase_Y1_wave1)
  fitY2_wave1 <- svyglm(Y2 ~ X + Z1 + Z2, family = quasibinomial, 
                        design = calibrated_twophase_Y2_wave1)
  
  # Get IFs
  inflB11_wave1 <- inf_fun_logit(fitY1_wave1)
  inflB12_wave1 <- inf_fun_logit(fitY2_wave1)
  infl_wave1 <- data.frame(as.numeric(rownames(inflB11_wave1)), 
                           inflB11_wave1[,"X"],
                           inflB12_wave1[,"X"])
  names(infl_wave1) <- c("id", "inflB11", "inflB12")
  phase2_wave1 <- dplyr::left_join(phase2_wave1, infl_wave1, by = "id")
  
  # Calculate phase 2 sampling weights
  phase2_wave1 <- phase2_wave1 |>
    dplyr::group_by(strata) |>
    dplyr::mutate(
      phase2_weight = dplyr::n() / sum(as.logical(sampled_wave1))
    ) |>
    dplyr::ungroup() |> as.data.frame()
  
  # Regress Latest IFs (computed above) on Phase 1 IFs
  resid_model_B11_wave1 <- lm(inflB11 ~ inflB11_phase1 +
                                inflB12_phase1 + inflBZ1Y1_phase1 +
                                inflBZ2Y1_phase1 + inflBZ1Y2_phase1 +
                                inflBZ2Y2_phase1 + strata,
                              data = phase2_wave1,
                              weights = phase2_weight,
                              na.action = na.exclude)
  resid_model_B12_wave1 <- lm(inflB12 ~ inflB11_phase1 +
                                inflB12_phase1 + inflBZ1Y1_phase1 +
                                inflBZ2Y1_phase1 + inflBZ1Y2_phase1 +
                                inflBZ2Y2_phase1 + strata,
                              data = phase2_wave1,
                              weights = phase2_weight,
                              na.action = na.exclude)
  
  # Get residuals
  phase2_wave1$residB11_wave1 <- resid(resid_model_B11_wave1)
  phase2_wave1$residB12_wave1 <- resid(resid_model_B12_wave1)
  
  
  # Calculate allocation for wave 2
  wave2_allocation <- allocate_wave(phase2_wave1,
                                              strata = "strata",
                                              y = c("residB11_wave1",
                                                       "residB12_wave1"),
                                              weights = my_weights,
                                              already_sampled = "sampled_wave1",
                                              nsample = n/4,
                                              method = "simple",
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
                                              inflB12_phase1 +
                                              inflBZ1Y1_phase1 +
                                              inflBZ2Y1_phase1 +
                                              inflBZ1Y2_phase1 +
                                              inflBZ2Y2_phase1 +
                                              strata,
                                            phase = 2,
                                            calfun = "raking")
  calibrated_twophase_Y2_wave2 <- calibrate(twophase_design_Y2,
                                            ~ inflB11_phase1 +
                                              inflB12_phase1 +
                                              inflBZ1Y1_phase1 +
                                              inflBZ2Y1_phase1 +
                                              inflBZ1Y2_phase1 +
                                              inflBZ2Y2_phase1 +
                                              strata,
                                            phase = 2,
                                            calfun = "raking")
  
  # Run models
  fitY1_wave2 <- svyglm(Y1 ~ X + Z1 + Z2 + Y2, family = quasibinomial,
                        design = calibrated_twophase_Y1_wave2)
  fitY2_wave2 <- svyglm(Y2 ~ X + Z1 + Z2, family = quasibinomial, 
                        design = calibrated_twophase_Y2_wave2)
  
  # Get IFs
  inflB11_wave2 <- inf_fun_logit(fitY1_wave2)
  inflB12_wave2 <- inf_fun_logit(fitY2_wave2)
  infl_wave2 <- data.frame(as.numeric(rownames(inflB11_wave2)), 
                           inflB11_wave2[,"X"],
                           inflB12_wave2[,"X"])
  names(infl_wave2) <- c("id", "inflB11", "inflB12")
  phase2_wave2 <- dplyr::left_join(phase2_wave2, infl_wave2, by = "id")
  
  # Calculate phase 2 sampling weights
  phase2_wave2 <- phase2_wave2 |>
    dplyr::group_by(strata) |>
    dplyr::mutate(
      phase2_weight = dplyr::n() / sum(as.logical(already_sampled))
    ) |>
    dplyr::ungroup() |> as.data.frame()
  
  # Regress Latest IFs (computed above) on Phase 1 IFs
  resid_model_B11_wave2 <- lm(inflB11 ~ inflB11_phase1 +
                                inflB12_phase1 + inflBZ1Y1_phase1 +
                                inflBZ2Y1_phase1 + inflBZ1Y2_phase1 +
                                inflBZ2Y2_phase1 + strata,
                              data = phase2_wave2,
                              weights = phase2_weight,
                              na.action = na.exclude)
  resid_model_B12_wave2 <- lm(inflB12 ~ inflB11_phase1 +
                                inflB12_phase1 + inflBZ1Y1_phase1 +
                                inflBZ2Y1_phase1 + inflBZ1Y2_phase1 +
                                inflBZ2Y2_phase1 + strata,
                              data = phase2_wave2,
                              weights = phase2_weight,
                              na.action = na.exclude)
  
  # Get residuals
  phase2_wave2$residB11_wave2 <- resid(resid_model_B11_wave2)
  phase2_wave2$residB12_wave2 <- resid(resid_model_B12_wave2)
  
  
  # Re-calculate allocations
  wave3_allocation <- allocate_wave(phase2_wave2,
                                              strata = "strata",
                                              y = c("residB11_wave2", 
                                                       "residB12_wave2"),
                                              weights = my_weights,
                                              already_sampled = "already_sampled",
                                              nsample = n/4,
                                              method = "simple",
                                              detailed = TRUE)
  
  
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
                                              inflB12_phase1 +
                                              inflBZ1Y1_phase1 +
                                              inflBZ2Y1_phase1 +
                                              inflBZ1Y2_phase1 +
                                              inflBZ2Y2_phase1 +
                                              strata,
                                            phase = 2,
                                            calfun = "raking")
  calibrated_twophase_Y2_wave3 <- calibrate(twophase_design_Y2,
                                            ~ inflB11_phase1 +
                                              inflB12_phase1 +
                                              inflBZ1Y1_phase1 +
                                              inflBZ2Y1_phase1 +
                                              inflBZ1Y2_phase1 +
                                              inflBZ2Y2_phase1 +
                                              strata,
                                            phase = 2,
                                            calfun = "raking")
  
  # Run models
  fitY1_wave3 <- svyglm(Y1 ~ X + Z1 + Z2 + Y2, family = quasibinomial,
                        design = calibrated_twophase_Y1_wave3)
  fitY2_wave3 <- svyglm(Y2 ~ X + Z1 + Z2, family = quasibinomial, 
                        design = calibrated_twophase_Y2_wave3)
  
  # Get IFs
  inflB11_wave3 <- inf_fun_logit(fitY1_wave3)
  inflB12_wave3 <- inf_fun_logit(fitY2_wave3)
  infl_wave3 <- data.frame(as.numeric(rownames(inflB11_wave3)), 
                           inflB11_wave3[,"X"],
                           inflB12_wave3[,"X"])
  names(infl_wave3) <- c("id", "inflB11", "inflB12")
  phase2_wave3 <- dplyr::left_join(phase2_wave3, infl_wave3, by = "id")
  
  # Calculate phase 2 sampling weights
  phase2_wave3 <- phase2_wave3 |>
    dplyr::group_by(strata) |>
    dplyr::mutate(
      phase2_weight = dplyr::n() / sum(as.logical(already_sampled))
    ) |>
    dplyr::ungroup() |> as.data.frame()
  
  # Regress Latest IFs (computed above) on Phase 1 IFs
  resid_model_B11_wave3 <- lm(inflB11 ~ inflB11_phase1 +
                                inflB12_phase1 + inflBZ1Y1_phase1 +
                                inflBZ2Y1_phase1 + inflBZ1Y2_phase1 +
                                inflBZ2Y2_phase1 + strata,
                              data = phase2_wave3,
                              weights = phase2_weight,
                              na.action = na.exclude)
  resid_model_B12_wave3 <- lm(inflB12 ~ inflB11_phase1 +
                                inflB12_phase1 + inflBZ1Y1_phase1 +
                                inflBZ2Y1_phase1 + inflBZ1Y2_phase1 +
                                inflBZ2Y2_phase1 + strata,
                              data = phase2_wave3,
                              weights = phase2_weight,
                              na.action = na.exclude)
  
  # Get residuals
  phase2_wave3$residB11_wave3 <- resid(resid_model_B11_wave3)
  phase2_wave3$residB12_wave3 <- resid(resid_model_B12_wave3)
  
  
  # Re-calculate allocation
  wave4_allocation <- allocate_wave(phase2_wave3,
                                              strata = "strata",
                                              y = c("residB11_wave3", 
                                                       "residB12_wave3"),
                                              weights = my_weights,
                                              already_sampled = "already_sampled",
                                              nsample = n/4,
                                              method = "simple",
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
  
  # Creating a survey object and calibrating the weights
  mydesign <- survey::twophase(
    id = list(~1, ~1), subset = ~as.logical(already_sampled),
    strata = list(NULL, ~strata), data = data,
  )
  infcalY1 <- survey::calibrate(mydesign, formula = ~ strata + inflB11_phase1 + 
                                  inflB12_phase1 + inflBZ1Y1_phase1 +
                                  inflBZ2Y1_phase1 + inflBZ1Y2_phase1 +
                                  inflBZ2Y2_phase1, 
                                phase = 2, calfun = "raking")
  infcalY2 <- survey::calibrate(mydesign, formula = ~ strata + inflB11_phase1 + 
                                  inflB12_phase1 + inflBZ1Y1_phase1 +
                                  inflBZ2Y1_phase1 + inflBZ1Y2_phase1 +
                                  inflBZ2Y2_phase1,
                                phase = 2, calfun = "raking")
  
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
