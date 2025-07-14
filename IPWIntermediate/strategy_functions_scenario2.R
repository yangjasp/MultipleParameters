#####
### Functions for Each Strategy
#####

####
## Strategy 1 - Case-control sampling. n/2 cases and n/2 controls
####

# @param n sample size
# @param data Phase 1 dataset with cols for Y1, Y2, X, Z1, Z2 and
# Y1_obs, Y2_obs, X_obs, Z1_obs, Z2_obs

Strategy1 <- function(n, simulations_df,
                      scenario, N = 10000){
  
  # Set up betas for scenario
  B0 <- simulations_df["B0", scenario]
  B1 <- simulations_df["B1", scenario]
  B2 <- simulations_df["B2", scenario]
  B3 <- simulations_df["B3", scenario]
  
  # Set up of X's
  corX1X2 <- simulations_df["corX1X2", scenario]
  prevZ <- simulations_df["prevZ", scenario]
  
  # Measurement error
  sensY <- simulations_df["sensY", scenario]
  specY <- simulations_df["specY", scenario]
  error_varX1 <- simulations_df["error_varX1", scenario]
  error_varX2 <- simulations_df["error_varX2", scenario]
  sensZ <- simulations_df["sensZ", scenario]
  specZ <- simulations_df["specZ", scenario]

  
  
  #####
  ## Generate Phase 1 Data
  #####
  
  id <- 1:N
  
  # Generate correlated covariates
  sigma <- matrix(c(1, corX1X2 , 0.1, corX1X2, 1,  0.25, 0.1, 
                    0.25, 1), nrow = 3)
  covs <- mvrnorm(N, mu = c(0,0,0), Sigma = sigma)
  X1 <- covs[,1]
  X2 <- covs[,2]
  Z <- ifelse(covs[,3] >= quantile(covs[,3], 1- prevZ), 1, 0)
  cor(cbind(X1, X2, Z)) 
  
  Y_probs <- exp(B0 + B1*X1 + B2*X2 + B3*Z)/(1 +  exp(B0 + B1*X1 + B2*X2 + B3*Z))
  Y <- rbinom(N, 1, Y_probs) # Compute realized value for Y2
  
  full_data <- data.frame(id, X1, X2, Z, Y)
  
  ## Test true model
  true_model <- glm(Y ~ X1 + X2 + Z, family = "binomial", data = full_data)
  
  # Find observed correlation
  cor_X1_X2 <- cor(X1, X2) 
  
  # Prevalence
  prev_Y <- table(Y)["1"]/length(Y)
  
  # Find true coefficients
  true_model <- glm(Y ~ X1 + X2 + Z, family = "binomial", data = full_data) 
  true_BX1 <- coef(true_model)["X1"]
  true_BX2 <- coef(true_model)["X2"]
  
  
  #####
  ## Generate error-prone Phase 1 data
  ######
  
  Z_prob_obs1 <- ifelse(full_data$Z > 0, sensZ, 1 - specZ)
  full_data$Z_obs <- rbinom(N, 1, Z_prob_obs1)
  XZY_obs <- mvrnorm(N, c(0,0,0,0),
                     matrix(c(error_varX1, 0.03, 0.02, 0.02,
                              0.03, error_varX2, 0.05, 0.4,
                              0.02, 0.05, 1, 0,
                              0.02, 0.4, 0, 1), nrow = 4))
  full_data$X1_obs <- full_data$X1 + XZY_obs[,1]
  full_data$X2_obs <- full_data$X2 + XZY_obs[,2]
  
  threshold_positiveY <- qnorm((sensY + 1)/2)
  threshold_negativeY <- qnorm((specY + 1)/2)
  threshold_positiveZ <- qnorm((sensZ + 1)/2)
  threshold_negativeZ <- qnorm((specZ + 1)/2)
  full_data$Y_obs <- with(full_data, ifelse((Y == 1 & abs(XZY_obs[, 4]) < threshold_positiveY) | 
                                              (Y == 0 & abs(XZY_obs[, 4]) > threshold_negativeY), 1, 0))
  
  full_data$Z_obs <- with(full_data, ifelse((Z == 1 & abs(XZY_obs[, 3]) < threshold_positiveZ) | 
                                               (Z == 0 & abs(XZY_obs[, 3]) > threshold_negativeZ),
                                             max(Z), min(Z)))
  
  #### Compute naive influence functions
  fit_phase1 <-  glm(Y_obs ~ X1_obs + X2_obs + Z_obs, 
                     family = "binomial", data = full_data)
  
  full_data$inflBX1_phase1 <- inf_fun_logit(fit_phase1)[,"X1_obs"]
  full_data$inflBX2_phase1 <- inf_fun_logit(fit_phase1)[,"X2_obs"]
  
  #### Divide into strata based on observed phase 1 data X
  full_data <- full_data |> 
    #dplyr::mutate(Y_strat = Y_obs, X2_strat = X2_obs, strata = interaction(Y_strat, X2_strat))
    dplyr::mutate(strata = Y_obs)
  
  data <- full_data
  
  #####
  #####
  ##### Phase 2
  
  ### Collect samples via case-control sampling
  design_CC <- data.frame("strata" = c(1,0), "n_to_sample" = c(n/2,n/2))
  
  data <- data |>
    optimall::sample_strata(strata = "strata", id = "id",
                            design_data = design_CC)
  
  # Generalized Raking
  twophase_design <- twophase(id = list(~1, ~1), 
                              strata = list(NULL, ~strata), 
                              subset = ~as.logical(sample_indicator), data = data)
  
  # Weights
  weight <- svyglm(Y ~ X1 + X2 + Z, family = quasibinomial, 
                   design = twophase_design)
  
  #########
  #### Raking
  #########
  
  # Creating a survey object and calibrating the weights
  mydesign <- survey::twophase(
    id = list(~1, ~1), subset = ~as.logical(sample_indicator),
    strata = list(NULL, ~strata), data = data
  )
  infcal <- survey::calibrate(mydesign, formula = ~ inflBX1_phase1 + inflBX2_phase1
                              + strata, phase = 2, calfun = "raking")
  
  # Fitting the outcome model:
  fit <- survey::svyglm(Y ~ X1 + X2 + Z, design = infcal, family = quasibinomial)
  
  # Check if truth is contained in interval
  confintIPWBX1 <- confint(weight)["X1",]
  confintIPWBX2 <- confint(weight)["X2",]
  confintGRBX1 <- confint(fit)["X1",]
  confintGRBX2 <- confint(fit)["X2",]
  coverIPWBX1 <- ifelse(B1 >= confintIPWBX1[1] & 
                          B1 <= confintIPWBX1[2],1,0)
  coverIPWBX2 <- ifelse(B2 >= confintIPWBX2[1] & 
                          B2 <= confintIPWBX2[2],1,0)
  coverGRBX1 <- ifelse(B1 >= confintGRBX1[1] &
                         B1 <= confintGRBX1[2],1,0)
  coverGRBX2 <- ifelse(B2 >= confintGRBX2[1] & 
                         B2 <= confintGRBX2[2],1,0)
  
  
  output <- c(fit$coefficients["X1"], fit$coefficients["X2"],
              coef(weight)["X1"], coef(weight)["X2"], SE(fit)["X1"],
              SE(fit)["X2"], SE(weight)["X1"], SE(weight)["X2"], 
              coverIPWBX1, coverIPWBX2, coverGRBX1, coverGRBX2, cor_X1_X2,
              prev_Y, prevZ, B1, B2)
  return(output)
}


#####
## Strategy 2- Simulteaneous Multi-wave sampling
#####

# This is looped part, so
Strategy2 <- function(n, simulations_df, scenario, N = 10000){
  
  # Set up betas for scenario
  B0 <- simulations_df["B0", scenario]
  B1 <- simulations_df["B1", scenario]
  B2 <- simulations_df["B2", scenario]
  B3 <- simulations_df["B3", scenario]
  
  # Set up of X's
  corX1X2 <- simulations_df["corX1X2", scenario]
  prevZ <- simulations_df["prevZ", scenario]
  
  # Measurement error
  sensY <- simulations_df["sensY", scenario]
  specY <- simulations_df["specY", scenario]
  error_varX1 <- simulations_df["error_varX1", scenario]
  error_varX2 <- simulations_df["error_varX2", scenario]
  sensZ <- simulations_df["sensZ", scenario]
  specZ <- simulations_df["specZ", scenario]
  
  
  
  #####
  ## Generate Phase 1 Data
  #####
  
  id <- 1:N
  
  # Generate correlated covariates
  sigma <- matrix(c(1, corX1X2 , 0.1, corX1X2, 1,  0.25, 0.1, 
                    0.25, 1), nrow = 3)
  covs <- mvrnorm(N, mu = c(0,0,0), Sigma = sigma)
  X1 <- covs[,1]
  X2 <- covs[,2]
  Z <- ifelse(covs[,3] >= quantile(covs[,3], 1- prevZ), 1, 0)
  cor(cbind(X1, X2, Z)) 
  
  Y_probs <- exp(B0 + B1*X1 + B2*X2 + B3*Z)/(1 +  exp(B0 + B1*X1 + B2*X2 + B3*Z))
  Y <- rbinom(N, 1, Y_probs) # Compute realized value for Y2
  
  full_data <- data.frame(id, X1, X2, Z, Y)
  
  ## Test true model
  true_model <- glm(Y ~ X1 + X2 + Z, family = "binomial", data = full_data)
  
  # Find observed correlation
  cor_X1_X2 <- cor(X1, X2) 
  
  # Prevalence
  prev_Y <- table(Y)["1"]/length(Y)
  
  # Find true coefficients
  true_model <- glm(Y ~ X1 + X2 + Z, family = "binomial", data = full_data) 
  true_BX1 <- coef(true_model)["X1"]
  true_BX2 <- coef(true_model)["X2"]
  
  
  #####
  ## Generate error-prone Phase 1 data
  ######
  
  Z_prob_obs1 <- ifelse(full_data$Z > 0, sensZ, 1 - specZ)
  full_data$Z_obs <- rbinom(N, 1, Z_prob_obs1)
  XZY_obs <- mvrnorm(N, c(0,0,0,0),
                     matrix(c(error_varX1, 0.03, 0.02, 0.02,
                              0.03, error_varX2, 0.05, 0.4,
                              0.02, 0.05, 1, 0,
                              0.02, 0.4, 0, 1), nrow = 4))
  full_data$X1_obs <- full_data$X1 + XZY_obs[,1]
  full_data$X2_obs <- full_data$X2 + XZY_obs[,2]
  
  threshold_positiveY <- qnorm((sensY + 1)/2)
  threshold_negativeY <- qnorm((specY + 1)/2)
  threshold_positiveZ <- qnorm((sensZ + 1)/2)
  threshold_negativeZ <- qnorm((specZ + 1)/2)
  full_data$Y_obs <- with(full_data, ifelse((Y == 1 & abs(XZY_obs[, 4]) < threshold_positiveY) | 
                                              (Y == 0 & abs(XZY_obs[, 4]) > threshold_negativeY), 1, 0))
  
  full_data$Z_obs <- with(full_data, ifelse((Z == 1 & abs(XZY_obs[, 3]) < threshold_positiveZ) | 
                                              (Z == 0 & abs(XZY_obs[, 3]) > threshold_negativeZ),
                                            max(Z), min(Z)))
  
  ####
  #### Divide into strata based on observed phase 1 data
  
  #### Compute naive influence functions
  fit_phase1 <-  glm(Y_obs ~ X1_obs + X2_obs + Z_obs, 
                     family = "binomial", data = full_data)
  full_data$inflBX1_phase1 <- inf_fun_logit(fit_phase1)[,"X1_obs"]
  full_data$inflBX2_phase1 <- inf_fun_logit(fit_phase1)[,"X2_obs"]
  
  #### Divide into strata based on observed phase 1 data X
  full_data <- full_data |>
    dplyr::mutate(Y_strat = Y_obs) |>
    split_strata(strata = "Y_strat",
                 split_var = "X1_obs",
                 type = "local quantile",
                 split_at = 0.5,
                 trunc = "X1") |>
    split_strata(strata = "new_strata",
                 split_var = "X2_obs",
                 type = "local quantile",
                 split_at = 0.5,
                 trunc = "X2")
    
  
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
  fit_phase1 <-  glm(Y_obs ~ X1_obs + X2_obs + Z_obs, 
                     family = "binomial", data = phase1_data)
  
  # Step 3: Calculate influence functions using Tong's function
  phase1_data$inflBX1 <- inf_fun_logit(fit_phase1)[,"X1_obs"]
  phase1_data$inflBX2 <- inf_fun_logit(fit_phase1)[,"X2_obs"]
  
  # Step 4: Determine optimum allocation for wave 1 with Wright algorithm
  allocation1 <- optimum_allocation(phase1_data,
                                    strata = "strata",
                                    y = "inflBX1",
                                    nsample = n/8, method = "Neyman")
  
  allocation2 <- optimum_allocation(phase1_data,
                                    strata = "strata",
                                    y = "inflBX2",
                                    nsample = n/8, method = "Neyman")
  
  # Step 5: Combine allocations for total allocation of n/4
  wave1_allocation <- dplyr::left_join(allocation1, allocation2, by = "strata") |>
    dplyr::mutate(stratum_size = stratum_size.x + stratum_size.y)
  
  # Also, remove inflBX1 and inflBX2 vars
  phase1_data <- subset(phase1_data, select = -c(inflBX1,
                                                 inflBX2))
  
  # First 1/4th according to wave1_allocation
  phase2_wave1 <- sample_strata(data = phase1_data,
                                strata = "strata",
                                id = "id",
                                design_data = wave1_allocation,
                                n_allocated = "stratum_size")
  names(phase2_wave1)[names(phase2_wave1) == "sample_indicator"] <- "sampled_wave1"
  
  # Sample
  phase2_wave1$X1 <- ifelse(phase2_wave1$sampled_wave1 == 1, full_data$X1 , NA)
  phase2_wave1$X2 <- ifelse(phase2_wave1$sampled_wave1 == 1, full_data$X2, NA)
  phase2_wave1$Z <- ifelse(phase2_wave1$sampled_wave1 == 1, full_data$Z , NA)
  
  #####
  ## Wave 2
  #####
  
  # Calculate influence functions using IPW
  twophase_design_wave1 <- twophase(id = list(~1, ~1), 
                                    strata = list(NULL, ~strata), 
                                    subset = ~as.logical(sampled_wave1), 
                                    data = phase2_wave1)
  fit_wave1 <- svyglm(Y ~ X1 + X2 + Z, family = quasibinomial, 
                      design = twophase_design_wave1)
  infl_wave1 <- inf_fun_logit(fit_wave1)
  infl_wave1 <- data.frame(as.numeric(rownames(infl_wave1)), 
                           infl_wave1[,"X1"],
                           infl_wave1[,"X2"])
  names(infl_wave1) <- c("id", "inflBX1", "inflBX2")
  phase2_wave1 <- dplyr::left_join(phase2_wave1, infl_wave1, by = "id")
  
  
  # Re-calculate allocations
  allocation1 <- allocate_wave(phase2_wave1,
                               strata = "strata",
                               y = "inflBX1", method = "iterative",
                               already_sampled = "sampled_wave1",
                               nsample = n/8, 
                               allocation_method = "Neyman")
  
  allocation2 <- allocate_wave(phase2_wave1,
                               strata = "strata",
                               y = "inflBX2",
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
  phase2_wave2$X1<- ifelse(phase2_wave2$sampled_wave1 == 1 |
                             phase2_wave2$sampled_wave2 == 1, full_data$X1 , NA)
  phase2_wave2$X2 <- ifelse(phase2_wave2$sampled_wave1 == 1|
                              phase2_wave2$sampled_wave2 == 1, full_data$X2, NA)
  phase2_wave2$Z <- ifelse(phase2_wave2$sampled_wave1 == 1|
                             phase2_wave2$sampled_wave2 == 1, full_data$Z , NA)
  
  # Also, remove inflBX1 and inflBX2 vars
  phase2_wave2 <- subset(phase2_wave2, select = -c(inflBX1,
                                                   inflBX2))
  
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
  fit_wave2 <- svyglm(Y ~ X1 + X2 + Z, family = quasibinomial, 
                      design = twophase_design_wave2)
  infl_wave2 <- inf_fun_logit(fit_wave2)
  infl_wave2 <- data.frame(as.numeric(rownames(infl_wave2)), 
                           infl_wave2[,"X1"],
                           infl_wave2[,"X2"])
  names(infl_wave2) <- c("id", "inflBX1", "inflBX2")
  phase2_wave2 <- dplyr::left_join(phase2_wave2, infl_wave2, by = "id")
  
  
  # Re-calculate allocations
  allocation1 <- allocate_wave(phase2_wave2,
                               strata = "strata",
                               y = "inflBX1",
                               method = "iterative",
                               nsample = n/8,
                               already_sampled = "already_sampled",
                               allocation_method = "Neyman")
  
  allocation2 <- allocate_wave(phase2_wave2,
                               strata = "strata",
                               y = "inflBX2",
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
  phase2_wave3$X1<- ifelse(phase2_wave3$sampled_wave1 == 1 |
                             phase2_wave3$sampled_wave2 == 1|
                             phase2_wave3$sampled_wave3 == 1, full_data$X1 , NA)
  phase2_wave3$X2 <- ifelse(phase2_wave3$sampled_wave1 == 1 |
                              phase2_wave3$sampled_wave2 == 1|
                              phase2_wave3$sampled_wave3 == 1, full_data$X2 , NA)
  phase2_wave3$Z <- ifelse(phase2_wave3$sampled_wave1 == 1|
                             phase2_wave3$sampled_wave2 == 1|
                             phase2_wave3$sampled_wave3 == 1, full_data$Z , NA)
  
  # Also, remove inflBX1 and inflBX2 vars
  phase2_wave3 <- subset(phase2_wave3, select = -c(inflBX1,
                                                   inflBX2))
  
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
  fit_wave3 <- svyglm(Y ~ X1 + X2 + Z, family = quasibinomial, 
                      design = twophase_design_wave3)
  infl_wave3 <- inf_fun_logit(fit_wave3)
  infl_wave3 <- data.frame(as.numeric(rownames(infl_wave3)), 
                           infl_wave3[,"X1"],
                           infl_wave3[,"X2"])
  names(infl_wave3) <- c("id", "inflBX1", "inflBX2")
  phase2_wave3 <- dplyr::left_join(phase2_wave3, infl_wave3, by = "id")
  
  
  # Re-calculate allocations
  allocation1 <- allocate_wave(phase2_wave3,
                               strata = "strata",
                               y = "inflBX1",
                               method = "iterative",
                               nsample = n/8, 
                               already_sampled = "already_sampled",
                               detailed = TRUE,
                               allocation_method = "Neyman")
  
  allocation2 <- allocate_wave(phase2_wave3,
                               strata = "strata",
                               y = "inflBX2",
                               method = "iterative",
                               nsample = n/8, 
                               already_sampled = "already_sampled",
                               detailed = TRUE,
                               allocation_method = "Neyman")
  
  # Set indicators for oversampling
  oversampled_X1 <- ifelse(all(allocation1$nsample_optimal == 
                                 allocation1$nsample_actual),
                           0, 1)
  oversampled_X2 <- ifelse(all(allocation2$nsample_optimal == 
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
  
  # Standardize observed covariates wave4
  phase2_wave4$X1<- ifelse(phase2_wave4$sampled_wave1 == 1 |
                             phase2_wave4$sampled_wave2 == 1|
                             phase2_wave4$sampled_wave3 == 1|
                             phase2_wave4$sampled_wave4 == 1, full_data$X1 , NA)
  phase2_wave4$X2<- ifelse(phase2_wave4$sampled_wave1 == 1 |
                             phase2_wave4$sampled_wave2 == 1|
                             phase2_wave4$sampled_wave3 == 1|
                             phase2_wave4$sampled_wave4 == 1, full_data$X2 , NA)
  phase2_wave4$Z <- ifelse(phase2_wave4$sampled_wave1 == 1|
                             phase2_wave4$sampled_wave2 == 1|
                             phase2_wave4$sampled_wave3 == 1|
                             phase2_wave4$sampled_wave4 == 1, full_data$Z , NA)
  
  # Remove inflBX1 and inflBX2 vars
  phase2_wave4 <- subset(phase2_wave4, select = -c(inflBX1,
                                                   inflBX2))
  
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
  weight <- svyglm(Y ~ X1 + X2 + Z, family = quasibinomial, 
                   design = twophase_design)
  
  #########
  #### Raking approach: Sentinel Vignette strategy
  #########
  data <- phase2_wave4
  
  # Creating a survey object and calibrating the weights
  mydesign <- survey::twophase(
    id = list(~1, ~1), subset = ~as.logical(already_sampled),
    strata = list(NULL, ~strata), data = data
  )
  infcal <- survey::calibrate(mydesign, formula = ~ inflBX1_phase1 + inflBX2_phase1
                              + strata, phase = 2, calfun = "raking")
  
  # Fitting the outcome model: conditional treatment effect of interest.
  fit <- survey::svyglm(Y ~ X1 + X2 + Z, design = infcal, family = quasibinomial)
  
  # Check if truth is contained in interval
  confintIPWBX1 <- confint(weight)["X1",]
  confintIPWBX2 <- confint(weight)["X2",]
  confintGRBX1 <- confint(fit)["X1",]
  confintGRBX2 <- confint(fit)["X2",]
  coverIPWBX1 <- ifelse(B1 >= confintIPWBX1[1] & B1 <= confintIPWBX1[2],1,0)
  coverIPWBX2 <- ifelse(B2 >= confintIPWBX2[1] & B2 <= confintIPWBX2[2],1,0)
  coverGRBX1 <- ifelse(B1 >= confintGRBX1[1] & B1 <= confintGRBX1[2],1,0)
  coverGRBX2 <- ifelse(B2 >= confintGRBX2[1] & B2 <= confintGRBX2[2],1,0)
  
  output <- c(fit$coefficients["X1"], fit$coefficients["X2"],
              coef(weight)["X1"], coef(weight)["X2"], SE(fit)["X1"],
              SE(fit)["X2"], SE(weight)["X1"], SE(weight)["X2"], 
              coverIPWBX1, coverIPWBX2, coverGRBX1, coverGRBX2,
              oversampled_X1, oversampled_X2, cor_X1_X2,
              prev_Y, prevZ, B1, B2)
  return(output)
}

#####
## Strategy 3 - Independent sequential
#####

# This is looped part, so Wave 1 allocation is required as input
Strategy3 <- function(n, simulations_df, scenario, N = 10000){
  
  # Set up betas for scenario
  B0 <- simulations_df["B0", scenario]
  B1 <- simulations_df["B1", scenario]
  B2 <- simulations_df["B2", scenario]
  B3 <- simulations_df["B3", scenario]
  
  # Set up of X's
  corX1X2 <- simulations_df["corX1X2", scenario]
  prevZ <- simulations_df["prevZ", scenario]
  
  # Measurement error
  sensY <- simulations_df["sensY", scenario]
  specY <- simulations_df["specY", scenario]
  error_varX1 <- simulations_df["error_varX1", scenario]
  error_varX2 <- simulations_df["error_varX2", scenario]
  sensZ <- simulations_df["sensZ", scenario]
  specZ <- simulations_df["specZ", scenario]
  
  
  
  #####
  ## Generate Phase 1 Data
  #####
  
  id <- 1:N
  
  # Generate correlated covariates
  sigma <- matrix(c(1, corX1X2 , 0.1, corX1X2, 1,  0.25, 0.1, 
                    0.25, 1), nrow = 3)
  covs <- mvrnorm(N, mu = c(0,0,0), Sigma = sigma)
  X1 <- covs[,1]
  X2 <- covs[,2]
  Z <- ifelse(covs[,3] >= quantile(covs[,3], 1- prevZ), 1, 0)
  cor(cbind(X1, X2, Z)) 
  
  Y_probs <- exp(B0 + B1*X1 + B2*X2 + B3*Z)/(1 +  exp(B0 + B1*X1 + B2*X2 + B3*Z))
  Y <- rbinom(N, 1, Y_probs) # Compute realized value for Y2
  
  full_data <- data.frame(id, X1, X2, Z, Y)
  
  ## Test true model
  true_model <- glm(Y ~ X1 + X2 + Z, family = "binomial", data = full_data)
  
  # Find observed correlation
  cor_X1_X2 <- cor(X1, X2) 
  
  # Prevalence
  prev_Y <- table(Y)["1"]/length(Y)
  
  # Find true coefficients
  true_model <- glm(Y ~ X1 + X2 + Z, family = "binomial", data = full_data) 
  true_BX1 <- coef(true_model)["X1"]
  true_BX2 <- coef(true_model)["X2"]
  
  
  #####
  ## Generate error-prone Phase 1 data
  ######
  
  Z_prob_obs1 <- ifelse(full_data$Z > 0, sensZ, 1 - specZ)
  full_data$Z_obs <- rbinom(N, 1, Z_prob_obs1)
  XZY_obs <- mvrnorm(N, c(0,0,0,0),
                     matrix(c(error_varX1, 0.03, 0.02, 0.02,
                              0.03, error_varX2, 0.05, 0.4,
                              0.02, 0.05, 1, 0,
                              0.02, 0.4, 0, 1), nrow = 4))
  full_data$X1_obs <- full_data$X1 + XZY_obs[,1]
  full_data$X2_obs <- full_data$X2 + XZY_obs[,2]
  
  threshold_positiveY <- qnorm((sensY + 1)/2)
  threshold_negativeY <- qnorm((specY + 1)/2)
  threshold_positiveZ <- qnorm((sensZ + 1)/2)
  threshold_negativeZ <- qnorm((specZ + 1)/2)
  full_data$Y_obs <- with(full_data, ifelse((Y == 1 & abs(XZY_obs[, 4]) < threshold_positiveY) | 
                                              (Y == 0 & abs(XZY_obs[, 4]) > threshold_negativeY), 1, 0))
  
  full_data$Z_obs <- with(full_data, ifelse((Z == 1 & abs(XZY_obs[, 3]) < threshold_positiveZ) | 
                                              (Z == 0 & abs(XZY_obs[, 3]) > threshold_negativeZ),
                                            max(Z), min(Z)))
  
  ####
  #### Divide into strata based on observed phase 1 data
  
  #### Compute naive influence functions
  fit_phase1 <-  glm(Y_obs ~ X1_obs + X2_obs + Z_obs, 
                     family = "binomial", data = full_data)
  full_data$inflBX1_phase1 <- inf_fun_logit(fit_phase1)[,"X1_obs"]
  full_data$inflBX2_phase1 <- inf_fun_logit(fit_phase1)[,"X2_obs"]
  
  #### Divide into strata based on observed phase 1 data X
  full_data <- full_data |>
    dplyr::mutate(Y_strat = Y_obs) |>
    split_strata(strata = "Y_strat",
                 split_var = "X1_obs",
                 type = "local quantile",
                 split_at = 0.5,
                 trunc = "X1") |>
    split_strata(strata = "new_strata",
                 split_var = "X2_obs",
                 type = "local quantile",
                 split_at = 0.5,
                 trunc = "X2")
  
  
  names(full_data)[names(full_data) == "new_strata"] <- "strata"
  
  data <- full_data
  
  #####
  #####
  ##### Phase 2
  
  ## Here, each wave chooses n/4 samples where n/8 are chosen to be optimal for
  ## BX1 and n/8 to be optimal for BX2.
  
  ####
  ## Wave 1
  ####
  
  # Step 1: Initialize Phase 1 data
  phase1_data <- full_data
  
  # Step 2: Logistic regression on Phase 1 data (Y1 only needed)
  fit_phase1 <-  glm(Y_obs ~ X1_obs + X2_obs + Z_obs, 
                     family = "binomial", data = phase1_data)
  
  # Step 3: Calculate influence functions using Tong's function
  phase1_data$inflBX1 <- inf_fun_logit(fit_phase1)[,"X1_obs"]
  
  # Step 4: Determine optimum allocation for wave 1 with Wright algorithm
  wave1_allocation <- optimum_allocation(phase1_data,
                                         strata = "strata",
                                         y = "inflBX1",
                                         nsample = n/4, method = "Neyman")
  
  # Step5: Now can remove inflBX1 column, 
  phase1_data <- subset(phase1_data, select = -c(inflBX1))
  
  
  # Sample First n/4 according to wave1_allocation
  phase2_wave1 <- sample_strata(data = phase1_data,
                                strata = "strata",
                                id = "id",
                                design_data = wave1_allocation,
                                n_allocated = "stratum_size")
  names(phase2_wave1)[names(phase2_wave1) == "sample_indicator"] <- "sampled_wave1"
  
  # Sample
  phase2_wave1$X1 <- ifelse(phase2_wave1$sampled_wave1 == 1, full_data$X1 , NA)
  phase2_wave1$X2 <- ifelse(phase2_wave1$sampled_wave1 == 1, full_data$X2, NA)
  phase2_wave1$Z <- ifelse(phase2_wave1$sampled_wave1 == 1, full_data$Z , NA)
  
  #####
  ## Wave 2
  #####
  
  # Calculate influence functions using IPW
  twophase_design_wave1 <- twophase(id = list(~1, ~1), 
                                    strata = list(NULL, ~strata), 
                                    subset = ~as.logical(sampled_wave1), 
                                    data = phase2_wave1)
  fit_wave1 <- svyglm(Y ~ X1 + X2 + Z, family = quasibinomial, 
                      design = twophase_design_wave1)
  infl_wave1 <- inf_fun_logit(fit_wave1)
  infl_wave1 <- data.frame(as.numeric(rownames(infl_wave1)), 
                           infl_wave1[,"X1"])
  names(infl_wave1) <- c("id", "inflBX1")
  phase2_wave1 <- dplyr::left_join(phase2_wave1, infl_wave1, by = "id")
  
  
  # Calculate allocation for Y1
  wave2_allocation <- allocate_wave(phase2_wave1,
                                    strata = "strata",
                                    y = "inflBX1",
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
  phase2_wave2$X1<- ifelse(phase2_wave2$sampled_wave1 == 1 |
                             phase2_wave2$sampled_wave2 == 1, full_data$X1 , NA)
  phase2_wave2$X2<- ifelse(phase2_wave2$sampled_wave1 == 1 |
                             phase2_wave2$sampled_wave2 == 1, full_data$X2 , NA)
  phase2_wave2$Z <- ifelse(phase2_wave2$sampled_wave1 == 1|
                             phase2_wave2$sampled_wave2 == 1, full_data$Z , NA)
  
  phase2_wave2 <- subset(phase2_wave2, select = -c(inflBX1))
  
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
  fit_wave2 <- svyglm(Y ~ X1 + X2 + Z, family = quasibinomial, 
                      design = twophase_design_wave2)
  infl_wave2 <- inf_fun_logit(fit_wave2)
  infl_wave2 <- data.frame(as.numeric(rownames(infl_wave2)), 
                           infl_wave2[,"X2"])
  names(infl_wave2) <- c("id", "inflBX2")
  phase2_wave2 <- dplyr::left_join(phase2_wave2, infl_wave2, by = "id")
  
  
  # Re-calculate allocation
  wave3_allocation <- allocate_wave(phase2_wave2,
                                    strata = "strata",
                                    y = "inflBX2",
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
  phase2_wave3$X1<- ifelse(phase2_wave3$sampled_wave1 == 1 |
                             phase2_wave3$sampled_wave2 == 1|
                             phase2_wave3$sampled_wave3 == 1, full_data$X1 , NA)
  phase2_wave3$X2<- ifelse(phase2_wave3$sampled_wave1 == 1 |
                             phase2_wave3$sampled_wave2 == 1|
                             phase2_wave3$sampled_wave3 == 1, full_data$X2 , NA)
  phase2_wave3$Z <- ifelse(phase2_wave3$sampled_wave1 == 1|
                             phase2_wave3$sampled_wave2 == 1|
                             phase2_wave3$sampled_wave3 == 1, full_data$Z , NA)
  
  phase2_wave3 <- subset(phase2_wave3, select = -c(inflBX2))
  
  
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
  fit_wave3 <- svyglm(Y ~ X1 + X2 + Z, family = quasibinomial, 
                      design = twophase_design_wave3)
  infl_wave3 <- inf_fun_logit(fit_wave3)
  infl_wave3 <- data.frame(as.numeric(rownames(infl_wave3)),
                           infl_wave3[,"X2"])
  names(infl_wave3) <- c("id", "inflBX2")
  phase2_wave3 <- dplyr::left_join(phase2_wave3, infl_wave3, by = "id")
  
  
  # Re-calculate allocation
  wave4_allocation <- allocate_wave(phase2_wave3,
                                    strata = "strata",
                                    y = "inflBX2",
                                    nsample = n/4, method = "iterative",
                                    already_sampled = "already_sampled",
                                    detailed = TRUE,
                                    allocation_method = "Neyman")
  
  # Check for oversampling
  oversampled_X2 <- ifelse(all(wave4_allocation$nsample_optimal == 
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
  
  # Standardize observed variables wave4
  phase2_wave4$X1<- ifelse(phase2_wave4$sampled_wave1 == 1 |
                             phase2_wave4$sampled_wave2 == 1|
                             phase2_wave4$sampled_wave3 == 1|
                             phase2_wave4$sampled_wave4 == 1, full_data$X1 , NA)
  phase2_wave4$X2<- ifelse(phase2_wave4$sampled_wave1 == 1 |
                             phase2_wave4$sampled_wave2 == 1|
                             phase2_wave4$sampled_wave3 == 1|
                             phase2_wave4$sampled_wave4 == 1, full_data$X2 , NA)
  phase2_wave4$Z <- ifelse(phase2_wave4$sampled_wave1 == 1|
                             phase2_wave4$sampled_wave2 == 1|
                             phase2_wave4$sampled_wave3 == 1|
                             phase2_wave4$sampled_wave4 == 1, full_data$Z , NA)
  
  phase2_wave4 <- subset(phase2_wave4, select = -c(inflBX2))
  
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
  weight <- svyglm(Y ~ X1 + X2 + Z, family = quasibinomial, 
                   design = twophase_design)
  
  #########
  #### Raking 
  #########
  data <- phase2_wave4
  
  # Creating a survey object and calibrating the weights
  mydesign <- survey::twophase(
    id = list(~1, ~1), subset = ~as.logical(already_sampled),
    strata = list(NULL, ~strata), data = data
  )
  infcal <- survey::calibrate(mydesign, formula = ~ inflBX1_phase1 + inflBX2_phase1
                              + strata, phase = 2, calfun = "raking")
  
  # Fitting the outcome model: conditional treatment effect of interest.
  fit <- survey::svyglm(Y ~ X1 + X2 + Z, design = infcal, family = quasibinomial)
  
  # Check if truth is contained in interval
  confintIPWBX1 <- confint(weight)["X1",]
  confintIPWBX2 <- confint(weight)["X2",]
  confintGRBX1 <- confint(fit)["X1",]
  confintGRBX2 <- confint(fit)["X2",]
  coverIPWBX1 <- ifelse(B1 >= confintIPWBX1[1] & B1 <= confintIPWBX1[2],1,0)
  coverIPWBX2 <- ifelse(B2 >= confintIPWBX2[1] & B2 <= confintIPWBX2[2],1,0)
  coverGRBX1 <- ifelse(B1 >= confintGRBX1[1] & B1 <= confintGRBX1[2],1,0)
  coverGRBX2 <- ifelse(B2 >= confintGRBX2[1] & B2 <= confintGRBX2[2],1,0)
  
  output <- c(fit$coefficients["X1"], fit$coefficients["X2"],
              coef(weight)["X1"], coef(weight)["X2"], SE(fit)["X1"],
              SE(fit)["X2"], SE(weight)["X1"], SE(weight)["X2"], 
              coverIPWBX1, coverIPWBX2, coverGRBX1, coverGRBX2, oversampled_X2,
              cor_X1_X2,
              prev_Y, prevZ, B1, B2)
  return(output)
}


#####
## Strategy 3.5 - Independent sequential with rare outcome first
#####

# This is looped part, so Wave 1 allocation is required as input
Strategy3.5 <- function(n, simulations_df, scenario, N = 10000){
  
  # Set up betas for scenario
  B0 <- simulations_df["B0", scenario]
  B1 <- simulations_df["B1", scenario]
  B2 <- simulations_df["B2", scenario]
  B3 <- simulations_df["B3", scenario]
  
  # Set up of X's
  corX1X2 <- simulations_df["corX1X2", scenario]
  prevZ <- simulations_df["prevZ", scenario]
  
  # Measurement error
  sensY <- simulations_df["sensY", scenario]
  specY <- simulations_df["specY", scenario]
  error_varX1 <- simulations_df["error_varX1", scenario]
  error_varX2 <- simulations_df["error_varX2", scenario]
  sensZ <- simulations_df["sensZ", scenario]
  specZ <- simulations_df["specZ", scenario]
  
  
  
  #####
  ## Generate Phase 1 Data
  #####
  
  id <- 1:N
  
  # Generate correlated covariates
  sigma <- matrix(c(1, corX1X2 , 0.1, corX1X2, 1,  0.25, 0.1, 
                    0.25, 1), nrow = 3)
  covs <- mvrnorm(N, mu = c(0,0,0), Sigma = sigma)
  X1 <- covs[,1]
  X2 <- covs[,2]
  Z <- ifelse(covs[,3] >= quantile(covs[,3], 1- prevZ), 1, 0)
  cor(cbind(X1, X2, Z)) 
  
  Y_probs <- exp(B0 + B1*X1 + B2*X2 + B3*Z)/(1 +  exp(B0 + B1*X1 + B2*X2 + B3*Z))
  Y <- rbinom(N, 1, Y_probs) # Compute realized value for Y2
  
  full_data <- data.frame(id, X1, X2, Z, Y)
  
  ## Test true model
  true_model <- glm(Y ~ X1 + X2 + Z, family = "binomial", data = full_data)
  
  # Find observed correlation
  cor_X1_X2 <- cor(X1, X2) 
  
  # Prevalence
  prev_Y <- table(Y)["1"]/length(Y)
  
  # Find true coefficients
  true_model <- glm(Y ~ X1 + X2 + Z, family = "binomial", data = full_data) 
  true_BX1 <- coef(true_model)["X1"]
  true_BX2 <- coef(true_model)["X2"]
  
  
  #####
  ## Generate error-prone Phase 1 data
  ######
  
  Z_prob_obs1 <- ifelse(full_data$Z > 0, sensZ, 1 - specZ)
  full_data$Z_obs <- rbinom(N, 1, Z_prob_obs1)
  XZY_obs <- mvrnorm(N, c(0,0,0,0),
                     matrix(c(error_varX1, 0.03, 0.02, 0.02,
                              0.03, error_varX2, 0.05, 0.4,
                              0.02, 0.05, 1, 0,
                              0.02, 0.4, 0, 1), nrow = 4))
  full_data$X1_obs <- full_data$X1 + XZY_obs[,1]
  full_data$X2_obs <- full_data$X2 + XZY_obs[,2]
  
  threshold_positiveY <- qnorm((sensY + 1)/2)
  threshold_negativeY <- qnorm((specY + 1)/2)
  threshold_positiveZ <- qnorm((sensZ + 1)/2)
  threshold_negativeZ <- qnorm((specZ + 1)/2)
  full_data$Y_obs <- with(full_data, ifelse((Y == 1 & abs(XZY_obs[, 4]) < threshold_positiveY) | 
                                              (Y == 0 & abs(XZY_obs[, 4]) > threshold_negativeY), 1, 0))
  
  full_data$Z_obs <- with(full_data, ifelse((Z == 1 & abs(XZY_obs[, 3]) < threshold_positiveZ) | 
                                              (Z == 0 & abs(XZY_obs[, 3]) > threshold_negativeZ),
                                            max(Z), min(Z)))
  
  ####
  #### Divide into strata based on observed phase 1 data
  
  #### Compute naive influence functions
  fit_phase1 <-  glm(Y_obs ~ X1_obs + X2_obs + Z_obs, 
                     family = "binomial", data = full_data)
  full_data$inflBX1_phase1 <- inf_fun_logit(fit_phase1)[,"X1_obs"]
  full_data$inflBX2_phase1 <- inf_fun_logit(fit_phase1)[,"X2_obs"]
  
  #### Divide into strata based on observed phase 1 data X
  full_data <- full_data |>
    dplyr::mutate(Y_strat = Y_obs) |>
    split_strata(strata = "Y_strat",
                 split_var = "X1_obs",
                 type = "local quantile",
                 split_at = 0.5,
                 trunc = "X1") |>
    split_strata(strata = "new_strata",
                 split_var = "X2_obs",
                 type = "local quantile",
                 split_at = 0.5,
                 trunc = "X2")
  
  
  names(full_data)[names(full_data) == "new_strata"] <- "strata"
  
  data <- full_data
  
  ####
  ## Phase 2
  ####
  
  ####
  ## Wave 1
  ####
  
  ## Here, 1st and 2nd waves (n/4 samples each) are optimized for B_X2 and 
  ## 3rd and 4th last n/4 are optimized for B_X1.
  
  # Step 1: Initialize Phase 1 data
  phase1_data <- full_data
  
  # Step 2: Logistic regression on Phase 1 data
  fit_phase1 <-  glm(Y_obs ~ X1_obs + X2_obs + Z_obs, 
                     family = "binomial", data = phase1_data)
  
  # Step 3: Calculate influence functions using Tong's function
  phase1_data$inflBX2 <- inf_fun_logit(fit_phase1)[,"X2_obs"]
  
  # Step 4: Determine optimum allocation for wave 1 with Wright algorithm,
  # which ensures two observations per stratum
  wave1_allocation <- optimum_allocation(phase1_data,
                                         strata = "strata",
                                         y = "inflBX2",
                                         nsample = n/4, method = "WrightII")
  
  # Step 6a: Now can remove inflBX1 column, 
  phase1_data <- subset(phase1_data, select = -c(inflBX2))
  
  
  # First n/4 according to wave1_allocation
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
  
  # Calculate influence functions using IPW
  twophase_design_wave1 <- twophase(id = list(~1, ~1), 
                                    strata = list(NULL, ~strata), 
                                    subset = ~as.logical(sampled_wave1), 
                                    data = phase2_wave1)
  fit_wave1 <- svyglm(Y ~ X1 + X2 + Z, family = quasibinomial, 
                      design = twophase_design_wave1)
  infl_wave1 <- inf_fun_logit(fit_wave1)
  infl_wave1 <- data.frame(as.numeric(rownames(infl_wave1)), 
                           infl_wave1[,"X2"])
  names(infl_wave1) <- c("id", "inflBX2")
  phase2_wave1 <- dplyr::left_join(phase2_wave1, infl_wave1, by = "id")
  
  
  # Calculate allocation for Y1
  wave2_allocation <- allocate_wave(phase2_wave1,
                                    strata = "strata",
                                    y = "inflBX2",
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
  phase2_wave2$X1<- ifelse(phase2_wave2$sampled_wave1 == 1 |
                             phase2_wave2$sampled_wave2 == 1, full_data$X1 , NA)
  phase2_wave2$X2<- ifelse(phase2_wave2$sampled_wave1 == 1 |
                             phase2_wave2$sampled_wave2 == 1, full_data$X2 , NA)
  phase2_wave2$Z <- ifelse(phase2_wave2$sampled_wave1 == 1|
                             phase2_wave2$sampled_wave2 == 1, full_data$Z , NA)
  
  phase2_wave2 <- subset(phase2_wave2, select = -c(inflBX2))
  
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
  fit_wave2 <- svyglm(Y ~ X1 + X2 + Z, family = quasibinomial, 
                      design = twophase_design_wave2)
  infl_wave2 <- inf_fun_logit(fit_wave2)
  infl_wave2 <- data.frame(as.numeric(rownames(infl_wave2)), 
                           infl_wave2[,"X1"])
  names(infl_wave2) <- c("id", "inflBX1")
  phase2_wave2 <- dplyr::left_join(phase2_wave2, infl_wave2, by = "id")
  
  
  # Re-calculate allocation
  wave3_allocation <- allocate_wave(phase2_wave2,
                                    strata = "strata",
                                    y = "inflBX1",
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
  phase2_wave3$X1<- ifelse(phase2_wave3$sampled_wave1 == 1 |
                             phase2_wave3$sampled_wave2 == 1|
                             phase2_wave3$sampled_wave3 == 1, full_data$X1 , NA)
  phase2_wave3$X2<- ifelse(phase2_wave3$sampled_wave1 == 1 |
                             phase2_wave3$sampled_wave2 == 1|
                             phase2_wave3$sampled_wave3 == 1, full_data$X2 , NA)
  phase2_wave3$Z <- ifelse(phase2_wave3$sampled_wave1 == 1|
                             phase2_wave3$sampled_wave2 == 1|
                             phase2_wave3$sampled_wave3 == 1, full_data$Z , NA)
  
  phase2_wave3 <- subset(phase2_wave3, select = -c(inflBX1))
  
  
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
  fit_wave3 <- svyglm(Y ~ X1 + X2 + Z, family = quasibinomial, 
                      design = twophase_design_wave3)
  infl_wave3 <- inf_fun_logit(fit_wave3)
  infl_wave3 <- data.frame(as.numeric(rownames(infl_wave3)),
                           infl_wave3[,"X1"])
  names(infl_wave3) <- c("id", "inflBX1")
  phase2_wave3 <- dplyr::left_join(phase2_wave3, infl_wave3, by = "id")
  
  
  # Re-calculate allocation
  wave4_allocation <- allocate_wave(phase2_wave3,
                                    strata = "strata",
                                    y = "inflBX1",
                                    nsample = n/4, method = "iterative",
                                    already_sampled = "already_sampled",
                                    detailed = TRUE,
                                    allocation_method = "Neyman")
  
  # Check for oversampling
  oversampled_X1 <- ifelse(all(wave4_allocation$nsample_optimal == 
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
  phase2_wave4$X1<- ifelse(phase2_wave4$sampled_wave1 == 1 |
                             phase2_wave4$sampled_wave2 == 1|
                             phase2_wave4$sampled_wave3 == 1|
                             phase2_wave4$sampled_wave4 == 1, full_data$X1 , NA)
  phase2_wave4$X2<- ifelse(phase2_wave4$sampled_wave1 == 1 |
                             phase2_wave4$sampled_wave2 == 1|
                             phase2_wave4$sampled_wave3 == 1|
                             phase2_wave4$sampled_wave4 == 1, full_data$X2 , NA)
  phase2_wave4$Z <- ifelse(phase2_wave4$sampled_wave1 == 1|
                             phase2_wave4$sampled_wave2 == 1|
                             phase2_wave4$sampled_wave3 == 1|
                             phase2_wave4$sampled_wave4 == 1, full_data$Z , NA)
  
  phase2_wave4 <- subset(phase2_wave4, select = -c(inflBX1))
  
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
  weight <- svyglm(Y ~ X1 + X2 + Z, family = quasibinomial, 
                   design = twophase_design)
  
  #########
  #### Raking 
  #########
  data <- phase2_wave4
  
  # Creating a survey object and calibrating the weights
  mydesign <- survey::twophase(
    id = list(~1, ~1), subset = ~as.logical(already_sampled),
    strata = list(NULL, ~strata), data = data
  )
  infcal <- survey::calibrate(mydesign, formula = ~ inflBX1_phase1 + inflBX2_phase1
                              + strata, phase = 2, calfun = "raking")
  
  # Fitting the outcome model: conditional treatment effect of interest.
  fit <- survey::svyglm(Y ~ X1 + X2 + Z, design = infcal, family = quasibinomial)
  
  # Check if truth is contained in interval
  confintIPWBX1 <- confint(weight)["X1",]
  confintIPWBX2 <- confint(weight)["X2",]
  confintGRBX1 <- confint(fit)["X1",]
  confintGRBX2 <- confint(fit)["X2",]
  coverIPWBX1 <- ifelse(B1 >= confintIPWBX1[1] & B1 <= confintIPWBX1[2],1,0)
  coverIPWBX2 <- ifelse(B2 >= confintIPWBX2[1] & B2 <= confintIPWBX2[2],1,0)
  coverGRBX1 <- ifelse(B1 >= confintGRBX1[1] & B1 <= confintGRBX1[2],1,0)
  coverGRBX2 <- ifelse(B2 >= confintGRBX2[1] & B2 <= confintGRBX2[2],1,0)
  
  output <- c(fit$coefficients["X1"], fit$coefficients["X2"],
              coef(weight)["X1"], coef(weight)["X2"], SE(fit)["X1"],
              SE(fit)["X2"], SE(weight)["X1"], SE(weight)["X2"], 
              coverIPWBX1, coverIPWBX2, coverGRBX1, coverGRBX2,
              oversampled_X1, cor_X1_X2,
              prev_Y, prevZ, B1, B2)
  return(output)
}


######
### Strategy 4 - A-optimal sampling
#####

# Step 5: Create function to estimate beta-hat at each iteration
Strategy4 <- function(n, simulations_df, scenario, N = 10000){
  
  # Set up betas for scenario
  B0 <- simulations_df["B0", scenario]
  B1 <- simulations_df["B1", scenario]
  B2 <- simulations_df["B2", scenario]
  B3 <- simulations_df["B3", scenario]
  
  # Set up of X's
  corX1X2 <- simulations_df["corX1X2", scenario]
  prevZ <- simulations_df["prevZ", scenario]
  
  # Measurement error
  sensY <- simulations_df["sensY", scenario]
  specY <- simulations_df["specY", scenario]
  error_varX1 <- simulations_df["error_varX1", scenario]
  error_varX2 <- simulations_df["error_varX2", scenario]
  sensZ <- simulations_df["sensZ", scenario]
  specZ <- simulations_df["specZ", scenario]
  
  
  
  #####
  ## Generate Phase 1 Data
  #####
  
  id <- 1:N
  
  # Generate correlated covariates
  sigma <- matrix(c(1, corX1X2 , 0.1, corX1X2, 1,  0.25, 0.1, 
                    0.25, 1), nrow = 3)
  covs <- mvrnorm(N, mu = c(0,0,0), Sigma = sigma)
  X1 <- covs[,1]
  X2 <- covs[,2]
  Z <- ifelse(covs[,3] >= quantile(covs[,3], 1- prevZ), 1, 0)
  cor(cbind(X1, X2, Z)) 
  
  Y_probs <- exp(B0 + B1*X1 + B2*X2 + B3*Z)/(1 +  exp(B0 + B1*X1 + B2*X2 + B3*Z))
  Y <- rbinom(N, 1, Y_probs) # Compute realized value for Y2
  
  full_data <- data.frame(id, X1, X2, Z, Y)
  
  ## Test true model
  true_model <- glm(Y ~ X1 + X2 + Z, family = "binomial", data = full_data)
  
  # Find observed correlation
  cor_X1_X2 <- cor(X1, X2) 
  
  # Prevalence
  prev_Y <- table(Y)["1"]/length(Y)
  
  # Find true coefficients
  true_model <- glm(Y ~ X1 + X2 + Z, family = "binomial", data = full_data) 
  true_BX1 <- coef(true_model)["X1"]
  true_BX2 <- coef(true_model)["X2"]
  
  
  #####
  ## Generate error-prone Phase 1 data
  ######
  
  Z_prob_obs1 <- ifelse(full_data$Z > 0, sensZ, 1 - specZ)
  full_data$Z_obs <- rbinom(N, 1, Z_prob_obs1)
  XZY_obs <- mvrnorm(N, c(0,0,0,0),
                     matrix(c(error_varX1, 0.03, 0.02, 0.02,
                              0.03, error_varX2, 0.05, 0.4,
                              0.02, 0.05, 1, 0,
                              0.02, 0.4, 0, 1), nrow = 4))
  full_data$X1_obs <- full_data$X1 + XZY_obs[,1]
  full_data$X2_obs <- full_data$X2 + XZY_obs[,2]
  
  threshold_positiveY <- qnorm((sensY + 1)/2)
  threshold_negativeY <- qnorm((specY + 1)/2)
  threshold_positiveZ <- qnorm((sensZ + 1)/2)
  threshold_negativeZ <- qnorm((specZ + 1)/2)
  full_data$Y_obs <- with(full_data, ifelse((Y == 1 & abs(XZY_obs[, 4]) < threshold_positiveY) | 
                                              (Y == 0 & abs(XZY_obs[, 4]) > threshold_negativeY), 1, 0))
  
  full_data$Z_obs <- with(full_data, ifelse((Z == 1 & abs(XZY_obs[, 3]) < threshold_positiveZ) | 
                                              (Z == 0 & abs(XZY_obs[, 3]) > threshold_negativeZ),
                                            max(Z), min(Z)))
  
  ####
  #### Divide into strata based on observed phase 1 data
  
  #### Compute naive influence functions
  fit_phase1 <-  glm(Y_obs ~ X1_obs + X2_obs + Z_obs, 
                     family = "binomial", data = full_data)
  full_data$inflBX1_phase1 <- inf_fun_logit(fit_phase1)[,"X1_obs"]
  full_data$inflBX2_phase1 <- inf_fun_logit(fit_phase1)[,"X2_obs"]
  
  #### Divide into strata based on observed phase 1 data X
  full_data <- full_data |>
    dplyr::mutate(Y_strat = Y_obs) |>
    split_strata(strata = "Y_strat",
                 split_var = "X1_obs",
                 type = "local quantile",
                 split_at = 0.5,
                 trunc = "X1") |>
    split_strata(strata = "new_strata",
                 split_var = "X2_obs",
                 type = "local quantile",
                 split_at = 0.5,
                 trunc = "X2")
  
  
  names(full_data)[names(full_data) == "new_strata"] <- "strata"
  
  data <- full_data
  
  ####
  ## Phase 2
  ####
  
  ####
  ## Wave 1
  ####
  
  # Step 1: Initialize phase1_data
  phase1_data <- full_data
  
  # Step 2: Logistic regression on Phase 1 data
  fit_phase1 <-  glm(Y_obs ~ X1_obs + X2_obs + Z_obs, 
                     family = "binomial", data = phase1_data)
  
  # Step 3: Calculate influence functions using Tong's function
  phase1_data$inflBX1 <- inf_fun_logit(fit_phase1)[,"X1_obs"]
  phase1_data$inflBX2 <- inf_fun_logit(fit_phase1)[,"X2_obs"]
  
  # Step 4: Determine optimum allocation for wave 1 with Wright algorithm
  wave1_allocation <- a_optimum_allocation(phase1_data,
                                           strata = "strata",
                                           nsample = n/4,
                                           vars = c("inflBX1", "inflBX2"),
                                           weights = c(0.5, 0.5),
                                           method = "Neyman")
  
  # Also, remove inflBX1 and inflBX2 vars
  phase1_data <- subset(phase1_data, select = -c(inflBX1,
                                                 inflBX2))
  
  # First n/4 according to wave1_allocation
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
  
  # Calculate influence functions using IPW
  twophase_design_wave1 <- twophase(id = list(~1, ~1), 
                                    strata = list(NULL, ~strata), 
                                    subset = ~as.logical(sampled_wave1), 
                                    data = phase2_wave1)
  fit_wave1 <- svyglm(Y ~ X1 + X2 + Z, family = quasibinomial, 
                      design = twophase_design_wave1)
  infl_wave1 <- inf_fun_logit(fit_wave1)
  infl_wave1 <- data.frame(as.numeric(rownames(infl_wave1)), 
                           infl_wave1[,"X1"],
                           infl_wave1[,"X2"])
  names(infl_wave1) <- c("id", "inflBX1", "inflBX2")
  phase2_wave1 <- dplyr::left_join(phase2_wave1, infl_wave1, by = "id")
  
  
  # Calculate allocation for wave 2
  wave2_allocation <- a_optimal_allocate_wave(phase2_wave1,
                                              strata = "strata",
                                              vars = c("inflBX1", "inflBX2"),
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
  
  # Sample
  phase2_wave2$X1<- ifelse(phase2_wave2$sampled_wave1 == 1 |
                             phase2_wave2$sampled_wave2 == 1, full_data$X1 , NA)
  phase2_wave2$X2<- ifelse(phase2_wave2$sampled_wave1 == 1 |
                             phase2_wave2$sampled_wave2 == 1, full_data$X2 , NA)
  phase2_wave2$Z <- ifelse(phase2_wave2$sampled_wave1 == 1|
                             phase2_wave2$sampled_wave2 == 1, full_data$Z , NA)
  
  phase2_wave2 <- subset(phase2_wave2, select = -c(inflBX1,
                                                   inflBX2))
  
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
  fit_wave2 <- svyglm(Y ~ X1 + X2 + Z, family = quasibinomial, 
                      design = twophase_design_wave2)
  infl_wave2 <- inf_fun_logit(fit_wave2)
  infl_wave2 <- data.frame(as.numeric(rownames(infl_wave2)), 
                           infl_wave2[,"X1"],
                           infl_wave2[,"X2"])
  names(infl_wave2) <- c("id", "inflBX1", "inflBX2")
  phase2_wave2 <- dplyr::left_join(phase2_wave2, infl_wave2, by = "id")
  
  
  # Re-calculate allocations
  wave3_allocation <- a_optimal_allocate_wave(phase2_wave2,
                                              strata = "strata",
                                              vars = c("inflBX1", "inflBX2"),
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
  phase2_wave3$X1<- ifelse(phase2_wave3$sampled_wave1 == 1 |
                             phase2_wave3$sampled_wave2 == 1|
                             phase2_wave3$sampled_wave3 == 1, full_data$X1 , NA)
  phase2_wave3$X2<- ifelse(phase2_wave3$sampled_wave1 == 1 |
                             phase2_wave3$sampled_wave2 == 1|
                             phase2_wave3$sampled_wave3 == 1, full_data$X2 , NA)
  phase2_wave3$Z <- ifelse(phase2_wave3$sampled_wave1 == 1|
                             phase2_wave3$sampled_wave2 == 1|
                             phase2_wave3$sampled_wave3 == 1, full_data$Z , NA)
  
  phase2_wave3 <- subset(phase2_wave3, select = -c(inflBX1, inflBX2))
  
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
  fit_wave3 <- svyglm(Y ~ X1 + X2 + Z, family = quasibinomial, 
                      design = twophase_design_wave3)
  fit_wave3 <- svyglm(Y ~ X1 + X2 + Z, family = quasibinomial, 
                      design = twophase_design_wave3)
  infl_wave3 <- inf_fun_logit(fit_wave3)
  infl_wave3 <- data.frame(as.numeric(rownames(infl_wave3)), 
                           infl_wave3[,"X1"],
                           infl_wave3[,"X2"])
  names(infl_wave3) <- c("id", "inflBX1", "inflBX2")
  phase2_wave3 <- dplyr::left_join(phase2_wave3, infl_wave3, by = "id")
  
  
  # Re-calculate allocation
  wave4_allocation <- a_optimal_allocate_wave(phase2_wave3,
                                              strata = "strata",
                                              vars = c("inflBX1", "inflBX2"),
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
  
  # Sample
  phase2_wave4$X1<- ifelse(phase2_wave4$sampled_wave1 == 1 |
                             phase2_wave4$sampled_wave2 == 1|
                             phase2_wave4$sampled_wave3 == 1|
                             phase2_wave4$sampled_wave4 == 1, full_data$X1 , NA)
  phase2_wave4$X2<- ifelse(phase2_wave4$sampled_wave1 == 1 |
                             phase2_wave4$sampled_wave2 == 1|
                             phase2_wave4$sampled_wave3 == 1|
                             phase2_wave4$sampled_wave4 == 1, full_data$X2 , NA)
  phase2_wave4$Z <- ifelse(phase2_wave4$sampled_wave1 == 1|
                             phase2_wave4$sampled_wave2 == 1|
                             phase2_wave4$sampled_wave3 == 1|
                             phase2_wave4$sampled_wave4 == 1, full_data$Z , NA)
  
  phase2_wave4 <- subset(phase2_wave4, select = -c(inflBX1,
                                                   inflBX2))
  
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
  weight <- svyglm(Y ~ X1 + X2 + Z, family = quasibinomial, 
                   design = twophase_design)
  
  #########
  #### Raking 
  #########
  data <- phase2_wave4
  
  # Creating a survey object and calibrating the weights
  mydesign <- survey::twophase(
    id = list(~1, ~1), subset = ~as.logical(already_sampled),
    strata = list(NULL, ~strata), data = data
  )
  infcal <- survey::calibrate(mydesign, formula = ~ inflBX1_phase1 + inflBX2_phase1
                              + strata, phase = 2, calfun = "raking")
  
  # Fitting the outcome model: conditional treatment effect of interest.
  fit <- survey::svyglm(Y ~ X1 + X2 + Z, design = infcal, family = quasibinomial)
  
  # Check if truth is contained in interval
  confintIPWBX1 <- confint(weight)["X1",]
  confintIPWBX2 <- confint(weight)["X2",]
  confintGRBX1 <- confint(fit)["X1",]
  confintGRBX2 <- confint(fit)["X2",]
  coverIPWBX1 <- ifelse(B1 >= confintIPWBX1[1] & B1 <= confintIPWBX1[2],1,0)
  coverIPWBX2 <- ifelse(B2 >= confintIPWBX2[1] & B2 <= confintIPWBX2[2],1,0)
  coverGRBX1 <- ifelse(B1 >= confintGRBX1[1] & B1 <= confintGRBX1[2],1,0)
  coverGRBX2 <- ifelse(B2 >= confintGRBX2[1] & B2 <= confintGRBX2[2],1,0)
  
  output <- c(fit$coefficients["X1"], fit$coefficients["X2"],
              coef(weight)["X1"], coef(weight)["X2"], SE(fit)["X1"],
              SE(fit)["X2"], SE(weight)["X1"], SE(weight)["X2"], 
              coverIPWBX1, coverIPWBX2, coverGRBX1, coverGRBX2,
              oversampled_A_optimal, cor_X1_X2,
              prev_Y, prevZ, B1, B2)
  return(output)
}


######
### Strategy 5 - True A-optimal sampling
#####

Strategy5 <- function(n, simulations_df, scenario, N = 10000){
  
  # Set up betas for scenario
  B0 <- simulations_df["B0", scenario]
  B1 <- simulations_df["B1", scenario]
  B2 <- simulations_df["B2", scenario]
  B3 <- simulations_df["B3", scenario]
  
  # Set up of X's
  corX1X2 <- simulations_df["corX1X2", scenario]
  prevZ <- simulations_df["prevZ", scenario]
  
  # Measurement error
  sensY <- simulations_df["sensY", scenario]
  specY <- simulations_df["specY", scenario]
  error_varX1 <- simulations_df["error_varX1", scenario]
  error_varX2 <- simulations_df["error_varX2", scenario]
  sensZ <- simulations_df["sensZ", scenario]
  specZ <- simulations_df["specZ", scenario]
  
  
  
  #####
  ## Generate Phase 1 Data
  #####
  
  id <- 1:N
  
  # Generate correlated covariates
  sigma <- matrix(c(1, corX1X2 , 0.1, corX1X2, 1,  0.25, 0.1, 
                    0.25, 1), nrow = 3)
  covs <- mvrnorm(N, mu = c(0,0,0), Sigma = sigma)
  X1 <- covs[,1]
  X2 <- covs[,2]
  Z <- ifelse(covs[,3] >= quantile(covs[,3], 1- prevZ), 1, 0)
  cor(cbind(X1, X2, Z)) 
  
  Y_probs <- exp(B0 + B1*X1 + B2*X2 + B3*Z)/(1 +  exp(B0 + B1*X1 + B2*X2 + B3*Z))
  Y <- rbinom(N, 1, Y_probs) # Compute realized value for Y2
  
  full_data <- data.frame(id, X1, X2, Z, Y)
  
  ## Test true model
  true_model <- glm(Y ~ X1 + X2 + Z, family = "binomial", data = full_data)
  
  # Find observed correlation
  cor_X1_X2 <- cor(X1, X2) 
  
  # Prevalence
  prev_Y <- table(Y)["1"]/length(Y)
  
  # Find true coefficients
  true_model <- glm(Y ~ X1 + X2 + Z, family = "binomial", data = full_data) 
  true_BX1 <- coef(true_model)["X1"]
  true_BX2 <- coef(true_model)["X2"]
  
  
  #####
  ## Generate error-prone Phase 1 data
  ######
  
  Z_prob_obs1 <- ifelse(full_data$Z > 0, sensZ, 1 - specZ)
  full_data$Z_obs <- rbinom(N, 1, Z_prob_obs1)
  XZY_obs <- mvrnorm(N, c(0,0,0,0),
                     matrix(c(error_varX1, 0.03, 0.02, 0.02,
                              0.03, error_varX2, 0.05, 0.4,
                              0.02, 0.05, 1, 0,
                              0.02, 0.4, 0, 1), nrow = 4))
  full_data$X1_obs <- full_data$X1 + XZY_obs[,1]
  full_data$X2_obs <- full_data$X2 + XZY_obs[,2]
  
  threshold_positiveY <- qnorm((sensY + 1)/2)
  threshold_negativeY <- qnorm((specY + 1)/2)
  threshold_positiveZ <- qnorm((sensZ + 1)/2)
  threshold_negativeZ <- qnorm((specZ + 1)/2)
  full_data$Y_obs <- with(full_data, ifelse((Y == 1 & abs(XZY_obs[, 4]) < threshold_positiveY) | 
                                              (Y == 0 & abs(XZY_obs[, 4]) > threshold_negativeY), 1, 0))
  
  full_data$Z_obs <- with(full_data, ifelse((Z == 1 & abs(XZY_obs[, 3]) < threshold_positiveZ) | 
                                              (Z == 0 & abs(XZY_obs[, 3]) > threshold_negativeZ),
                                            max(Z), min(Z)))
  
  ####
  #### Divide into strata based on observed phase 1 data
  
  #### Compute naive influence functions
  fit_phase1 <-  glm(Y_obs ~ X1_obs + X2_obs + Z_obs, 
                     family = "binomial", data = full_data)
  full_data$inflBX1_phase1 <- inf_fun_logit(fit_phase1)[,"X1_obs"]
  full_data$inflBX2_phase1 <- inf_fun_logit(fit_phase1)[,"X2_obs"]
  
  #### True model
  true_modelY <-  glm(Y ~ X1 + X2 + Z, 
                      family = "binomial", data = full_data)
  
  #### Divide into strata based on observed phase 1 data X
  full_data <- full_data |>
    dplyr::mutate(Y_strat = Y_obs) |>
    split_strata(strata = "Y_strat",
                 split_var = "X1_obs",
                 type = "local quantile",
                 split_at = 0.5,
                 trunc = "X1") |>
    split_strata(strata = "new_strata",
                 split_var = "X2_obs",
                 type = "local quantile",
                 split_at = 0.5,
                 trunc = "X2")
  
  
  names(full_data)[names(full_data) == "new_strata"] <- "strata"
  
  data <- full_data
  
  ####
  ## Phase 2
  ####
  
  # Step 1: data
  data <- full_data
  
  # Step 2: Find true influence functions
  data$inflBX1_true <- inf_fun_logit(true_modelY)[,"X1"]
  data$inflBX2_true <- inf_fun_logit(true_modelY)[,"X2"]
  
  # Step 3: Determine optimum allocation forphase 2
  allocation <- a_optimum_allocation(data,
                                           strata = "strata",
                                           nsample = n,
                                           vars = c("inflBX1_true", "inflBX2_true"),
                                           weights = c(0.5, 0.5),
                                           method = "Neyman")
  
  # Sample according to wave1_allocation
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
  weight <- svyglm(Y ~ X1 + X2 + Z, family = quasibinomial, 
                   design = twophase_design)
  
  #########
  #### Raking
  #########
  
  # Creating a survey object and calibrating the weights
  mydesign <- survey::twophase(
    id = list(~1, ~1), subset = ~as.logical(sample_indicator),
    strata = list(NULL, ~strata), data = data
  )
  infcal <- survey::calibrate(mydesign, formula = ~ inflBX1_phase1 + inflBX2_phase1
                              + strata, phase = 2, calfun = "raking")
  
  # Fitting the outcome model:
  fit <- survey::svyglm(Y ~ X1 + X2 + Z, design = infcal, family = quasibinomial)
  
  # Check if truth is contained in interval
  confintIPWBX1 <- confint(weight)["X1",]
  confintIPWBX2 <- confint(weight)["X2",]
  confintGRBX1 <- confint(fit)["X1",]
  confintGRBX2 <- confint(fit)["X2",]
  coverIPWBX1 <- ifelse(B1 >= confintIPWBX1[1] & 
                          B1 <= confintIPWBX1[2],1,0)
  coverIPWBX2 <- ifelse(B2 >= confintIPWBX2[1] & 
                          B2 <= confintIPWBX2[2],1,0)
  coverGRBX1 <- ifelse(B1 >= confintGRBX1[1] &
                         B1 <= confintGRBX1[2],1,0)
  coverGRBX2 <- ifelse(B2 >= confintGRBX2[1] & 
                         B2 <= confintGRBX2[2],1,0)
  
  
  output <- c(fit$coefficients["X1"], fit$coefficients["X2"],
              coef(weight)["X1"], coef(weight)["X2"], SE(fit)["X1"],
              SE(fit)["X2"], SE(weight)["X1"], SE(weight)["X2"], 
              coverIPWBX1, coverIPWBX2, coverGRBX1, coverGRBX2, cor_X1_X2,
              prev_Y, prevZ, B1, B2)
  return(output)
 
}
