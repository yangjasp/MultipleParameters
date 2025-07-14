######
### Functions for Multiple Outcomes Simulations
######

# 1. Influence function logistic regression (with help from Tong Chen)

# @param fit: A logistic model fit

inf_fun_logit <- function(fit){
  dm <- model.matrix(fit)
  Ihat <- (t(dm) %*% (dm * fit$fitted.values * (1 - fit$fitted.values))) /nrow(dm)
  infl <- (dm * resid(fit, type = "response")) %*% solve(Ihat)
}


# 2. A-optimality
###
### A-optimality function - same output as optimum_allocation. Given 2 inputs
a_optimum_allocation <- function(data, strata = "strata", 
                                 nsample, vars, weights,
                                 method){
  P <- length(vars)
  N_k <- table(data[,strata]) # before removing NAs
  data <- data[!is.na(data[,vars[1]]),]
  var_list <- list()
  a_var_list <- list()
  for (i in seq_along(vars)){
    variances <- tapply(data[,vars[i]], data[,strata], var)
    var_list[[i]] <- variances
    a_var_list[[i]] <- weights[i]*variances
  }
  sums <- rowSums(t(dplyr::bind_rows(a_var_list))) # Gives sum per stratum
  df <- data.frame("strata" = names(N_k), "N_k" = as.vector(N_k),
                   "sqrt_a_var_sum" = sqrt(sums))
  
  # Now run Neyman on this df
  output <- optimum_allocation(df, strata = "strata", sd_h = "sqrt_a_var_sum",
                               N_h = "N_k", nsample = nsample, 
                               method = method)
  
  return(output)
}

# 3. allocate_wave() with Neyman solution rather than Wright for speed (adapted
# from 'optimal' function)
allocate_wave_Neyman <- function(data,
                                 strata,
                                 y, already_sampled,
                                 nsample,
                                 method = c("iterative","simple"),
                                 detailed = FALSE,
                                 allocation_method = "Neyman") {
  key <- stratum_size <- wave1_size <- npop <- difference <-
    nsample_prior <- n_to_sample <- nsample_actual <-
    nsample_optimal <- sd <- NULL # bind global vars as necessary
  if (is.matrix(data)) {
    data <- data.frame(data)
  }
  if (is.data.frame(data) == FALSE) {
    stop("Input data must be a dataframe or matrix with named columns.")
  }
  if (all(strata %in% names(data)) == FALSE) {
    stop("'strata' must be a character string or vector of
    strings matching column names of data.")
  }
  if (y %in% names(data) == FALSE) {
    stop("'y' must be a character string matching a column name of data.")
  }
  if (already_sampled %in% names(data) == FALSE) {
    stop("'already_sampled' must be a character string matching a column name of
           data.")
  }
  if (inherits(detailed, "logical") == FALSE) {
    stop("'detailed' must be a logical value.")
  }
  if (length(table(data[, already_sampled])) != 2) {
    stop("'already_sampled' must be a character string matching a column in
         'data' that has a binary indicator for whether each unit
         was already sampled. If no units have been sampled yet,
         use 'optimum_allocation'.")
  }
  if (("Y" %in% data[, already_sampled] == FALSE & 1 %in%
       data[, already_sampled] == FALSE) | anyNA(data[, already_sampled])) {
    stop("'already_sampled' column must contain '1' (numeric) or 'Y'
         (character) as indicators that a unit was sampled in a
         previous wave and cannot contain NAs. If no units have
         been sample, use 'optimum_allocation.")
  }
  if (nsample + sum(data[, already_sampled] == "Y") +
      sum(data[, already_sampled] == 1) > length(data[, y])) {
    stop("Total sample size across waves, taken as nsampled in
         already_sampled + nsample, is larger than the population size.")
  }
  method <- match.arg(method)
  # Find the total sample size and optimally allocate that
  nsampled <- sum(data[, already_sampled] == "Y" | data[, already_sampled] == 1)
  output1 <- optimall::optimum_allocation(
    data = data,
    strata = strata,
    y = y,
    nsample = nsample + nsampled,
    allow.na = TRUE,
    method = allocation_method
  )
  # Optimal for total sample size
  
  # Create groups from strata argument and determine the prior
  # sample size for each
  y <- enquo(y)
  strata <- enquo(strata)
  key_q <- enquo(already_sampled)
  wave1_df <- data %>%
    dplyr::select(!!strata, !!y, !!key_q)
  group <- interaction(dplyr::select(wave1_df, !!strata))
  wave1_df <- cbind(group, wave1_df)
  wave1_df <- dplyr::select(wave1_df, 1, !!y, !!key_q)
  # Only columns of interest
  names(wave1_df) <- c("group", "y", "key")
  wave1_summary <- wave1_df %>%
    dplyr::group_by(group) %>%
    dplyr::summarize(wave1_size = sum(key == 1 | key == "Y"))
  
  names(output1)[1] <- "group"
  comp_df <- dplyr::inner_join(output1, wave1_summary, by = "group")
  comp_df <- dplyr::mutate(comp_df,
                           difference = stratum_size - wave1_size,
                           n_avail = npop - wave1_size
  )
  
  # For the simple case in which no strata have been oversampled
  if (all(comp_df$difference >= 0)) {
    comp_df <- comp_df %>%
      dplyr::rename(
        nsample_optimal = stratum_size,
        nsample_prior = wave1_size,
        n_to_sample = difference
      ) %>%
      dplyr::mutate(nsample_actual = nsample_prior + n_to_sample)
    if (detailed == FALSE) {
      comp_df <- comp_df %>%
        dplyr::select(
          "strata" = group, npop, nsample_actual,
          nsample_prior, n_to_sample
        )
    } else if (detailed == TRUE) {
      comp_df <- comp_df %>%
        dplyr::select(
          "strata" = group, npop, nsample_optimal,
          nsample_actual, nsample_prior,
          n_to_sample, sd
        )
    }
    return(comp_df)
  }
  
  # If some Strata have been oversampled. Basic, non-iterative method.
  if (any(comp_df$difference < 0) & method == "simple") {
    temp <- dplyr::filter(comp_df, difference <= 0)
    n_oversampled <- -sum(temp$difference)
    closed_groups <- (temp$group)
    nsampled_in_closed_groups <- sum(temp$wave1_size)
    
    open_groups <- dplyr::filter(comp_df, difference > 0)$group
    open_df <- wave1_df %>%
      dplyr::filter(group %in% open_groups)
    open_output <- optimall::optimum_allocation(
      data = open_df,
      strata = "group",
      y = "y",
      nsample = nsample + nsampled - nsampled_in_closed_groups,
      allow.na = TRUE,
      method = allocation_method
    )
    names(open_output)[1] <- "group"
    open_output <- dplyr::inner_join(open_output, wave1_summary, by = "group")
    open_output <- dplyr::mutate(
      open_output,
      difference = stratum_size - wave1_size,
      n_avail = npop - wave1_size
    )
    
    open_output <- open_output %>%
      dplyr::rename(
        nsample_optimal = stratum_size,
        nsample_prior = wave1_size
      ) %>%
      dplyr::mutate(
        n_to_sample = difference,
        nsample_actual = nsample_prior + n_to_sample
      ) %>%
      dplyr::select(
        "strata" = group,
        npop,
        nsample_actual,
        nsample_prior,
        n_to_sample
      )
    
    closed_output <- temp %>%
      dplyr::rename(
        nsample_optimal = stratum_size,
        nsample_prior = wave1_size
      ) %>%
      dplyr::mutate(
        n_to_sample = 0,
        nsample_actual = nsample_prior
      ) %>%
      dplyr::select(
        "strata" = group,
        npop,
        nsample_actual,
        nsample_prior,
        n_to_sample
      )
    
    output_df <- rbind(closed_output, open_output)
    if (detailed == TRUE) {
      output_df <- dplyr::inner_join(
        output_df,
        dplyr::select(output1,
                      "nsample_optimal" = stratum_size,
                      sd,
                      "strata" = group
        ),
        by = "strata"
      )
      output_df <- dplyr::select(
        output_df, strata, npop, nsample_optimal,
        nsample_actual, nsample_prior, n_to_sample, sd
      )
    }
    output_df <- dplyr::arrange(output_df, strata)
    if (any(output_df$n_to_sample < 0)) {
      warning("The simple method yielded strata with negative
              n_to_sample values due to many groups being
              oversampled in prior waves. Switching to
              method = 'iterative'.")
      did_simple_work <- FALSE
      method <- "iterative"
      rm(output_df, closed_output, open_output, closed_groups, open_groups)
    } else {
      return(output_df)
    }
  }
  # Now, iterative method
  
  if (any(comp_df$difference < 0) & method == "iterative") {
    closed_groups_df <- data.frame()
    
    while (any(comp_df$difference < 0)) {
      # Find most oversampled group. Add that group to the closed strata.
      closed_groups_df <- rbind(
        closed_groups_df,
        dplyr::filter(
          comp_df,
          difference ==
            min(difference)
        )
      )
      nsampled_in_closed_groups <- sum(closed_groups_df$wave1_size)
      closed_groups <- (closed_groups_df$group)
      
      # Filter comp_df, remove the smallest group
      open_groups_names <- dplyr::filter(
        comp_df,
        difference !=
          min(difference)
      )$group
      open_df <- wave1_df %>%
        dplyr::filter(group %in% open_groups_names)
      
      # Run optimal allocation on this filtered df of open groups
      outputn <- optimall::optimum_allocation(
        data = open_df, strata = "group", y = "y",
        nsample = nsample + nsampled - nsampled_in_closed_groups,
        allow.na = TRUE, method = allocation_method
      )
      
      # Re-join with (cleaned) input data to  get new differences
      names(outputn)[1] <- "group"
      comp_df <- dplyr::inner_join(outputn, wave1_summary, by = "group")
      comp_df <- dplyr::mutate(comp_df,
                               difference = stratum_size - wave1_size,
                               n_avail = npop - wave1_size
      )
    }
    open_output <- comp_df %>%
      dplyr::rename(
        nsample_optimal = stratum_size,
        nsample_prior = wave1_size
      ) %>%
      dplyr::mutate(
        n_to_sample = difference,
        nsample_actual = nsample_prior + n_to_sample
      ) %>%
      dplyr::select(
        "strata" = group, npop, nsample_actual, nsample_prior,
        n_to_sample
      )
    
    closed_output <- closed_groups_df %>%
      dplyr::rename(
        nsample_optimal = stratum_size,
        nsample_prior = wave1_size
      ) %>%
      dplyr::mutate(
        n_to_sample = 0,
        nsample_actual = nsample_prior
      ) %>%
      dplyr::select(
        "strata" = group, npop, nsample_actual, nsample_prior,
        n_to_sample
      )
    output_df <- rbind(closed_output, open_output)
    if (detailed == TRUE) {
      output_df <- dplyr::inner_join(
        output_df,
        dplyr::select(output1,
                      "nsample_optimal" = stratum_size,
                      sd,
                      "strata" = group
        ),
        by = "strata"
      )
      output_df <- dplyr::select(
        output_df, strata, npop, nsample_optimal,
        nsample_actual, nsample_prior, n_to_sample, sd
      )
    }
    output_df <- dplyr::arrange(output_df, strata)
    return(output_df)
  }
}

# 4. allocate_wave() with a-optimality (adapted from 'optimall' function)
a_optimal_allocate_wave <- function(data,
                                    strata, vars, weights,
                                    already_sampled,
                                    nsample, 
                                    method = c("iterative","simple"),
                                    detailed = FALSE,
                                    allocation_method = "Neyman") {
  key <- stratum_size <- wave1_size <- npop <- difference <-
    nsample_prior <- n_to_sample <- nsample_actual <-
    nsample_optimal <- sd <- NULL # bind global vars as necessary
  if (is.matrix(data)) {
    data <- data.frame(data)
  }
  if (is.data.frame(data) == FALSE) {
    stop("Input data must be a dataframe or matrix with named columns.")
  }
  if (all(strata %in% names(data)) == FALSE) {
    stop("'strata' must be a character string or vector of
    strings matching column names of data.")
  }
  if (any(vars %in% names(data) == FALSE)) {
    stop("'y' must be a character string matching a column name of data.")
  }
  if (already_sampled %in% names(data) == FALSE) {
    stop("'already_sampled' must be a character string matching a column name of
           data.")
  }
  if (inherits(detailed, "logical") == FALSE) {
    stop("'detailed' must be a logical value.")
  }
  if (length(table(data[, already_sampled])) != 2) {
    stop("'already_sampled' must be a character string matching a column in
         'data' that has a binary indicator for whether each unit
         was already sampled. If no units have been sampled yet,
         use 'optimum_allocation'.")
  }
  if (("Y" %in% data[, already_sampled] == FALSE & 1 %in%
       data[, already_sampled] == FALSE) | anyNA(data[, already_sampled])) {
    stop("'already_sampled' column must contain '1' (numeric) or 'Y'
         (character) as indicators that a unit was sampled in a
         previous wave and cannot contain NAs. If no units have
         been sample, use 'optimum_allocation.")
  }
  if (nsample + sum(data[, already_sampled] == "Y") +
      sum(data[, already_sampled] == 1) > length(data[, vars[1]])) {
    stop("Total sample size across waves, taken as nsampled in
         already_sampled + nsample, is larger than the population size.")
  }
  method <- match.arg(method)
  # Find the total sample size and optimally allocate that
  nsampled <- sum(data[, already_sampled] == "Y" | data[, already_sampled] == 1)
  output1 <- a_optimum_allocation(
    data = data,
    strata = strata,
    vars = vars,
    weights = weights,
    nsample = nsample + nsampled,
    method = allocation_method
  )
  # Optimal for total sample size
  
  # Create groups from strata argument and determine the prior
  # sample size for each
  raw_vars <- vars
  vars <- enquo(vars)
  strata <- enquo(strata)
  key_q <- enquo(already_sampled)
  wave1_df <- data %>%
    dplyr::select(!!strata, !!vars, !!key_q)
  group <- interaction(dplyr::select(wave1_df, !!strata))
  wave1_df <- cbind(group, wave1_df)
  wave1_df <- dplyr::select(wave1_df, 1, !!vars, !!key_q)
  # Only columns of interest
  names(wave1_df) <- c("group", raw_vars, "key")
  wave1_summary <- wave1_df %>%
    dplyr::group_by(group) %>%
    dplyr::summarize(wave1_size = sum(key == 1 | key == "Y"))
  
  names(output1)[1] <- "group"
  comp_df <- dplyr::inner_join(output1, wave1_summary, by = "group")
  comp_df <- dplyr::mutate(comp_df,
                           difference = stratum_size - wave1_size,
                           n_avail = npop - wave1_size
  )
  
  # For the simple case in which no strata have been oversampled
  if (all(comp_df$difference >= 0)) {
    comp_df <- comp_df %>%
      dplyr::rename(
        nsample_optimal = stratum_size,
        nsample_prior = wave1_size,
        n_to_sample = difference
      ) %>%
      dplyr::mutate(nsample_actual = nsample_prior + n_to_sample)
    if (detailed == FALSE) {
      comp_df <- comp_df %>%
        dplyr::select(
          "strata" = group, npop, nsample_actual,
          nsample_prior, n_to_sample
        )
    } else if (detailed == TRUE) {
      comp_df <- comp_df %>%
        dplyr::select(
          "strata" = group, npop, nsample_optimal,
          nsample_actual, nsample_prior,
          n_to_sample, sd
        )
    }
    return(comp_df)
  }
  
  # If some Strata have been oversampled. Basic, non-iterative method.
  if (any(comp_df$difference < 0) & method == "simple") {
    temp <- dplyr::filter(comp_df, difference <= 0)
    n_oversampled <- -sum(temp$difference)
    closed_groups <- (temp$group)
    nsampled_in_closed_groups <- sum(temp$wave1_size)
    
    open_groups <- dplyr::filter(comp_df, difference > 0)$group
    open_df <- wave1_df %>%
      dplyr::filter(group %in% open_groups)
    
    open_df$group <- factor(open_df$group, levels = open_groups) # remove closed
    # groups as factors
    
    open_output <- a_optimum_allocation(
      data = open_df,
      strata = "group",
      vars = raw_vars,
      weights = weights,
      nsample = nsample + nsampled - nsampled_in_closed_groups,
      method = allocation_method
    )
    names(open_output)[1] <- "group"
    open_output <- dplyr::inner_join(open_output, wave1_summary, by = "group")
    open_output <- dplyr::mutate(
      open_output,
      difference = stratum_size - wave1_size,
      n_avail = npop - wave1_size
    )
    
    open_output <- open_output %>%
      dplyr::rename(
        nsample_optimal = stratum_size,
        nsample_prior = wave1_size
      ) %>%
      dplyr::mutate(
        n_to_sample = difference,
        nsample_actual = nsample_prior + n_to_sample
      ) %>%
      dplyr::select(
        "strata" = group,
        npop,
        nsample_actual,
        nsample_prior,
        n_to_sample
      )
    
    closed_output <- temp %>%
      dplyr::rename(
        nsample_optimal = stratum_size,
        nsample_prior = wave1_size
      ) %>%
      dplyr::mutate(
        n_to_sample = 0,
        nsample_actual = nsample_prior
      ) %>%
      dplyr::select(
        "strata" = group,
        npop,
        nsample_actual,
        nsample_prior,
        n_to_sample
      )
    
    output_df <- rbind(closed_output, open_output)
    if (detailed == TRUE) {
      output_df <- dplyr::inner_join(
        output_df,
        dplyr::select(output1,
                      "nsample_optimal" = stratum_size,
                      sd,
                      "strata" = group
        ),
        by = "strata"
      )
      output_df <- dplyr::select(
        output_df, strata, npop, nsample_optimal,
        nsample_actual, nsample_prior, n_to_sample, sd
      )
    }
    output_df <- dplyr::arrange(output_df, strata)
    if (any(output_df$n_to_sample < 0)) {
      warning("The simple method yielded strata with negative
              n_to_sample values due to many groups being
              oversampled in prior waves. Switching to
              method = 'iterative'.")
      did_simple_work <- FALSE
      method <- "iterative"
      rm(output_df, closed_output, open_output, closed_groups, open_groups)
    } else {
      return(output_df)
    }
  }
  # Now, iterative method
  
  if (any(comp_df$difference < 0) & method == "iterative") {
    closed_groups_df <- data.frame()
    
    while (any(comp_df$difference < 0)) {
      # Find most oversampled group. Add that group to the closed strata.
      closed_groups_df <- rbind(
        closed_groups_df,
        dplyr::filter(
          comp_df,
          difference ==
            min(difference)
        )
      )
      nsampled_in_closed_groups <- sum(closed_groups_df$wave1_size)
      closed_groups <- (closed_groups_df$group)
      
      # Filter comp_df, remove the smallest group
      open_groups_names <- dplyr::filter(
        comp_df,
        difference !=
          min(difference)
      )$group
      open_df <- wave1_df %>%
        dplyr::filter(group %in% open_groups_names)
      
      
      open_df$group <- factor(open_df$group, levels = open_groups_names) # remove closed
      # groups as factors
      
      # Run optimal allocation on this filtered df of open groups
      outputn <-  a_optimum_allocation(
        data = open_df,
        strata = "group",
        vars = raw_vars,
        weights = weights,
        nsample = nsample + nsampled - nsampled_in_closed_groups,
        method = allocation_method
      )
      
      # Re-join with (cleaned) input data to  get new differences
      names(outputn)[1] <- "group"
      comp_df <- dplyr::inner_join(outputn, wave1_summary, by = "group")
      comp_df <- dplyr::mutate(comp_df,
                               difference = stratum_size - wave1_size,
                               n_avail = npop - wave1_size
      )
    }
    open_output <- comp_df %>%
      dplyr::rename(
        nsample_optimal = stratum_size,
        nsample_prior = wave1_size
      ) %>%
      dplyr::mutate(
        n_to_sample = difference,
        nsample_actual = nsample_prior + n_to_sample
      ) %>%
      dplyr::select(
        "strata" = group, npop, nsample_actual, nsample_prior,
        n_to_sample
      )
    
    closed_output <- closed_groups_df %>%
      dplyr::rename(
        nsample_optimal = stratum_size,
        nsample_prior = wave1_size
      ) %>%
      dplyr::mutate(
        n_to_sample = 0,
        nsample_actual = nsample_prior
      ) %>%
      dplyr::select(
        "strata" = group, npop, nsample_actual, nsample_prior,
        n_to_sample
      )
    output_df <- rbind(closed_output, open_output)
    if (detailed == TRUE) {
      output_df <- dplyr::inner_join(
        output_df,
        dplyr::select(output1,
                      "nsample_optimal" = stratum_size,
                      sd,
                      "strata" = group
        ),
        by = "strata"
      )
      output_df <- dplyr::select(
        output_df, strata, npop, nsample_optimal,
        nsample_actual, nsample_prior, n_to_sample, sd
      )
    }
    output_df <- dplyr::arrange(output_df, strata)
    return(output_df)
  }
}
