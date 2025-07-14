####
## Scenario 2
####
library(ggplot2)
library(patchwork)

ggplot(phase2_wave4, aes(x = inflBX2_phase1, y = inflBX2, color = strata)) +
  geom_point() +
  theme_minimal()

ggplot(phase2_wave4, aes(x = strata, y = residBX2_wave3, color = strata)) +
  geom_point() +
  theme_minimal()


# Find the common y-axis range
ymin <- min(c(phase2_wave4$residBX1_wave3, phase2_wave4$residBX2_wave3), na.rm = TRUE)
ymax <- max(c(phase2_wave4$residBX1_wave3, phase2_wave4$residBX2_wave3), na.rm = TRUE)

# Now build the plots using coord_cartesian() with the same limits
p1 <- ggplot(phase2_wave4, aes(x = strata, y = residBX2_wave3, color = strata)) +
  geom_boxplot() +
  theme_minimal() +
  ggtitle("residBX2_wave3") +
  coord_cartesian(ylim = c(ymin, ymax))

p2 <- ggplot(phase2_wave4, aes(x = strata, y = residBX1_wave3, color = strata)) +
  geom_boxplot() +
  theme_minimal() +
  ggtitle("residBX1_wave3") +
  coord_cartesian(ylim = c(ymin, ymax))

# Combine the plots side by side
p1 + p2

# Try
phase2_wave4 %>%
  group_by(strata) %>%
  summarise(
    sd_residBX2_wave3 = sd(residBX2_wave3, na.rm = TRUE),
    sd_residBX1_wave3 = sd(residBX1_wave3, na.rm = TRUE),
    size = n()
  ) %>%
  dplyr::mutate(allocation_fracX2 = (size* sd_residBX2_wave3)/
                  sum(size * sd_residBX2_wave3),
                allocation_fracX1 = (size* sd_residBX1_wave3)/
                  sum(size * sd_residBX1_wave3),
                sd_residBX2_wave3_std = sd_residBX2_wave3/max(sd_residBX2_wave3[8]),
                sd_residBX1_wave3_std = sd_residBX1_wave3/max(sd_residBX1_wave3))

# Or try with true IFs in strat 5
data %>%
  group_by(strata) %>%
  summarise(
    sd_residBX2 = sd(residBX2, na.rm = TRUE),
    sd_residBX1 = sd(residBX1, na.rm = TRUE),
    size = n()
  ) %>%
  dplyr::mutate(allocation_fracX2 = (size* sd_residBX2)/
                  sum(size * sd_residBX2),
                allocation_fracX1 = (size* sd_residBX1)/
                  sum(size * sd_residBX1))

ggplot(data, aes(x = strata, y = residBX1, color = X1)) +
  geom_jitter(width = 0.2, alpha = 0.7) +
  labs(x = "Strata", y = "Residual of BX1", color = "Y") +
  theme_minimal()

ggplot(data, aes(x = residBX1, y = residBX2, color = interaction(as.factor(Y), as.factor(Y_obs)))) +
  geom_point() +
  labs(color = "Y") +
  theme_minimal()

ggplot(data, aes(x = residBX1, y = residBX2, color = as.factor(strata))) +
  geom_point() +
  labs(color = "Y") +
  theme_minimal()


#####
## Scenario 1
#####

# Find the common y-axis range
ymin <- min(c(phase2_wave4$residB11_wave3, phase2_wave4$residB12_wave3), na.rm = TRUE)
ymax <- max(c(phase2_wave4$residB11_wave3, phase2_wave4$residB12_wave3), na.rm = TRUE)

# Now build the plots using coord_cartesian() with the same limits
p3 <- ggplot(phase2_wave4, aes(x = strata, y = residB11_wave3, color = strata)) +
  geom_boxplot() +
  theme_minimal() +
  ggtitle("residBX2_wave3") +
  coord_cartesian(ylim = c(ymin, ymax))

p4 <- ggplot(phase2_wave4, aes(x = strata, y = residB12_wave3, color = strata)) +
  geom_boxplot() +
  theme_minimal() +
  ggtitle("residBX1_wave3") +
  coord_cartesian(ylim = c(ymin, ymax))

# Combine the plots side by side
p3 + p4

# Try
phase2_wave4 %>%
  group_by(strata) %>%
  summarise(
    sd_residBX2_wave3 = sd(residBX2_wave3, na.rm = TRUE),
    sd_residBX1_wave3 = sd(residBX1_wave3, na.rm = TRUE)
  )

#####
### Data example
#####


# Find the common y-axis range
ymin <- min(c(phase2_wave4$resid_ADE_wave3, phase2_wave4$resid_Death_wave3), na.rm = TRUE)
ymax <- max(c(phase2_wave4$resid_ADE_wave3, phase2_wave4$resid_Death_wave3), na.rm = TRUE)

# Now build the plots using coord_cartesian() with the same limits
p5 <- ggplot(phase2_wave4, aes(x = strata, y = resid_ADE_wave3, color = strata)) +
  geom_boxplot() +
  theme_minimal() +
  ggtitle("resid_ADE_wave3") +
  coord_cartesian(ylim = c(ymin, ymax))

p6 <- ggplot(phase2_wave4, aes(x = strata, y = resid_Death_wave3, color = strata)) +
  geom_boxplot() +
  theme_minimal() +
  ggtitle("resid_Death_wave3") +
  coord_cartesian(ylim = c(ymin, ymax))

# Combine the plots side by side
p5 + p6

# Try
phase2_wave4 %>%
  group_by(strata) %>%
  summarise(
    sd_resid_ADE_wave3 = sd(resid_ADE_wave3, na.rm = TRUE),
    sd_resid_Death_wave3 = sd(resid_Death_wave3, na.rm = TRUE)
  )

# Now build the plots using coord_cartesian() with the same limits
p7 <- ggplot(phase2_wave4, aes(x = strata, y = inflB_ADE_CDE4_wave3, color = strata)) +
  geom_boxplot() +
  theme_minimal() +
  ggtitle("infl_ADE_wave3") +
  coord_cartesian(ylim = c(ymin, ymax))

p8 <- ggplot(phase2_wave4, aes(x = strata, y = inflB_Death_CDE4_wave3, color = strata)) +
  geom_boxplot() +
  theme_minimal() +
  ggtitle("infl_Death_wave3") +
  coord_cartesian(ylim = c(ymin, ymax))

# Combine the plots side by side
p7 + p8

phase2_wave4 %>%
  group_by(strata) %>%
  summarise(
    sd_infl_ADE_wave3 = sd(inflB_ADE_CDE4_wave3, na.rm = TRUE),
    sd_infl_Death_wave3 = sd(inflB_Death_CDE4_wave3, na.rm = TRUE)
  )


######
### Scenario 3
######
#### Compute true influence functions
fitY1_true <-  glm(Y1 ~ X1 + X2 + Z + Y2, 
                     family = "binomial", data = full_data)
fitY2_true <-  glm(Y2 ~ X1 + X2 + Z, 
                     family = "binomial", data = full_data)

full_data$inflB11_true <- inf_fun_logit(fitY1_true)[,"X1"]
full_data$inflB12_true <- inf_fun_logit(fitY1_true)[,"X2"]
full_data$inflB21_true <- inf_fun_logit(fitY2_true)[,"X1"]
full_data$inflB22_true <- inf_fun_logit(fitY2_true)[,"X2"]

#### Residuals
# Regress Latest IFs (computed above) on Phase 1 IFs
resid_model_B11 <- lm(full_data$inflB11_true ~ full_data$inflB11_phase1,
                            na.action = na.exclude)
resid_model_B21 <- lm(full_data$inflB21_true ~ full_data$inflB21_phase1,
                            na.action = na.exclude)
resid_model_B12 <- lm(full_data$inflB12_true ~ full_data$inflB12_phase1,
                            na.action = na.exclude)
resid_model_B22 <- lm(full_data$inflB22_true ~ full_data$inflB22_phase1,
                            na.action = na.exclude)

# Get residuals
full_data$residB11 <- resid(resid_model_B11)
full_data$residB21 <- resid(resid_model_B21)
full_data$residB12 <- resid(resid_model_B12)
full_data$residB22 <- resid(resid_model_B22)

phase2_wave4 %>%
  group_by(strata) %>%
  summarise(
    sd_residB11 = sd(residB11, na.rm = TRUE),
    sd_residB21 = sd(residB21, na.rm = TRUE),
    sd_residB12 = sd(residB12, na.rm = TRUE),
    sd_residB22 = sd(residB22, na.rm = TRUE)
  )
