####
## Set Simulation Conditions
####

simulations_df1 <- data.frame(row.names = c("B01", "B11", "B21", "B31", "B41",
                                            "B02", "B12", "B22", "B32",
                                            "sensY1", "specY1",
                                            "sensY2", "specY2","error_varX",
                                            "sensZ1", "specZ1","error_varZ2",
                                            "data_seed"),
                              "Scenario1A" = c(-1.5, 0.4, 0, 0.3, 0, 
                                               -0.5, 0.2, 0.5, 0, 
                                               0.95, 0.99, 0.90, 0.95,
                                               0.15, 0.9, 0.95, 0.1,
                                               1),
                              "Scenario1B" = c(-1.5, 0.4, 0, 0.3, 0, 
                                               -0.5, 0.2, 0.5, 0, 
                                               0.85, 0.9, 0.8, 0.85, 
                                               0.5, 0.9, 0.95, 0.1,
                                               1),
                              "Scenario1C" = c(-3.1, 0.4, 1.0, 0.7, 1.9, 
                                               -0.8, 0.2, 1.3, 0.8, 
                                               0.95, 0.99, 0.90, 0.95,
                                               0.15, 0.9, 0.95, 0.1,
                                               2501),
                              "Scenario1D" = c(-3.1, 0.4, 1.0, 0.7, 1.9, 
                                               -0.8, 0.2, 1.3, 0.8,
                                               0.85, 0.9, 0.8, 0.85,
                                               0.5, 0.9, 0.95, 0.1,
                                               2501))

simulations_df2 <- data.frame(row.names = c("B0", "B1", "B2", "B3",
                                            "sensY", "specY", "error_varX1",
                                            "error_varX2",
                                            "prevZ",
                                            "sensZ", "specZ", "corX1X2",
                                            "data_seed"),
                              "Scenario2A" = c(-2.1, 0.3, 0.7, 0.7, 0.9, 0.95,
                                               0.15, 0.5, 0.25, 0.9, 0.95, 0.05,
                                               5001),
                              "Scenario2B" = c(-2.1, 0.3, 0.7, 0.7, 0.85, 0.9,
                                               0.4, 0.6, 0.25, 0.9, 0.95, 0.05,
                                               5001),
                              "Scenario2C" = c(-2.1, 0.3, 0.7, 0.7, 0.9, 0.95,
                                               0.15, 0.5, 0.25, 0.9, 0.95, 0.65,
                                               7501),
                              "Scenario2D" = c(-2.1, 0.3, 0.7, 0.7, 0.85, 0.9,
                                               0.4, 0.6, 0.25, 0.9, 0.95, 0.65,
                                               7501))

simulations_df3 <- data.frame(row.names = c("B01", "B11", "B21", "B31", "B41",
                                            "B02", "B12", "B22", "B32",
                                            "sensY1", "specY1",
                                            "sensY2", "specY2",
                                            "error_varX1",
                                            "error_varX2", "prevZ",
                                            "sensZ", "specZ", "corX1X2",
                                            "data_seed"),
                              "Scenario3A" = c(-1.5, 0.4, 0.6, 0.3, 0.3, 
                                               -2.1, 0.3, 0.7, 0.7, 
                                               0.95, 0.99, 0.9, 0.95,
                                               0.15, 0.5, 0.4,
                                               0.9, 0.95, 0.3, 
                                               10001),
                              "Scenario3B" =  c(-1.5, 0.4, 0.6, 0.3, 0.3, 
                                                -2.1, 0.3, 0.7, 0.7, 
                                                0.85, 0.9, 0.8, 0.85,
                                                0.4, 0.6, 0.4,
                                                0.9, 0.95, 0.3, 
                                                10001))
                             
  