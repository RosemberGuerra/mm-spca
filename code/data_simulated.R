############################################
# Data Generation for Simulation Exercise  #
# Author: RI Guerra Urzola                 #
# Date: 11-13-2024                         #
############################################

# Clear the workspace to avoid conflicts with previous variables or functions
rm(list = ls())

# Load external functions for data generation
source('datageneration.R')

# Set the seed for reproducibility of random number generation
set.seed(123)

############################################
# Simulation Parameters                    #
############################################
# S: Number of simulations to run
# n: Number of observations in each dataset
# p: Number of variables to test (small, medium, large settings)
S <- 1000
n <- 100
p <- c(20, 100, 200)

# Create a grid of parameter combinations for simulations
param_grid <- expand.grid(n = n, p = p, s = 1:S)

############################################
# Data Distribution Parameters             #
############################################
# mu: Mean of the normal distribution
# sigma: Standard deviation of the normal distribution
mu <- 0
sigma <- 1

############################################
# Progress Bar Setup                       #
############################################
# Initialize a progress bar to track simulation progress
pb <- txtProgressBar(min = 0, max = nrow(param_grid), style = 3)

############################################
# Simulation Loop                          #
############################################
# For each combination of parameters, generate and save a dataset
for (i in 1:nrow(param_grid)) {
  # Update the progress bar
  setTxtProgressBar(pb, i)
  
  # Extract current parameter settings
  n <- param_grid$n[i]
  p <- param_grid$p[i]
  s <- param_grid$s[i]
  
  # Generate a data matrix with normal distribution
  X <- generate_ordered_eigen_data(n, p, mu, sigma)
  
  # Scale the data matrix to have zero mean and unit variance
  X <- scale(X)
  
  # Save the generated data matrix to a file
  write.table(X, 
              file = paste0('../data/simdata_', n, '_', p, '_', s, '.txt'), 
              row.names = FALSE, 
              col.names = FALSE)
}

# Close the progress bar after all simulations are completed
close(pb)
