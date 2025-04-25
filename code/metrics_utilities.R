#' @title Metrics and Utility Functions
#' @author Rosember Guerra Urzola
#' @description This script contains utility functions and performance metrics for evaluating sparse PCA results.
#' @date 2024-04-25
 

### Performance Metrics and Utility Functions ######

## Compare sparse PCA patterns to PCA patterns
pca_diff <- function(w, X) {
  nz <- which(w != 0)
  X <- X[, nz]
  w <- w[nz]
  w_pca <- prcomp(X)$rotation[, 1]
  norm(abs(w) - abs(w_pca), type = "2")
}

## Define penalty functions for objective calculations
Penalty_functions <- list(
  l1 = function(x, lambda) lambda * abs(x),
  l0 = function(x, lambda) lambda * (x != 0),
  scad = function(x, lambda) {
    a <- 3.7
    ifelse(abs(x) <= lambda, 
           lambda * sum(abs(x)), 
           ifelse(abs(x) <= a * lambda, 
                  (2 * lambda * a * abs(x) - x^2 - lambda^2) / (2 * (a - 1)), 
                  0.5 * (a + 1) * lambda^2))
  }
)

## Objective function for sparse PCA
space_objective <- function(X, w, alpha, penalty = 'l1') {
  if (!(penalty %in% names(Penalty_functions))) {
    stop("Invalid penalty type.")
  }
  penaltyfun <- Penalty_functions[[penalty]]
  norm(X %*% w, type = "F") - sum(penaltyfun(w, alpha))
}

## Objective function difference with PCA
obj_pca_diff <- function(X, w, alpha = 0, penalty = 'l1') {
  nz <- which(w != 0)
  w_pca <- rep(0, ncol(X))
  w_pca[nz] <- prcomp(X[, nz])$rotation[, 1]
  space_objective(X, w, alpha, penalty) - space_objective(X, w_pca, alpha, penalty)
}

## Variance calculation
variance <- function(X, w) {
  w_pca <- prcomp(X)
  (t(w) %*% t(X) %*% X %*% w) / w_pca$sdev[1]^2
}

## Adjusted variance calculation
adj_variance <- function(X, w) {
  indexw <- which(w != 0)
  X_adj <- X[, indexw]
  w_adj <- prcomp(X_adj)$rotation[, 1]
  w_pca <- prcomp(X)
  (t(w_adj) %*% t(X_adj) %*% X_adj %*% w_adj) / w_pca$sdev[1]^2
}

## Cardinality (number of non-zero elements)
cardinality <- function(w) {
  sum(w != 0)
}

## Permutation Test for Sparse PCA Results
permutation_test <- function(values_A, values_B, n_permutations = 10000, alternative = "two.sided") {
  observed_diff <- mean(values_A) - mean(values_B)
  combined <- c(values_A, values_B)
  n_A <- length(values_A)
  permuted_diffs <- replicate(n_permutations, {
    permuted <- sample(combined)
    mean(permuted[1:n_A]) - mean(permuted[(n_A + 1):length(combined)])
  })
  p_value <- switch(alternative,
                    "two.sided" = mean(abs(permuted_diffs) >= abs(observed_diff)),
                    "greater" = mean(permuted_diffs >= observed_diff),
                    "less" = mean(permuted_diffs <= observed_diff),
                    stop("Invalid alternative argument."))
  list(observed_diff = observed_diff, p_value = p_value)
}
