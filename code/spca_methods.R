#' @title Sparse PCA Methods
#' @description Functions for performing Sparse PCA using methods.
#' @date 2025-04-25
#' @author Rosember Guerra Urzola


### Penalty Functions ###

Penalty_functions_sol <- list(
  l1 = function(x, lambda) {
    sign(x) * pmax(abs(x) - lambda, 0)
  },
  l0 = function(x, lambda) {
    x * (abs(x) > lambda)
  },
  scad = function(x, lambda) {
    a <- 3.7
    ifelse(abs(x) <= 2 * lambda, 
           sign(x) * pmax(abs(x) - lambda, 0), 
           ifelse(abs(x) <= a * lambda, 
                  sign(x) * ((a - 1) * abs(x) - a * lambda) / (a - 2), 
                  x * (abs(x) > a * lambda)))
  }
)

### Alternating Sparse PCA Function#######

alt_spca <- function(X, w0 = NULL, alpha, penalty = 'l1', maxiter = 1000, tol = 1e-6) {
  # Perform sparse PCA using alternating optimization
  # Parameters:
  # - X: Data matrix
  # - w0: Initial weights (optional)
  # - alpha: Regularization parameter
  # - penalty: Type of penalty ('l1', 'l0', 'scad')
  # - maxiter: Maximum number of iterations
  # - tol: Convergence tolerance

  # Ensure valid penalty type
  if (!(penalty %in% names(Penalty_functions_sol))) {
    stop("Invalid penalty type. Available options are: 'l1', 'l0', 'scad'")
  }
  if (alpha < 0) stop("alpha must be non-negative.")
  if (alpha == 0) return(list(w = prcomp(X)$rotation[, 1], iter = 0, time = 0))
  
  penaltyfun <- Penalty_functions_sol[[penalty]]
  if (is.null(w0)) w0 <- rnorm(ncol(X))
  
  iter <- 0
  runningtime <- system.time({
    while (iter < maxiter) {
      # z-step: Project and normalize
      z <- X %*% w0
      z <- z / norm(z, type = "2")
      
      # w-step: Apply penalty and normalize
      w <- penaltyfun(t(X) %*% z, alpha)
      if (all(w == 0)) break
      w <- w / norm(w, type = "2")
      
      if (norm(w - w0, type = "2") < tol) break
      
      w0 <- w
      iter <- iter + 1
    }
  })
  return(list(w = w, iter = iter, time = runningtime[[3]]))
}

### Multivariate Sparse PCA with Deflation   #########

alt_spca_multi <- function(X, W0 = NULL, K, alpha, penalty = 'l1') {
  # Perform multivariate sparse PCA with deflation
  if (K > ncol(X)) stop("K cannot exceed the number of columns in X.")
  if (is.null(W0)) W0 <- matrix(rnorm(ncol(X) * K), ncol(X), K)
  if (length(alpha) == 1) alpha <- rep(alpha, K)
  else if (length(alpha) != K) stop("alpha must be a scalar or a vector of length K.")
  
  W <- matrix(0, ncol(X), K)
  for (k in 1:K) {
    W[, k] <- alt_spca(X, w0 = W0[, k], alpha = alpha[k], penalty = penalty)$w
    X <- X - X %*% W[, k] %*% t(W[, k])
  }
  return(W)
}

