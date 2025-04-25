# Data set generations for numerical experiments
# RI Guerra Urzola
# 30/10/2024

normal_data = function(n, p, mu, sigma){
  # n: number of samples
  # p: number of features
  # mu: mean of the normal distribution
  # sigma: standard deviation of the normal distribution
  X = matrix(rnorm(n*p, mu, sigma), n, p)
  return(X)
}

generate_ordered_eigen_data <- function(n, p, mu = 0, sigma = 1, eigenvalues = NULL) {
  # n: number of samples
  # p: number of features
  # mu: mean of the normal distribution
  # sigma: standard deviation of the normal distribution
  # eigenvalues: optional vector of decreasing eigenvalues
  
  if (is.null(eigenvalues)) {
    # Generate a decreasing sequence of positive eigenvalues
    eigenvalues <- sort(runif(p, 1, 10), decreasing = TRUE)
  }
  
  # Generate a random orthonormal matrix (eigenvectors)
  Q <- qr.Q(qr(matrix(rnorm(p^2), p, p)))
  
  # Construct the covariance matrix with controlled eigenvalues
  Sigma <- Q %*% diag(eigenvalues) %*% t(Q)
  
  # Generate multivariate normal data
  library(MASS)
  X <- mvrnorm(n, rep(mu, p), Sigma)
  
  return(X)
}
