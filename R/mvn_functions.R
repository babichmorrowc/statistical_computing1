# Functions for generating samples from a multivariate normal ------------------
# Function using 2 for loops:
mvn_generator_2_for_loops <- function(dimension, mu, Sigma, n_samples) {
  if(!matrixcalc::is.positive.definite(Sigma)) {
    stop("Sigma is not positive definite!")
  }
  
  # Find the eigendecomposition of Sigma inverse
  ev_Sigma_inv <- eigen(solve(Sigma))
  U <- ev_Sigma_inv$vectors
  eigenvalues <- ev_Sigma_inv$values
  D <- diag(eigenvalues)
  if(!all.equal(U %*% D %*% solve(U), solve(Sigma))) {
    stop("Something has gone wrong in your eigendecomposition!")
  }
  
  x_samples <- matrix(0, nrow = dimension, ncol = n_samples)
  for (i in 1:n_samples) {
    y_i <- matrix(0, nrow = dimension, ncol = 1)
    for (j in 1:dimension) {
      y_i[j] <- rnorm(1, mean = 0, sd = sqrt(1/eigenvalues[j]))
    }
    x_i <- U %*% y_i + mu
    x_samples[,i] <- x_i
  }
  
  return(x_samples)
  
}

# Vectorized function
mvn_generator_vector <- function(dimension, mu, Sigma, n_samples) {
  if(!matrixcalc::is.positive.definite(Sigma)) {
    stop("Sigma is not positive definite!")
  }
  
  # Find the eigendecomposition of Sigma inverse
  ev_Sigma_inv <- eigen(solve(Sigma))
  U <- ev_Sigma_inv$vectors
  eigenvalues <- ev_Sigma_inv$values
  D <- diag(eigenvalues)
  if(!all.equal(U %*% D %*% solve(U), solve(Sigma))) {
    stop("Something has gone wrong in your eigendecomposition!")
  }
  
  y_samples <- matrix(rnorm(n = dimension*n_samples,
                            mean = rep(0, dimension*n_samples),
                            sd = rep(sqrt(1/eigenvalues), n_samples)),
                      nrow = dimension,
                      ncol = n_samples)
  x_samples <- U %*% y_samples + matrix(rep(mu, n_samples),
                                        nrow = dimension,
                                        ncol = n_samples)
  
  return(x_samples)
  
}