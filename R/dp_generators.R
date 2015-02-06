#' Sample for the Dirichlet Process
#' 
#' Returns a sample (non-unique) from the dirichlet process with a specified
#' concentration parameter (alpha) and a base distribution
#' 
#' @param n The number of samples to sample
#' @param alpha The concentration parameter for the dirichlet process
#' @param base The base distribution (commonly refered to as H) for the dirichlet process
#' @keywords dirichlet process
#' @export
#' @examples
#' base <- function(n) rpois(n, 20)
#' samp <- dp_sample(10, 5, base)
#'
dp_sample <- function(n, alpha, base) {
  n <- 100
  beta <- rbeta(n, 1, alpha) 
  betas <- c(1, cumprod(1 - beta))[1:n]
  pis <- betas * beta
  y <- base(n) 
  theta <- sample(y, prob = pis, replace = TRUE) 
  return (sample(theta, n, replace=TRUE))
}

#' A R function that returns a generator for the dirichlet process
#' 
#' Returns a generator for a dirichlet process (DP) of a specified alpha and base distribution. 
#' 
#' @param alpha The concentration parameter for the dirichlet process
#' @param base The base distribution (commonly refered to as H) for the dirichlet process
#' @keywords dirichlet process
#' @export
#' @examples
#' base <- function(n) rnorm(n, mean=5, sd=2)
#' dp <- dp_generate(2, base)
#' dp(5)
#
dp_generator <- function( alpha, base ) {
  n <- 100
  beta <- rbeta(n, 1, alpha) 
  betas <- c(1, cumprod(1 - beta))[1:n]
  pis <- betas * beta
  y <- base(n) 
  theta <- sample(y, prob = pis, replace = TRUE) 
  
  return (function(num) {
    return (sample(theta, num, replace = TRUE))
  })
}


#' Sample for the Pitman-Yor Process
#' 
#' Returns a sample (non-unique) from the pitman-yor process with a specified
#' concentration parameter (alpha), tuning paramter (d) and a base distribution
#' 
#' @param n The number of samples to sample
#' @param alpha The concentration parameter for the PY process
#' @param d The tuning parameter for the PY process
#' @param base The base distribution (commonly refered to as H) for the PY process
#' @keywords pitman-yor process
#' @export
#' @examples
#' base <- function(n) rpois(n, 20)
#' samp <- py_sample(10, 5, 2, base)
#'
py_sample <- function(n, alpha, d, base) {
  n <- 100
  beta <- rbeta(n, 1 - d, alpha + n*d )
  betas <- c(1, cumprod(1 - beta))[1:n]
  pis <- betas * beta
  y <- base(n)
  theta <- sample(y, prob = pis, replace = TRUE)
  return (sample(theta, num, replace = TRUE))
}

#' A generator for the Pitman-Yor Process in R
#' 
#' Returns a generator for a Pitman-Yor process, given a specified $d$, $alpha$, and a base distribution
#' 
#' @param alpha The concentration parameter 
#' @param d The tuning paramter for the Pitman-Yor. If d=1 we have the Dirichlet Process
#' @param base The base distribution for the Pitman-Yor process
#' @keywords pitman-yor, process
#' @export
#' @examples
#' base <- function(n) rnorm(n, mean=5, sd=2)
#' py <- py_generate(3, 2, base)
#' py(5)
#' 
py_generator <- function( alpha, d, base ) {
  n <- 100
  beta <- rbeta(n, 1 - d, alpha + n*d )
  betas <- c(1, cumprod(1 - beta))[1:n]
  pis <- betas * beta
  y <- base(n)
  theta <- sample(y, prob = pis, replace = TRUE)
  
  return (function(num) {
    return (sample(theta, num, replace = TRUE))
  })
}
