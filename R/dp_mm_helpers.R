# sample and mean and precision for a component using the normal-gamma distribution
# given the specified mean, inverse variance, alpha (shape) and beta (scale)
rnorm_gam <- function(mu, prec, alpha, beta) {
  p <- rgamma(1, alpha, beta)
  m <- rnorm(1, mu, 1 / (p * prec))
  return (c(mean = m, prec = p))
}


# draws one sample from the base distribution, conditional on the observation x
# obj$mu1 <- 0.1
# obj$mu2, var1, var2
sample_base_conditional <- function(x, mu, kappa, alpha, beta) {
  mu_new <- (mu + x) / 2
  kap_new <- kappa
  alpha_new <- alpha + 1/2
  beta_new <- beta
  return (rnorm_gam(mu_new, kap_new, alpha_new, beta_new))
}
