# Represent the dirichlet process mixture model as an S3 object
# to allow plotting and summary capabaility
dp_mm <- structure(list(), class = "dp_mm")

#' Dirichlet Process Mixture Modeling with Gibbs Sampling
#' 
#' Implements a Dirichlet Process Mixture Model with Gibbs Sampling. The base distribution for the Dirichlet 
#' Process is a Normal-Gamma distribution, which results in a predictive posterior of the t-distribution. After 
#' Gibbs Sampling, this function returns a cluster parameters for each of the clusters, as the Dirichlet Process Mixture Model
#' is a non-parametric clustering algorithm. The alpha (concentration parameter for the DP) is resampled after each iteration
#' using the soft-max trick, and each of the cluster parameters are resampled at each iteration as well. 
#' 
#' @param data The data (a data.frame) to be passed in to the DP mixture model
#' @param alpha_0 The alpha value (concentration parameter) for the DP
#' @param num_iters the number of iterations to run the gibbs sampler
#' @param burnin The number of iterations to 'burn' for the Markov Chain
#' @param useSM Use the softmax trick to configure the 'best' alpha value. Default is FALSE
#' @export
#' @keywords dirichlet process, mixture model, clustering
#' @examples
#' # Makeup some normal mixture model
#' x <- c(rnorm(50, 3, 1), rnorm(50, 6, .5), rnorm(50, 9, 1.5))
#' # hist(x, 10) 
#' output <- dp_mm(x, alpha_0=5,num_iters=100) 
#' plot(output)
#' summary(output)
#'
dp_mm <- function(data,alpha_0,num_iters=100, burnin=50,useSM=FALSE ...) {
  if (!is.numeric(data) | alpha_0 == 0) stop("Mixture model not properly initialized")
  # If our data is normal, we can put a normal prior on the mean 
  # and an independent Γ prior on the reciprocal of the variance
  
  # set our hyperparameters based on the mean and variance of the data
  mu <- mean(data)
  kappa <- 1 / var(data) 
  alpha <- 1
  beta <- 1
  
  # DP data members, initialized to maximum size
  n_points <- length(data)
  Phi <- matrix(NA, ncol = 2, nrow = n_points, dimnames = list(c(), c("mean", "var")))
  N <- rep(0, n_points) # num/cluster
  C <- rep(0, n_points) # labels for each obs x_i
  
  # Sample component parameters from the polya urn process
  # according to the current proportions and the base distribution
  j <- 2
  Phi[1, ] <- rnorm_gam(mu, kappa, alpha, beta)
  N[1] <- 1
  C[1] <- 1
  for (i in 2:n_points) {
    # make a new cluster, or sample from empirical
    if (runif(1) < (alpha_0 / (alpha_0 + n_points))) {
      Phi[j, ] <- rnorm_gam(mu, kappa, alpha, beta)
      N[j] <- 1
      C[i] <- j
      j <- j + 1
    } else {
      p <- N / sum(N)
      ci <- sum( runif(1) > cumsum( p )) + 1
      N[ci] <- N[ci] + 1
      C[i] <- ci
    }
  }
  
  # Algorithm 2 from Neal (2000), a two-step gibbs sampling method
  # First resample all of the indicators, then sample component parameters
  # separating these leads to faster mixing rates
  # Note: we don't do the full bayesian update of the base parameters for the mean and variance
  
  # make some variables
  posterior <- NULL

  # set up graph stuff
  panels <- round(sqrt(round(num_iters / 20)))
  par(mfrow=c(panels, panels))

  # Gibbs Sampler! 
  for (i in 1:num_iters) {
    # every so often, print the mixture model, and save the relative info
    if (i > burnin & i %% 10 == 0) {
      # plot
      hist(data, main=paste("Iteration:", i))
      min_max <- range(data)
      x_range <- seq(min_max[1], min_max[2], 0.1)
      y <- rep(0, length(x_range))
      for (i in which(N != 0)) {
        y <- y + N[i] * dnorm(x_range, Phi[i,"mean"], Phi[i, "var"])
      }
      lines(x_range, y, col="red"), ...)

      # save relevant info
      post <- list()
      post$Phi <- Phi[rowSums(1*is.na(Phi)) == 0,] 
      post$N <- N[N > 0]
      post$alpha <- alpha
      post$k <- length(unique(C))
      posterior <- list(posterior, post)
    }
    # Step 1: resample components (i.e the Phi values) 
    # In this step we first remove the jth data point from the model, 
    # and resample by either sampling from the base distribution conditional on x[j], 
    # or drawing from the empirical distribution. 
    for (j in 1:n_points) {
      # remove the jth data point from our model
      clust_index <- C[j]
      N[clust_index] <- N[clust_index] - 1
      if (N[clust_index] == 0) { Phi[clust_index, ] <- NA }
      
      # and now we add it back in 
      if (runif(1) < alpha_0/(alpha_0 + n_points)) {
        new_idx <- min(which(is.na(Phi)*1 == 1))
        Phi[new_idx, ] <- sample_base_conditional(data[j], mu, kappa, alpha, beta)
        N[new_idx] <- N[new_idx] + 1
        C[j] <- new_idx
      } else {                                    # draw from empirical (rich get richer)
        non_na <- which(N != 0)
        likelihood <- rep(1, length(N))
        likelihood[non_na] <- dnorm(data[j], Phi[non_na, "mean"], Phi[non_na, "var"])
                                                  # find likelihoods based on density of each Phi
        p <- N * likelihood
        p <- p / sum(p)
        new_idx <- sum(runif(1) > cumsum(p)) + 1  # sample a component
        N[new_idx] <- N[new_idx] + 1 	            # increment the cluster
        C[j] <- new_idx 				                  # record the cluster assignment
      }
    }	
    
    
    # Step 2: Update Component Parameters
    # Here we find the mean of the observations w/in each cluster
    # then we update the parameters for the normal-gamma dist. 
    # (only interested in the non-empty clusters)
    clusters <- unique(C)#non_empty <- which(apply(!is.na(Phi), 1, any) == TRUE)
    for (clstr in clusters) {
        xbar <- mean(data[C==clstr])              
        num <- N[clstr] 				                	# number of points in the cluster
        prec <- Phi[clstr, "var"]                 
                                                  # update prior, resample w/ new priors
        Phi[clstr, ] <- rnorm_gam(mu = (xbar*num*prec+mu*kappa)/(num*prec+kappa),
                                  prec = num*prec + kappa, 
                                  alpha = alpha + 1/2,
                                  beta = beta) 
    }
  
    
    # Resample the alpha for the DP prior, conditional on the current distribution of clusters
    # Use a soft-max trick to control how likely we were to sample the ML alpha
    if (useSM) {
      k <- length(N[N != 0])
      incr <- n_points / 500
      A <- seq(incr, n_points, incr)
      p <- rep(1, length(A))
      for (i in 1:length(A)) {
        a <- A[i]
        p[i] <- (k - 3/2)*log(a) - (1/(2*a)) + log(gamma(a)) - log(gamma(n_points+1))
      }
      beta = 5
      bp <- rexp(1) ^ (beta*p)
      pp <- bp / sum(bp)
      ai <- sum(runif(1) > cumsum(pp)) + 1
      alpha_0 <- A[ai]
    }
  }
  last_iter <- posterior[[length(posterior)]]
  last_iter$data <- data
  obj <- structure(last_iter, class = "dp_mm")
  return (obj)
}




#' Plot the Dirichlet Process Mixture Model output returned from the 'dp_mm' function
#' 
#' Returns a ggplot of the dirichlet process mixture model, with each of the clusters colored accordingly
#' 
#' @param obj The output of the 'dp_mm' function
#' @keywords plot, dpmm
#' @export
#' @examples
#' # Make up some data
#' data <- c(rnorm(50, 4, 1), rnorm(50, 8, .5))
#' output <- dp_mm(data, 3, num_iters=100)
#' plot(obj, col="red", lty=3)
#' 
plot.dp_mm <- function(obj, ...) {
  hist(obj$data)
  min_max <- range(obj$data)
  x_range <- seq(min_max[1], min_max[2], 0.1)
  y <- rep(0, length(x_range))
  for (i in which(N != 0)) {
        y <- y + N[i] * dnorm(x_range, Phi[i,"mean"], Phi[i, "var"])
  }
  lines(x_range, y, ...)
}

#' Get a summary for the dirichlet process mixture model
#' 
#' This function returns a summary of the dirichlet process mixture model. The summary includes the cluster parameters
#' at each iteration, the cluster labels, and more information about the mixture model. 
#' 
#' @param obj The output of the 'dp_mm' function. 
#' @keywords summary, dpmm
#' @export
#' @examples
#' # Make up some data
#' data <- c(rnorm(50, 4, 1), rnorm(50, 8, .5))
#' output <- dp_mm(data, 3, num_iters=100)
#' summary(output)
#'
summary.dp_mm <- function(obj) {
  non_na <- which(N != 0)
  summ <- data.frame(Phi=obj$Phi, N=obj$N)
  attr(summ, "num_clusters") <- obj$k
  attr(summ, "alpha") <- obj$alpha

  return (as.data.frame(Phi=obj$Phi, N=obj$N))
}


