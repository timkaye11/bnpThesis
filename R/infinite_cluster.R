# A S3 object to represent the infinite clustering mixture model. Allows for plotting and summary functions
infinite_mm <- structure(list(), class = "infinite_mm")

#' Infinite dimensional clustering for R
#' 
#' Performs infinite dimensional clustering, as specified by Algorithm 1 in 
#'  "___" by Jordan, Blei, et al. Unlike the typical k-means
#' clustering algorithm when the value of *k* is specified, this algorithm clusters the observations according to a dispersion
#' parameter (alpha). 
#'
#'@param data The data to be clustered by the algorithm 
#'@param alpha The dispersion parameter for the clustering algorithm
#'@param distance The distance metric used in the clustering algorithm. The default is set to Euclidean distance, 
#'  but other options include: manhattan, mahalanobis, and cosine distance
#'@export
#'@keywords infinite, clustering, dp
#'@examples
#'output <- infinite_cluster(data, alpha=2, distance="Euclidian")
infinite_cluster <- function(data, alpha, distance="Euclidian")  {
  # define the different distance functions
  euc_distance <- function(x, y) { sqrt(sum(x-y)^2) }
  man_distance <- function(x, y) { sqrt(sum(x-y)) }
  cos_distance <- function(x, y) { cos(x, y) }
  
  # establish the distance function 
  dist_func <- NULL
  switch(distance, 
         Euclidian = (dist_func <- euc_distance),
         Manhattan = (dist_func <- man_distance),
         Cosine = (dist_func <- cos_distance)
  )
  
  # initialize variables to save later
  n 	    <- nrow(data)
  labels  <- rep(1, n) 
  k 	 	<- 1
  iters   <- 0
  centers <- matrix( apply(data[, ], 2, mean) , ncol = ncol(data) )
  colnames(centers) <- colnames(data)
  
  # loop until convergence
  while (iters < maxiters) {
    iters <- iters + 1
    # assign points to labels 
    for (i in 1:n) {
      distances <- apply( centers, 1, distance, data[i, ] )  
      # if no 'centers' are the the observation, make a new 'center'
      if ( min(distances) > alpha ) {
        k <- k + 1
        labels[i] <- k
        centers <- rbind( centers, data[i, ] )
      }else {
        labels[i] <- which( distances == min(distances) )			
      }
    }
    
    # re-calculate cluster means after each iteration
    for (i in 1:k) {
      clust <- labels == i
      if (sum(clust) > 1) {
        centers[i, ] <- apply( data[ labels == i, ], 2, mean )
      }
    } 
  }
  values = list( labels = labels, 
                 centers = centers, 
                 k = k,
                 data = data,
                 distance = distance)
  return (values)
}

#' Plot the output of the infinite mixture model 
#' 
#' Plot Plot PLot 
#' @param output The S3 object returned from the infinite_cluster function
#' @export
#' @keywords plot
#' 
plot.infinite_mm <- function(output) {
  plot(output$data[,1], output$data[,2], col = output$labels)
}

#' Classify a new observation with the infinite mm
#'
#' Picks the closest cluster for a new observation, or a set of observations
#'
#' @param new_obs A vector of observations to classify
#' @param model The output from the infinite_mm function (S3 object)
#' @export
#' @keywords predict, infinite_mm
#' 
predict.infinite_mm <- function(new_obs, model) {
  centers <- model$centers
  if (length(new.obs) != ncol(centers) ) {
    stop("New Observation is of wrong dimension")
  }
  # get predictions for all of the observations
  lapply(new_obs, function(x) { 
    distances <- apply( centers, 1, model$distance, x)
    return (which(distances == min(distances)))  
  })
}