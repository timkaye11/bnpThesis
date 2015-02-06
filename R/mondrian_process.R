#' A R function for a one-dimensional mondrian process
#' 
#' The mondrian process (http://danroy.org/papers/RoyTeh-NIPS-2009.pdf) is a point process
#' that is used in the Mondrian forest algorithm. It creates partitions based on a budget value, 
#' and a specified interval
#' 
#' @param budget The budget value for the mondrian process. Higher values lead to more splits
#' @param interval The interval to create the splits on. The default is set to [0,1]
#' @keywords mondrian process, process
#' @export 
#' @examples
#' mondrian <- mondrian_1d(5)
#' mondrian
#'
mondrian_1d <- function( budget, interval = c(0, 1) ) {
  # the range for the interval 
  interval_range <- diff(interval)
  
  # the cost function, which is reduced from the budget
  cost <- rexp(1, interval_range)
  
  # update the budget
  new_budget <- budget - cost
  # the base case for the recursion
  if (new_budget < 0) {
    return (interval)
  }
  
  split_point <- runif(1, interval[1], interval[2])
  interval_left <- c(interval[1], split_point)
  interval_right <- c(split_point, interval[2])
  # recurse on the left and right sides w/ the new budget
  return (rbind(mondrian_1d(new_budget, interval_left),
                mondrian_1d(new_budget, interval_right)))
}

#' Multi-dimensional Mondrian Process
#' 
#' Returns the split locations for all dimensions in a multi dimensional mondrian process, 
#' from a specified domain matrix and a budget. More information about this can be found in 
#' (http://danroy.org/papers/RoyTeh-NIPS-2009.pdf). 
#' 
#' @param budget The budget for the Mondrian Process
#' @param domain A matrix/dataframe of intervals for the mondrian process. The default is set to [0,1][0,1]
#' @export
#' @keywords mondrian process
#' @examples
#' splits <- mondrian_2d(5)
#' 
mondrian_2d <- function( budget, domain = matrix(c(0,0,1,1), nrow =2)) {
  # cumulative sum of the interval sizes
  interval <- c(0, cumsum(apply(domain, 1, diff)))
  
  # calculate the cost based on the largest interval
  cost <- rexp(1, interval[length(interval)])
  if (cost > budget) {
    return (domain)
  }
  else {	
    new_budget <- budget - cost
    cut <- runif(1, 0, 1) * interval[length(interval)]
    d <- bisect(interval, cut) 
    split_size <- cut - interval[d]
    
    # handle the splits differently depending on the location of the bisection
    if (d == 1) {
      piece_left  <- rbind( c(domain[d, 1], domain[d, 1] + split_size), domain[d+1, ])
      piece_right <- rbind( c(domain[d, 1] + split_size, domain[d, 2]), domain[d+1, ])
    }
    else {
      piece_left  <- rbind( domain[1, ], c(domain[d, 1], domain[d, 1] + split_size))
      piece_right <- rbind( domain[1, ], c(domain[d, 1] + split_size, domain[d, 2]))
    }
    # recursively complete the mondrian process w/ new budget
    return (list(mondrian_2d(new_budget, piece_left),
                 mondrian_2d(new_budget, piece_right)))
  }
}

#' Plots of Mondrian Process, for visualization purposes
#' 
#' A R function that plots of the Mondrian process for a specified budget value and interval ranges. 
#' Currently only supports the one dimensional and two dimensional cases
#'
#' @param budget The budget for the mondrian process
#' @param interval The interval for the mondrian process
#' @param pretty An option to add color and text to the visualization
#' @export
#' @keywords plot, mondrian
#' @examples
#' plot_mondrian(5) # single dimension
#' plot_mondrian(5, interval=matrix(c(0,0,1,1), nrow=2)) # two dimensions
#' 
plot_mondrian <- function(budget, interval=c(0,1), pretty=TRUE...) {
  # create a new canvas
  plot.new()
  
  # handle the multivariate and univariate cases separately
  if (is.matrix(interval) | is.data.frame(interval)) {
    splits <- mondrian_2d(budget, interval)
    locations <- seq(2, ncol(splits), 2)
    start <- 1
    
    # time to draw some rectangles!
    for (loc in locations) {
      split <- splits[, start:loc]
      start <- loc + 1
      
      # the sizes for the rectangle to be drawn
      width <- split[1,2] - split[1,1]
      height <- split[2,2] - split[2,1]
      area <- width * height
      
      # dimensions for the rectangle to be drawn
      x.lo <- split[1,1]
      x.hi <- split[1,2]
      y.lo <- split[2,1]
      y.hi <- split[2,2]
      
      # add color and text boxes to the rectangle if specified! 
      if (pretty) {
        # randomly pick a color. Color has no meaning, just for aesthetic purposes
        rect(x.lo, y.lo, x.hi, y.hi, col = color[sample(1:20, 1)])
      
        # write down the area on the rectangle (again for visualization purposes)
        if (width > height) {
          text(mean(c(x.lo, x.hi)), mean(c(y.lo, y.hi)), paste("Area: ", round(area, 3)), cex = 2 * area , ... ) 
        } else {
          text(mean(c(x.lo, x.hi)), mean(c(y.lo, y.hi)), paste("Area: ", round(area, 3)), cex = 2 * area, srt = 90, ... )
        }
      } 
      # in case pretty is set to FALSE
      else {
        rect(x.lo, y.lo, x.hi, y.hi)
      }
    }
  }
  # make a different plot for the one dimensional case
  else {
    splits <- mondrian_1d(budget, interval)
    axis(1, at=interval)
    rug( unique(unlist(splits[,2])), lwd=2, ...)
  }
}
