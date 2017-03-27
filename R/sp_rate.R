#' Simple Projection Model with Growth Rate
#'
#' A simple projection model with a constnat growth rate for the population.
#'
#' @param n A numeric value for the number of projection steps to run the model.
#' @param NO A numeric value for the inital total population size.
#' @param r A numeric value for the assumed constant future demographic growth rate.
#'
#' @return A vector of length n + 1 containing the intial total population size (N0) and the projected total population size when assuming a future growth rate of r
#' @export
#'
#' @examples
sp_rate <- function(n = NULL, N0 = NULL, r = NULL){
  NN <- rep(NA, times = n + 1)
  NN[1] <- N0
  for(i in 1:n){
    NN[i+1] <- NN[i] * (1 + r)
  }
  return(NN)
}

