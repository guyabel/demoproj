#' Simple Projection Model Based on Crude Demographic Rates and Net Migration Counts.
#'
#' A simple projection model based on crude demographic rates and net migration counts,
#' \deqn{ N_{t+1} = \left( N_t + \frac{M_{t,t+1}}{2} \right) \times (1 + b_{t,t+1} - d_{t,t+1}) + \frac{M_{t,t+1}}{2} }
#'
#' @param n Numeric value for the number of projection periods to run the model.
#' @param N0 Numeric value for the initial total population size.
#' @param b,d,M Numeric value for the crude birth and death rates and net migration counts in each projection period.
#'
#' If \code{sp_net0} is used a single value is required.
#'
#' If \code{sp_net} is used a vector of period specific rates (and net migration counts) is required. If a single rate values (and net migration count)  are passed to \code{sp_net} a vector based on constant assumptions in all future rates (and net migration counts) will be constructed.
#' @return A vector of length \code{n + 1} containing the initial total population size (\code{N0}) and the projected total population sizes given \code{b}, \code{d} and \code{M}.
#'
#' \code{sp_net0} produces population projections based stricly on constant future rates and counts.
#'
#' \code{sp_net} produces population projections based non-constant future rates and counts.
#'
#' @export
#'
#' @examples
#' # constant future assumptions
#' sp_net0(n = 20, N0 = 100, b = 0.04, d = 0.02, M = 3)
#'
#' # non-constant future assumptions
#' bb <- rnorm(n = 20, mean = 0.04, sd = 0.001)
#' MM <- rnorm(n = 20, mean = 3)
#' sp_net(n = 20, N0 = 100, b = bb, d = 0.02, M = MM)
sp_net0 <- function(n = NULL, N0 = NULL, b = NULL, d = NULL, M = NULL){
  NN <- rep(NA, times = n + 1)
  NN[1] <- N0
  for(i in 1:n){
    NN[i+1] <- (NN[i] + M/2) * (1 + b - d) + M/2
  }
  return(NN)
}
#' @export
#' @rdname sp_net0
sp_net <- function(n = NULL, N0 = NULL, b = NULL, d = NULL, M = NULL){
  if(length(b) == 1)
    b <- rep(x = b, times = n)
  if(length(d) == 1)
    d <- rep(x = d, times = n)
  if(length(M) == 1)
    M <- rep(x = M, times = n)

  if(length(b) != n)
    stop("b vector needs to be the length of the n argument.")
  if(length(d) != n)
    stop("d vector needs to be the length of the n argument.")
  if(length(M) != n)
    stop("M vector needs to be the length of the n argument.")

  NN <- rep(NA, times = n + 1)
  NN[1] <- N0
  for(i in 1:n){
    NN[i+1] <- (NN[i] + M[i]/2) * (1 + b[i] - d[i]) + M[i]/2
  }
  return(NN)
}
