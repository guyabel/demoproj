#' Simple Projection Model Based on Population Growth Rates
#'
#' A simple projection model based on  population growth rates,
#' \deqn{ N_{t} = N_{t} \times (1 + r_{t,t+1}) }
#'
#' @param n Numeric value for the number of projection periods to run the model.
#' @param N0 Numeric value for the initial total population size.
#' @param r Numeric value for the population growth rate in each projection period.
#'
#' If \code{sp_rate0} is used a single value is required.
#'
#' If \code{sp_rate} is used a vector of period specific rates is required. If a single rate value is passed to \code{sp_rate} a vector based on constant assumptions in all future rates will be constructed.
#' @return A vector of length \code{n + 1} containing the initial total population size (\code{N0}) and the projected total population sizes given \code{r}.
#'
#' \code{sp_rate0} produces population projections based strictly on constant future rates.
#'
#' \code{sp_rate} produces population projections based non-constant future rates.
#'
#' @export
#'
#' @examples
#' # constant future assumptions
#' sp_rate0(n = 20, N0 = 100, r = 0.05)
#'
#' # non-constant future assumptions
#' rr <- rnorm(n = 20, mean = 0.05, sd = 0.01)
#' sp_rate(n = 20, N0 = 100, r = rr)
sp_rate0 <- function(n = NULL, N0 = NULL, r = NULL){
  NN <- rep(NA, times = n + 1)
  NN[1] <- N0
  for(i in 1:n){
    NN[i+1] <- NN[i] * (1 + r)
  }
  return(NN)
}
#' @export
#' @rdname sp_rate0
sp_rate <- function(n = NULL, N0 = NULL, r = NULL){
  if(length(r) == 1)
    r <- rep(x = r, times = n)

  if(length(r) != n)
    stop("r vector needs to be the length of the n argument.")

  NN <- rep(NA, times = n + 1)
  NN[1] <- N0
  for(i in 1:n){
    NN[i+1] <- NN[i] * (1 + r[i])
  }
  return(NN)
}
