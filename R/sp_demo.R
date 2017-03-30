#' Simple Projection Model Based on Crude Demographic Rates
#'
#' A simple projection model based on crude demographic rates,
#' \deqn{ N_{t+1} =  N_t + \times (1 + b_{t,t+1} - d_{t,t+1} + m_{t,t+1}) }
#'
#' @param n Numeric value for the number of projection periods to run the model.
#' @param N0 Numeric value for the initial total population size.
#' @param b,d,m Numeric value for the crude birth, death and net migration rate in each projection period.
#'
#' If \code{sp_demo0} is used a single value is required.
#'
#' If \code{sp_demo} is used a vector of period specific rates is required. If a single rate values are passed to \code{sp_demo} a vector based on constant assumptions in all future rates will be constructed.
#' @return A vector of length \code{n + 1} containing the intial total population size (\code{N0}) and the projected total population sizes given \code{b}, \code{d} and \code{m}.
#'
#' \code{sp_demo0} produces population projections based strictly on constant future rates.
#'
#' \code{sp_demo} produces population projections based non-constant future rates.
#'
#' @export
#'
#' @examples
#' # constant future assumptions
#' sp_demo0(n = 20, N0 = 100, b = 0.04, d = 0.02, m = 0.03)
#'
#' # non-constant future assumptions
#' dd <- rnorm(n = 20, mean = 0.02, sd = 0.001)
#' mm <- rnorm(n = 20, mean = 0.03, sd = 0.01)
#' sp_demo(n = 20, N0 = 100, b = 0.04, d = dd, m = mm)
sp_demo0 <- function(n = NULL, N0 = NULL, b = NULL, d = NULL, m = NULL){
  NN <- rep(NA, times = n + 1)
  NN[1] <- N0
  for(i in 1:n){
    NN[i+1] <- NN[i] * (1 + b - d + m)
  }
  return(NN)
}
#' @export
#' @rdname sp_demo0
sp_demo <- function(n = NULL, N0 = NULL, b = NULL, d = NULL, m = NULL){
  if(length(b) == 1)
    b <- rep(x = b, times = n)
  if(length(d) == 1)
    d <- rep(x = d, times = n)
  if(length(m) == 1)
    M <- rep(x = m, times = n)

  if(length(b) != n)
    stop("b vector needs to be the length of the n argument.")
  if(length(d) != n)
    stop("d vector needs to be the length of the n argument.")
  if(length(m) != n)
    stop("m vector needs to be the length of the n argument.")

  NN <- rep(NA, times = n + 1)
  NN[1] <- N0
  for(i in 1:n){
    NN[i+1] <- NN[i] * (1 + b[i] - d[i] + m[i])
  }
  return(NN)
}
