#' Simple Projection Model Based on Crude Demographic Rates and Immigration Counts
#'
#' A simple projection model based on crude fertility, survival and emigration rates and immigration counts,
#' \deqn{ N_{t+1} = \left( N_t + \frac{I_{t,t+1}}{2} \right) \times (1 + b_{t,t+1} - d_{t,t+1} - e_{t,t+1}) + \frac{I_{t,t+1}}{2} }
#'
#' @param n Numeric value for the number of projection periods to run the model.
#' @param N0 Numeric value for the initial total population size.
#' @param b,d,e,I Numeric value for the crude birth, death and emigration rate and immigration count in each projection period.
#'
#' If \code{sp_open0} is used a single value is required.
#'
#' If \code{sp_open} is used a vector of period specific rates (and immigration counts) is required. If a single rate values (and immigration count)  are passed to \code{sp_open} a vector based on constant assumptions in all future rates (and immigration counts) will be constructed.
#' @return A vector of length \code{n + 1} containing the initial total population size (\code{N0}) and the projected total population sizes given \code{b}, \code{d}, \code{e} and \code{I}.
#'
#' \code{sp_open0} produces population projections based strictly on constant future rates and counts.
#'
#' \code{sp_open} produces population projections based non-constant future rates and counts.
#'
#' @export
#'
#' @examples
#' # constant future assumptions
#' sp_open0(n = 20, N0 = 100, b = 0.04, d = 0.02, e = 0.01, I = 4)
#'
#' # non-constant future assumptions
#' bb <- rnorm(n = 20, mean = 0.04, sd = 0.001)
#' ee <- rnorm(n = 20, mean = 0.01, sd = 0.001)
#' sp_open(n = 20, N0 = 100, b = bb, d = 0.02, e = ee, I = 4)
sp_open0 <- function(n = NULL, N0 = NULL, b = NULL, d = NULL, e = NULL, I = NULL){
  NN <- rep(NA, times = n + 1)
  NN[1] <- N0
  for(i in 1:n){
    NN[i+1] <- (NN[i] + I/2) * (1 + b - d - e) + I/2
  }
  return(NN)
}
#' @export
#' @rdname sp_open0
sp_open <- function(n = NULL, N0 = NULL, b = NULL, d = NULL, e = NULL, I = NULL){
  if(length(b) == 1)
    b <- rep(x = b, times = n)
  if(length(d) == 1)
    d <- rep(x = d, times = n)
  if(length(e) == 1)
    e <- rep(x = e, times = n)
  if(length(I) == 1)
    I <- rep(x = I, times = n)

  if(length(b) != n)
    stop("b vector needs to be the length of the n argument.")
  if(length(d) != n)
    stop("d vector needs to be the length of the n argument.")
  if(length(e) != n)
    stop("e vector needs to be the length of the n argument.")
  if(length(I) != n)
    stop("e vector needs to be the length of the n argument.")

  NN <- rep(NA, times = n + 1)
  NN[1] <- N0
  for(i in 1:n){
    NN[i+1] <- (NN[i] + I[i]/2) * (1 + b[i] - d[i] - e[i]) + I[i]/2
  }
  return(NN)
}
