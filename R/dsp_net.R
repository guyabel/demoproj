#' Title
#'
#' @param n
#' @param N0
#' @param b
#' @param d
#' @param M
#'
#' @return
#' @export
#'
#' @examples
dsp_net <- function(n = NULL, N0 = NULL, b = NULL, d = NULL, M = NULL){
  NN <- rep(NA, times = n + 1)
  NN[1] <- N0
  for(i in 1:n){
    NN[i+1] <- (NN[i] + M[i]/2) * (1 + b[i] - d[i]) + M[i]/2
  }
  return(NN)
}
