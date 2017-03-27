#' Title
#'
#' @param n
#' @param N0
#' @param b
#' @param d
#' @param m
#'
#' @return
#' @export
#'
#' @examples
sp_demo <- function(n = NULL, N0 = NULL, b = NULL, d = NULL, m = NULL){
  NN <- rep(NA, times = n + 1)
  NN[1] <- N0
  for(i in 1:n){
    NN[i+1] <- NN[i] * (1 + b - d + m)
  }
  return(NN)
}
