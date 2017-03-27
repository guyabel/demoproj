#' Title
#'
#' @param n
#' @param N0
#' @param b
#' @param d
#' @param e
#' @param I
#'
#' @return
#' @export
#'
#' @examples
sp_open <- function(n = NULL, N0 = NULL, b = NULL, d = NULL, e = NULL, I = NULL){
  NN <- rep(NA, times = n + 1)
  NN[1] <- N0
  for(i in 1:n){
    NN[i+1] <- (NN[i] + I/2) * (1 + b - d - e) + I/2
  }
  return(NN)
}
