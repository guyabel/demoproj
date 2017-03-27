#' Title
#'
#' @param n
#' @param N0
#' @param r
#'
#' @return
#' @export
#'
#' @examples
dsp_rate <- function(n = NULL, N0 = NULL, r = NULL){
  NN <- rep(NA, times = n + 1)
  NN[1] <- N0
  for(i in 1:n){
    NN[i+1] <- NN[i] * (1 + r[i])
  }
  return(NN)
}
