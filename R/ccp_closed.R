#' Cohort Component Population Projection Model based on a Closed Population.
#'
#' A cohort component projection model based on a closed population,
#' \deqn{ \mathbf N(t+1) =  \mathbf L[t, t+1] \mathbf N(t)  }
#' where the Leslie matrix, eqn{\mathbf L}, is created given user defined age specific fertility and survivorship rates.
#'
#' @param n Numeric value for the number of projection steps.
#' @param x Vector containing a character string of age group labels.
#' @param p Numeric value for step size of the population projection.
#' @param Nx_f,Nx_m Vectors containing numeric values of the initial female and male population size in each age group (\code{x}).
#' @param sx_f,sx_m,fx Vectors containing numeric values of the age specific female and male survival and fertility rates.
#'
#' If \code{ccp_closed0} is used a single vector is required.
#'
#' If \code{ccp_closed} is used a matrix of period specific rates is required, where rows represent age groups (females in top rows, males in bottom rows) and columns future periods. If a single vector is passed to \code{ccp_closed} a matrix based on constant assumptions in all future rates will be constructed.
#' @param sn,sex_ratio Numeric value of the survivorship of new-born babies from birth to the end of the interval and the sex ratio at birth of new-born babies.
#'
#' If \code{ccp_closed0} is used a single value is required.
#'
#' If \code{ccp_closed} is used a vector of period specific values is required. If a single value is passed to \code{ccp_closed} a vector based on constant assumptions in all future rates will be constructed.
#' @param tidy_output Logical value to indicate if projection output should be in a tidy data format (\code{TRUE}, the default) or as a \code{matrix} where rows represent age and gender groups and columns the projection year.
#' @param age_lab,gender_lab Vector containing a character string of age and gender group labels. Only used if projection output is in a tidy data format. See \code{tidy_pp} for more details.
#' @param ... Additional arguments passed to \code{\link{tidy_pp}} to build a tidy data frame (if chosen from \code{tidy_output}).
#'
#' @return Projected populations by age and gender for \code{n} future steps, given the age specific fertility and survivorship rates. Depending on the \code{tidy_output} value the projections will be returned as either a matrix or a tibble. Both versions contain the initial population sizes given in \code{Nx}.
#'
#' \code{ccp_closed0} produces population projections based strictly on constant future rates.
#'
#' \code{ccp_closed} produces population projections based non-constant future rates.
#'
#' @export
#'
#' @examples
#' df0 <- sweden1993
#'
#' # matrix output
#' ccp_closed0(n = 5, x = df0$x, p = 5, Nx_f = df0$Nx_f, Nx_m = df0$Nx_f,
#'             sx_f = df0$sx_f, sx_m = df0$sx_f,
#'             fx = df0$fx, sn_f = df0$Lx_f[1]/(5*100000), sn_m = df0$Lx_m[1]/(5*100000),
#'             tidy_output = FALSE)
#'
#' # tidy data frame output
#' ccp_closed0(n = 5, x = df0$x, p = 5, Nx = df0$Nx_f,
#'              sx = df0$sx_f, fx = df0$fx, sn = df0$Lx_f[1]/(5*100000),
#'              year0 = 1993, age_lab = df0$age)
#'
#' # setting up non-constant future age specific fertility rates
#' ff <- matrix(df0$fx, nrow = length(df0$fx), ncol = 5)
#' ff <- sweep(ff, 2, seq(from = 1, to = 1.5, length = 5), "*")
#' # tfr increase
#' 5 * colSums(ff)
#' # run projection with increasing fx, sx remains constant
#' ccp_closed(n = 5, x = df0$x, p = 5, Nx_f = df0$Nx_f, Nx_m = df0$Nx_f,
#'             sx_f = df0$sx_f, sx_m = df0$sx_f,
#'             fx = ff, sn_f = df0$Lx_f[1]/(5*100000), sn_m = df0$Lx_m[1]/(5*100000)),
#'             tidy_output = FALSE)
ccp_closed0 <- function(n = NULL, x = df0$x, p = NULL, Nx_f= NULL, Nx_m= NULL,
                        sx_f = NULL, sx_m = NULL,
                        fx = NULL, sn_f = NULL, sn_m = NULL, sex_ratio = 1/(1 + 1.05),
                        tidy_output = TRUE, age_lab = x, gender_lab = c("Female", "Male"), ...){
  require(dplyr)
  xx <- length(x)

  if(length(sx_f) != xx | length(sx_m) != xx)
    stop("sx_f and sx_m must be of the same length as x")
  if(length(fx) != xx)
    stop("fx must be of the same length as x")

  #template block matrices
  L_f <- L_m <- B_m <- Z <- matrix(0, nrow = xx, ncol = xx)
  #females
  L_f[2:xx, 1:(xx-1)] <- diag(sx_f[-xx])
  L_f[xx, xx] <- sx_f[xx]
  L_f[1, 1:xx] <- p * sn_f * sex_ratio * 0.5 * (fx + lead(fx) * sx_f)
  #males (surviving)
  L_m[2:xx, 1:(xx-1)] <- diag(sx_m[-xx])
  L_m[xx, xx] <- sx_m[xx]
  #males (births)
  B_m[1, 1:xx] <- p * sn_m * (1 - sex_ratio) * 0.5 * (fx + lead(fx) * sx_f)
  #bring the blocks together
  L1 <- cbind(L_f, Z)
  L2 <- cbind(B_m, L_m)
  L <- rbind(L1, L2)
  L[is.na(L)] <- 0

  WW <- matrix(NA, nrow = xx*2, ncol = n+1)
  WW[,1] <- c(Nx_f, Nx_m)
  for(i in 1:n){
    WW[,i+1] <- L %*% WW[,i]
  }

  if(tidy_output == TRUE)
    WW <- tidy_pp(proj_mat = WW, steps = p, age_lab = age_lab, gender_lab = gender_lab, ...)
  return(WW)
}
#' @export
#' @rdname ccp_closed0
ccp_closed <- function(n = NULL, x = df0$x, p = NULL, Nx_f= NULL, Nx_m= NULL,
                      sx_f = NULL, sx_m = NULL,
                      fx = NULL, sn_f = NULL, sn_m = NULL, sex_ratio = 1/(1 + 1.05),
                      tidy_output = TRUE, age_lab = x, gender_lab = c("Female", "Male"), ...){
  require(dplyr)
  xx <- length(x)

  if(!is.matrix(sx_f))
    sx_f <- matrix(sx_f, nrow = xx, ncol = n)
  if(!is.matrix(sx_m))
    sx_m <- matrix(sx_m, nrow = xx, ncol = n)
  if(!is.matrix(fx))
    fx <- matrix(fx, nrow = xx, ncol = n)
  if(length(sn_f) == 1)
    sn_f <- rep(x = sn_f, times = n)
  if(length(sn_m) == 1)
    sn_m <- rep(x = sn_m, times = n)
  if(length(sex_ratio) == 1)
    sex_ratio <- rep(x = sex_ratio, times = n)

  if(nrow(sx_f) != xx | ncol(sx_f) != n)
    stop("sx_f must have the same number of rows as x and the same number of columns as n")
  if(nrow(sx_m) != xx | ncol(sx_m) != n)
    stop("sx_m must have the same number of rows as x and the same number of columns as n")
  if(nrow(fx) != xx | ncol(fx) != n)
    stop("fx must have the same number of rows as x and the same number of columns as n")
  if(length(sn_f) != n)
    stop("sn_f must have the same values as n")
  if(length(sn_m) != n)
    stop("sn_m must have the same values as n")
  if(length(sex_ratio) != n)
    stop("sex_ratio must have the same values as n")

  WW <- matrix(NA, nrow = xx*2, ncol = n+1)
  WW[,1] <- c(Nx_f, Nx_m)
  for(i in 1:n){
    L_f <- L_m <- B_m <- Z <- matrix(0, nrow = xx, ncol = xx)
    L_f[2:xx, 1:(xx-1)] <- diag(sx_f[-xx, i])
    L_f[xx, xx] <- sx_f[xx, i]
    L_f[1, 1:xx] <- p * sn_f[i] * sex_ratio[i] * 0.5 * (fx[, i] + lead(fx[, i]) * sx_f[, i])
    L_m[2:xx, 1:(xx-1)] <- diag(sx_m[-xx, i])
    L_m[xx, xx] <- sx_m[xx]
    B_m[1, 1:xx] <- p * sn_m[i] * (1 - sex_ratio[i]) * 0.5 * (fx[, i] + lead(fx[, i]) * sx_f[, i])
    L1 <- cbind(L_f, Z)
    L2 <- cbind(B_m, L_m)
    L <- rbind(L1, L2)
    L[is.na(L)] <- 0

    WW[,i+1] <- L %*% WW[,i]
  }

  if(tidy_output == TRUE)
    WW <- tidy_pp(proj_mat = WW, steps = p, age_lab = age_lab, gender_lab = gender_lab, ...)
  return(WW)
}
