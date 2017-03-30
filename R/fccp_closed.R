#' Cohort Component Population Projection Model based on a Closed Female Population.
#'
#' A cohort component projection model based on a closed female population,
#' \deqn{ \mathbf N(t+1) =  \mathbf L[t, t+1] \mathbf N(t)  }
#' where the Leslie matrix, eqn{\mathbf L}, is created given user defined age specific fertility and survivorship rates.
#'
#' @param n Numeric value for the number of projection steps.
#' @param x Vector containing a character string of age group labels.
#' @param p Numeric value for step size of the population projection.
#' @param Nx Vector containing numeric values of the initial female population size in each age group (\code{x}).
#' @param sx,fx Vectors containing numeric values of the age specific female survival and fertility rates.
#'
#' If \code{fccp_closed0} is used a single vector is required.
#'
#' If \code{fccp_closed} is used a matrix of period specific rates is required, where rows represent age groups and columns future periods. If a single vector is passed to \code{fccp_closed} a matrix based on constant assumptions in all future rates will be constructed.
#' @param sn,sex_ratio Numeric value of the survivorship of new-born female babies from birth to the end of the interval and the sex ratio at birth of new-born babies.
#'
#' If \code{fccp_closed0} is used a single value is required.
#'
#' If \code{fccp_closed} is used a vector of period specific values is required. If a single value is passed to \code{fccp_closed} a vector based on constant assumptions in all future rates will be constructed.
#' @param tidy_output Logical value to indicate if projection output should be in a tidy data format (\code{TRUE}, the default) or as a \code{matrix} where rows represent age and gender groups and columns the projection year.
#' @param age_lab,gender_lab Vector containing a character string of age and gender group labels. Only used if projection output is in a tidy data format. See \code{tidy_pp} for more details.
#' @param ... Additional arguments passed to \code{\link{tidy_pp}} to build a tidy data frame (if chosen from \code{tidy_output}).
#'
#' @return Projected populations by age and gender for \code{n} future steps, given the age specific fertility and survivorship rates. Depending on the \code{tidy_output} value the projections will be returned as either a matrix or a tibble. Both versions contain the initial population sizes given in \code{Nx}.
#'
#' \code{fccp_closed0} produces population projections based strictly on constant future rates.
#'
#' \code{fccp_closed} produces population projections based non-constant future rates.
#'
#' @export
#'
#' @examples
#' df0 <- sweden1993
#'
#' # matrix output
#' fccp_closed0(n = 5, x = df0$x, p = 5, Nx = df0$Nx_f,
#'              sx = df0$sx_f, fx = df0$fx, sn = df0$Lx_f[1]/(5*100000),
#'              tidy_output = FALSE)
#'
#' # tidy data frame output
#' fccp_closed0(n = 5, x = df0$x, p = 5, Nx = df0$Nx_f,
#'              sx = df0$sx_f, fx = df0$fx, sn = df0$Lx_f[1]/(5*100000),
#'              year0 = 1993, age_lab = df0$age)
#'
#' # setting up non-constant future age specific fertility rates
#' ff <- matrix(df0$fx, nrow = length(df0$fx), ncol = 5)
#' ff <- sweep(ff, 2, seq(from = 1, to = 1.5, length = 5), "*")
#' # tfr increase
#' 5 * colSums(ff)
#' # run projection with increasing fx, sx remains constant
#' fccp_closed(n = 5, x = df0$x, p = 5, Nx = df0$Nx_f,
#'             sx = df0$sx_f, fx = ff, sn = df0$Lx_f[1]/(5*100000),
#'             tidy_output = FALSE)
fccp_closed0 <- function(n = NULL, x = NULL, p = NULL, Nx = NULL,
                         sx = NULL,
                         fx = NULL, sn = NULL, sex_ratio = 1/(1 + 1.05),
                         tidy_output = TRUE, age_lab = x, gender_lab = "Female", ...){
  require(dplyr)
  xx <- length(x)

  if(length(sx) != xx)
    stop("sx must be of the same length as x")
  if(length(fx) != xx)
    stop("fx must be of the same length as x")

  L <- matrix(0, nrow = xx, ncol = xx)
  L[2:xx, 1:(xx-1)] <- diag(sx[-xx])
  L[xx, xx] <- sx[xx]
  L[1, 1:xx] <- p * sn * sex_ratio * 0.5 * (fx + dplyr::lead(fx) * sx)
  L[is.na(L)] <- 0

  WW <- matrix(NA, nrow = xx, ncol = n + 1)
  WW[,1] <- Nx
  for(i in 1:n){
    WW[,i+1] <- L %*% WW[,i]
  }

  if(tidy_output == TRUE)
    WW <- tidy_pp(proj_mat = WW, steps = p, age_lab = age_lab, gender_lab = gender_lab, ...)
  return(WW)
}
#' @export
#' @rdname fccp_closed0
fccp_closed <- function(n = NULL, x = NULL, p = NULL, Nx = NULL,
                        sx = NULL,
                        fx = NULL, sn = NULL, sex_ratio = 1/(1 + 1.05),
                        tidy_output = TRUE, age_lab = x, gender_lab = "Female", ...){
  xx <- length(x)

  if(!is.matrix(sx))
    sx <- matrix(sx, nrow = xx, ncol = n)
  if(!is.matrix(fx))
    fx <- matrix(fx, nrow = xx, ncol = n)
  if(length(sn) == 1)
    sn <- rep(x = sn, times = n)
  if(length(sex_ratio) == 1)
    sex_ratio <- rep(x = sex_ratio, times = n)

  if(nrow(sx) != xx | ncol(sx) != n)
    stop("sx must have the same number of rows as x and the same number of columns as n")
  if(nrow(fx) != xx | ncol(fx) != n)
    stop("fx must have the same number of rows as x and the same number of columns as n")
  if(length(sn) != n)
    stop("sn must have the same values as n")
  if(length(sex_ratio) != n)
    stop("sex_ratio must have the same values as n")

  WW <- matrix(NA, nrow = xx, ncol = n + 1)
  WW[,1] <- Nx
  for(i in 1:n){
    L <- matrix(0, nrow = xx, ncol = xx)
    L[2:xx, 1:(xx-1)] <- diag(sx[-xx, i])
    L[xx, xx] <- sx[xx, i]
    L[1, 1:xx] <- p * sn[i] * sex_ratio[i] * 0.5 * (fx[,i] + dplyr::lead(fx[,i]) * sx[,i])
    L[is.na(L)] <- 0

    WW[,i+1] <- L %*% WW[,i]
  }

  if(tidy_output == TRUE)
    WW <- tidy_pp(proj_mat = WW, steps = p, age_lab = age_lab, gender_lab = gender_lab, ...)
  return(WW)
}
