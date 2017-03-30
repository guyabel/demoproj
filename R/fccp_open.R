#' Cohort Component Population Projection Model based on a Female Population with Emigration and Immigration.
#'
#' A cohort component projection model based on a closed female population,
#' \deqn{ \mathbf N(t+5) =  \mathbf L[t,t+5] \left( \mathbf N(t) + \frac{1}{2} \mathbf I[t,t+5] \right) + \frac{1}{2} \mathbf I[t,t+5]    }
#' where the Leslie matrix, eqn{\mathbf L}, is created given user defined age specific fertility and survivorship rates.
#'
#' @param n Numeric value for the number of projection steps.
#' @param x Vector containing a character string of age group labels.
#' @param p Numeric value for step size of the population projection.
#' @param Nx Vector containing numeric values of the initial female population size in each age group (\code{x}).
#' @param sx,fx,ex,Ix Vectors containing numeric values of the age specific female survival, fertility and emigration rates and immigration counts.
#'
#' If \code{fccp_open0} is used a single vector is required.
#'
#' If \code{fccp_open} is used a matrix of period specific rates (or counts for immigration) is required, where rows represent age groups and columns future periods. If a single vector is passed to \code{fccp_open} a matrix based on constant assumptions in all future rates and immigration counts will be constructed.
#' @param sn,sex_ratio Numeric value of the survivorship of new-born female babies from birth to the end of the interval and the sex ratio at birth of new-born babies.
#'
#' If \code{fccp_open0} is used a single value is required.
#'
#' If \code{fccp_open} is used a vector of period specific values is required. If a single value is passed to \code{fccp_open} a vector based on constant assumptions in all future rates will be constructed.
#' @param tidy_output Logical value to indicate if projection output should be in a tidy data format (\code{TRUE}, the default) or as a \code{matrix} where rows represent age and gender groups and columns the projection year.
#' @param age_lab,gender_lab Vector containing a character string of age and gender group labels. Only used if projection output is in a tidy data format. See \code{tidy_pp} for more details.
#' @param ... Additional arguments passed to \code{\link{tidy_pp}} to build a tidy data frame (if chosen from \code{tidy_output}).
#'
#' @return Projected populations by age and gender for \code{n} future steps, given the age specific fertility, survivorship and emigration rates and immigration counts. Depending on the \code{tidy_output} value the projections will be returned as either a matrix or a tibble. Both versions contain the initial population sizes given in \code{Nx}.
#'
#' \code{fccp_open0} produces population projections based strictly on constant future rates.
#'
#' \code{fccp_open} produces population projections based non-constant future rates.
#'
#' @export
#'
#' @examples
#' data(sweden1993)
#' df0 <- sweden1993 %>%
#'     #immigration from estimated count and net age pattern
#'     mutate(Ix_f = 138000 * Mx_f/sum(Mx_f),
#'     Ix_f = round(Ix_f),
#'     #emigration as remainder
#'     Ex_f = Ix_f - Mx_f ,
#'     #emigration rate for projection model
#'     ex_f = Ex_f/Nx_f)
#'
#' # matrix output
#' fccp_open0(n = 5, x = df0$x, p = 5, Nx = df0$Nx_f,
#'            sx = df0$sx_f,
#'            fx = df0$fx, sn = df0$Lx_f[1]/(5*100000),
#'            ex = df0$ex_f, Ix = df0$Ix_f,
#'            tidy_output = FALSE)
#'
#' # tidy data frame output
# fccp_open0(n = 5, x = df0$x, p = 5, Nx = df0$Nx_f,
#            sx = df0$sx_f,
#            fx = df0$fx, sn = df0$Lx_f[1]/(5*100000),
#            ex = df0$ex_f, Ix = df0$Ix_f,
#            year0 = 1993, age_lab = df0$age)
#'
#' # setting up non-constant future immigrant counts
#' II <- matrix(df0$Ix_f, nrow = length(df0$Ix_f), ncol = 5)
#' II <- sweep(II, 2, seq(from = 1, to = 1.5, length = 5), "*")
#' # immigration increase
#' colSums(II)
#' # run projection with increasing immigration, fx, sx and ex remains constant
#' fccp_open(n = 5, x = df0$x, p = 5, Nx = df0$Nx_f,
#'           sx = df0$sx_f,
#'           fx = df0$fx, sn = df0$Lx_f[1]/(5*100000),
#'           ex = df0$ex_f, Ix = II,
#'           tidy_output = FALSE)
fccp_open0 <- function(n = NULL, x = df0$x, p = NULL, Nx = NULL,
                       sx = NULL,
                       fx = NULL, sn = NULL, sex_ratio = 1/(1 + 1.05),
                       ex = NULL, Ix = NULL,
                       tidy_output = TRUE, age_lab = x, gender_lab = "Female", ...){
  require(dplyr)
  xx <- length(x)

  if(length(sx) != xx)
    stop("sx must be of the same length as x")
  if(length(fx) != xx)
    stop("fx must be of the same length as x")
  if(length(ex) != xx)
    stop("ex must be of the same length as x")
  if(length(ex) != xx)
    stop("Ix must be of the same length as x")

  L <- matrix(0, nrow = xx, ncol = xx)
  L[2:xx, 1:(xx-1)] <- diag(sx[-xx] - ex[-xx])
  L[xx, xx] <- sx[xx] - ex[xx]
  L[1, 1:xx] <- p * sn * sex_ratio * 0.5 * (fx + lead(fx) * (sx - ex))
  L[is.na(L)] <- 0

  WW <- matrix(NA, nrow = xx, ncol = n + 1)
  WW[,1] <- Nx
  for(i in 1:n){
    WW[,i+1] <- L %*% (WW[,i] + Ix/2) + Ix/2
  }

  if(tidy_output == TRUE)
    WW <- tidy_pp(proj_mat = WW, steps = p, age_lab = age_lab, gender_lab = gender_lab, ...)
  return(WW)
}
#' @export
#' @rdname fccp_open0
fccp_open <- function(n = NULL, x = df0$x, p = NULL, Nx = NULL,
                      sx = NULL,
                      fx = NULL, sn = NULL, sex_ratio = 1/(1 + 1.05),
                      ex = NULL, Ix = NULL,
                      tidy_output = TRUE, age_lab = x, gender_lab = "Female", ...){
  require(dplyr)
  xx <- length(x)

  if(!is.matrix(sx))
    sx <- matrix(sx, nrow = xx, ncol = n)
  if(!is.matrix(fx))
    fx <- matrix(fx, nrow = xx, ncol = n)
  if(!is.matrix(ex))
    ex <- matrix(ex, nrow = xx, ncol = n)
  if(!is.matrix(Ix))
    Ix <- matrix(Ix, nrow = xx, ncol = n)
  if(length(sn) == 1)
    sn <- rep(x = sn, times = n)
  if(length(sex_ratio) == 1)
    sex_ratio <- rep(x = sex_ratio, times = n)

  if(nrow(sx) != xx | ncol(sx) != n)
    stop("sx must have the same number of rows as x and the same number of columns as n")
  if(nrow(fx) != xx | ncol(fx) != n)
    stop("fx must have the same number of rows as x and the same number of columns as n")
  if(nrow(ex) != xx | ncol(ex) != n)
    stop("ex must have the same number of rows as x and the same number of columns as n")
  if(nrow(Ix) != xx | ncol(Ix) != n)
    stop("Ix must have the same number of rows as x and the same number of columns as n")
  if(length(sn) != n)
    stop("sn must have the same values as n")
  if(length(sex_ratio) != n)
    stop("sex_ratio must have the same values as n")

  WW <- matrix(NA, nrow = xx, ncol = n + 1)
  WW[,1] <- Nx
  for(i in 1:n){
    L <- matrix(0, nrow = xx, ncol = xx)
    L[2:xx, 1:(xx-1)] <- diag(sx[-xx, i] - ex[-xx, i])
    L[xx, xx] <- sx[xx, i] - ex[xx, i]
    L[1, 1:xx] <- p * sn[i] * sex_ratio[i] * 0.5 * (fx[,i] + lead(fx[,i]) * (sx[,i] - ex[,i]) )
    L[is.na(L)] <- 0

    WW[,i+1] <- L %*% (WW[,i] + Ix[,i]/2) + Ix[,i]/2
  }

  if(tidy_output == TRUE)
    WW <- tidy_pp(proj_mat = WW, steps = p, age_lab = age_lab, gender_lab = gender_lab, ...)
  return(WW)
}
