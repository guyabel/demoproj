#' Cohort Component Population Projection Model based on a Population with Net Migration.
#'
#' A cohort component projection model based on a closed population,
#' \deqn{ N(t+5) =  L[t, t+1] \left( N(t) + \frac{1}{2} M[t, t+1] \right) + \frac{1}{2} M[t, t+1]   }
#' where the Leslie matrix, eqn{L}, is created given user defined age specific fertility and survivorship rates.
#'
#' @param n Numeric value for the number of projection steps.
#' @param x Vector containing a character string of age group labels.
#' @param p Numeric value for step size of the population projection.
#' @param Nx_f,Nx_m Vectors containing numeric values of the initial female and male population size in each age group (\code{x}).
#' @param sx_f,sx_m,fx,Mx_f,Mx_m Vectors containing numeric values of the age specific female and male survival and fertility rates and female and male net migration counts.
#'
#' If \code{ccp_net0} is used a single vector is required.
#'
#' If \code{ccp_net} is used a matrix of period specific rates (or counts for net migration) is required, where rows represent age groups (females in top rows, males in bottom rows) and columns future periods. If a single vector is passed to \code{ccp_net} a matrix based on constant assumptions in all future rates and net migration counts will be constructed.
#' @param sn_f,sn_m,sex_ratio Numeric value of the female and male survivorship of new-born babies from birth to the end of the interval and the sex ratio at birth of new-born babies.
#'
#' If \code{ccp_net0} is used a single value is required.
#'
#' If \code{ccp_net} is used a vector of period specific values is required. If a single value is passed to \code{ccp_net} a vector based on constant assumptions in all future rates will be constructed.
#' @param tidy_output Logical value to indicate if projection output should be in a tidy data format (\code{TRUE}, the default) or as a \code{matrix} where rows represent age and gender groups and columns the projection year.
#' @param age_lab,gender_lab Vector containing a character string of age and gender group labels. Only used if projection output is in a tidy data format. See \code{tidy_pp} for more details.
#' @param ... Additional arguments passed to \code{\link{tidy_pp}} to build a tidy data frame (if chosen from \code{tidy_output}).
#'
#' @return Projected populations by age and gender for \code{n} future steps, given the fertility and survivorship rates and net migration counts. Depending on the \code{tidy_output} value the projections will be returned as either a matrix or a tibble. Both versions contain the initial population sizes given in \code{Nx}.
#'
#' \code{ccp_net0} produces population projections based strictly on constant future rates.
#'
#' \code{ccp_net} produces population projections based non-constant future rates.
#'
#' @export
#'
#' @examples
#' library(dplyr)
#' df0 <- sweden1993 %>%
#'    # set male net migration to the same as female
#'    mutate(Mx_m = Mx_f)
#'
#' # matrix output
#' ccp_net0(n = 5, x = df0$x, p = 5,Nx_f = df0$Nx_f, Nx_m = df0$Nx_m,
#'          sx_f = df0$sx_f, sx_m = df0$sx_m,
#'          fx = df0$fx, sn_f = df0$Lx_f[1]/(5*100000), sn_m = df0$Lx_m[1]/(5*100000),
#'          Mx_f = df0$Mx_f, Mx_m = df0$Mx_m,
#'          tidy_output = FALSE)
#'
#' # tidy data frame output
#' ccp_net0(n = 5, x = df0$x, p = 5,Nx_f = df0$Nx_f, Nx_m = df0$Nx_m,
#'          sx_f = df0$sx_f, sx_m = df0$sx_m,
#'          fx = df0$fx, sn_f = df0$Lx_f[1]/(5*100000), sn_m = df0$Lx_m[1]/(5*100000),
#'          Mx_f = df0$Mx_f, Mx_m = df0$Mx_m,
#'          year0 = 1993, age_lab = df0$age)
#'
#' # setting up non-constant future male net migrant counts
#' MM <- matrix(df0$Mx_m, nrow = length(df0$Mx_m), ncol = 5)
#' MM <- sweep(MM, 2, seq(from = 1, to = 1.5, length = 5), "*")
#' # net migration increase
#' colSums(MM)
#' # run projection with increasing net migration, fx and sx remains constant
#' ccp_net(n = 5, x = df0$x, p = 5, Nx_f = df0$Nx_f, Nx_m = df0$Nx_m,
#'         sx_f = df0$sx_f, sx_m = df0$sx_m,
#'         fx = df0$fx, sn_f = df0$Lx_f[1]/(5*100000), sn_m = df0$Lx_m[1]/(5*100000),
#'         Mx_f = df0$Mx_f, Mx_m = MM,
#'         tidy_output = FALSE)
ccp_net0 <- function(n = NULL, x = NULL, p = NULL, Nx_f= NULL, Nx_m= NULL,
                     sx_f = NULL, sx_m = NULL,
                     fx = NULL, sn_f = NULL, sn_m = NULL, sex_ratio = 1/(1 + 1.05),
                     Mx_f = NULL, Mx_m = NULL,
                     tidy_output = TRUE, age_lab = x, gender_lab = c("Female", "Male"), ...){
  xx <- length(x)

  if(length(sx_f) != xx | length(sx_m) != xx)
    stop("sx_f and sx_m must be of the same length as x")
  if(length(fx) != xx)
    stop("fx must be of the same length as x")
  if(length(Mx_f) != xx | length(Mx_m) != xx)
    stop("Mx_f and Mx_m must be of the same length as x")

  #template block matrices
  L_f <- L_m <- B_m <- Z <- matrix(0, nrow = xx, ncol = xx)
  #females
  L_f[2:xx, 1:(xx-1)] <- diag(sx_f[-xx])
  L_f[xx, xx] <- sx_f[xx]
  L_f[1, 1:xx] <- p * sn_f * sex_ratio * 0.5 * (fx + dplyr::lead(fx) * sx_f)
  #males (surviving)
  L_m[2:xx, 1:(xx-1)] <- diag(sx_m[-xx])
  L_m[xx, xx] <- sx_m[xx]
  #males (births)
  B_m[1, 1:xx] <- p * sn_m * (1 - sex_ratio) * 0.5 * (fx + dplyr::lead(fx) * sx_f)
  #bring the blocks together
  L1 <- cbind(L_f, Z)
  L2 <- cbind(B_m, L_m)
  L <- rbind(L1, L2)
  L[is.na(L)] <- 0

  WW <- matrix(NA, nrow = xx*2, ncol = n+1)
  WW[,1] <- c(Nx_f, Nx_m)
  Mx <- c(Mx_f, Mx_m)
  for(i in 1:n){
    WW[,i+1] <- L %*% (WW[,i] + Mx/2) + Mx/2
  }

  if(tidy_output == TRUE)
    WW <- tidy_pp(proj_mat = WW, steps = p, age_lab = age_lab, gender_lab = gender_lab, ...)
  return(WW)
}
#' @export
#' @rdname ccp_net0
ccp_net <- function(n = NULL, x = NULL, p = NULL, Nx_f= NULL, Nx_m= NULL,
                     sx_f = NULL, sx_m = NULL,
                     fx = NULL, sn_f = NULL, sn_m = NULL, sex_ratio = 1/(1 + 1.05),
                     Mx_f = NULL, Mx_m = NULL,
                     tidy_output = TRUE, age_lab = x, gender_lab = c("Female", "Male"), ...){

  xx <- length(x)

  if(!is.matrix(sx_f))
    sx_f <- matrix(sx_f, nrow = xx, ncol = n)
  if(!is.matrix(sx_m))
    sx_m <- matrix(sx_m, nrow = xx, ncol = n)
  if(!is.matrix(fx))
    fx <- matrix(fx, nrow = xx, ncol = n)
  if(!is.matrix(Mx_f))
    Mx_f <- matrix(Mx_f, nrow = xx, ncol = n)
  if(!is.matrix(Mx_m))
    Mx_m <- matrix(Mx_m, nrow = xx, ncol = n)
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
  if(nrow(Mx_f) != xx | ncol(Mx_f) != n)
    stop("Mx_f must have the same number of rows as x and the same number of columns as n")
  if(nrow(Mx_m) != xx | ncol(Mx_m) != n)
    stop("Mx_m must have the same number of rows as x and the same number of columns as n")
  if(length(sn_f) != n)
    stop("sn_f must have the same values as n")
  if(length(sn_m) != n)
    stop("sn_m must have the same values as n")
  if(length(sex_ratio) != n)
    stop("sex_ratio must have the same values as n")

  WW <- matrix(NA, nrow = xx*2, ncol = n+1)
  WW[,1] <- c(Nx_f, Nx_m)
  Mx <- rbind(Mx_f, Mx_m)
  for(i in 1:n){
    L_f <- L_m <- B_m <- Z <- matrix(0, nrow = xx, ncol = xx)
    L_f[2:xx, 1:(xx-1)] <- diag(sx_f[-xx, i])
    L_f[xx, xx] <- sx_f[xx, i]
    L_f[1, 1:xx] <- p * sn_f[i] * sex_ratio[i] * 0.5 * (fx[, i] + dplyr::lead(fx[, i]) * sx_f[, i])
    L_m[2:xx, 1:(xx-1)] <- diag(sx_m[-xx, i])
    L_m[xx, xx] <- sx_m[xx, i]
    B_m[1, 1:xx] <- p * sn_m[i] * (1 - sex_ratio[i]) * 0.5 * (fx[, i] + dplyr::lead(fx[, i]) * sx_f[, i])
    L1 <- cbind(L_f, Z)
    L2 <- cbind(B_m, L_m)
    L <- rbind(L1, L2)
    L[is.na(L)] <- 0

    WW[,i+1] <- L %*% (WW[,i] + Mx[, i]/2) + Mx[, i]/2
  }

  if(tidy_output == TRUE)
    WW <- tidy_pp(proj_mat = WW, steps = p, age_lab = age_lab, gender_lab = gender_lab, ...)
  return(WW)
}
