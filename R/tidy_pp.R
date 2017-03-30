#' Create a Tidy Data Frame from Population Projection Matrix Objects
#'
#' Gathers a population projection matrix and converts it into a tidy data frame with labels for projection year, age and gender groups.
#'
#' @param proj_mat Matrix of population projection results. Columns represent projection years, rows represent age and gender groups. If composed of both male and female groups, female population projections should be in the top half of the matrix.
#' @param year0 Numeric value for the base year of the population projection.
#' @param steps Numeric value for step size of the population projection.
#' @param age_lab Vector containing a character string of age group labels. The length of \code{age_lab} multiplied by the length of \code{gender_lab} must be equal to the number of rows in \code{proj_mat}
#' @param gender_lab Vector containing a character string of gender group labels. The length of \code{age_lab} multiplied by the length of \code{gender_lab} must be equal to the number of rows in \code{proj_mat}
#' @param year_col Character string for the name of the year labels column.
#' @param age_col Character string for the name of the age labels column.
#' @param gender_col Character string for the name of the gender labels column.
#' @param pop_col Character string for the name of the projected population counts column.
#' @param denom Numeric value of a demominator to divide the population size of each age-gender group by.
#'
#' @return A tibble containing four columns. The first three contain information on the year, age and gender of the projection population. The fourth contain information in the projected population sizes or each year-age-gender combination.
#'
#' @export
#' @importFrom magrittr %>%
#' @importFrom tidyr gather
#'
#' @examples
#' df0 <- sweden1993
#'
#' # matrix output
#' pp <- fccp_closed0(n = 5, x = df0$x, p = 5, Nx = df0$Nx_f,
#'                    sx = df0$sx_f, fx = df0$fx, sn = df0$Lx_f[1]/(5*100000),
#'                    tidy_output = FALSE)
#' pp
#'
#' tidy_pp(proj_mat = pp, year0 = 1993, steps = 5,
#'         age_lab = df0$age, gender_lab = "Female")
tidy_pp <- function(proj_mat = NULL, year0 = 0, steps = NULL,
                    age_lab = NULL, gender_lab = c("Female", "Male"),
                    year_col = "Year", age_col = "Age", gender_col = "Gender", pop_col = "Population",
                    denom = 1000){
  g <- length(gender_lab)
  a <- length(age_lab)
  if(nrow(proj_mat) != a*g)
    stop("Number of rows of proj_mat must match the length of age_lab * length of gender_lab")

  df0 <- proj_mat %>%
    dplyr::tbl_df() %>%
    #to allow for one or two genders in proj_mat
    dplyr::mutate(age = rep(age_lab, times = g),
           gender = rep(gender_lab, each = a) ) %>%
    tidyr::gather(key = year, value = pop, -age, -gender) %>%
    #format Year and Population depending on the argument inputs
    dplyr::mutate(year = gsub(pattern = "V", replacement = "", x = year),
           year = as.numeric(year),
           year = year0 + (year*steps - steps),
           pop = pop/denom) %>%
    dplyr::select(year, gender, age, pop) %>%
    stats::setNames(c(year_col, gender_col, age_col, pop_col))
  return(df0)
}
