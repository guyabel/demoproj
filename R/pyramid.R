#' Population Pyramid Using ggplot2
#'
#' Draws a population pyramid using ggplot2.
#'
#' @param data data frame containing columns on age, gender and population size.
#' @param m Upper limit on the horizontal axis
#' @param gender_lab Vector of length two containing a character string of gender group labels in \code{data}. Female label first, male second. If left as \code{NULL}, the function will guess labels from the unique values in the \code{gender_col}.
#' @param age_col Character string for the name of the age column in \code{data}.
#' @param gender_col Character string for the name of the gender labels column in \code{data}.
#' @param pop_col Character string for the name of the projected population counts column in \code{data}.
#'
#' @return A ggplot bar charts forming a classic population pyramid.
#' @export
#'
#' @examples
#' library(dplyr)
#' library(tidyr)
#'
#' # data into tidy format for plotting
#' df0 <- sweden1993 %>%
#'   select(age, Nx_m, Nx_f) %>%
#'   gather(key = gender, value = pop, -age) %>%
#'   mutate(gender = gsub("Nx_", "", gender))
#' df0
#'
#' # plot
#' pyramid(data = df0,
#'         gender_lab = c("f","m"),
#'         age_col = "age", gender_col = "gender", pop_col = "pop")
pyramid <- function(data = NULL, m = NULL,
                    gender_lab = NULL,
                    age_col = "Age", gender_col = "Gender", pop_col = "Population"){

  Population <- Age <- Gender <- NULL
  df0 <- dplyr::rename_(data, Age = age_col, Gender = gender_col, Population = pop_col)
  if(is.null(m))
    m <- max(df0$Population)
  if(is.null(gender_lab))
    gender_lab <- unique(df0$Gender)
  if(length(gender_lab)>2)
    stop("Maximum of two gender labels permitted.")
  g <- ggplot2::ggplot(data = df0,
                       mapping = ggplot2::aes(x = Age, y = Population, fill = Gender)) +
    ggplot2::geom_bar(data = dplyr::filter(df0, Gender == gender_lab[1]),
                      stat = "identity") +
    ggplot2::geom_bar(data = dplyr::filter(df0, Gender == gender_lab[2]),
                      stat = "identity",
                      position = "identity",
                      mapping = ggplot2::aes(y = -Population)) +
    #set limits based -m and m
    ggplot2::scale_y_continuous(labels = abs, limits = c(-m, m)) +
    ggplot2::coord_flip() +
    ggplot2::labs(x = "Age", y = "Populaton", fill = "Gender")
  return(g)
}

