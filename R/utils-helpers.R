#' Calculate ranges for each row in a table
#'
#' This function is useful for preparing the objective functions of the RAM and
#' BAM models. It calculates the range (\code{max{x_ij, j = 1, ..., n} -
#' min{x_ij, j = 1, ..., n}}), lower-sided range (\code{x_io -
#' min{x_ij, j = 1, ..., n}}), and upper-sided range
#' (\code{max{y_rj, j = 1, ..., n} - y_ro}), where \code{i = 1, ..., m} and
#' \code{r = 1, ..., s} are the input and output indices respectively, and
#' \code{n} is the number of DMUs.
#' @param x A numeric matrix or data frame.
#' @param range_type String vector. Choices are \code{"standard"},
#' \code{"lower_sided"}, and \code{"upper_sided"}.
#'
#' @return The selected \code{range_type} for each row of the matrix.
#'
#' @examples
#' x <- matrix(c(2, 12, 2, 8, 5, 5, 10, 4, 10, 6, 3, 13),
#'   ncol = 6,
#'   byrow = TRUE
#' )
#' calc_range(x, "standard")
#' calc_range(x, "lower_sided")
#' calc_range(x, "upper_sided")
calc_range <- function(x, range_type = c("standard", "lower_sided",
                                         "upper_sided")) {

  if (range_type == "standard") {
    range_value <- apply(x, 1, function(x) {
      max(x) - min(x)
    })
  } else if (range_type == "lower_sided") {
    range_value <- x - apply(x, 1, min)
  } else if (range_type == "upper_sided") {
    range_value <- apply(x, 1, max) - x
  }

  return(range_value)
}
