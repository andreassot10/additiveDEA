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

prep_constraints_coefficients <- function(x, y, orientation, rts, model,
                                          m = NULL, s = NULL, n = NULL) {
  stopifnot(
    'Orientation must be "in", "out" or "in_out"' =
      orientation %in% c("in", "out", "in_out")
  )

  if (is.null(m)) {
    m <- nrow(x)
  }
  if (is.null(s)) {
    s <- nrow(y)
  }
  if (is.null(n)) {
    n <- ncol(x)
  }

  constraints_x <- cbind(x, diag(m), matrix(0, m, s))
  constraints_y <- cbind(y, matrix(0, s, m), -diag(s))
  constraints_x_y <- rbind(constraints_x, constraints_y)

  if (rts != "crs") {
    constraints_rts <- c(rep(1, n), rep(0, m + s))
  } else {
    constraints_rts <- NULL
  }

  constraint_coefficients <- rbind(constraints_x_y, constraints_rts)

  if (model == "sbm") {
    # In SBM, the constraints coefficients depend on x0 and y0. In any additive
    # model these coefficients are actually the RHS constants. We can thus
    # retrieve the RHS from any additive model and use them as constraint
    # coefficients for the SBM model
    variable_constraint_coefficients <- prep_rhs(x, y,
      orientation = "in_out",
      rts, model = "additive"
    )

    # Add coefficient for tau variable (denominator constraint) and correct sign
    # of coefficients of tau in the input, output and RTS (if applicable)
    # constraints
    variable_constraint_coefficients <- lapply(
      variable_constraint_coefficients,
      function(x) {
        # The constraints for the inputs, outputs and RTS are multiplied by tau,
        # and brought over to the LHS. They should all be negative then
        x <- -x

        # Denominator constraint
        c(1, x)
      }
    )

    # Output slack weights
    wt_y <- turner::matrix_to_blocks(
      1 / (s * zero_to_inf(y)),
      rep(1, n),
      byrow = FALSE
    )

    # Lambda coefficients (all zero) and slack coefficients (first constraint)
    first_constraint_coefficients <- lapply(
      wt_y,
      function(x) {
        c(rep(0, n), rep(0, m), x)
      }
    )

    # Add coefficient for tau variable
    constraint_coefficients <- purrr::map2(
      function(.x, .y) {
        # When the model is oriented, the tau variable is unnecessary. Instead
        # of adding/removing constraints based on the orientation, it is easier
        # to maintain the dimensions of the constraints table for "in_out" and
        # modify the constraint coefficients accordingly. Here, when the model
        # is input- or output-oriented, the first constraint just needs to state
        # that tau equals 1, thus all other coefficients should be zero.
        if (orientation != "in_out") {
          .x <- 0 * .x
        }
        aux <- rbind(.x, constraint_coefficients)
        cbind(.y, aux)
      },
      .x = first_constraint_coefficients,
      .y = variable_constraint_coefficients
    )
  }

  return(constraint_coefficients)
}

prep_directions <- function(orientation, rts, model, m = NULL, s = NULL) {
  directions <- rep("==", m + s)

  if (rts == "vrs") {
    directions <- c(directions, "==")
  } else if (rts == "drs") {
    directions <- c(directions, "<=")
  } else if (rts == "irs") {
    directions <- c(directions, ">=")
  }

  if (model == "sbm") {
    directions <- c("==", directions)
  }

  return(directions)
}

prep_objective_function <- function(x, y, orientation, fixed,
                                    model = c(
                                      "sbm", "additive",
                                      "bam", "lov-past",
                                      "mip", "ram"
                                    ),
                                    m = NULL, s = NULL, n = NULL) {
  if (is.null(m)) {
    m <- nrow(x)
  }
  if (is.null(s)) {
    s <- nrow(y)
  }
  if (is.null(n)) {
    n <- ncol(x)
  }

  if (model == "sbm") {
    if (orientation != "out") {
      wt_x <- -1 / (m * zero_to_inf(x))
      wt_y <- 0 * y
    } else {
      wt_x <- 0 * x
      wt_y <- 1 / (s * zero_to_inf(y))
    }
  } else if (model == "additive") {
    wt_x <- matrix(1, m, n)
    wt_y <- matrix(1, s, n)
  } else if (model == "bam") {
    lower_range_x <- calc_range(x, range_type = "lower_sided")
    upper_range_y <- calc_range(y, range_type = "upper_sided")

    wt_x <- matrix(1 / ((m + s) * zero_to_inf(lower_range_x)), m, n)
    wt_y <- matrix(1 / ((m + s) * zero_to_inf(upper_range_y)), s, n)
  } else if (model == "lov-past") {
    sd_x <- apply(x, 1, sd)
    sd_y <- apply(y, 1, sd)

    wt_x <- 1 / zero_to_inf(sd_x)
    wt_y <- 1 / zero_to_inf(sd_y)
  } else if (model == "mip") {
    wt_x <- 1 / zero_to_inf(x)
    wt_y <- 1 / zero_to_inf(y)
  } else if (model == "ram") {
    range_x <- calc_range(x, range_type = "standard")
    range_y <- calc_range(y, range_type = "standard")

    wt_x <- matrix(1 / ((m + s) * zero_to_inf(range_x)), m, n)
    wt_y <- matrix(1 / ((m + s) * zero_to_inf(range_y)), s, n)
  }

  if (orientation == "in") {
    wt_y <- wt_y * 0
  } else if (orientation == "out") {
    wt_x <- wt_x * 0
  }

  wt <- rbind(wt_x, wt_y)
  wt[fixed, ] <- 0

  objective_function_weights_list <- turner::matrix_to_blocks(
    wt,
    rep(1, ncol(wt)),
    byrow = FALSE
  )

  # Add lambda coefficients
  objective_function_weights_list <- lapply(
    objective_function_weights_list,
    function(x) {
      c(rep(0, n), x)
    }
  )

  if (model == "sbm") {
    objective_function_weights_list <- lapply(
      objective_function_weights_list,
      function(x) {
        c(1, x)
      }
    )
  }

  return(objective_function_weights_list)
}

prep_rhs <- function(x, y, orientation, rts, model,
                     m = NULL, s = NULL, n = NULL) {
  if (is.null(m)) {
    m <- nrow(x)
  }
  if (is.null(s)) {
    s <- nrow(y)
  }
  if (is.null(n)) {
    n <- ncol(x)
  }

  if (model != "sbm") {
    rhs_list <- turner::matrix_to_blocks(
      rbind(x, y),
      rep(1, n),
      byrow = FALSE
    )
  } else {
    rhs_list <- lapply(
      1:n,
      function(z) {
        c(1, rep(0, m), rep(0, s))
      }
    )
  }

  if (rts != "crs") {
    rhs_list <- lapply(
      rhs_list,
      function(x) {
        if (model == "sbm") {
          c(x, 0)
        } else {
          c(x, 1)
        }
      }
    )
  }

  return(rhs_list)
}

#' Convert zero values to \code{Inf}
#'
#' This function is useful when slacks in the objective function are divided by
#' zeroes, as a zero in the denominator would result in an infeasible linear
#' program. Converting zeroes to \code{Inf} will result in \code{1 / Inf = 0},
#' thus making the linear program feasible.
#'
#' @param x A numeric vector, matrix or data frame.
#'
#' @return A modified version of the the original numeric vector, matrix or data
#' frame, though with any zeroes now converted to \code{Inf}.
#' @export
#'
#' @examples
zero_to_inf <- function(x) {
  x[x == 0] <- Inf
  return(x)
}

zero_to_penalty <- function(x, denominator = 10 ^ -6) {
  x[x == 0] <- 1 / denominator
  return(x)
}
