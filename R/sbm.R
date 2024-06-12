sbm <- function(base, noutput, fixed = NULL, rts = "vrs", orientation = "in_out",
                bound = NULL,
                whichDMUs = NULL, print.status = FALSE) {
  s <- noutput
  m <- ncol(base) - s
  n <- nrow(base)

  base.X <- t(base[, 1:m])
  base.Y <- t(base[, (m + 1):(s + m)])

  f.obj <- prep_objective_function(base.X, base.Y, orientation, fixed,
    model = "sbm", m, s, n
  )

  f.rhs <- prep_rhs(base.X, base.Y, orientation, rts, model = "sbm", m, s, n)

  f.con <- prep_constraints_coefficients(base.X, base.Y, orientation, rts,
    model = "sbm",
    m, s, n
  )

  f.dir <- prep_directions(orientation, rts, model = "sbm", m, s)

  linear_programmes <- mapply(
    function(x, y, z) {
      lpSolve::lp(
        ifelse(orientation == "out", "max", "min"),
        objective.in = x,
        const.mat = z,
        const.dir = f.dir,
        const.rhs = y
      )
    },
    x = f.obj,
    y = f.rhs,
    z = f.con,
    SIMPLIFY = FALSE
  )

  efficiency_scores <- sapply(
    linear_programmes,
    function(x) {
      ifelse(orientation == "out", 1 / x$objval, x$objval)
    }
  )

  tau_values <- sapply(
    linear_programmes,
    function(x) {
      x$solution[1]
    }
  )

  lambdas_and_slacks <- sapply(
    linear_programmes,
    function(x) {
      x$solution[-1] # Exclude tau
    }
  )

  lambdas_and_slacks <- sweep(lambdas_and_slacks, 2, tau_values, "/")

  optimal_solution <- cbind(efficiency_scores, t(lambdas_and_slacks))

  colnames(optimal_solution) <- c(
    "eff",
    paste("lambda_", 1:n, sep = ""),
    paste("slack_x", 1:m, sep = ""),
    paste("slack_y", 1:s, sep = "")
  )

  return(optimal_solution)
}
