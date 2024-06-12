gem <- function(base, noutput, fixed = NULL, rts = "vrs", orientation = "in_out",
                bound = NULL,
                model = c("additive", "bam", "lov-past", "mip", "ram"),
                whichDMUs = NULL, print.status = FALSE) {
  model <- tolower(model)
  s <- noutput
  m <- ncol(base) - s
  n <- nrow(base)

  base.X <- t(base[, 1:m])
  base.Y <- t(base[, (m + 1):(s + m)])

  f.obj <- prep_objective_function(
    base.X, base.Y, orientation, fixed,
    model
  )

  f.rhs <- prep_rhs(base.X, base.Y, orientation, rts, model, m, s, n)

  f.con <- prep_constraints_coefficients(
    base.X, base.Y, orientation, rts, model,
    m, s, n
  )

  f.dir <- prep_directions(orientation, rts, model, m, s)

  linear_programmes <- mapply(
    function(x, y) {
      lpSolve::lp(
        "max",
        objective.in = x,
        const.mat = f.con,
        const.dir = f.dir,
        const.rhs = y
      )
    },
    x = f.obj,
    y = f.rhs,
    SIMPLIFY = FALSE
  )

  inefficiency_scores <- sapply(
    linear_programmes,
    function(x) {
      x$objval
    }
  )

  lambdas_and_slacks <- sapply(
    linear_programmes,
    function(x) {
      x$solution
    }
  )

  optimal_solution <- cbind(inefficiency_scores, t(lambdas_and_slacks))

  colnames(optimal_solution) <- c(
    "eff",
    paste("lambda_", 1:n, sep = ""),
    paste("slack_x", 1:m, sep = ""),
    paste("slack_y", 1:s, sep = "")
  )

  return(optimal_solution)
}
