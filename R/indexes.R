indexes <- function(base, model, measure = c('gdf', 'brwz', 'med', 'sbm')) {

  slacks <- model %>%
    as.data.frame() %>%
    dplyr::select(dplyr::starts_with('slack'))

  slacks_x <- slacks %>%
    dplyr::select(dplyr::starts_with('slack_x'))

  slacks_y <- slacks %>%
    dplyr::select(dplyr::starts_with('slack_y'))

  m <- ncol(slacks_x)
  s <- ncol(slacks_y)

  inputs <- base[1:m]
  outputs <- base[(m + 1):ncol(base)]

  targets_x <- inputs - slacks_x
  targets_y <- outputs + slacks_y

  ratios_x <- targets_x / inputs
  ratios_y <- outputs / targets_y

  brwz <- rowMeans(ratios_x) * rowMeans(ratios_y)

  gdf <- apply(ratios_x, 1, DescTools::Gmean) /
    apply(ratios_y ^ -1, 1, DescTools::Gmean)

  med <- 1 - (rowSums(slacks_x / inputs) +
                rowSums(slacks_y / targets_y)) / ncol(base)

  sbm <- rowMeans(ratios_x) / rowMeans(ratios_y ^ -1)
}
