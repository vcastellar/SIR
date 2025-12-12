#' Calcular métricas de crecimiento
#'
#' Calcula el porcentaje de incremento diario y los días estimados para
#' duplicar los casos usando medias móviles.
#'
#' @param infected Vector numérico con casos infectados acumulados.
#' @return Una lista con `inc_pct` (porcentaje de incremento diario) y
#'   `doubling_days` (estimación de días para duplicar casos).
#' @examples
#' growth <- compute_growth_metrics(c(100, 120, 150, 200, 260))
#' growth$doubling_days
#' @importFrom zoo rollapply
compute_growth_metrics <- function(infected) {
  inc_pct <- rollapply(infected, width = 2, FUN = function(x) (x[2] - x[1]) / x[1], align = "left")
  alpha <- rollapply(inc_pct, width = 4, FUN = mean)
  doubling_days <- log(2) / log(1 + alpha)
  doubling_days[is.infinite(doubling_days)] <- NA
  list(inc_pct = inc_pct, doubling_days = doubling_days)
}
