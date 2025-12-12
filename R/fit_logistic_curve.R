#' Ajustar una curva logística
#'
#' Ajusta una curva logística simple a los datos observados de casos
#' infectados mediante mínimos cuadrados.
#'
#' @param infected Vector con los casos infectados observados.
#' @param dates Vector de fechas correspondientes a `infected`.
#' @return Una lista con `parameters` (parámetros optimizados) y `fit`
#'   (valores ajustados para cada fecha).
#' @examples
#' ajuste <- fit_logistic_curve(c(100, 120, 150), dates = as.Date("2020-03-01") + 0:2)
#' ajuste$parameters
#' @importFrom stats optim
fit_logistic_curve <- function(infected, dates) {
  Day <- seq_along(infected)
  logistic <- function(parms, x) {
    parms[3] / (1 + exp(-(x - parms[2]) / parms[1]))
  }

  rss <- function(parms) {
    out <- logistic(parms, x = Day)
    sum((infected - out) ^ 2)
  }

  opt <- optim(c(1, 1, max(infected) * 2), rss, method = "L-BFGS-B")
  fit <- logistic(opt$par, 1:length(dates))
  list(parameters = opt$par, fit = fit)
}
