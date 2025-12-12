#' Graficar los resultados del modelo
#'
#' Genera las gráficas de las trayectorias SIR ajustadas y la curva
#' logística comparada con los datos observados.
#'
#' @param fit_sir Lista devuelta por [fit_sir_model()].
#' @param fit_logistic Lista devuelta por [fit_logistic_curve()].
#' @param infected_dates Fechas correspondientes a los casos observados.
#' @param infected Vector con los casos infectados observados.
#' @param population Tamaño total de la población susceptible.
#' @return Invisiblemente, la lista `fit_sir` con la columna de fechas añadida.
#' @examples
#' # Se necesitan resultados previos de ajuste para generar las gráficas
#' @importFrom graphics plot points lines legend grid
plot_results <- function(fit_sir, fit_logistic, infected_dates, infected, population) {
  fechas_modelo <- seq.Date(from = infected_dates[1], length.out = nrow(fit_sir$fit), by = "day")
  fit_sir$fit$time <- fechas_modelo

  plot(fit_sir$fit$time, fit_sir$fit$S, type = "l", lwd = 2, ylim = c(0, population))
  lines(fit_sir$fit$time, fit_sir$fit$I, lwd = 2, col = "red")
  lines(fit_sir$fit$time, fit_sir$fit$R, lwd = 2, col = "green")
  points(infected_dates, infected, pch = 19)
  grid()
  legend("left", c("Susceptibles", "Infectados", "Recobrados"), lty = 1, lwd = 2, col = 1:3, inset = 0.05)

  plot(fit_sir$fit$time, fit_sir$fit$I, type = "l", lwd = 2, col = "red", pch = 19)
  points(infected_dates, infected, pch = 19)
  grid()

  plot(infected_dates, fit_logistic$fit, type = "l", lwd = 2, col = "blue",
       main = "Curva logística ajustada", ylab = "Casos", xlab = "Fecha")
  points(infected_dates, infected, pch = 19, col = "red")
  grid()

  invisible(fit_sir)
}
