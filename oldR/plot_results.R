#
# plot_results <- function(fit_sir, nfected_dates, infected, population) {
#   fechas_modelo <- seq.Date(from = infected_dates[1], length.out = nrow(fit_sir$fit), by = "day")
#   fit_sir$fit$time <- fechas_modelo
#
#   plot(fit_sir$fit$time, fit_sir$fit$S, type = "l", lwd = 2, ylim = c(0, population))
#   lines(fit_sir$fit$time, fit_sir$fit$I, lwd = 2, col = "red")
#   lines(fit_sir$fit$time, fit_sir$fit$R, lwd = 2, col = "green")
#   points(infected_dates, infected, pch = 19)
#   grid()
#   legend("left", c("Susceptibles", "Infectados", "Recobrados"), lty = 1, lwd = 2, col = 1:3, inset = 0.05)
# }
