#' Ajustar el modelo SIR por mínimos cuadrados
#'
#' Estima los parámetros `beta` y `gamma` del modelo SIR ajustando la serie de
#' infectados observados mediante mínimos cuadrados y devuelve la trayectoria
#' proyectada.
#'
#' @param infected Vector con los casos infectados observados.
#' @param population Tamaño total de la población susceptible.
#' @param horizon_days Días adicionales a proyectar más allá de los datos
#'   observados.
#' @return Una lista con `parameters` (valores estimados de `beta` y `gamma`) y
#'   `fit` (data.frame con la trayectoria simulada de `S`, `I` y `R`).
#' @examples
#' resultado <- fit_sir_model(c(100, 120, 150), population = 1000000, horizon_days = 30)
#' resultado$parameters
#' @importFrom deSolve ode
#' @importFrom stats optim
fit_sir_model <- function(infected, population, horizon_days) {
  Day <- seq_along(infected)
  init <- c(S = population - infected[1], I = infected[1], R = 0)

  rss <- function(parameters) {
    names(parameters) <- c("beta", "gamma")
    out <- ode(y = init, times = Day, func = function(time, state, parameters) {
      SIR(time, state, parameters, population)
    }, parms = parameters)
    fit <- out[, 3]
    mean((infected - fit) ^ 2)
  }

  opt <- optim(c(0.5, 0.5), rss, method = "Nelder-Mead", control = list(maxit = 500, factr = 1e4))
  opt_par <- setNames(opt$par, c("beta", "gamma"))

  t <- 1:(horizon_days + length(infected))
  fit <- data.frame(ode(y = init, times = t, func = function(time, state, parameters) {
    SIR(time, state, parameters, population)
  }, parms = opt_par))

  list(parameters = opt_par, fit = fit)
}
