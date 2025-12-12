#' Sistema de ecuaciones SIR
#'
#' Define el sistema de ecuaciones diferenciales ordinarias para el modelo SIR
#' clásico.
#'
#' @param time Tiempo actual de la simulación.
#' @param state Vector con los valores de `S`, `I` y `R`.
#' @param parameters Vector con los parámetros `beta` y `gamma`.
#' @param N Tamaño total de la población.
#' @return Lista con las derivadas de `S`, `I` y `R` para usar con
#'   [deSolve::ode()].
#' @examples
#' SIR(time = 0, state = c(S = 999, I = 1, R = 0), parameters = c(beta = 0.4, gamma = 0.1), N = 1000)
SIR <- function(time, state, parameters, N) {
  par <- as.list(c(state, parameters))
  with(par, {
    dS <- -beta / N * I * S
    dI <- beta / N * I * S - gamma * I
    dR <- gamma * I
    list(c(dS, dI, dR))
  })
}
