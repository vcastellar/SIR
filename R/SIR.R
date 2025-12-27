#' SIR model with cumulative infections and incidence output
#' @name sir_c
#' @description
#' Defines the right-hand side (derivatives) of an SIR compartmental model
#' extended with an auxiliary cumulative-infections state \code{C}. This function
#' is intended to be passed as the \code{func} argument to \code{deSolve::ode()}.
#'
#' @details
#' ## State variables
#' The state vector \code{state} must contain (at least) the following named
#' components:
#' \describe{
#'   \item{S}{Number of susceptible individuals at time \eqn{t}.}
#'   \item{I}{Number of infectious (active) individuals at time \eqn{t}.}
#'   \item{R}{Number of removed/recovered individuals at time \eqn{t}.}
#'   \item{C}{Cumulative number of infections up to time \eqn{t}.}
#' }
#'
#' ## Parameters
#' The parameter vector \code{parms} must contain:
#' \describe{
#'   \item{beta}{Transmission rate (per day).}
#'   \item{gamma}{Recovery/removal rate (per day).}
#' }
#'
#' ## Model equations
#' The total population is computed as \eqn{N = S + I + R}. The incidence
#' (new infections per day) is:
#' \deqn{\lambda(t) = \beta \frac{S(t) I(t)}{N}}
#'
#' and the ODE system is:
#' \deqn{
#' \begin{aligned}
#' \frac{dS}{dt} &= -\lambda(t) \\
#' \frac{dI}{dt} &= \lambda(t) - \gamma I(t) \\
#' \frac{dR}{dt} &= \gamma I(t) \\
#' \frac{dC}{dt} &= \lambda(t) \\
#' \end{aligned}
#' }
#' The function returns a list with the derivatives and an additional named output
#' \code{incidence} equal to \eqn{\lambda(t)}, which will appear as an extra column
#' in the output of \code{deSolve::ode()}.
#'
#' @param time Numeric scalar. Current time \eqn{t} at which derivatives are
#'   evaluated (provided by \code{deSolve::ode()}).
#' @param state Named numeric vector of current state values. Must include
#'   \code{S}, \code{I}, \code{R}, and \code{C}.
#' @param parms Named numeric vector of model parameters. Must include
#'   \code{beta} and \code{gamma}.
#'
#' @return A list with two elements:
#' \describe{
#'   \item{1}{A numeric vector of derivatives \code{c(dS, dI, dR, dC)} (in this order).}
#'   \item{incidence}{A numeric scalar with the incidence \eqn{\lambda(t)} at \code{time}.}
#' }
#'
#' @examples
#' library(deSolve)
#'
#' parms <- c(beta = 0.30, gamma = 0.10)
#' init <- c(S = 999, I = 1, R = 0, C = 1)
#' times <- 0:160
#'
#' out <- deSolve::ode(y = init, times = times, func = sir_c, parms = parms, method = "lsoda")
#' head(out)
#'
#' # Plot infectious prevalence I(t)
#' plot(out[, "time"], out[, "I"], type = "l", xlab = "Days", ylab = "I(t)")
#'
#' # Plot incidence lambda(t)
#' plot(out[, "time"], out[, "incidence"], type = "l", xlab = "Days", ylab = "Incidence")
#'
#' @export
sir <- function(time, state, parms) {
  with(as.list(c(state, parms)), {
    N <- S + I + R
    lambda <- beta * S * I / N          # incidencia (nuevas infecciones / día)
    dS <- -lambda
    dI <-  lambda - gamma * I
    dR <-  gamma * I
    dC <-  lambda                        # acumulado de infecciones
    list(c(dS, dI, dR, dC), incidence = lambda)
  })
}


#' SIRS model with cumulative infections and incidence output
#' @export
sirs <- function(time, state, parms) {
  with(as.list(c(state, parms)), {
    N <- S + I + R
    lambda <- beta * S * I / N

    dS <- -lambda + omega * R
    dI <-  lambda - gamma * I
    dR <-  gamma * I - omega * R
    dC <-  lambda

    list(c(dS, dI, dR, dC), incidence = lambda)
  })
}
