
#-------------------------------------------------------------------------------
# modelo SIR
#-------------------------------------------------------------------------------
#' @include constructor.R
#' @keywords internal
#' @noRd
sir_rhs <- function(time, state, parms) {
  with(as.list(c(state, parms)), {

    N <- S + I + R

    lambda <- beta * S * I / N

    dS <- -lambda
    dI <-  lambda - gamma * I
    dR <-  gamma * I


    list(
      c(dS, dI, dR),

      ## flujos instantáneos
      incidence = lambda,
      recovery  = gamma * I
    )
  })
}



#' SIR epidemic model
#'
#' @name SIR_MODEL
#' @description
#' An \code{epi_model} object representing a deterministic **SIR**
#' (Susceptible–Infectious–Recovered) compartmental epidemic model.
#'
#' The model describes the spread of an infection in a closed population where
#' susceptible individuals become infectious and subsequently recover with
#' permanent immunity.
#'
#' @details
#' ## State variables
#' The model is defined in terms of the following state variables:
#' \describe{
#'   \item{S(t)}{Number of susceptible individuals at time \eqn{t}.}
#'   \item{I(t)}{Number of infectious (actively infected) individuals at time \eqn{t}.}
#'   \item{R(t)}{Number of removed/recovered individuals at time \eqn{t}.}
#' }
#'
#' The total population size is conserved:
#' \deqn{N = S(t) + I(t) + R(t).}
#'
#' ## Model variables
#' The SIR model declares the following variables:
#' \describe{
#'   \item{\code{"S"}}{Susceptible population size.}
#'   \item{\code{"I"}}{Infectious population size.}
#'   \item{\code{"R"}}{Recovered (removed) population size.}
#'   \item{\code{"incidence"}}{Instantaneous rate of new infections
#'     \eqn{\lambda(t)} returned by the model's right-hand side.}
#'   \item{\code{"recovery"}}{Instantaneous recovery flow
#'     \eqn{\gamma I(t)} returned by the model's right-hand side.}
#' }
#'
#' All declared variables may be used as observables in generic utilities
#' and summary methods built around \code{epi_model} objects.
#'
#' ## Parameters
#' The SIR model depends on the following parameters:
#' \describe{
#'   \item{beta}{Transmission rate (per day).}
#'   \item{gamma}{Recovery/removal rate (per day).}
#' }
#'
#' ## Model equations
#' New infections occur at rate
#' \deqn{\lambda(t) = \beta \frac{S(t)\, I(t)}{N}.}
#'
#' The system of ordinary differential equations is:
#' \deqn{
#' \begin{aligned}
#' \frac{dS}{dt} &= -\lambda(t), \\
#' \frac{dI}{dt} &= \lambda(t) - \gamma I(t), \\
#' \frac{dR}{dt} &= \gamma I(t). \\
#' \end{aligned}
#' }
#'
#' ## Usage
#' This predefined model object is intended to be used with generic utilities
#' such as \code{\link{simulate_epi}}, \code{\link{plot.sim_epi}}, and
#' \code{\link{summary.sim_epi}} that operate on \code{epi_model} objects.
#'
#' @format
#' An object of class \code{"epi_model"}.
#'
#' @examples
#' ## Simulate a SIR epidemic
#' sim <- simulate_epi(
#'   model = SIR_MODEL,
#'   times = 0:200,
#'   parms = c(beta = 0.3, gamma = 0.1),
#'   init  = c(S = 1e6, I = 10, R = 0)
#' )
#'
#' plot(sim)
#'
#' ## Plot incidence
#' plot(sim, what = "incidence")
#'
#'
#' @seealso
#' \code{\link{simulate_epi}},
#' \code{\link{epi_model}}
#'
#' @export
SIR_MODEL <- epi_model(
  name = "SIR",
  rhs = sir_rhs,
  par_names = c("beta", "gamma"),

  lower = c(beta = 1e-8, gamma = 1e-8),
  upper = c(beta = 2,    gamma = 1),

  defaults = c(beta = 0.3, gamma = 0.1),

  init = c(S  = 1e6, I  = 10, R  = 0),

  states = c("S", "I", "R"),
  flows = c("incidence", "recovery")
)
