
#' Right-hand side for SIR model with vital dynamics
#'
#' @description
#' Internal function defining the system of differential equations
#' for \code{SIR_VITAL_MODEL}.
#'
#' @param time Time variable (required by ODE solver).
#' @param state Named numeric vector of state variables.
#' @param parms Named numeric vector of model parameters.
#' @param ... Additional arguments passed by the solver.
#'
#' @return
#' A list containing derivatives and declared flows.
#'
#' @keywords internal
sir_vital_rhs <- function(time, state, parms) {
  with(as.list(c(state, parms)), {

    N <- S + I + R
    infection <- beta * S * I / N

    births  <- mu * N
    recovery <- gamma * I

    death_S <- mu * S
    death_I <- mu * I
    death_R <- mu * R

    dS <- births - infection - death_S
    dI <- infection - recovery - death_I
    dR <- recovery - death_R

    list(
      c(dS, dI, dR),
      births    = births,
      infection = infection,
      recovery  = recovery,
      death_S   = death_S,
      death_I   = death_I,
      death_R   = death_R
    )
  })
}

#-------------------------------------------------------------------------------
# SIR model with vital dynamics and constant population
#-------------------------------------------------------------------------------
#' SIR epidemic model with vital dynamics and constant population
#'
#' @name SIR_VITAL_MODEL
#' @description
#' An \code{epi_model} object representing a deterministic **SIR**
#' (Susceptible–Infectious–Recovered) epidemic model with demographic turnover
#' and **strictly constant population size**.
#'
#' Natural deaths occur in all compartments at a constant per-capita rate, and
#' births exactly balance deaths at each instant, ensuring that the total
#' population size remains constant over time.
#'
#' @details
#' ## State variables
#' The model is defined in terms of the following state variables:
#' \describe{
#'   \item{S(t)}{Number of susceptible individuals at time \eqn{t}.}
#'   \item{I(t)}{Number of infectious individuals at time \eqn{t}.}
#'   \item{R(t)}{Number of recovered individuals at time \eqn{t}.}
#' }
#'
#' The total population size is conserved:
#' \deqn{N = S(t) + I(t) + R(t).}
#'
#' ## Flows
#' The following epidemiological flows are declared:
#' \describe{
#'   \item{births}{Recruitment of new susceptible individuals exactly balancing
#'     natural deaths.}
#'   \item{infection}{Incidence of new infections.}
#'   \item{recovery}{Recovery of infectious individuals.}
#'   \item{death_S}{Natural deaths among susceptibles.}
#'   \item{death_I}{Natural deaths among infectious.}
#'   \item{death_R}{Natural deaths among recovered.}
#' }
#'
#' ## Parameters
#' The model depends on the following parameters:
#' \describe{
#'   \item{beta}{Transmission rate.}
#'   \item{gamma}{Recovery rate.}
#'   \item{mu}{Natural mortality rate (per capita).}
#' }
#'
#' ## Model equations
#' New infections occur at rate
#' \deqn{\lambda(t) = \beta \frac{S(t)\, I(t)}{N}.}
#'
#' The system of ordinary differential equations is:
#' \deqn{
#' \begin{aligned}
#' \frac{dS}{dt} &= \mu N - \lambda(t) - \mu S, \\
#' \frac{dI}{dt} &= \lambda(t) - \gamma I - \mu I, \\
#' \frac{dR}{dt} &= \gamma I - \mu R.
#' \end{aligned}
#' }
#'
#' @format
#' An object of class \code{"epi_model"}.
#'
#' @examples
#' ## Simulate a SIR epidemic with constant population size
#' sim <- simulate_epi(
#'   model = SIR_VITAL_MODEL,
#'   times = 0:300,
#'   parms = c(
#'     beta  = 0.4,
#'     gamma = 0.1,
#'     mu    = 0.01
#'   ),
#'   init = c(S = 0.99, I = 0.01, R = 0)
#' )
#'
#' plot(sim)
#'
#' ## Plot incidence
#' plot(sim, what = "infection")
#'
#' @seealso
#' \code{\link{simulate_epi}},
#' \code{\link{epi_model}},
#' \code{\link{SIR_MODEL}}
#'
#' @export


SIR_VITAL_MODEL <- epi_model(
  name = "SIR with vital dynamics",
  rhs  = sir_vital_rhs,
  par_names   = c("beta", "gamma", "mu"),

  states = c("S", "I", "R"),
  flows  = c(
    "births",
    "infection",
    "recovery",
    "death_S",
    "death_I",
    "death_R"
  )
)
