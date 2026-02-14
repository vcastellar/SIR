
#-------------------------------------------------------------------------------
# modelo SEIRS
#-------------------------------------------------------------------------------
#' @include constructor.R
#' @keywords internal
#' @noRd
seirs_rhs <- function(time, state, parms) {
  with(as.list(c(state, parms)), {

    N <- S + E + I + R

    lambda <- beta * S * I / N

    dS <- -lambda + omega * R
    dE <-  lambda - sigma * E
    dI <-  sigma * E - gamma * I
    dR <-  gamma * I - omega * R

    list(
      c(dS, dE, dI, dR), incidence = sigma * E
    )
  })
}


#' SEIRS epidemic model with latent period and waning immunity
#'
#' @name SEIRS_MODEL
#' @description
#' An \code{epi_model} object representing a deterministic **SEIRS**
#' (Susceptible–Exposed–Infectious–Recovered–Susceptible) compartmental epidemic
#' model with a latent (exposed) period and waning immunity.
#'
#' The model describes the spread of an infection in a closed population where
#' susceptible individuals become infected at rate \eqn{\lambda(t)} and enter
#' the exposed compartment \code{E}. Exposed individuals progress to the
#' infectious compartment at rate \code{sigma}. Infectious individuals recover
#' at rate \code{gamma}, and recovered individuals lose immunity at rate
#' \code{omega}, returning to the susceptible compartment.
#'
#' @details
#' ## State variables
#' The model is defined in terms of the following state variables:
#' \describe{
#'   \item{S(t)}{Number of susceptible individuals at time \eqn{t}.}
#'   \item{E(t)}{Number of exposed (infected but not yet infectious) individuals at time \eqn{t}.}
#'   \item{I(t)}{Number of infectious (actively infected) individuals at time \eqn{t}.}
#'   \item{R(t)}{Number of recovered individuals with temporary immunity at time \eqn{t}.}
#' }
#'
#' The total population size is conserved:
#' \deqn{N = S(t) + E(t) + I(t) + R(t).}
#'
#' ## Model variables
#' The SEIRS model declares the following variables:
#' \describe{
#'   \item{\code{"S"}}{Susceptible population size.}
#'   \item{\code{"E"}}{Exposed (latent) population size.}
#'   \item{\code{"I"}}{Infectious population size.}
#'   \item{\code{"R"}}{Recovered (temporarily immune) population size.}
#'   \item{\code{"incidence"}}{Rate of progression from \code{E} to \code{I},
#'     \eqn{\sigma E(t)}, representing the instantaneous incidence of new
#'     infectious cases returned by the model's right-hand side.}
#' }
#'
#' All declared variables may be used as observables in generic utilities
#' and summary methods built around \code{epi_model} objects.
#'
#' ## Parameters
#' The SEIRS model depends on the following parameters:
#' \describe{
#'   \item{beta}{Transmission rate (per day).}
#'   \item{sigma}{Rate of progression from exposed to infectious (per day);
#'     \eqn{1/\sigma} is the mean latent period.}
#'   \item{gamma}{Recovery/removal rate from infectious to recovered (per day);
#'     \eqn{1/\gamma} is the mean infectious period.}
#'   \item{omega}{Rate of waning immunity from \code{R} back to \code{S} (per day);
#'     \eqn{1/\omega} is the mean immunity duration.}
#' }
#'
#' ## Model equations
#' New infections occur at rate
#' \deqn{\lambda(t) = \beta \frac{S(t)\, I(t)}{N}.}
#'
#' Case incidence (entries into \code{I}) occurs at rate
#' \deqn{\text{incidence}(t) = \sigma E(t).}
#'
#' The system of ordinary differential equations is:
#' \deqn{
#' \begin{aligned}
#' \frac{dS}{dt} &= -\lambda(t) + \omega R(t), \\
#' \frac{dE}{dt} &= \lambda(t) - \sigma E(t), \\
#' \frac{dI}{dt} &= \sigma E(t) - \gamma I(t), \\
#' \frac{dR}{dt} &= \gamma I(t) - \omega R(t).
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
#' ## Simulate a SEIRS epidemic
#' sim <- simulate_epi(
#'   model = SEIRS_MODEL,
#'   times = 0:300,
#'   parms = c(beta = 0.3, sigma = 0.2, gamma = 0.14, omega = 0.01),
#'   init  = c(S = 1e6, E = 0, I = 20, R = 0)
#' )
#'
#' plot(sim)
#'
#' ## Plot observed incidence
#' plot(sim, what = "incidence")
#'
#'
#' @seealso
#' \code{\link{simulate_epi}},
#' \code{\link{epi_model}}
#'
#' @export
SEIRS_MODEL <- epi_model(
  name = "SEIRS",
  rhs = seirs_rhs,
  par_names = c("beta", "sigma", "gamma", "omega"),
  lower = c(beta = 1e-8, sigma = 1e-8, gamma = 1e-8, omega = 1e-8),
  upper = c(beta = 5,    sigma = 2,    gamma = 2,    omega = 1),
  defaults = c(beta = 0.3, sigma = 0.2, gamma = 0.14, omega = 0.01),
  init = c("S" = 1e6, "I" = 10, "R" = 0, "E" = 0),
  states = c("S", "E", "I", "R"),
  flows = c("incidence")
)
