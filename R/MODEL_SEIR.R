#-------------------------------------------------------------------------------
# modelo SEIR
#-------------------------------------------------------------------------------
#' @include constructor.R
#' @keywords internal
#' @noRd
seir_rhs <- function(time, state, parms) {
  with(as.list(c(state, parms)), {

    N <- S + E + I + R

    lambda <- beta * S * I / N

    dS <- -lambda
    dE <-  lambda - sigma * E
    dI <-  sigma * E - gamma * I
    dR <-  gamma * I

    list(c(dS, dE, dI, dR), incidence = sigma * E)
  })
}

#' SEIR epidemic model with latent (exposed) period
#'
#' @name SEIR_MODEL
#' @description
#' An \code{epi_model} object representing a deterministic **SEIR**
#' (Susceptible–Exposed–Infectious–Recovered) compartmental epidemic model
#' with a latent (exposed) period.
#'
#' The model describes the spread of an infection in a closed population where
#' susceptible individuals become infected at rate \eqn{\lambda(t)} and enter
#' the exposed compartment \code{E}. Exposed individuals progress to the
#' infectious compartment at rate \code{sigma} and subsequently recover with
#' permanent immunity.
#'
#' @details
#' ## State variables
#' The model is defined in terms of the following state variables:
#' \describe{
#'   \item{S(t)}{Number of susceptible individuals at time \eqn{t}.}
#'   \item{E(t)}{Number of exposed (infected but not yet infectious) individuals at time \eqn{t}.}
#'   \item{I(t)}{Number of infectious (actively infected) individuals at time \eqn{t}.}
#'   \item{R(t)}{Number of recovered (immune) individuals at time \eqn{t}.}
#' }
#'
#' The total population size is conserved:
#' \deqn{N = S(t) + E(t) + I(t) + R(t).}
#'
#' ## Model outputs
#' The SEIR model declares the following outputs:
#' \describe{
#'   \item{\code{"S"}}{Susceptible population size.}
#'   \item{\code{"E"}}{Exposed (latent) population size.}
#'   \item{\code{"I"}}{Infectious population size.}
#'   \item{\code{"R"}}{Recovered (immune) population size.}
#'   \item{\code{"incidence"}}{Rate of progression from \code{E} to \code{I},
#'     \eqn{\sigma E(t)}, representing the instantaneous incidence of new
#'     infectious cases returned by the model's right-hand side.}
#' }
#'
#' All declared outputs may be used as observables in generic utilities such as
#' \code{\link{fit_epi_model}} via the \code{target} argument.
#'
#' ## Parameters
#' The SEIR model depends on the following parameters:
#' \describe{
#'   \item{beta}{Transmission rate (per day).}
#'   \item{sigma}{Rate of progression from exposed to infectious (per day);
#'     \eqn{1/\sigma} is the mean latent period.}
#'   \item{gamma}{Recovery/removal rate from infectious to recovered (per day);
#'     \eqn{1/\gamma} is the mean infectious period.}
#' }
#'
#' ## Model equations
#' New infections occur at rate
#' \deqn{\lambda(t) = \beta \frac{S(t)\, I(t)}{N}.}
#'
#' Progression from exposed to infectious occurs at rate
#' \deqn{\text{incidence}(t) = \sigma E(t).}
#'
#' The system of ordinary differential equations is:
#' \deqn{
#' \begin{aligned}
#' \frac{dS}{dt} &= -\lambda(t), \\
#' \frac{dE}{dt} &= \lambda(t) - \sigma E(t), \\
#' \frac{dI}{dt} &= \sigma E(t) - \gamma I(t), \\
#' \frac{dR}{dt} &= \gamma I(t).
#' \end{aligned}
#' }
#'
#' ## Usage
#' This predefined model object is intended to be used with generic utilities
#' such as \code{\link{simulate_epi}}, \code{\link{fit_epi_model}}, and
#' \code{\link{predict.fit_epi_model}} that operate on \code{epi_model} objects.
#'
#' @format
#' An object of class \code{"epi_model"}.
#'
#' @examples
#' ## Simulate a SEIR epidemic
#' sim <- simulate_epi(
#'   model = SEIR_MODEL,
#'   times = 0:200,
#'   parms = c(beta = 0.3, sigma = 0.2, gamma = 0.14),
#'   init  = c(S = 1e6, E = 5, I = 10, R = 0),
#'   obs   = "poisson"
#' )
#'
#' plot(sim)
#'
#' ## Plot observed incidence
#' plot(sim, what = "incidence")
#'
#' ## Fit the model to observed incidence
#' fit_inc <- fit_epi_model(
#'   x = sim$incidence$inc,
#'   model = SEIR_MODEL,
#'   init = SEIR_MODEL$init,
#'   target = "incidence"
#' )
#'
#' fit_inc
#'
#' @seealso
#' \code{\link{simulate_epi}},
#' \code{\link{fit_epi_model}},
#' \code{\link{new_epi_model}}
#'
#' @export
SEIR_MODEL <- epi_model(
  name = "SEIR",
  rhs = seir_rhs,
  par_names = c("beta", "sigma", "gamma"),
  lower = c(beta = 1e-8, sigma = 1e-8, gamma = 1e-8),
  upper = c(beta = 5,    sigma = 2,    gamma = 2),
  defaults = c(beta = 0.3, sigma = 0.2, gamma = 0.14),
  init = c("S" = 1e6, "E" = 5, "I" = 10, "R" = 0),
  states = c("S", "E", "I", "R"),
  flows  = c("incidence"),
  roles = list(
    susceptible = "S",
    exposed     = "E",
    infectious  = "I",
    recovered   = "R",
    incidence   = "incidence"
  )
)
