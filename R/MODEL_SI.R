
#-------------------------------------------------------------------------------
# modelo SI
#-------------------------------------------------------------------------------
#' @include constructor.R
#' @keywords internal
#' @noRd
si_rhs <- function(time, state, parms) {
  with(as.list(c(state, parms)), {
    N <- S + I
    lambda <- beta * S * I / N
    dS <- -lambda
    dI <-  lambda
    list(c(dS, dI), incidence = lambda)
  })
}


#' SI epidemic model
#'
#' @name SI_MODEL
#' @description
#' An \code{epi_model} object representing a deterministic **SI**
#' (Susceptible–Infectious) compartmental epidemic model.
#'
#' The model describes the spread of an infection in a closed population where
#' individuals move irreversibly from the susceptible compartment \code{S} to the
#' infectious compartment \code{I}. No recovery or removal process is included.
#'
#' @details
#' ## State variables
#' The model is defined in terms of the following state variables:
#' \describe{
#'   \item{S(t)}{Number of susceptible individuals at time \eqn{t}.}
#'   \item{I(t)}{Number of infectious individuals at time \eqn{t}.}
#' }
#'
#' The total population size is conserved:
#' \deqn{N = S(t) + I(t).}
#'
#' ## Model outputs
#' The SI model declares the following outputs:
#' \describe{
#'   \item{\code{"S"}}{Susceptible population size.}
#'   \item{\code{"I"}}{Infectious population size.}
#'   \item{\code{"incidence"}}{Instantaneous rate of new infections
#'     \eqn{\lambda(t)} returned by the model's right-hand side.}
#' }
#'
#' All declared outputs may be used as observables in generic utilities such as
#' \code{\link{fit_epi_model}} via the \code{target} argument.
#'
#' ## Parameters
#' The SI model depends on a single parameter:
#' \describe{
#'   \item{beta}{Transmission rate (per day).}
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
#' \frac{dI}{dt} &= \lambda(t).
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
#' ## Simulate an SI epidemic
#' sim <- simulate_epi(
#'   model = SI_MODEL,
#'   times = 0:100,
#'   parms = c(beta = 0.4),
#'   obs = "negbin",
#'   init = SI_MODEL$init
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
#' \code{\link{fit_epi_model}},
#' \code{\link{new_epi_model}}
#'
#' @export

SI_MODEL <- epi_model(
  name = "SI",
  rhs = si_rhs,
  state_names = c("S", "I"),
  par_names = c("beta"),
  states = c("S", "I"),
  flows = c("incidence"),
  defaults = c(beta = 0.3),
  init = c(S = 999999, I = 1),
  roles = list(
    susceptible = "S",
    infectious  = "I",
    incidence   = "incidence"
  )
)
