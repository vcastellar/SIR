#-------------------------------------------------------------------------------
# modelo SIR
#------------------------------------------------------------------------------

#' @keywords internal
#' @noRd
sir_rhs <- function(time, state, parms) {
  with(as.list(c(state, parms)), {
    N <- S + I + R
    lambda <- beta * S * I / N
    dS <- -lambda
    dI <-  lambda - gamma * I
    dR <-  gamma * I
    list(c(dS, dI, dR), incidence = lambda)
  })
}

#' #' @keywords internal
#' #' @noRd
#' make_init_sir <- function(N, I0 = 10, R0 = 0) {
#'   c(S = N - I0 - R0, I = I0, R = R0)
#' }


#' SIR epidemic model with cumulative infections
#' @name SIR_MODEL
#' @description
#' An \code{epi_model} object representing a deterministic **SIR**
#' (Susceptible–Infectious–Recovered) compartmental epidemic model.
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
#' The total population size is given by
#' \deqn{N = S(t) + I(t) + R(t)}
#' which is conserved by the model dynamics.
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
#' \frac{dR}{dt} &= \gamma I(t), \\
#' \end{aligned}
#' }
#'
#'
#' ## Usage
#' This predefined model object is intended to be used with generic utilities
#' such as \code{\link{simulate_epi}} and model-fitting functions that operate
#' on \code{epi_model} objects.
#'
#' @format
#' An object of class \code{"epi_model"}.
#'
#' @examples
#' ## Simulate a SIR epidemic without an observation model
#' sim <- simulate_epi(
#'   model = SIR_MODEL,
#'   n_days = 200,
#'   parms = c(beta = 0.3, gamma = 0.1),
#'   init_args = list(N = 1e6, I0 = 20, R0 = 0),
#'   obs = "none"
#' )
#'
#' plot(sim$states$time, sim$states$I, type = "l",
#'      xlab = "Days", ylab = "I(t)",
#'      main = "SIR model: infectious individuals")
#'
#' @seealso
#' \code{\link{simulate_epi}}, \code{\link{new_epi_model}}
#'
#' @export
SIR_MODEL <- new_epi_model(
  name = "SIR",
  rhs = sir_rhs,
  state_names = c("S", "I", "R"),
  par_names = c("beta", "gamma"),
  lower = c(beta = 1e-8, gamma = 1e-8),
  upper = c(beta = 2,    gamma = 1),
  defaults = c(beta = 0.3, gamma = 0.1),
  init = c("S" = 1e6, "I" = 10, "R" = 0),
  output = list(incidence_col = "incidence")
)


#-------------------------------------------------------------------------------
# modelo SIRS
#-------------------------------------------------------------------------------
#' @keywords internal
#' @noRd
sirs_rhs <- function(time, state, parms) {
  with(as.list(c(state, parms)), {
    N <- S + I + R
    lambda <- beta * S * I / N
    dS <- -lambda + omega * R
    dI <-  lambda - gamma * I
    dR <-  gamma * I - omega * R
    list(c(dS, dI, dR), incidence = lambda)
  })
}

#' #' @keywords internal
#' #' @noRd
#' make_init_sirs <- function(N, I0 = 10, R0 = 0) {
#'   c(S = N - I0 - R0, I = I0, R = R0)
#' }

#' SIRS epidemic model with waning immunity
#' @name SIRS_MODEL
#' @description
#' An \code{epi_model} object representing a deterministic **SIRS**
#' (Susceptible–Infectious–Recovered–Susceptible) compartmental epidemic model
#' with waning immunity. Individuals who recover from infection lose immunity
#' at rate \code{omega} and return to the susceptible compartment.
#'
#'
#' @details
#' ## State variables
#' The model is defined in terms of the following state variables:
#' \describe{
#'   \item{S(t)}{Number of susceptible individuals at time \eqn{t}.}
#'   \item{I(t)}{Number of infectious (actively infected) individuals at time \eqn{t}.}
#'   \item{R(t)}{Number of recovered (temporarily immune) individuals at time \eqn{t}.}
#' }
#'
#' The total population size is given by
#' \deqn{N = S(t) + I(t) + R(t),}
#' which is conserved by the model dynamics.
#'
#' ## Parameters
#' The SIRS model depends on the following parameters:
#' \describe{
#'   \item{beta}{Transmission rate (per day).}
#'   \item{gamma}{Recovery/removal rate (per day).}
#'   \item{omega}{Rate of waning immunity from \code{R} back to \code{S} (per day).}
#' }
#'
#' ## Model equations
#' New infections occur at rate
#' \deqn{\lambda(t) = \beta \frac{S(t)\, I(t)}{N}.}
#'
#' The system of ordinary differential equations is:
#' \deqn{
#' \begin{aligned}
#' \frac{dS}{dt} &= -\lambda(t) + \omega R(t), \\
#' \frac{dI}{dt} &= \lambda(t) - \gamma I(t), \\
#' \frac{dR}{dt} &= \gamma I(t) - \omega R(t). \\
#' \end{aligned}
#' }
#'
#'
#' ## Usage
#' This predefined model object is intended to be used with generic utilities
#' such as \code{\link{simulate_epi}} and model-fitting functions that operate
#' on \code{epi_model} objects.
#'
#' @format
#' An object of class \code{"epi_model"}.
#'
#' @examples
#' ## Simulate a SIRS epidemic without an observation model
#' sim <- simulate_epi(
#'   model = SIRS_MODEL,
#'   n_days = 200,
#'   parms = c(beta = 0.3, gamma = 0.1, omega = 0.02),
#'   init = list(N = 1e6, I0 = 20, R0 = 0),
#'   obs = "none"
#' )
#'
#' plot(sim$states$time, sim$states$I, type = "l",
#'      xlab = "Days", ylab = "I(t)",
#'      main = "SIRS model: infectious individuals")
#'
#' @seealso
#' \code{\link{simulate_epi}}, \code{\link{new_epi_model}}
#'
#' @export
SIRS_MODEL <- new_epi_model(
  name = "SIRS",
  rhs = sirs_rhs,
  state_names = c("S", "I", "R"),
  par_names = c("beta", "gamma", "omega"),
  lower = c(beta = 1e-8, gamma = 1e-8, omega = 1e-8),
  upper = c(beta = 2,    gamma = 1,    omega = 1),
  defaults = c(beta = 0.3, gamma = 0.1, omega = 0.02),
  init = c("S" = 1e6, "I" = 20, "R" = 0)
)




#-------------------------------------------------------------------------------
# modelo SEIR
#-------------------------------------------------------------------------------
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

#' #' @keywords internal
#' #' @noRd
#' make_init_seir <- function(N, I0 = 10, R0 = 0, E0 = 0) {
#'   c(
#'     S = N - E0 - I0 - R0,
#'     E = E0,
#'     I = I0,
#'     R = R0
#'   )
#' }

#' SEIR epidemic model with latent (exposed) period
#' @name SEIR_MODEL
#' @description
#' An \code{epi_model} object representing a deterministic **SEIR**
#' (Susceptible–Exposed–Infectious–Recovered) compartmental epidemic model
#' with a latent period. Individuals become infected at rate \eqn{\lambda(t)}
#' and enter the exposed compartment \code{E}. Exposed individuals progress
#' to the infectious compartment at rate \code{sigma}.
#'
#' The model is extended with an auxiliary state variable \code{C(t)} that
#' tracks the cumulative number of **cases** over time (defined as transitions
#' from \code{E} to \code{I}), and it provides the instantaneous incidence
#' as an additional model output.
#'
#' @details
#' ## State variables
#' The model is defined in terms of the following state variables:
#' \describe{
#'   \item{S(t)}{Number of susceptible individuals at time \eqn{t}.}
#'   \item{E(t)}{Number of exposed (infected but not yet infectious) individuals at time \eqn{t}.}
#'   \item{I(t)}{Number of infectious (actively infected) individuals at time \eqn{t}.}
#'   \item{R(t)}{Number of recovered (immune) individuals at time \eqn{t}.}
#'   \item{C(t)}{Cumulative number of **cases** up to time \eqn{t} (entries into \code{I}).}
#' }
#'
#' The total population size is given by
#' \deqn{N = S(t) + E(t) + I(t) + R(t),}
#' which is conserved by the model dynamics.
#'
#' ## Parameters
#' The SEIR model depends on the following parameters:
#' \describe{
#'   \item{beta}{Transmission rate (per day).}
#'   \item{sigma}{Rate of progression from \code{E} to \code{I} (per day), so \eqn{1/\sigma} is the mean latent period.}
#'   \item{gamma}{Recovery/removal rate from \code{I} to \code{R} (per day), so \eqn{1/\gamma} is the mean infectious period.}
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
#' \frac{dR}{dt} &= \gamma I(t), \\
#' \frac{dC}{dt} &= \sigma E(t).
#' \end{aligned}
#' }
#'
#' The auxiliary state \eqn{C(t)} provides a continuous-time analogue of cumulative
#' **case** incidence (entries into \code{I}). This choice aligns the model with
#' common daily case-count time series.
#'
#' ## Usage
#' This predefined model object is intended to be used with generic utilities
#' such as \code{\link{simulate_epi}} and model-fitting functions that operate
#' on \code{epi_model} objects.
#'
#' @format
#' An object of class \code{"epi_model"}.
#'
#' @examples
#' ## Simulate a SEIR epidemic without an observation model
#' sim <- simulate_epi(
#'   model = SEIR_MODEL,
#'   n_days = 200,
#'   parms = c(beta = 0.3, sigma = 0.2, gamma = 0.14),
#'   init_args = list(N = 1e6, E0 = 0, I0 = 20, R0 = 0),
#'   obs = "poisson"
#' )
#'
#' plot(sim)
#'
#' @seealso
#' \code{\link{simulate_epi}}, \code{\link{new_epi_model}}
#'
#' @export
SEIR_MODEL <- new_epi_model(
  name = "SEIR",
  rhs = seir_rhs,
  state_names = c("S", "E", "I", "R"),
  par_names = c("beta", "sigma", "gamma"),
  lower = c(beta = 1e-8, sigma = 1e-8, gamma = 1e-8),
  upper = c(beta = 5,    sigma = 2,    gamma = 2),
  defaults = c(beta = 0.3, sigma = 0.2, gamma = 0.14),
  init = c("S" = 1e6, "E" = 5, "I" = 10, "R" = 0),
  output = list(
    incidence_col = "incidence",
    incidence_desc = "entries into I (cases)"
  )
)



#-------------------------------------------------------------------------------
# modelo SEIRS
#-------------------------------------------------------------------------------
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
      c(dS, dE, dI, dR),
      incidence = sigma * E
    )
  })
}


#' #' @keywords internal
#' #' @noRd
#' make_init_seirs <- function(N, I0 = 10, R0 = 0, E0 = 0) {
#'   c(
#'     S = N - E0 - I0 - R0,
#'     E = E0,
#'     I = I0,
#'     R = R0
#'   )
#' }

#' SEIRS epidemic model with latent period and waning immunity
#' @name SEIRS_MODEL
#' @description
#' An \code{epi_model} object representing a deterministic **SEIRS**
#' (Susceptible–Exposed–Infectious–Recovered–Susceptible) compartmental epidemic
#' model with a latent period and waning immunity. Individuals become infected
#' at rate \eqn{\lambda(t)} and enter the exposed compartment \code{E}. Exposed
#' individuals progress to \code{I} at rate \code{sigma}. Recovered individuals
#' lose immunity at rate \code{omega} and return to \code{S}.
#'
#' @details
#' ## State variables
#' \describe{
#'   \item{S(t)}{Number of susceptible individuals at time \eqn{t}.}
#'   \item{E(t)}{Number of exposed (infected but not yet infectious) individuals at time \eqn{t}.}
#'   \item{I(t)}{Number of infectious individuals at time \eqn{t}.}
#'   \item{R(t)}{Number of recovered (temporarily immune) individuals at time \eqn{t}.}
#' }
#'
#' The total population size is
#' \deqn{N = S(t) + E(t) + I(t) + R(t),}
#' which is conserved by the model dynamics.
#'
#' ## Parameters
#' \describe{
#'   \item{beta}{Transmission rate (per day).}
#'   \item{sigma}{Rate of progression from \code{E} to \code{I} (per day); \eqn{1/\sigma} is the mean latent period.}
#'   \item{gamma}{Recovery/removal rate from \code{I} to \code{R} (per day); \eqn{1/\gamma} is the mean infectious period.}
#'   \item{omega}{Rate of waning immunity from \code{R} back to \code{S} (per day); \eqn{1/\omega} is the mean immunity duration.}
#' }
#'
#' ## Model equations
#' New infections occur at rate
#' \deqn{\lambda(t) = \beta \frac{S(t)\, I(t)}{N}.}
#'
#' Case incidence (entries into \code{I}) occurs at rate
#' \deqn{\text{incidence}(t) = \sigma E(t).}
#'
#' The ODE system is:
#' \deqn{
#' \begin{aligned}
#' \frac{dS}{dt} &= -\lambda(t) + \omega R(t), \\
#' \frac{dE}{dt} &= \lambda(t) - \sigma E(t), \\
#' \frac{dI}{dt} &= \sigma E(t) - \gamma I(t), \\
#' \frac{dR}{dt} &= \gamma I(t) - \omega R(t). \\
#' \end{aligned}
#' }
#'
#' @format An object of class \code{"epi_model"}.
#'
#' @examples
#' ## Simulate a SEIRS epidemic
#' sim <- simulate_epi(
#'   model = SEIRS_MODEL,
#'   n_days = 300,
#'   parms = SEIRS_MODEL$default,
#'   init_args = list(N = 1e6, E0 = 0, I0 = 20, R0 = 0),
#'   obs = "poisson"
#' )
#' plot(sim)
#'
#' @seealso
#' \code{\link{simulate_epi}}, \code{\link{new_epi_model}}
#'
#' @export
SEIRS_MODEL <- new_epi_model(
  name = "SEIRS",
  rhs = seirs_rhs,
  state_names = c("S", "E", "I", "R"),
  par_names = c("beta", "sigma", "gamma", "omega"),
  lower = c(beta = 1e-8, sigma = 1e-8, gamma = 1e-8, omega = 1e-8),
  upper = c(beta = 5,    sigma = 2,    gamma = 2,    omega = 1),
  defaults = c(beta = 0.3, sigma = 0.2, gamma = 0.14, omega = 0.01),
  init = c("S" = 1e6, "I" = 10, "R" = 0, "E" = 0),
  output = list(
    incidence_col = "incidence",
    incidence_desc = "entries into I (cases)"
  )
)




