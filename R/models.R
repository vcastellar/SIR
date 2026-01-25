#-------------------------------------------------------------------------------
# modelo SI
#-------------------------------------------------------------------------------

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
#' ## Fit the model to observed incidence
#' fit_inc <- fit_epi_model(
#'   x = sim$incidence$inc,
#'   model = SI_MODEL,
#'   init = SI_MODEL$init,
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

SI_MODEL <- new_epi_model(
  name = "SI",
  rhs = si_rhs,
  state_names = c("S", "I"),
  par_names = c("beta"),
  outputs = c("S", "I", "incidence"),
  defaults = c(beta = 0.3),
  init = c(S = 999999, I = 1)
)

#-------------------------------------------------------------------------------
# modelo SIR
#-------------------------------------------------------------------------------

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
#' ## Model outputs
#' The SIR model declares the following outputs:
#' \describe{
#'   \item{\code{"S"}}{Susceptible population size.}
#'   \item{\code{"I"}}{Infectious population size.}
#'   \item{\code{"R"}}{Recovered (removed) population size.}
#'   \item{\code{"incidence"}}{Instantaneous rate of new infections
#'     \eqn{\lambda(t)} returned by the model's right-hand side.}
#' }
#'
#' All declared outputs may be used as observables in generic utilities such as
#' \code{\link{fit_epi_model}} via the \code{target} argument.
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
#' such as \code{\link{simulate_epi}}, \code{\link{fit_epi_model}}, and
#' \code{\link{predict.fit_epi_model}} that operate on \code{epi_model} objects.
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
#' ## Plot observed incidence (if an observation model is used)
#' plot(sim, what = "incidence")
#'
#' ## Fit the model to observed incidence
#' fit_inc <- fit_epi_model(
#'   x = sim$incidence$inc,
#'   model = SIR_MODEL,
#'   init = SIR_MODEL$init,
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
SIR_MODEL <- new_epi_model(
  name = "SIR",
  rhs = sir_rhs,
  state_names = c("S", "I", "R"),
  par_names = c("beta", "gamma"),
  lower = c(beta = 1e-8, gamma = 1e-8),
  upper = c(beta = 2,    gamma = 1),
  defaults = c(beta = 0.3, gamma = 0.1),
  init = c("S" = 1e6, "I" = 10, "R" = 0),
  outputs = c("S", "I", "R", "incidence")
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


#' SIRS epidemic model with waning immunity
#'
#' @name SIRS_MODEL
#' @description
#' An \code{epi_model} object representing a deterministic **SIRS**
#' (Susceptible–Infectious–Recovered–Susceptible) compartmental epidemic model
#' with waning immunity.
#'
#' The model describes the spread of an infection in a closed population where
#' susceptible individuals become infectious, subsequently recover, and may
#' lose immunity over time, returning to the susceptible compartment at rate
#' \code{omega}.
#'
#' @details
#' ## State variables
#' The model is defined in terms of the following state variables:
#' \describe{
#'   \item{S(t)}{Number of susceptible individuals at time \eqn{t}.}
#'   \item{I(t)}{Number of infectious (actively infected) individuals at time \eqn{t}.}
#'   \item{R(t)}{Number of recovered individuals with temporary immunity at time \eqn{t}.}
#' }
#'
#' The total population size is conserved:
#' \deqn{N = S(t) + I(t) + R(t).}
#'
#' ## Model outputs
#' The SIRS model declares the following outputs:
#' \describe{
#'   \item{\code{"S"}}{Susceptible population size.}
#'   \item{\code{"I"}}{Infectious population size.}
#'   \item{\code{"R"}}{Recovered (temporarily immune) population size.}
#'   \item{\code{"incidence"}}{Instantaneous rate of new infections
#'     \eqn{\lambda(t)} returned by the model's right-hand side.}
#' }
#'
#' All declared outputs may be used as observables in generic utilities such as
#' \code{\link{fit_epi_model}} via the \code{target} argument.
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
#' ## Usage
#' This predefined model object is intended to be used with generic utilities
#' such as \code{\link{simulate_epi}}, \code{\link{fit_epi_model}}, and
#' \code{\link{predict.fit_epi_model}} that operate on \code{epi_model} objects.
#'
#' @format
#' An object of class \code{"epi_model"}.
#'
#' @examples
#' ## Simulate a SIRS epidemic
#' sim <- simulate_epi(
#'   model = SIRS_MODEL,
#'   times = 0:200,
#'   parms = c(beta = 0.3, gamma = 0.1, omega = 0.02),
#'   init  = c(S = 1e6, I = 20, R = 0)
#' )
#'
#' plot(sim)
#'
#' ## Plot observed incidence (if an observation model is used)
#' plot(sim, what = "incidence")
#'
#' ## Fit the model to observed incidence
#' fit_inc <- fit_epi_model(
#'   x = sim$incidence$inc,
#'   model = SIRS_MODEL,
#'   init = SIRS_MODEL$init,
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

SIRS_MODEL <- new_epi_model(
  name = "SIRS",
  rhs = sirs_rhs,
  state_names = c("S", "I", "R"),
  par_names = c("beta", "gamma", "omega"),
  lower = c(beta = 1e-8, gamma = 1e-8, omega = 1e-8),
  upper = c(beta = 2,    gamma = 1,    omega = 1),
  defaults = c(beta = 0.3, gamma = 0.1, omega = 0.02),
  init = c("S" = 1e6, "I" = 20, "R" = 0),
  outputs = c("S", "I", "R", "incidence")
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
SEIR_MODEL <- new_epi_model(
  name = "SEIR",
  rhs = seir_rhs,
  state_names = c("S", "E", "I", "R"),
  par_names = c("beta", "sigma", "gamma"),
  lower = c(beta = 1e-8, sigma = 1e-8, gamma = 1e-8),
  upper = c(beta = 5,    sigma = 2,    gamma = 2),
  defaults = c(beta = 0.3, sigma = 0.2, gamma = 0.14),
  init = c("S" = 1e6, "E" = 5, "I" = 10, "R" = 0),
  outputs = c("S", "E", "I", "R", "incidence")
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
#' ## Model outputs
#' The SEIRS model declares the following outputs:
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
#' All declared outputs may be used as observables in generic utilities such as
#' \code{\link{fit_epi_model}} via the \code{target} argument.
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
#' such as \code{\link{simulate_epi}}, \code{\link{fit_epi_model}}, and
#' \code{\link{predict.fit_epi_model}} that operate on \code{epi_model} objects.
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
#'   init  = c(S = 1e6, E = 0, I = 20, R = 0),
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
#'   model = SEIRS_MODEL,
#'   init = SEIRS_MODEL$init,
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
SEIRS_MODEL <- new_epi_model(
  name = "SEIRS",
  rhs = seirs_rhs,
  state_names = c("S", "E", "I", "R"),
  par_names = c("beta", "sigma", "gamma", "omega"),
  lower = c(beta = 1e-8, sigma = 1e-8, gamma = 1e-8, omega = 1e-8),
  upper = c(beta = 5,    sigma = 2,    gamma = 2,    omega = 1),
  defaults = c(beta = 0.3, sigma = 0.2, gamma = 0.14, omega = 0.01),
  init = c("S" = 1e6, "I" = 10, "R" = 0, "E" = 0),
  outputs = c("S", "E", "I", "R", "incidence")
)




