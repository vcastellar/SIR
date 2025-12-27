#' Simulate an epidemic model (SIR or SIRS) and generate observed incidence
#' @name simulate_epi
#' @description
#' Simulates a deterministic compartmental epidemic model solved as an ODE system
#' (via \code{deSolve::ode()}) and produces both the latent epidemic trajectories
#' and a simple observation process for reported cases.
#'
#' The function supports \strong{SIR} and \strong{SIRS} dynamics, selected via
#' the \code{model} argument. The model is extended with an auxiliary compartment
#' \eqn{C(t)} that tracks the cumulative number of \emph{infection events} through
#' the incidence \eqn{\lambda(t)}.
#'
#' @details
#' ## Models
#' Let \eqn{N = S(t) + I(t) + R(t)} and define the incidence (new infections per day):
#' \deqn{\lambda(t) = \beta \frac{S(t) I(t)}{N}}
#'
#' ### SIR (\code{model = "sir"})
#' \deqn{
#' \begin{aligned}
#' \frac{dS}{dt} &= -\lambda(t) \\
#' \frac{dI}{dt} &= \lambda(t) - \gamma I(t) \\
#' \frac{dR}{dt} &= \gamma I(t)
#' \end{aligned}
#' }
#'
#' ### SIRS (\code{model = "sirs"})
#' In SIRS, recovered individuals lose immunity at rate \eqn{\omega} and return to
#' the susceptible compartment:
#' \deqn{
#' \begin{aligned}
#' \frac{dS}{dt} &= -\lambda(t) + \omega R(t) \\
#' \frac{dI}{dt} &= \lambda(t) - \gamma I(t) \\
#' \frac{dR}{dt} &= \gamma I(t) - \omega R(t)
#' \end{aligned}
#' }
#'
#' ## Cumulative infections
#' To provide a cumulative infection curve, the model includes:
#' \deqn{\frac{dC}{dt} = \lambda(t)}
#' so that \eqn{C(t)} accumulates infection events over time (continuous-time
#' analogue of cumulative incidence). Under \code{model = "sirs"}, reinfections are
#' possible, so \eqn{C(t)} may exceed \eqn{N}.
#'
#' ## Observation model
#' The function converts the latent incidence \eqn{\lambda(t)} into expected reported
#' counts using a reporting fraction \eqn{\rho}:
#' \deqn{\mu(t) = \rho \lambda(t)}
#' and then samples:
#' * \code{obs = "poisson"}: \eqn{Y(t) \sim \mathrm{Poisson}(\mu(t))}
#' * \code{obs = "negbin"}: \eqn{Y(t) \sim \mathrm{NegBin}(\mu(t), size)}
#'
#' The observed cumulative cases are computed as \code{cumsum(inc_obs)}.
#'
#' @param n_days Integer. Number of days to simulate (output will include day 0 to
#'   day \code{n_days}).
#' @param model Character. Epidemic model to simulate. One of \code{"sir"} or
#'   \code{"sirs"}.
#' @param N Numeric. Total population size.
#' @param beta Numeric. Transmission rate (per day).
#' @param gamma Numeric. Recovery/removal rate (per day).
#' @param omega Numeric. Loss-of-immunity rate (per day). Only used when
#'   \code{model = "sirs"}. Setting \eqn{\omega = 0} recovers SIR-like behavior.
#' @param I0 Numeric. Initial number of infectious individuals at day 0.
#' @param R0 Numeric. Initial number of removed/recovered individuals at day 0.
#' @param rho Numeric in \eqn{[0, 1]}. Reporting fraction mapping true incidence to
#'   expected observed cases.
#' @param obs Character. Observation model for reported incidence. Either
#'   \code{"poisson"} or \code{"negbin"}.
#' @param size Numeric. Dispersion (size) parameter for the negative binomial
#'   observation model. Larger values imply less overdispersion.
#' @param seed Optional integer. If provided, sets the random seed for reproducible
#'   observation draws.
#'
#' @return A named list with:
#' \describe{
#'   \item{params}{A list with the simulation parameters used (including \code{model}).}
#'   \item{states}{A data frame with columns \code{time}, \code{S}, \code{I},
#'     \code{R}, and \code{C} (cumulative infection events).}
#'   \item{incidence_true}{A data frame with columns \code{time} and \code{inc}
#'     (latent incidence \eqn{\lambda(t)}).}
#'   \item{incidence_obs}{A data frame with columns \code{time} and \code{inc}
#'     (reported incidence sampled from the observation model).}
#'   \item{cumulative_obs}{A data frame with columns \code{time} and
#'     \code{cases_cum} (cumulative reported cases).}
#' }
#'
#' @examples
#' # SIR simulation
#' sim_sir <- simulate_epi(
#'   n_days = 300, model = "sir",
#'   N = 1e6, beta = 0.30, gamma = 0.10,
#'   I0 = 20, rho = 0.25, obs = "negbin", size = 15, seed = 1
#' )
#' plot(sim_sir)
#' summary(sim_sir)
#' print(sim_sir)
#'
#' # SIRS simulation (average immunity duration ~ 180 days)
#' sim_sirs <- simulate_epi(
#'   n_days = 600, model = "sirs",
#'   N = 1e6, beta = 0.30, gamma = 0.10, omega = 1/180,
#'   I0 = 20, rho = 0.25, obs = "poisson", seed = 1
#' )
#' plot(sim_sirs)
#' summary(sim_sirs)
#' print(sim_sirs)
#' @export
simulate_epi <- function(
    n_days = 200,
    model = c("sir", "sirs"),
    N = 1e6,
    beta = 0.35,
    gamma = 0.10,
    omega = 0.0,        # solo se usa en SIRS (tasa pérdida inmunidad, 1/día)
    I0 = 10,
    R0 = 0,
    rho = 1,            # tasa de notificación (0-1)
    obs = c("negbin", "poisson"),
    size = 20,          # dispersión para NegBin
    seed = NULL
) {
  model <- match.arg(model)
  obs   <- match.arg(obs)
  if (!is.null(seed)) set.seed(seed)

  # RHS (sistema ODE) según modelo
  rhs <- switch(
    model,
    sir = function(time, state, parms) {
      with(as.list(c(state, parms)), {
        Ntot <- S + I + R
        lambda <- beta * S * I / Ntot
        dS <- -lambda
        dI <-  lambda - gamma * I
        dR <-  gamma * I
        dC <-  lambda
        list(c(dS, dI, dR, dC), incidence = lambda)
      })
    },
    sirs = function(time, state, parms) {
      with(as.list(c(state, parms)), {
        Ntot <- S + I + R
        lambda <- beta * S * I / Ntot
        dS <- -lambda + omega * R
        dI <-  lambda - gamma * I
        dR <-  gamma * I - omega * R
        dC <-  lambda
        list(c(dS, dI, dR, dC), incidence = lambda)
      })
    }
  )

  # Estado inicial
  init <- c(
    S = N - I0 - R0,
    I = I0,
    R = R0,
    C = I0 + R0
  )

  times <- 0:n_days

  # Parámetros según modelo
  parms <- if (model == "sir") {
    c(beta = beta, gamma = gamma)
  } else {
    c(beta = beta, gamma = gamma, omega = omega)
  }

  out <- deSolve::ode(
    y = init,
    times = times,
    func = rhs,
    parms = parms,
    method = "lsoda"
  )
  out <- as.data.frame(out)

  # Incidencia "verdadera"
  inc_true <- out$incidence

  # Observación: reportados ~ rho * incidencia
  mu <- pmax(rho * inc_true, 0)

  inc_obs <- switch(
    obs,
    poisson = stats::rpois(length(mu), lambda = mu),
    negbin  = stats::rnbinom(length(mu), mu = mu, size = size)
  )

  cum_obs <- cumsum(inc_obs)

  res <- list(
    params = list(
      model = model, N = N, beta = beta, gamma = gamma, omega = omega,
      I0 = I0, R0 = R0, rho = rho, obs = obs, size = size
    ),
    states = out[, c("time", "S", "I", "R", "C")],
    incidence_true = data.frame(time = out$time, inc = inc_true),
    incidence_obs  = data.frame(time = out$time, inc = inc_obs),
    cumulative_obs = data.frame(time = out$time, cases_cum = cum_obs)
  )

  class(res) <- "sim_epi"
  return(res)
}


