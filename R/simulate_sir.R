#' Simulate an SIR epidemic and generate observed incidence
#' @name simulate_sir
#' @description
#' Simulates a deterministic SIR model (solved as an ODE system) and produces both
#' the latent epidemic trajectories and a simple observation process for reported
#' cases. The model is extended with an auxiliary compartment \eqn{C(t)} that
#' tracks the cumulative number of infections via the incidence (force of infection).
#'
#' @details
#' ## Model
#' The SIR dynamics are:
#' \deqn{
#' \begin{aligned}
#' \frac{dS}{dt} &= -\lambda(t) \\
#' \frac{dI}{dt} &= \lambda(t) - \gamma I(t) \\
#' \frac{dR}{dt} &= \gamma I(t)
#' \end{aligned}
#' }
#' where the incidence \eqn{\lambda(t)} (new infections per day) is defined as:
#' \deqn{\lambda(t) = \beta \frac{S(t) I(t)}{N}}
#' with \eqn{N = S(t) + I(t) + R(t)}.
#'
#' To provide a cumulative infection curve, the model includes:
#' \deqn{\frac{dC}{dt} = \lambda(t)}
#' so that \eqn{C(t)} accumulates infections over time (continuous-time analogue
#' of cumulative incidence).
#'
#' ## Observation model
#' The function converts the latent incidence \eqn{\lambda(t)} into reported counts
#' using a reporting fraction \eqn{\rho}:
#' \deqn{\mu(t) = \rho \lambda(t)}
#' and then samples:
#' * \code{obs = "poisson"}: \eqn{Y(t) \sim \mathrm{Poisson}(\mu(t))}
#' * \code{obs = "negbin"}: \eqn{Y(t) \sim \mathrm{NegBin}(\mu(t), size)}
#'
#' The observed cumulative cases are computed as \code{cumsum(inc_obs)}.
#'
#' @param n_days Integer. Number of days to simulate (output will include day 0 to
#'   day \code{n_days}).
#' @param N Numeric. Total population size.
#' @param beta Numeric. Transmission rate (per day).
#' @param gamma Numeric. Recovery/removal rate (per day).
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
#'   \item{params}{A list with the simulation parameters used.}
#'   \item{states}{A data frame with columns \code{time}, \code{S}, \code{I},
#'     \code{R}, and \code{C} (cumulative infections).}
#'   \item{incidence_true}{A data frame with columns \code{time} and \code{inc}
#'     (latent incidence \eqn{\lambda(t)}).}
#'   \item{incidence_obs}{A data frame with columns \code{time} and \code{inc}
#'     (reported incidence sampled from the observation model).}
#'   \item{cumulative_obs}{A data frame with columns \code{time} and
#'     \code{cases_cum} (cumulative reported cases).}
#' }
#'
#' @examples
#' sim <- simulate_sir(
#'   n_days = 300,
#'   N = 1e6,
#'   beta = 0.30,
#'   gamma = 0.10,
#'   I0 = 20,
#'   rho = 0.25,
#'   obs = "negbin",
#'   size = 15,
#'   seed = 1
#' )
#'
#' head(sim$incidence_obs)
#' head(sim$cumulative_obs)
#'
#' sim$states$I
#'
#' plot(sim$states$time, sim$states$I, type = "l",
#'      xlab = "Días", ylab = "I(t) (infectados activos)", main = "SIR: infectados activos")
#'
#' plot(sim$incidence_obs$time, sim$incidence_obs$inc, type = "h",
#'      xlab = "Días", ylab = "Incidencia observada", main = "Incidencia observada (conteos)")
#'
#' @export
simulate_sir <- function(n_days = 200, N = 1e6, beta = 0.35, gamma = 0.10,
    I0 = 10, R0 = 0,
    rho = 1,          # tasa de notificación (0-1)
    obs = c("negbin", "poisson"),
    size = 20,          # dispersión para NegBin (más grande = menos overdisp)
    seed = NULL
) {
  obs <- match.arg(obs)
  if (!is.null(seed)) set.seed(seed)

  # Estado inicial
  init <- c(S = N - I0 - R0,
            I = I0,
            R = R0,
            C = I0 + R0) # C: acumulado de infectados (aprox)
  times <- 0:n_days


  out <- deSolve::ode(
    y = init,
    times = times,
    func = sir,
    parms = c(beta = beta, gamma = gamma),
    method = "lsoda"
  )
  out <- as.data.frame(out)

  # Incidencia "verdadera" por día (en la rejilla diaria)
  inc_true <- out$incidence

  # Observación: conteos reportados ~ rho * incidencia
  mu <- pmax(rho * inc_true, 0)

  inc_obs <- switch(
    obs,
    poisson = rpois(length(mu), lambda = mu),
    negbin  = rnbinom(length(mu), mu = mu, size = size)
  )

  # Acumulado observado a partir de incidencia observada
  cum_obs <- cumsum(inc_obs)

  list(
    params = list(N = N, beta = beta, gamma = gamma, I0 = I0, rho = rho, obs = obs, size = size),
    states = out[, c("time", "S", "I", "R", "C")],
    incidence_true = data.frame(time = out$time, inc = inc_true),
    incidence_obs  = data.frame(time = out$time, inc = inc_obs),
    cumulative_obs = data.frame(time = out$time, cases_cum = cum_obs)
  )
}


