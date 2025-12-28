#' Simulate an epidemic model defined by an `epi_model` object
#' @name simulate_epi
#' @description
#' Simulates a deterministic compartmental epidemic model (e.g. SIR, SIRS, or any
#' user-defined extension) specified via an \code{epi_model} object. The model is
#' solved as an ODE system using \code{deSolve::ode()}, and can optionally include
#' an observation process to generate reported incidence counts.
#'
#' The function is model-agnostic: all model-specific information (state variables,
#' parameters, ODE right-hand side, bounds, and output conventions) is read from
#' the supplied \code{epi_model}.
#'
#' @details
#' ## Model structure
#' The \code{model} argument must be an object of class \code{"epi_model"}, typically
#' created with \code{new_epi_model()}. At minimum, the model must define:
#' \itemize{
#'   \item a right-hand side ODE function (\code{model$rhs});
#'   \item the required state variables (\code{model$state_names});
#'   \item the model parameters (\code{model$par_names}).
#' }
#'
#' Optionally, the model may also define:
#' \itemize{
#'   \item default parameter values (\code{model$defaults});
#'   \item parameter bounds (\code{model$lower}, \code{model$upper});
#'   \item a helper function to construct initial conditions
#'     (\code{model$make_init});
#'   \item standard output column names for incidence and cumulative infections
#'     (\code{model$output}).
#' }
#'
#' ## Time grid
#' By default, the model is simulated on a daily grid from day 0 to
#' \code{n_days}. Alternatively, a custom numeric vector of times can be supplied
#' via \code{times}.
#'
#' ## Observation model
#' The latent incidence produced by the ODE model can be converted into reported
#' incidence counts using a simple observation model:
#' \describe{
#'   \item{\code{obs = "poisson"}}{Reported counts are drawn from a Poisson
#'     distribution with mean \eqn{\mu(t) = \rho \lambda(t)}.}
#'   \item{\code{obs = "negbin"}}{Reported counts are drawn from a negative binomial
#'     distribution with mean \eqn{\mu(t) = \rho \lambda(t)} and dispersion
#'     parameter \code{size}.}
#'   \item{\code{obs = "none"}}{No observation process is applied; observed
#'     incidence and cumulative counts are returned as \code{NA}.}
#' }
#'
#' @param model An object of class \code{"epi_model"} defining the epidemic model
#'   to simulate (e.g. SIR, SIRS).
#' @param n_days Integer. Number of days to simulate. Ignored if \code{times} is
#'   provided.
#' @param parms Named numeric vector of model parameters. Must match
#'   \code{model$par_names}. Any missing parameters are taken from
#'   \code{model$defaults}, if available.
#' @param init Named numeric vector with initial state values. If \code{NULL}
#'   (default), the function uses \code{model$make_init} together with
#'   \code{init_args}.
#' @param init_args Named list of arguments passed to \code{model$make_init} to
#'   construct the initial state (e.g. \code{N}, \code{I0}, \code{R0}).
#' @param times Optional numeric vector of time points at which to solve the ODE
#'   system. If supplied, it overrides \code{n_days}.
#' @param rho Numeric in \eqn{[0,1]}. Reporting fraction mapping true incidence to
#'   expected observed incidence.
#' @param obs Character string specifying the observation model. One of
#'   \code{"poisson"}, \code{"negbin"}, or \code{"none"}.
#' @param size Numeric. Dispersion (size) parameter for the negative binomial
#'   observation model. Larger values imply less overdispersion.
#' @param seed Optional integer. If provided, sets the random seed for reproducible
#'   observation draws.
#' @param method Character string. Integration method passed to
#'   \code{deSolve::ode()} (default: \code{"lsoda"}).
#'
#' @return A named list with the following components:
#' \describe{
#'   \item{model}{Name of the epidemic model.}
#'   \item{params}{List of model parameters used in the simulation.}
#'   \item{states}{Data frame with columns \code{time} and the model state
#'     variables (e.g. \code{S}, \code{I}, \code{R}, \code{C}).}
#'   \item{incidence_true}{Data frame with columns \code{time} and \code{inc}
#'     containing the latent model incidence.}
#'   \item{incidence_obs}{Data frame with columns \code{time} and \code{inc}
#'     containing the observed (reported) incidence counts.}
#'   \item{cumulative_obs}{Data frame with columns \code{time} and
#'     \code{cases_cum} containing cumulative observed cases.}
#' }
#'
#' @examples
#' ## SIR simulation
#' sim <- simulate_epi(
#'   model = SIR_MODEL,
#'   n_days = 200,
#'   parms = c(beta = 0.30, gamma = 0.10),
#'   init_args = list(N = 1e6, I0 = 20, R0 = 0),
#'   rho = 0.3,
#'   obs = "poisson",
#'   seed = 1
#' )
#'
#' plot(sim)
#'
#' ## SIRS simulation
#' sim2 <- simulate_epi(
#'   model = SIRS_MODEL,
#'   n_days = 200,
#'   parms = c(beta = 0.30, gamma = 0.10, omega = 0.02),
#'   init_args = list(N = 1e6, I0 = 20, R0 = 0),
#'   obs = "none"
#' )
#' plot(sim2)
#'
#' @seealso
#' \code{\link{new_epi_model}}, \code{\link{deSolve::ode}}
#'
#' @export

simulate_epi <- function(model,
                         n_days = 200,
                         parms = NULL,
                         init = NULL,
                         init_args = list(N = 1e6, I0 = 10, R0 = 0),
                         times = NULL,
                         rho = 1,
                         obs = c("negbin", "poisson", "none"),
                         size = 20,
                         seed = NULL,
                         method = "lsoda") {
  stopifnot(inherits(model, "epi_model"))

  obs <- match.arg(obs)
  if (!is.null(seed)) set.seed(seed)

  # 1) Tiempos
  if (is.null(times)) times <- 0:n_days

  # 2) Parámetros: defaults + user override, ordenados por par_names
  theta <- model$defaults
  if (is.null(theta)) theta <- setNames(rep(NA_real_, length(model$par_names)), model$par_names)

  if (!is.null(parms)) {
    stopifnot(is.numeric(parms))
    if (is.null(names(parms))) stop("`parms` debe ser un vector *nombrado* con parámetros del modelo.")
    unknown <- setdiff(names(parms), model$par_names)
    if (length(unknown) > 0) stop("Parámetros desconocidos: ", paste(unknown, collapse = ", "))
    theta[names(parms)] <- parms
  }

  if (any(is.na(theta))) {
    miss <- names(theta)[is.na(theta)]
    stop("Faltan parámetros: ", paste(miss, collapse = ", "),
         ". Pásalos en `parms` o define `defaults` en el modelo.")
  }
  theta <- theta[model$par_names]

  # 3) Bounds (si están)
  if (!is.null(model$lower)) {
    bad <- theta < model$lower[model$par_names]
    if (any(bad)) stop("Parámetros por debajo del lower: ", paste(names(theta)[bad], collapse = ", "))
  }
  if (!is.null(model$upper)) {
    bad <- theta > model$upper[model$par_names]
    if (any(bad)) stop("Parámetros por encima del upper: ", paste(names(theta)[bad], collapse = ", "))
  }

  # 4) Estado inicial
  if (is.null(init)) {
    if (!is.function(model$make_init)) {
      stop("No se proporcionó `init` y el modelo no tiene `make_init`.")
    }
    init <- do.call(model$make_init, init_args)
  }

  if (is.null(names(init))) stop("`init` debe ser un vector *nombrado*.")
  miss_state <- setdiff(model$state_names, names(init))
  if (length(miss_state) > 0) stop("Faltan estados en `init`: ", paste(miss_state, collapse = ", "))
  init <- init[model$state_names]

  # 5) Integración ODE
  out <- deSolve::ode(
    y = init,
    times = times,
    func = model$rhs,
    parms = theta,
    method = method
  )
  out <- as.data.frame(out)

  # 6) Extraer incidencia y acumulado de forma estándar (según el modelo)
  inc_col <- model$output$incidence_col %||% "incidence"
  cum_col <- model$output$cumulative_col %||% "C"

  if (!inc_col %in% names(out)) stop("La salida ODE no contiene la columna de incidencia '", inc_col, "'.")
  if (!cum_col %in% names(out)) warning("La salida ODE no contiene la columna acumulada '", cum_col, "'.")

  inc_true <- out[[inc_col]]

  # 7) Modelo de observación (opcional)
  if (obs == "none") {
    inc_obs <- rep(NA_integer_, length(inc_true))
  } else {
    mu <- pmax(rho * inc_true, 0)
    inc_obs <- switch(
      obs,
      poisson = stats::rpois(length(mu), lambda = mu),
      negbin  = stats::rnbinom(length(mu), mu = mu, size = size)
    )
  }

  cum_obs <- if (obs == "none") rep(NA_real_, length(inc_true)) else cumsum(inc_obs)

  res <- list(
    model = model$name,
    params = as.list(theta),
    states = out[, c("time", model$state_names), drop = FALSE],
    incidence_true = data.frame(time = out$time, inc = inc_true),
    incidence_obs  = data.frame(time = out$time, inc = inc_obs),
    cumulative_obs = data.frame(time = out$time, cases_cum = cum_obs)
  )


  class(res) <- "sim_epi"
  return(res)
}

# helper: operador "si NULL entonces"
`%||%` <- function(a, b) if (is.null(a)) b else a

