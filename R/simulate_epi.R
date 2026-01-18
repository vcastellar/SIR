#' Simulate an epidemic model defined by an \code{epi_model} object
#'
#' @name simulate_epi
#'
#' @description
#' Simulates a deterministic compartmental epidemic model specified via an
#' \code{\link{epi_model}} object. The model is solved as a system of ordinary
#' differential equations (ODEs) using \code{\link[deSolve]{ode}}. Optionally,
#' a stochastic observation process can be applied to a model-defined
#' incidence process to generate reported incidence counts.
#'
#' The function is fully model-agnostic: all model-specific information
#' (state variables, parameters, default values, ODE right-hand side, and
#' optional incidence definition) is read directly from the supplied
#' \code{epi_model}.
#'
#' @details
#' ## Model structure
#' The \code{model} argument must be an object of class \code{"epi_model"},
#' typically created with \code{\link{new_epi_model}}. At minimum, the model
#' must define:
#' \itemize{
#'   \item an ODE right-hand side function (\code{model$rhs});
#'   \item the state variables (\code{model$state_names});
#'   \item the model parameters (\code{model$par_names}).
#' }
#'
#' Optionally, the model may define a latent incidence process by returning
#' a variable named \code{"incidence"} as an additional output of the ODE
#' right-hand side function. If no such variable is defined, no incidence
#' time series is produced.
#'
#' ## Time grid
#' By default, the model is simulated on a regular daily grid from day 0 to
#' \code{n_days}. Alternatively, a custom numeric vector of time points can be
#' supplied via \code{times}, in which case \code{n_days} is ignored.
#'
#' ## Initial conditions
#' The initial state \code{init} must be supplied explicitly as a named numeric
#' vector whose names match \code{model$state_names}. No automatic construction
#' of initial conditions is performed.
#'
#' ## Observation model
#' If the model defines a latent incidence process, it can optionally be
#' converted into reported incidence counts using a simple observation model:
#' \describe{
#'   \item{\code{obs = "poisson"}}{Reported counts are drawn from a Poisson
#'     distribution with mean equal to the latent incidence.}
#'   \item{\code{obs = "negbin"}}{Reported counts are drawn from a negative binomial
#'     distribution with mean equal to the latent incidence and dispersion
#'     parameter \code{size}.}
#'   \item{\code{obs = "none"}}{No stochastic observation model is applied; the
#'     observed incidence is taken to be equal to the latent incidence.}
#' }
#'
#' If the model does not define a latent incidence process, the arguments
#' \code{obs} and \code{size} are ignored and no incidence outputs are produced.
#'
#' ## Numerical integration
#' The ODE system is solved using \code{\link[deSolve]{ode}}. Additional arguments
#' passed via \code{...} are forwarded directly to \code{deSolve::ode()}, allowing
#' full control over the numerical integration. In particular, the integration
#' method can be specified using the \code{method} argument.
#'
#' If \code{method} is not supplied, \code{deSolve::ode()} uses its default
#' integration method (typically \code{"lsoda"}).
#'
#' @param model An object of class \code{"epi_model"} defining the epidemic model
#'   to simulate.
#' @param n_days Integer. Number of days to simulate. Ignored if \code{times} is
#'   provided.
#' @param parms Named numeric vector of model parameters. Names must match
#'   \code{model$par_names}. Any missing parameters are taken from
#'   \code{model$defaults}, if available.
#' @param init Named numeric vector giving the initial values of the state
#'   variables. Names must exactly match \code{model$state_names}.
#' @param times Optional numeric vector of time points at which to solve the ODE
#'   system. If supplied, it overrides \code{n_days}.
#' @param obs Character string specifying the observation model. One of
#'   \code{"none"}, \code{"poisson"}, or \code{"negbin"}.
#' @param size Numeric. Dispersion (size) parameter for the negative binomial
#'   observation model.
#' @param seed Optional integer. If provided, sets the random seed for reproducible
#'   observation draws.
#' @param ... Additional arguments passed directly to
#'   \code{\link[deSolve]{ode}}, such as \code{method}, \code{rtol}, or \code{atol}.
#'
#' @return
#' An object of class \code{"sim_epi"}, containing:
#' \describe{
#'   \item{model}{The \code{epi_model} object used to generate the simulation.}
#'   \item{params}{A list of model parameter values used in the simulation.}
#'   \item{states}{A data frame with columns \code{time} and the model state
#'     variables, containing the simulated state trajectories.}
#'   \item{incidence}{A data frame with columns \code{time} and \code{inc}
#'     containing the observed (reported) incidence counts, or \code{NULL} if
#'     no incidence process is defined by the model.}
#'   \item{incidence_cum}{A data frame with columns \code{time} and
#'     \code{cases_cum} containing cumulative observed cases, or \code{NULL}.}
#' }
#'
#' @examples
#' ## ------------------------------------------------------------------
#' ## Example 1: SIR model with model-defined incidence
#' ## ------------------------------------------------------------------
#'
#' # In the built-in SIR model, the RHS defines a latent incidence
#' # (typically beta * S * I / N), returned as an extra output of the ODE.
#'
#' sim <- simulate_epi(
#'   model = SIR_MODEL,
#'   n_days = 200,
#'   parms = c(beta = 0.30, gamma = 0.10),
#'   init  = c(S = 1e6 - 10, I = 10, R = 0, C = 10),
#'   seed  = 1
#' )
#'
#' # Plot state trajectories defined by the model
#' plot(sim)
#'
#' # Plot observed incidence (requires the model to define incidence)
#' plot(sim, what = "incidence")
#'
#'
#' ## ------------------------------------------------------------------
#' ## Example 2: Using a different numerical integration method
#' ## ------------------------------------------------------------------
#'
#' # Additional arguments passed via `...` are forwarded to deSolve::ode().
#' # Here we explicitly select a fixed-step Runge--Kutta solver.
#'
#' sim_rk4 <- simulate_epi(
#'   model = SIR_MODEL,
#'   n_days = 200,
#'   parms = c(beta = 0.30, gamma = 0.10),
#'   init  = c(S = 1e6 - 10, I = 10, R = 0, C = 10),
#'   method = "rk4"
#' )
#'
#' plot(sim_rk4)
#'
#'
#' ## ------------------------------------------------------------------
#' ## Example 3: Model without an incidence definition
#' ## ------------------------------------------------------------------
#'
#' # If the model RHS does not return a variable named "incidence",
#' # no incidence time series is produced.
#'
#' # For such models, sim$incidence and sim$incidence_cum are NULL,
#' # and plotting incidence will result in an error.
#'
#' sim_no_inc <- simulate_epi(
#'   model = SI_MODEL,
#'   n_days = 100,
#'   parms = c(beta = 0.25),
#'   init  = c(S = 999, I = 1),
#'   obs   = "none"
#' )
#'
#' plot(sim_no_inc)              # valid: plots model states
#' # plot(sim_no_inc, what = "incidence")  # error: no incidence defined
#'
#'
#' ## ------------------------------------------------------------------
#' ## Example 4: Stochastic observation model
#' ## ------------------------------------------------------------------
#'
#' # When a model defines incidence, a stochastic observation process
#' # can be applied to generate reported case counts.
#'
#' sim_obs <- simulate_epi(
#'   model = SIR_MODEL,
#'   n_days = 150,
#'   parms = c(beta = 0.35, gamma = 0.12),
#'   init  = c(S = 1e5 - 5, I = 5, R = 0, C = 5),
#'   obs   = "negbin",
#'   size  = 20,
#'   seed  = 123
#' )
#'
#' plot(sim_obs, what = "incidence")

#'
#' @seealso
#' \code{\link{new_epi_model}}, \code{\link[deSolve]{ode}},
#' \code{\link{plot.sim_epi}}, \code{\link{print.sim_epi}}
#'
#' @export
simulate_epi <- function(model,
                         n_days = 200,
                         parms = NULL,
                         init = NULL,
                         times = NULL,
                         obs = c("none", "poisson", "negbin"),
                         size = 20,
                         seed = NULL,
                         ...) {

  stopifnot(inherits(model, "epi_model"))
  obs <- match.arg(obs)
  if (!is.null(seed)) set.seed(seed)

  # 1) Tiempos
  if (is.null(times)) times <- 0:n_days

  # 2) Parámetros: defaults + override
  theta <- model$defaults
  if (is.null(theta)) {
    theta <- setNames(rep(NA_real_, length(model$par_names)), model$par_names)
  }

  if (!is.null(parms)) {
    stopifnot(is.numeric(parms), !is.null(names(parms)))
    unknown <- setdiff(names(parms), model$par_names)
    if (length(unknown) > 0) {
      stop("Unknown parameters: ", paste(unknown, collapse = ", "))
    }
    theta[names(parms)] <- parms
  }

  if (any(is.na(theta))) {
    stop("Missing parameters: ",
         paste(names(theta)[is.na(theta)], collapse = ", "))
  }
  theta <- theta[model$par_names]

  # 3) Estado inicial
  if (is.null(init)) {
    stop("Provide `init`.")
  }
  init <- unlist(init)

  missing_states <- setdiff(model$state_names, names(init))
  if (length(missing_states) > 0) {
    stop("Missing initial states: ",
         paste(missing_states, collapse = ", "))
  }
  init <- init[model$state_names]

  # 4) Integración ODE
  out <- deSolve::ode(
    y = init,
    times = times,
    func = model$rhs,
    parms = theta,
    ...
  )
  out <- as.data.frame(out)

  # 5) Incidencia latente definida por el modelo (si existe)
  if ("incidence" %in% names(out)) {
    inc_true <- out$incidence[-1]
    time_inc <- out$time[-1]
  } else {
    inc_true <- NULL
    time_inc <- NULL
  }

  # 6) Modelo de observación
  if (is.null(inc_true)) {

    inc_obs <- NULL
    inc_cum <- NULL

  } else if (obs == "none") {

    inc_obs <- inc_true
    inc_cum <- cumsum(inc_obs)

  } else {

    inc_obs <- switch(
      obs,
      poisson = stats::rpois(length(inc_true), lambda = inc_true),
      negbin  = stats::rnbinom(length(inc_true), mu = inc_true, size = size)
    )
    inc_cum <- cumsum(inc_obs)
  }

  # 7) Resultado
  res <- list(
    model = model,
    params = as.list(theta),
    states = out[, c("time", model$state_names), drop = FALSE],
    incidence = if (is.null(inc_obs)) NULL
    else data.frame(time = time_inc, inc = inc_obs),
    incidence_cum = if (is.null(inc_cum)) NULL
    else data.frame(time = time_inc, cases_cum = inc_cum)
  )

  class(res) <- "sim_epi"
  res
}
