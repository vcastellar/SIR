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
#' The model is simulated at the time points specified by the numeric vector
#' \code{times}. This vector must be strictly increasing and of length >= 2.
#' The first value typically corresponds to the initial time (e.g. \code{0}).
#'
#' This design follows the interface of \code{\link[deSolve]{ode}} and provides
#' full control over the temporal resolution of the simulation, including
#' irregular or non-integer time grids.
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
#' @param times Numeric vector of time points at which to solve the ODE system.
#'   Must be strictly increasing and of length >= 2. The first value typically
#'   corresponds to the initial time (e.g. \code{0}).
#' @param time_unit Character string specifying the unit of time associated
#'   with \code{times}. One of \code{"days"}, \code{"weeks"}, \code{"month"},
#'   or \code{"year"}. This is used only for printing and plotting labels and
#'   does not affect the numerical simulation.
#'   The optional \code{time_unit} argument can be used to label the time axis
#'   in printed summaries and plots, but does not affect the numerical solution.
#' @param parms Named numeric vector of model parameters. Names must match
#'   \code{model$par_names}. Any missing parameters are taken from
#'   \code{model$defaults}, if available.
#' @param init Named numeric vector giving the initial values of the state
#'   variables. Names must exactly match \code{model$state_names}.
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
#' sim <- simulate_epi(
#'   model = SIR_MODEL,
#'   times = 0:200,
#'   time_unit = "week",
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
#' sim_rk4 <- simulate_epi(
#'   model = SIR_MODEL,
#'   times = 0:200,
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
#' sim_no_inc <- simulate_epi(
#'   model = SI_MODEL,
#'   times = 0:100,
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
#' sim_obs <- simulate_epi(
#'   model = SIR_MODEL,
#'   times = 0:150,
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
                         times,
                         time_unit = c("days", "weeks", "month", "year"),
                         parms = NULL,
                         init = NULL,
                         obs = c("none", "poisson", "negbin"),
                         size = 20,
                         seed = NULL,
                         ...) {

  stopifnot(inherits(model, "epi_model"))
  obs <- match.arg(obs)
  time_unit <- match.arg(time_unit)
  if (!is.null(seed)) set.seed(seed)

  # 1) Tiempos
  if (missing(times)) {
    stop("Provide `times` as a strictly increasing numeric vector.")
  }

  if (!is.numeric(times) || length(times) < 2 || anyNA(times)) {
    stop("`times` must be a numeric vector of length >= 2 with no missing values.")
  }

  if (any(diff(times) <= 0)) {
    stop("`times` must be strictly increasing.")
  }

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
    else data.frame(time = time_inc, cases_cum = inc_cum),
    time_unit = time_unit
  )

  class(res) <- "sim_epi"
  res
}


