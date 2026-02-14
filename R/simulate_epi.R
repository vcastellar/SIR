#' Simulate an epidemic model defined by an \code{epi_model} object
#'
#' @name simulate_epi
#'
#' @description
#' Simulates a deterministic compartmental epidemic model specified via an
#' \code{\link{epi_model}} object. The model is solved as a system of ordinary
#' differential equations (ODEs) using \code{\link[deSolve]{ode}}.
#'
#' The function is fully model-agnostic: all model-specific information
#' (state variables, parameters, default values, ODE right-hand side, and
#' optional flow definitions) is read directly from the supplied
#' \code{epi_model}.
#'
#' @details
#' ## Model structure
#' The \code{model} argument must be an object of class \code{"epi_model"},
#' typically created with \code{\link{epi_model}}. At minimum, the model
#' must define:
#' \itemize{
#'   \item an ODE right-hand side function (\code{model$rhs});
#'   \item the state variables (\code{model$states});
#'   \item the model parameters (\code{model$par_names}).
#' }
#'
#' Optionally, the model may declare one or more flows via \code{model$flows}.
#' If present, these flows must be returned by the ODE right-hand side function
#' as additional named flows and will be extracted into the simulation result.
#'
#' No assumptions are made about the semantic meaning of states or flows
#' (e.g. incidence, prevalence, cumulative counts).
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
#' ## Parameters
#' Model parameters are provided via the named numeric vector \code{parms}.
#' Names must match \code{model$par_names}. If the model defines default
#' parameter values in \code{model$defaults}, these are used as a baseline and
#' may be partially overridden by \code{parms}. All parameters must be fully
#' specified either via \code{parms} or \code{model$defaults}.
#'
#' ## Initial conditions
#' The initial state \code{init} must be supplied as a named numeric vector
#' whose names exactly match \code{model$states}. If \code{init} is not provided,
#' the model must define default initial conditions in \code{model$init}.
#'
#' No automatic construction of initial conditions is performed.
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
#' @param parms Named numeric vector of model parameters. Names must match
#'   \code{model$par_names}. Any missing parameters are taken from
#'   \code{model$defaults}, if available.
#' @param init Named numeric vector giving the initial values of the state
#'   variables. Names must exactly match \code{model$states}. If not provided,
#'   \code{model$init} must be defined.
#' @param seed Optional integer. If provided, sets the random seed. This argument
#'   is included for consistency with other functions but does not affect the
#'   deterministic ODE solution.
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
#'   \item{flows}{A data frame with columns \code{time} and the declared model
#'     flows, or \code{NULL} if the model defines no flows.}
#'   \item{time_unit}{Character string giving the time unit associated with
#'     the simulation.}
#' }
#'
#' @examples
#' ## ------------------------------------------------------------------
#' ## Example 1: SIR model
#' ## ------------------------------------------------------------------
#'
#' sim <- simulate_epi(
#'   model = SIR_MODEL,
#'   times = 0:200
#' )
#'
#' plot(sim)
#'
#'
#' ## ------------------------------------------------------------------
#' ## Example 2: Using a different numerical integration method
#' ## ------------------------------------------------------------------
#'
#' sim <- simulate_epi(
#'   model = SIR_MODEL,
#'   times = 0:200,
#'   parms = c(beta = 0.30, gamma = 0.10),
#'   init  = c(S = 1e6 - 10, I = 10, R = 0),
#'   method = "rk4"
#' )
#'
#' plot(sim)
#'
#'
#' ## ------------------------------------------------------------------
#' ## Example 3: Model without flows
#' ## ------------------------------------------------------------------
#'
#' sim <- simulate_epi(
#'   model = SI_MODEL,
#'   times = 30:100,
#'   parms = c(beta = 0.25),
#'   init  = c(S = 999, I = 1)
#' )
#'
#' plot(sim)
#'
#' @seealso
#' \code{\link{epi_model}}, \code{\link[deSolve]{ode}},
#' \code{\link{plot.sim_epi}}, \code{\link{print.sim_epi}}
#'
#' @export
simulate_epi <- function(model,
                         times,
                         time_unit = "days",
                         parms = NULL,
                         init = NULL,
                         seed = NULL,
                         ...) {

  stopifnot(inherits(model, "epi_model"))

  time_unit <- match.arg(
    time_unit,
    choices = c("days", "weeks", "month", "year")
  )

  if (!is.null(seed)) set.seed(seed)

  ## -------------------------------------------------------------------------
  ## 1) Tiempos
  ## -------------------------------------------------------------------------
  if (!is.numeric(times) || length(times) < 2 || anyNA(times)) {
    stop("`times` must be a numeric vector of length >= 2 with no missing values.")
  }
  if (any(diff(times) <= 0)) {
    stop("`times` must be strictly increasing.")
  }

  ## -------------------------------------------------------------------------
  ## 2) Parámetros
  ## -------------------------------------------------------------------------
  theta <- model$defaults

  if (is.null(theta)) {
    theta <- setNames(
      rep(NA_real_, length(model$par_names)),
      model$par_names
    )
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
    stop(
      "Missing parameters: ",
      paste(names(theta)[is.na(theta)], collapse = ", ")
    )
  }

  theta <- theta[model$par_names]

  ## -------------------------------------------------------------------------
  ## 3) Condiciones iniciales
  ## -------------------------------------------------------------------------
  if (is.null(init)) {
    if (is.null(model$init)) {
      stop("Provide `init` or define `model$init`.")
    }
    init <- model$init
  }

  init <- unlist(init)

  missing_states <- setdiff(model$states, names(init))
  if (length(missing_states) > 0) {
    stop(
      "Missing initial states: ",
      paste(missing_states, collapse = ", ")
    )
  }

  init <- init[model$states]

  ## -------------------------------------------------------------------------
  ## 4) Integración ODE
  ## -------------------------------------------------------------------------
  out <- deSolve::ode(
    y     = init,
    times = times,
    func  = model$rhs,
    parms = theta,
    ...
  )

  out <- as.data.frame(out)

  ## -------------------------------------------------------------------------
  ## 5) Extraer estados y flujos declarados
  ## -------------------------------------------------------------------------
  states_df <- out[, c("time", model$states), drop = FALSE]

  flows_df <- NULL
  if (!is.null(model$flows)) {
    missing_flows <- setdiff(model$flows, names(out))
    if (length(missing_flows) > 0) {
      stop(
        "Model RHS did not return declared flows: ",
        paste(missing_flows, collapse = ", ")
      )
    }
    flows_df <- out[, c("time", model$flows), drop = FALSE]
  }

  ## -------------------------------------------------------------------------
  ## 6) Resultado
  ## -------------------------------------------------------------------------
  res <- list(
    model     = model,
    params    = as.list(theta),
    states    = states_df,
    flows     = flows_df,
    time_unit = time_unit
  )

  class(res) <- "sim_epi"
  res
}
