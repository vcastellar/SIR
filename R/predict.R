#' Predict epidemic dynamics from a fitted epidemic model
#'
#' @name predict.fit_epi_model
#' @method predict fit_epi_model
#'
#' @description
#' Computes forward predictions from a fitted epidemic model by solving the
#' underlying ODE system using the estimated parameters.
#'
#' The prediction is deterministic and represents the expected epidemic
#' trajectory conditional on the fitted parameters.
#'
#' @param object An object of class \code{"fit_epi_model"}.
#' @param n_days Integer. Number of days to predict ahead.
#' @param init Named numeric vector of initial states (e.g. \code{S, I, R, C}).
#'   If \code{NULL}, initial conditions are constructed using \code{init_args}.
#' @param init_args Named list passed to \code{model$make_init()} (e.g.
#'   \code{list(N = 1e6, I0 = 10, R0 = 0)}).
#' @param times Optional numeric vector of time points. Overrides \code{n_days}.
#' @param type Character. One of \code{"states"}, \code{"incidence"}, or \code{"both"}.
#' @param method Integration method passed to \code{deSolve::ode()}.
#' @param ... Currently unused.
#'
#' @return
#' A list containing predicted epidemic quantities:
#' \describe{
#'   \item{states}{Predicted compartment trajectories.}
#'   \item{incidence}{Predicted incidence time series.}
#' }
#'
#' @examples
#' \dontrun{
#' sim <- simulate_epi(
#'   model = SIRS_MODEL,
#'   n_days = 200,
#'   parms = c(beta = 0.30, gamma = 0.10, omega = 0.02),
#'   init_args = list(N = 1e6, I0 = 20, R0 = 0),
#'   obs = "poisson"
#' )
#' plot(sim)
#' x <- sim$incidence_obs$inc
#' plot(x, type = "l", xlab = "Day", ylab = "Incidence")
#' fit <- fit_epi_model(x, model = SIRS_MODEL, init = list(I = 6, N = 1e6))
#' ini0 <- fit$ini0
#' out <- deSolve::ode(
#'   y = ini0,
#'   times = 1:fit$x_len,
#'   func = SIR_MODEL$rhs,
#'   parms = fit$par,
#'   method = "lsoda"
#' )
#' fin <- tail(out, n = 1)
#' iniF <- c(S = fin[2], I = fin[3], R = fin[4], C = fin[5])
#' pred <- predict(
#'   object = fit,
#'   n_days = 300,
#'   init = iniF,
#'   type = "both"
#' )
#'
#' plot(pred$states$time, pred$states$I, type = "l", ylim = c(0, max(c(pred$states$I, pred$states$IR, pred$states$S))))
#' lines(pred$states$R, col = "red", lty = 2)
#' lines(pred$states$S, col = "blue", lty = 2)
#' }
#'
#' @export
predict.fit_epi_model <- function(object,
                                  n_days,
                                  init = NULL,
                                  init_args = NULL,
                                  times = NULL,
                                  type = c("both", "states", "incidence"),
                                  method = "lsoda",
                                  ...) {

  stopifnot(inherits(object, "fit_epi_model"))
  type <- match.arg(type)

  model <- object$model
  parms <- object$par

  # 1) Time grid
  if (is.null(times)) {
    stopifnot(length(n_days) == 1, n_days >= 1)
    times <- 0:n_days
  }

  # 2) Initial state
  if (is.null(init)) {
    if (is.null(init_args)) {
      stop("Either `init` or `init_args` must be provided.")
    }
    if (!is.function(model$make_init)) {
      stop("Model does not define `make_init()`.")
    }
    init <- do.call(model$make_init, init_args)
  }

  if (is.null(names(init))) {
    stop("`init` must be a named numeric vector.")
  }

  init <- init[model$state_names]

  # 3) ODE integration
  out <- deSolve::ode(
    y = init,
    times = times,
    func = model$rhs,
    parms = parms,
    method = "lsoda"
  )

  out <- as.data.frame(out)

  # 4) Extract outputs
  inc_col <- model$output$incidence_col %||% "incidence"

  res_states <- out[, c("time", model$state_names), drop = FALSE]
  res_inc <- data.frame(
    time = out$time,
    inc = out[[inc_col]]
  )

  # 5) Return according to type
  switch(
    type,
    states = res_states,
    incidence = res_inc,
    both = list(states = res_states, incidence = res_inc)
  )
}

