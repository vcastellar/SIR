new_predict_epi_model <- function(fit,
                                  times,
                                  init,
                                  states,
                                  outputs,
                                  ode_control) {

  stopifnot(inherits(fit, "fit_epi_model"))

  structure(
    list(
      fit         = fit,
      times       = times,
      init        = init,        # estado final del fit
      states      = states,      # SOLO predicción
      outputs     = outputs,
      ode_control = ode_control
    ),
    class = "predict_epi_model"
  )
}




#' Predict epidemic dynamics from a fitted epidemic model
#'
#' @name predict.fit_epi_model
#' @method predict fit_epi_model
#'
#' @description
#' Computes forward predictions from a fitted deterministic epidemic model by
#' solving the underlying system of ordinary differential equations (ODEs) using
#' the parameter estimates obtained during model fitting.
#'
#' The prediction is *conditional on the fitted model*: the returned object
#' contains both the original \code{"fit_epi_model"} object and the simulated
#' trajectories obtained by propagating the model forward in time from a given
#' initial state.
#'
#' @param object An object of class \code{"fit_epi_model"} representing a fitted
#'   epidemic model.
#' @param times Numeric vector of time points at which to compute predictions.
#'   Must be strictly increasing and of length at least two.
#' @param init Named numeric vector specifying the initial state of the model
#'   at the first prediction time. If \code{NULL} (default), the initial state
#'   used during model fitting (\code{object$ini0}) is reused.
#' @param ... Additional arguments passed directly to
#'   \code{\link[deSolve]{ode}} to control numerical integration (e.g.
#'   \code{method}, \code{rtol}, \code{atol}). These arguments override any ODE
#'   control settings stored in the fitted object.
#'
#' @return
#' An object of class \code{"predict_epi_model"} containing:
#' \itemize{
#'   \item the fitted epidemic model (\code{fit});
#'   \item the time grid used for prediction;
#'   \item the initial state used for prediction;
#'   \item the predicted state trajectories;
#'   \item the predicted model outputs;
#'   \item the numerical integration settings used.
#' }
#'
#' This object preserves full traceability between the fitted model and the
#' resulting predictions.
#'
#' @examples
#' \dontrun{
#' ## ------------------------------------------------------------
#' ## Fit a model and generate forward predictions
#' ## ------------------------------------------------------------
#' sim <- simulate_epi(
#'   model = SEIRS_MODEL,
#'   times = 0:200,
#'   parms = SEIRS_MODEL$defaults,
#'   init  = SEIRS_MODEL$init,
#'   obs   = "negbin"
#' )
#' plot(sim)
#' ## Fit the model to observed incidence
#' fit <- fit_epi_model(
#'   x     = sim$incidence$inc,
#'   model = SEIRS_MODEL,
#'   target = "incidence",
#'   init  = SEIRS_MODEL$init
#' )
#'
#' ## Use the final simulated state as the prediction initial condition
#' init_pred <- tail(sim$states, n = 1)[, -1]
#'
#' ## Generate predictions on a future time grid
#' pred <- predict(
#'   object = fit,
#'   times  = 201:800
#' )
#'
#' plot(pred)
#'
#' ## Inspect predicted trajectories
#' plot(pred$states$time, pred$states$I, type = "l")
#' }
#'
#' @seealso
#' \code{\link{fit_epi_model}}, \code{\link{simulate_epi}}
#'
#' @export
predict.fit_epi_model <- function(object,
                                  times,
                                  init = NULL,
                                  method = "lsoda",
                                  ...) {

  stopifnot(inherits(object, "fit_epi_model"))

  model <- object$model
  parms <- object$par
  states <- model$state_names

  ## ------------------------------------------------------------
  ## 1) Reconstruir trayectoria del fit
  ## ------------------------------------------------------------
  times_fit <- seq_along(object$x) - 1

  out_fit <- deSolve::ode(
    y      = object$ini0,
    times  = times_fit,
    func   = model$rhs,
    parms  = parms,
    method = "lsoda"
  )
  out_fit <- as.data.frame(out_fit)

  last_state <- out_fit[nrow(out_fit), states, drop = FALSE]

  ## ------------------------------------------------------------
  ## 2) Validar / fijar init
  ## ------------------------------------------------------------
  if (!is.null(init)) {
    init <- as.numeric(init[states])
    if (any(abs(init - as.numeric(last_state)) > 1e-8)) {
      stop(
        "`init` does not match the final fitted state.\n",
        "Predictions must continue the fitted trajectory.",
        call. = FALSE
      )
    }
  }

  init_pred <- as.numeric(last_state)
  names(init_pred) <- states

  ## ------------------------------------------------------------
  ## 3) Integrar predicción
  ## ------------------------------------------------------------
  out_pred <- deSolve::ode(
    y      = init_pred,
    times  = times,
    func   = model$rhs,
    parms  = parms,
    method = method,
    ...
  )
  out_pred <- as.data.frame(out_pred)

  ## ------------------------------------------------------------
  ## 4) Construir objeto
  ## ------------------------------------------------------------
  new_predict_epi_model(
    fit         = object,
    times       = times,
    init        = init_pred,
    states      = out_pred[, c("time", states), drop = FALSE],
    outputs     = model$outputs,
    ode_control = list(method = method)
  )
}



#' Plot fitted and predicted state trajectories
#'
#' @name plot.predict_epi_model
#' @method plot predict_epi_model
#'
#' @param x Object of class "predict_epi_model".
#' @param lwd_fit Line width for fitted trajectory.
#' @param lwd_pred Line width for predicted trajectory.
#' @param lty_pred Line type for predicted trajectory.
#' @param col Optional vector of colors for states.
#' @param ... Additional graphical parameters.
#'
#' @return Invisibly returns x.
#'
#' @export
plot.predict_epi_model <- function(x,
                                   lwd = 2,
                                   col = NULL,
                                   ...) {

  stopifnot(inherits(x, "predict_epi_model"))

  fit   <- x$fit
  model <- fit$model
  states <- model$state_names
  n_st <- length(states)

  if (is.null(col)) col <- seq_len(n_st)

  ## ------------------------------------------------------------
  ## 1) Reconstruir fit
  ## ------------------------------------------------------------
  times_fit <- seq_along(fit$x) - 1

  out_fit <- deSolve::ode(
    y      = fit$ini0,
    times  = times_fit,
    func   = model$rhs,
    parms  = fit$par,
    method = "lsoda"
  )
  out_fit <- as.data.frame(out_fit)

  ## ------------------------------------------------------------
  ## 2) Concatenar (SIN duplicar frontera)
  ## ------------------------------------------------------------
  out_pred <- x$states

  time_all <- c(
    out_fit$time,
    out_pred$time[-1]
  )

  states_all <- rbind(
    out_fit[, states, drop = FALSE],
    out_pred[-1, states, drop = FALSE]
  )

  ## ------------------------------------------------------------
  ## 3) Plot
  ## ------------------------------------------------------------
  matplot(
    time_all,
    states_all,
    type = "l",
    lwd  = lwd,
    col  = col,
    xlab = "Time",
    ylab = "State values",
    main = paste("Fitted + predicted trajectories –", model$name),
    ...
  )

  abline(
    v   = max(out_fit$time),
    lty = 3,
    col = "grey50"
  )

  legend(
    "topright",
    legend = states,
    col = col,
    lty = 1,
    lwd = lwd,
    bty = "n"
  )

  invisible(x)
}
