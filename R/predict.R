#' Predict epidemic dynamics from a fitted epidemic model
#'
#' @name predict.fit_epi_model
#' @method predict fit_epi_model
#'
#' @description
#' Computes forward predictions from a fitted epidemic model by solving the
#' underlying ODE system using the estimated parameters.
#'
#' @param object An object of class \code{"fit_epi_model"}.
#' @param n_days Integer. Number of days to predict ahead.
#' @param init Named numeric vector of initial states. If \code{NULL}, the initial
#'   state is reconstructed using \code{object$ini0}.
#' @param times Optional numeric vector of time points. Overrides \code{n_days}.
#'   Defaults to \code{0:(n_days - 1)} when \code{NULL}.
#' @param method Integration method passed to \code{deSolve::ode()}.
#' @param ... Currently unused.
#'
#' @return
#' An object of class \code{"epi_model_predict"}.
#' @examples
#' \dontrun{
#' sim <- simulate_epi(
#'   model = SIR_MODEL,
#'   n_days = 200,
#'   parms = c(beta = 0.30, gamma = 0.10),
#'   init = list(S = 1e6, I = 20, R = 0),
#'   obs = "poisson"
#' )
#' plot(sim)
#' inc_obs <- sim$incidence$inc
#' plot(seq_along(inc_obs) - 1, inc_obs,
#'      type = "l", xlab = "Day", ylab = "Incidence")
#' fit <- fit_epi_model(inc_obs,
#'                      loss = "logrmse",
#'                      model = SIRS_MODEL,
#'                      init = list(I0 = 6, N = 1e6))
#' init <- tail(sim$states, n = 1)[, -1]
#' pred <- predict(
#'   object = fit,
#'   n_days = 1000,
#'   init = init
#' )
#' plot(pred)
#' plot(sim)
#' plot(pred$states$time, pred$states$I)
#' plot(pred$states$I, col = "red", lty = 2)
#' lines(pred$states$S, col = "blue", lty = 2)
#' }
#'
#' @export
predict.fit_epi_model <- function(object,
                                  n_days,
                                  init = NULL,
                                  times = NULL,
                                  method = "lsoda",
                                  ...) {

  stopifnot(inherits(object, "fit_epi_model"))

  model <- object$model
  parms <- object$par
  x     <- object$x

  # 1) time grids
  if (is.null(times)) {
    times_pred <- 0:(n_days - 1)
  } else {
    times_pred <- times
  }

  times_obs <- 0:(length(x) - 1)

  # 2) initial state
  if (is.null(init)) {
    init <- object$ini
  }

  init <- init[model$state_names]
  if (any(is.na(init))) {
    stop("`init` must contain all model state variables.")
  }

  init <- as.numeric(init[model$state_names])
  names(init) <- model$state_names


  if (any(!is.finite(init))) {
    stop("`init` contains non-finite values.")
  }

  if (any(init < 0)) {
    stop("`init` must be non-negative.")
  }

  # 3) integrate ODE
  out <- deSolve::ode(
    y = init,
    times = times_pred,
    func = model$rhs,
    parms = parms,
    method = method
  )
  out <- as.data.frame(out)

  inc_col <- model$output$incidence_col %||% "incidence"
  if (!inc_col %in% names(out)) {
    stop("ODE output does not contain incidence column '", inc_col, "'.")
  }

  res <- list(
    model      = model,
    par        = parms,
    x          = x,
    times_obs  = times_obs,
    times_pred = out$time,
    incidence  = out[[inc_col]],
    states     = out[, c("time", model$state_names), drop = FALSE]
  )

  class(res) <- "epi_model_predict"
  res
}


#' Plot predictions from a fitted epidemic model
#'
#' @name plot.epi_model_predict
#'
#' @description
#' Plot method for objects of class \code{"epi_model_predict"}.
#' The function displays the observed incidence data used for fitting together
#' with the model-predicted incidence trajectory.
#'
#' @param x An object of class \code{"epi_model_predict"}.
#' @param type_obs Plot type for observed data (default: \code{"h"}).
#' @param col_obs Color for observed incidence.
#' @param col_pred Color for predicted incidence.
#' @param lwd_pred Line width for predicted incidence.
#' @param ... Further graphical parameters passed to \code{\link{plot}}.
#'
#' @return
#' Invisibly returns the input object \code{x}.
#'
#' @export
plot.epi_model_predict <- function(x,
                                   type_obs = "h",
                                   col_obs = "black",
                                   col_pred = "red",
                                   lwd_pred = 2,
                                   ...) {

  stopifnot(inherits(x, "epi_model_predict"))

  xlim <- range(c(x$times_obs, x$times_pred))

  # observed
  plot(x$times_obs, x$x,
       type = type_obs,
       col = col_obs,
       xlab = "Time",
       ylab = "Incidence",
       main = paste("Observed vs predicted incidence –", x$model$name),
       xlim = xlim)

  # predicted
  incidencia <- diff(x$states$C)
  lines(x$times_pred, x$incidence,
        col = col_pred,
        lwd = lwd_pred)

  legend("topright",
         legend = c("Observed", "Predicted"),
         col    = c(col_obs, col_pred),
         lty    = c(1, 1),
         lwd    = c(1, lwd_pred),
         bty    = "n")

  invisible(x)
}
