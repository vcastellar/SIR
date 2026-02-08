#' Plot a simulated epidemic
#'
#' @name plot.sim_epi
#'
#' @description
#' Plot method for objects of class \code{"sim_epi"} as returned by
#' \code{\link{simulate_epi}}. The function provides a flexible visualization
#' interface for simulated epidemic trajectories by allowing the user to plot
#' either all model states, the observed incidence time series, or a single
#' model-defined output.
#'
#' Unlike earlier versions, this method does not assume a specific
#' compartmental structure (e.g. SIR or SEIR). All available quantities that can
#' be plotted are taken directly from the outputs declared in the underlying
#' \code{\link{epi_model}} used to generate the simulation.
#'
#' @details
#' The behavior of the function depends on the value of the \code{what} argument:
#'
#' \itemize{
#'   \item{\code{what = "states"}}{Plots the time evolution of all state variables
#'     defined in \code{x$model$states}, using the simulated trajectories
#'     stored in \code{x$states}. The time variable is taken from the
#'     \code{time} column of \code{x$states} and is not considered a state.}
#'
#'   \item{\code{what = "incidence"}}{Plots the observed incidence counts generated
#'     by the observation model and stored in \code{x$incidence}. If no observation
#'     model was specified during simulation (i.e. \code{obs = "none"}), this
#'     option results in an error.}
#'
#'   \item{\code{what = <output>}}{Plots a single model-defined output, where
#'     \code{<output>} is the name of one of the quantities declared in
#'     \code{x$model$outputs} (e.g. \code{"I"}, \code{"S"}, \code{"R"},
#'     or \code{"incidence"}). The corresponding trajectory is extracted from
#'     \code{x$states}.}
#' }
#'
#' If available, the time unit stored in the \code{sim_epi} object
#' (i.e. \code{x$time_unit}) is used to label the time axis. This affects only the
#' plot labels and does not change the numerical values of time.
#'
#' @param x Object of class \code{"sim_epi"} as returned by
#'   \code{\link{simulate_epi}}.
#' @param what Character string specifying what to plot. One of:
#'   \code{"states"} (all state variables),
#'   \code{"incidence"} (observed incidence counts),
#'   or the name of a single model output declared in
#'   \code{x$model$outputs} (e.g. \code{"I"}, \code{"S"}).
#' @param scale Character string specifying the scale for state plots.
#'   One of \code{"auto"}, \code{"full"}, \code{"small"}, or \code{"log"}.
#'   This argument is only used when \code{what = "states"}.
#' @param ... Additional graphical parameters passed to base plotting
#'   functions such as \code{\link{plot}} and \code{\link{matplot}}.
#'
#' @return
#' Invisibly returns the input object \code{x}.
#'
#' @seealso
#' \code{\link{simulate_epi}}, \code{\link{summary.sim_epi}},
#' \code{\link{print.sim_epi}}
#'
#' @examples
#' sird_rhs <- function(time, state, parms) {
#'   with(as.list(c(state, parms)), {
#'     N <- S + I + R
#'     lambda <- beta * S * I / N
#'
#'     dS <- -lambda
#'     dI <-  lambda - gamma * I - mu * I
#'     dR <-  gamma * I
#'     dD <-  mu * I
#'
#'     list(
#'       c(dS, dI, dR, dD),
#'       incidence = lambda
#'     )
#'   })
#' }
#'
#' SIRD_MODEL <- epi_model(
#'   name      = "SIRD",
#'   rhs       = sird_rhs,
#'   states    = c("S", "I", "R", "D"),
#'   flows     = c("incidence"),
#'   par_names = c("beta", "gamma", "mu"),
#'   roles     = list(
#'     susceptible = "S",
#'     infectious  = "I",
#'     recovered   = "R",
#'     deceased    = "D"
#'   )
#' )
#'
#' # Plot observed incidence (requires an observation model)
#' simulate_epi(model = SIRD_MODEL,
#'              times = 0:200)
#'
#' # Plot a single state variable
#' plot(sim, what = "I")
#'
#' # Plot all flows
#' plot(sim, what = "flows")
#'
#'
#' # If the model defines Rt as a flow, it can be plotted directly
#' plot(sim, what = "Rt")
#'
#' @export
plot.sim_epi <- function(x,
                         what = "states",
                         scale = c("auto", "full", "small", "log"),
                         ...) {

  scale <- match.arg(scale)

  if (!inherits(x, "sim_epi")) {
    stop("Object must be of class 'sim_epi'.")
  }

  if (is.null(x$model) || !inherits(x$model, "epi_model")) {
    stop("Invalid 'sim_epi' object: missing epi_model in x$model.")
  }

  if (is.null(x$states) || !"time" %in% names(x$states)) {
    stop("Invalid 'sim_epi' object: missing $states with column 'time'.")
  }

  model   <- x$model
  st      <- x$states
  states  <- model$states
  outputs <- model$outputs

  ## Etiqueta del eje X según time_unit
  unit <- x$time_unit
  xlab <- if (is.null(unit) || !nzchar(unit)) {
    "Time"
  } else {
    paste0("Time (", unit, ")")
  }

  ## ===========================================================================
  ## CASE 1: all states
  ## ===========================================================================
  if (identical(what, "states")) {

    missing <- setdiff(states, names(st))
    if (length(missing) > 0) {
      stop(
        "States missing in x$states: ",
        paste(missing, collapse = ", ")
      )
    }

    y <- st[, states, drop = FALSE]

    if (scale == "auto") {
      max_vals <- apply(y, 2, max, na.rm = TRUE)
      ratio <- max(max_vals) / max(min(max_vals[max_vals > 0]), 1)
      scale <- if (is.finite(ratio) && ratio > 100) "small" else "full"
    }

    if (scale == "small") {
      ylim <- range(y, finite = TRUE)
      ylab <- "States (zoomed)"
      log_arg <- ""
    } else if (scale == "log") {
      y[y <= 0] <- NA_real_
      ylim <- range(y, finite = TRUE)
      ylab <- "States (log scale)"
      log_arg <- "y"
    } else {
      ylim <- c(0, max(y, na.rm = TRUE))
      ylab <- "State values"
      log_arg <- ""
    }

    matplot(
      st$time, y,
      type = "l",
      lty = 1,
      lwd = 2,
      xlab = xlab,
      ylab = ylab,
      main = paste("Simulation:", model$name),
      ylim = ylim,
      log = log_arg,
      ...
    )

    legend(
      "topright",
      legend = states,
      col = seq_along(states),
      lty = 1,
      lwd = 2,
      bty = "n"
    )

    ## ===========================================================================
    ## CASE 2: all flows
    ## ===========================================================================
  } else if (identical(what, "flows")) {

    fl <- x$flows

    if (is.null(fl) || !"time" %in% names(fl)) {
      stop("Invalid 'sim_epi' object: missing $flows with column 'time'.")
    }

    vars <- setdiff(names(fl), "time")

    if (length(vars) == 0) {
      stop("No flows to plot.")
    }

    y <- fl[, vars, drop = FALSE]

    matplot(
      fl$time, y,
      type = "l",
      lty = 1,
      lwd = 2,
      xlab = xlab,
      ylab = "Flows",
      main = paste("Flows:", model$name),
      ...
    )

    legend(
      "topright",
      legend = vars,
      col = seq_along(vars),
      lty = 1,
      lwd = 2,
      bty = "n"
    )

    ## ===========================================================================
    ## CASE 3: single state
    ## ===========================================================================
  } else if (is.character(what) &&
             length(what) == 1 &&
             what %in% states) {

    if (!what %in% names(st)) {
      stop(
        "State '", what, "' not found in x$states."
      )
    }

    y <- st[[what]]

    plot(
      st$time, y,
      type = "l",
      lwd = 2,
      xlab = xlab,
      ylab = what,
      main = paste("Simulation:", model$name, "-", what),
      ...
    )

    ## ===========================================================================
    ## CASE 4: single flow
    ## ===========================================================================
  } else if (is.character(what) &&
             length(what) == 1 &&
             !is.null(x$flows) &&
             what %in% names(x$flows)) {

    fl <- x$flows

    plot(
      fl$time, fl[[what]],
      type = "l",
      lwd = 2,
      xlab = xlab,
      ylab = what,
      main = paste("Flow:", model$name, "-", what),
      ...
    )

  } else {
    stop(
      "`what` must be 'states', 'flows', or one of: ",
      paste(c(states, setdiff(names(x$flows), "time")), collapse = ", "),
      call. = FALSE
    )
  }

  invisible(x)
}





#' Summarize a simulated epidemic
#'
#' @name summary.sim_epi
#'
#' @description
#' Summary method for objects of class \code{"sim_epi"} as returned by
#' \code{\link{simulate_epi}}. The function computes and returns a small set of
#' epidemiologically meaningful summary quantities derived from the simulation.
#'
#' @details
#' The summary includes:
#' \describe{
#'   \item{model}{The epidemic model simulated (e.g. \code{"sir"} or \code{"sirs"}).}
#'   \item{R0}{The basic reproduction number, computed as
#'     \eqn{R_0 = \beta / \gamma}.}
#'   \item{peak_I}{The maximum number of infectious individuals observed during
#'     the simulation.}
#'   \item{total_infections}{The total number of infection events, computed as
#'     the maximum value of the cumulative infection variable \code{C(t)}.}
#' }
#'
#' This method is automatically dispatched when calling \code{summary()} on an
#' object of class \code{"sim_epi"}.
#'
#' @param object Object of class \code{"sim_epi"} as returned by
#'   \code{\link{simulate_epi}}.
#' @param ... Currently unused; included for compatibility with generic
#'   \code{\link{summary}}.
#'
#' @return
#' A named list with summary statistics describing the simulated epidemic.
#'
#' @seealso
#' \code{\link{simulate_epi}}, \code{\link{plot.sim_epi}}
#'
#' @examples
#' sim <- simulate_epi(n_days = 300, model = SIRS_MODEL, omega = 1/180, seed = 1)
#'
#' summary(sim)
#'
#' @export
summary.sim_epi <- function(object, ...) {

  st  <- object$states
  roles <- object$model$roles

  res <- list(
    model = object$model
  )

  if ("infectious" %in% names(roles)) {
    inf <- get_role(object, "infectious")
    res$peak_I <- max(inf, na.rm = TRUE)
    res$time_peak_I <- st$time[which.max(inf)]
  } else if ("I" %in% names(st)) {
    res$peak_I <- max(st$I, na.rm = TRUE)
    res$time_peak_I <- st$time[which.max(st$I)]
  }

  if ("incidence" %in% names(roles)) {
    inc <- get_role(object, "incidence")
    res$total_infections <- sum(inc, na.rm = TRUE)
  }

  res
}



#' Print a simulated epidemic
#'
#' @name print.sim_epi
#'
#' @description
#' Print method for objects of class \code{"sim_epi"} as returned by
#' \code{\link{simulate_epi}}. The function provides a concise, human-readable
#' summary of a simulated epidemic, including the underlying epidemic model,
#' the simulation time horizon, the parameter values used, and selected outcome
#' metrics derived from the simulation results.
#'
#' Unlike earlier versions, this method does not assume a specific compartmental
#' structure (such as SIR or SEIR). Model-specific information is obtained directly
#' from the \code{\link{epi_model}} object stored in \code{x$model}.
#'
#' @details
#' This method is automatically dispatched when an object of class
#' \code{"sim_epi"} is printed. It is intended to give a quick overview of the
#' simulation without displaying the full internal structure of the object.
#'
#' If available, the time unit stored in the \code{sim_epi} object
#' (i.e. \code{x$time_unit}) is used to label the simulation horizon and
#' reported event times. This affects only the printed labels and does not
#' change the numerical values of time.
#'
#' The printed output typically includes:
#' \itemize{
#'   \item the name of the epidemic model,
#'   \item the simulated time horizon,
#'   \item the parameter values used in the simulation,
#'   \item basic outcome summaries (e.g. peak prevalence or total infections),
#'     when such quantities are available.
#' }
#'
#' @param x Object of class \code{"sim_epi"} as returned by
#'   \code{\link{simulate_epi}}.
#' @param ... Additional arguments ignored by this method; included for
#'   compatibility with the generic \code{\link{print}} function.
#'
#' @return
#' Invisibly returns the input object \code{x}.
#'
#' @seealso
#' \code{\link{simulate_epi}}, \code{\link{summary.sim_epi}},
#' \code{\link{plot.sim_epi}}
#'
#' @examples
#' sim <- simulate_epi(
#'   model = SIR_MODEL,
#'   times = 0:200,
#'   time_unit = "days",
#'   parms = c(beta = 0.30, gamma = 0.10),
#'   init  = c(S = 999990, I = 10, R = 0, C = 10),
#'   seed  = 1
#' )
#'
#' sim
#'
#' @export
print.sim_epi <- function(x, ...) {

  stopifnot(inherits(x, "sim_epi"))

  model  <- x$model
  params <- x$params
  states <- x$states

  # Etiqueta de unidad de tiempo
  unit <- x$time_unit
  unit_lbl <- if (is.null(unit) || !nzchar(unit)) {
    "time units"
  } else {
    unit
  }

  cat("Epidemic simulation\n")
  cat("-------------------\n")

  cat("Model:            ", model$name, "\n", sep = "")
  cat("Time horizon:     ",
      max(states$time, na.rm = TRUE), " ", unit_lbl, "\n", sep = "")

  if ("N" %in% names(params)) {
    cat("Population (N):   ",
        format(params$N, scientific = TRUE), "\n", sep = "")
  }

  cat("\nParameters\n")
  for (nm in names(params)) {
    cat("  ", nm, ": ", params[[nm]], "\n", sep = "")
  }

  cat("\nOutcomes\n")

  ## ---------------------------------------------------------------------------
  ## Peak incidence (flow)
  ## ---------------------------------------------------------------------------
  if ("incidence" %in% names(x$model$roles)) {
    inc <- get_role(x, "incidence")

    t_inc <- x$flows$time
    peak_inc <- max(inc, na.rm = TRUE)

    if (is.finite(peak_inc)) {
      time_inc <- t_inc[which.max(inc)]
      cat("  Peak incidence:       ", round(peak_inc), "\n", sep = "")
      cat("  Time of inc. peak:    ",
          time_inc, " ", unit_lbl, "\n", sep = "")
    }
  }

  ## ---------------------------------------------------------------------------
  ## Peak infectious (state I)
  ## ---------------------------------------------------------------------------
  if ("infectious" %in% names(x$model$roles)) {
    inf <- get_role(x, "infectious")
    peak_I <- max(inf, na.rm = TRUE)

    if (is.finite(peak_I)) {
      time_I <- states$time[which.max(inf)]
      cat("  Peak infectious:      ", round(peak_I), "\n", sep = "")
      cat("  Time of I peak:       ",
          time_I, " ", unit_lbl, "\n", sep = "")
    }
  } else if ("I" %in% names(states)) {
    peak_I <- max(states$I, na.rm = TRUE)

    if (is.finite(peak_I)) {
      time_I <- states$time[which.max(states$I)]
      cat("  Peak infectious:      ", round(peak_I), "\n", sep = "")
      cat("  Time of I peak:       ",
          time_I, " ", unit_lbl, "\n", sep = "")
    }
  }

  if (!is.null(x$incidence_cum)) {
    total_C <- max(x$incidence_cum$cases_cum, na.rm = TRUE)
    cat("  Total infections:     ", round(total_C), "\n", sep = "")
  }

  invisible(x)
}
