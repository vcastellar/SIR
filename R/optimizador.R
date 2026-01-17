#' Fitted epidemic model (trajectory matching)
#'
#' @name fit_epi_model-class
#'
#' @description
#' Objects of class \code{"fit_epi_model"} represent the result of fitting a
#' deterministic epidemic model to incidence data by minimizing a
#' trajectory-matching loss (RMSE or log-RMSE).
#'
#' No probabilistic observation model is assumed. The fitted parameters should
#' be interpreted as providing the best deterministic approximation to the
#' observed incidence trajectory under the chosen loss.
#'
#' @section Contents:
#' A \code{"fit_epi_model"} object is a named list with components:
#' \describe{
#'   \item{model}{The fitted \code{epi_model} object.}
#'   \item{par}{Named numeric vector of fitted parameters.}
#'   \item{loss}{Loss minimized (\code{"rmse"} or \code{"logrmse"}).}
#'   \item{value}{Final loss value.}
#'   \item{optim}{Object returned by \code{stats::optim()}.}
#'   \item{convergence}{Convergence code from \code{optim()}.}
#'   \item{message}{Optional convergence message.}
#'   \item{init}{Initial condition arguments passed by the user.}
#'   \item{ini0}{Initial state vector used internally.}
#'   \item{x}{Observed incidence data.}
#' }
#'
#' @seealso
#' \code{\link{fit_epi_model}}, \code{\link{predict.fit_epi_model}}
#'
new_fit_epi_model <- function(model, par, optim, loss, init, ini0, x) {
  structure(
    list(
      model = model,
      par = par,
      loss = loss,
      value = optim$value,
      optim = optim,
      convergence = optim$convergence,
      message = optim$message,
      init = init,
      ini0 = ini0,
      x = x
    ),
    class = "fit_epi_model"
  )
}


objective_logrmse <- function(model, x, ini0, times, eps = 1e-8) {

  inc_col <- model$output$incidence_col %||% "incidence"

  function(log_theta) {

    theta <- exp(log_theta)
    names(theta) <- model$par_names

    out <- try(
      deSolve::ode(
        y = ini0,
        times = times,
        func  = model$rhs,
        parms = theta,
        method = "lsoda"
      ),
      silent = TRUE
    )
    if (inherits(out, "try-error")) return(1e30)

    mu <- as.data.frame(out)[[inc_col]]

    if (any(!is.finite(mu)) || any(mu < 0)) return(1e30)

    val <- sqrt(mean((log1p(x) - log1p(mu))^2))

    if (!is.finite(val)) 1e30 else val
  }
}







#' Fit an epidemic model by trajectory matching
#'
#' @name fit_epi_model
#'
#' @description
#' Fits a deterministic epidemic model to incidence data by minimizing a
#' trajectory-matching loss (RMSE or log-RMSE) between observed and simulated
#' incidence.
#'
#' @param x Numeric vector of observed incidence.
#' @param model An \code{epi_model} object.
#' @param loss Loss function: \code{"logrmse"} (default) or \code{"rmse"}.
#' @param init Named list passed to \code{model$make_init()}.
#' @param n_starts Number of random multi-start initializations.
#' @param control Control list passed to \code{optim()}.
#' @param seed Optional random seed.
#'
#' @return An object of class \code{"fit_epi_model"}.
#' @examples
#' \dontrun{
#' ## Simulate a SIRS epidemic and fit the model to observed incidence
#' sim <- simulate_epi(
#'   model = SIR_MODEL,
#'   n_days = 200,
#'   parms = SIR_MODEL$defaults,
#'   init = list(S = 1e6, I = 20, R = 0),
#'   seed = 22
#'  )
#' plot(sim)
#' x <- sim$incidence$inc
#' lines(x)
#' ## Fit with log-RMSE (no likelihood, size ignored)
#' fit_rmse <- fit_epi_model(
#'   x = x,
#'   model = SIR_MODEL,
#'   loss = "logrmse",
#'   init = list(S = 1e6, I= 6, R = 0)
#'   )
#' print(fit_rmse)
#'}
#' @export
fit_epi_model <- function(x,
                           model = SIR_MODEL,
                           loss = c("logrmse", "rmse"),
                           init = NULL,
                           n_starts = 100,
                           control = list(maxit = 500),
                           seed = NULL,
                           eps = 1e-6,
                           verbose = TRUE) {

  stopifnot(inherits(model, "epi_model"))
  loss <- match.arg(loss)

  if (is.null(init)) {
    stop("Model does not define init.")
  }

  ## --- initial state ----------------------------------------------------------
  ini0 <- unlist(init)
  ini0 <- ini0[model$state_names]

  times <- seq_along(x) - 1

  ## --- objective --------------------------------------------------------------
  fn <- objective_logrmse(
    model = model,
    x     = x,
    ini0  = ini0,
    eps   = eps,
    times = times
  )

  ## --- bounds (log-scale) -----------------------------------------------------
  lower <- log(model$lower[model$par_names])
  upper <- log(model$upper[model$par_names])

  if (!is.null(seed)) set.seed(seed)

  ## ========================================================================== ##
  ## 1) CHEAP MULTI-START
  ## ========================================================================== ##

  control_ms <- list(
    maxit = min(50L, control$maxit %||% 50L)
  )

  best <- list(value = Inf, par = NULL)

  if (verbose) {
    message("Multi-start optimization\n")
    message("------------------------\n")
  }

  for (i in seq_len(n_starts)) {

    par0 <- runif(length(lower), min = lower, max = upper)

    opt_i <- optim(
      par = par0,
      fn  = fn,
      method = "L-BFGS-B",
      lower = lower,
      upper = upper,
      control = control_ms
    )

    if (is.finite(opt_i$value) && opt_i$value < best$value) {
      best <- list(value = opt_i$value, par = opt_i$par)
      improved <- TRUE
    } else {
      improved <- FALSE
    }

    if (verbose) {
      message(
        sprintf(
          "  start %3d/%d | value = %.4g%s\n",
          i, n_starts, opt_i$value,
          if (improved) "  *best*" else ""
        )
      )
    }
  }

  if (is.null(best$par)) {
    stop("Multi-start phase failed: no valid optimization result.")
  }

  ## ========================================================================== ##
  ## 2) FINAL OPTIMIZATION
  ## ========================================================================== ##

  if (verbose) {
    message("\nFinal optimization\n")
    message("------------------\n")
    message(sprintf("Starting value: %.4g\n", best$value))
  }

  opt <- optim(
    par = best$par,
    fn  = fn,
    method = "L-BFGS-B",
    lower = lower,
    upper = upper,
    control = control
  )

  if (verbose) {
    message(sprintf("Final value: %.4g\n", opt$value))
    message("Convergence: ", opt$convergence, "\n", sep = "")
    if (!is.null(opt$message) && nzchar(opt$message)) {
      message("Message: ", opt$message, "\n", sep = "")
    }
  }

  par_hat <- setNames(exp(opt$par), model$par_names)

  ## --- build fitted object ----------------------------------------------------
  new_fit_epi_model(
    model = model,
    par = par_hat,
    optim = opt,
    loss = loss,
    init = init,
    ini0 = ini0,
    x = x
  )
}

#' Print a fitted epidemic model
#'
#' @name print.fit_epi_model
#' @description
#' Print method for objects of class \code{"fit_epi_model"}.
#'
#' The function displays a concise, human-readable summary of the result of
#' fitting an epidemic model to incidence data. The printed output includes
#' the model name, the loss function used for fitting, the final objective
#' value, convergence information from the optimizer, and the estimated
#' model parameters.
#'
#' This method is automatically called when an object of class
#' \code{"fit_epi_model"} is printed at the console.
#'
#' @param x An object of class \code{"fit_epi_model"} as returned by
#'   \code{\link{fit_epi_model}}.
#' @param ... Further arguments (ignored).
#'
#' @details
#' The reported objective value corresponds to the value of the loss function
#' evaluated at the estimated parameter vector. Convergence information and
#' optional messages are taken directly from the underlying optimization
#' routine.
#'
#' The printed parameter estimates represent point estimates obtained by
#' minimizing the chosen loss function; no uncertainty quantification is
#' provided by this method.
#'
#' @return
#' Invisibly returns the input object \code{x}.
#'
#' @examples
#' \dontrun{
#' fit <- fit_epi_model(
#'   x = incidence_data,
#'   model = SIR_MODEL
#' )
#'
#' fit
#' }
#'
#' @export
print.fit_epi_model <- function(x, ...) {
  cat("<fit_epi_model>\n")
  cat("  Model: ", x$model$name, "\n", sep = "")
  cat("  Loss:  ", x$loss, "\n", sep = "")
  cat("  Value: ", signif(x$value, 4), "\n", sep = "")
  cat("  Convergence: ", x$convergence, "\n", sep = "")
  if (!is.null(x$message) && nzchar(x$message)) {
    cat("  Message: ", x$message, "\n", sep = "")
  }
  cat("  Parameters:\n")
  for (nm in names(x$par)) {
    cat("    ", nm, ": ", signif(x$par[[nm]], 4), "\n", sep = "")
  }
  invisible(x)
}
