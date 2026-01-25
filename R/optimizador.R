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
  #'   \item{target}{Model output used for fitting (e.g. \code{"incidence"}, \code{"I"}).}
  #'   \item{value}{Final loss value.}
  #'   \item{optim}{Object returned by \code{stats::optim()}.}
  #'   \item{convergence}{Convergence code from \code{optim()}.}
  #'   \item{message}{Optional convergence message.}
  #'   \item{init}{Initial condition arguments passed by the user.}
  #'   \item{ini0}{Initial state vector used internally.}
  #'   \item{x}{Observed data used for fitting.}
  #' }
#'
#' @seealso
#' \code{\link{fit_epi_model}}, \code{\link{predict.fit_epi_model}}
#'
new_fit_epi_model <- function(model,
                              par,
                              optim,
                              loss,
                              target,
                              init,
                              ini0,
                              x,
                              ode_control = NULL) {

  stopifnot(
    inherits(model, "epi_model"),
    is.character(target),
    length(target) == 1,
    target %in% model$outputs
  )

  structure(
    list(
      model       = model,
      par         = par,
      loss        = loss,
      target      = target,
      value       = optim$value,
      optim       = optim,
      convergence = optim$convergence,
      message     = optim$message,
      init        = init,
      ini0        = ini0,
      x           = x,
      ode_control = ode_control
    ),
    class = "fit_epi_model"
  )
}


objective_logrmse <- function(model,
                              x,
                              ini0,
                              times,
                              target = "incidence",
                              eps = 1e-8,
                              ode_control = list()) {

  stopifnot(is.character(target), length(target) == 1)

  function(log_theta) {

    ## 0) parámetros
    theta <- exp(log_theta)
    names(theta) <- model$par_names

    ## 1) integrar ODE
    out <- try(
      do.call(
        deSolve::ode,
        c(
          list(
            y     = ini0,
            times = times,
            func  = model$rhs,
            parms = theta
          ),
          ode_control
        )
      ),
      silent = TRUE
    )

    if (inherits(out, "try-error")) return(1e30)

    out <- as.data.frame(out)

    ## 2) extraer observable
    if (!target %in% names(out)) return(1e30)
    mu <- out[[target]]

    ## 3) validaciones duras
    if (length(mu) != length(x)) return(1e30)
    if (any(!is.finite(mu)) || any(mu < 0)) return(1e30)

    ## 4) loss
    val <- sqrt(mean((log1p(x + eps) - log1p(mu + eps))^2))
    if (!is.finite(val)) return(1e30)

    val
  }
}

#' Fit an epidemic model by trajectory matching
#'
#' @name fit_epi_model
#'
#' @description
#' Fits a deterministic compartmental epidemic model to observed data by
#' minimizing a trajectory-matching loss (RMSE or log-RMSE) between the observed
#' time series and the corresponding model-predicted trajectory.
#'
#' The observable used for fitting is selected via the \code{target} argument
#' and must correspond to one of the outputs declared by the underlying
#' \code{\link{epi_model}} (e.g. \code{"incidence"}, \code{"I"}, \code{"S"}).
#'
#' No probabilistic observation model is assumed. Parameter estimates should be
#' interpreted as providing the best deterministic approximation to the observed
#' trajectory under the chosen loss.
#'
#' @param x Numeric vector of observed data used for fitting.
#' @param model An object of class \code{"epi_model"} defining the epidemic model.
#' @param loss Loss function to minimize. One of \code{"logrmse"} (default)
#'   or \code{"rmse"}.
#' @param init Named list or vector specifying the initial state of the model.
#'   Names must match \code{model$state_names}.
#' @param target Character string specifying which model output is matched to
#'   \code{x}. Must be one of \code{model$outputs} (e.g. \code{"incidence"},
#'   \code{"I"}, \code{"S"}).
#' @param n_starts Integer. Number of random multi-start initializations used to
#'   reduce sensitivity to local minima.
#' @param control Control list passed to \code{\link{optim}} for the final
#'   optimization step (e.g. \code{maxit}).
#' @param seed Optional integer. If provided, sets the random seed for
#'   reproducible multi-start initialization.
#' @param eps Small positive constant used for numerical stability in the
#'   log-RMSE loss.
#' @param verbose Logical. If \code{TRUE}, progress information from the
#'   multi-start and final optimization phases is printed.
#' @param optim_method Character string specifying the optimization algorithm
#'   passed to \code{\link{optim}}. Defaults to \code{"L-BFGS-B"}.
#'   Methods that do not support box constraints (e.g. \code{"Nelder-Mead"})
#'   are handled automatically.
#' @param ... Additional arguments passed directly to \code{\link[deSolve]{ode}}
#'   for numerical integration of the ODE system (e.g. \code{method},
#'   \code{rtol}, \code{atol}).
#'
#' @return
#' An object of class \code{"fit_epi_model"} containing the fitted parameters,
#' optimization diagnostics, the selected fitting target, and the ODE control
#' settings used during fitting.
#'
#' @examples
#' \dontrun{
#' ## ------------------------------------------------------------
#' ## Example 1: Fit incidence using log-RMSE
#' ## ------------------------------------------------------------
#' sim <- simulate_epi(
#'   model = SIR_MODEL,
#'   times = 0:200,
#'   parms = SIR_MODEL$defaults,
#'   init  = list(S = 1e6, I = 20, R = 0),
#'   seed  = 22
#' )
#'
#' x_inc <- sim$incidence$inc
#'
#' fit_inc <- fit_epi_model(
#'   x = x_inc,
#'   model = SIR_MODEL,
#'   target = "incidence",
#'   loss = "logrmse",
#'   init = list(S = 1e6, I = 6, R = 0)
#' )
#'
#' print(fit_inc)
#'
#' ## ------------------------------------------------------------
#' ## Example 2: Fit a state variable using a different optimizer
#' ## ------------------------------------------------------------
#' fit_I <- fit_epi_model(
#'   x = sim$states$I,
#'   model = SIR_MODEL,
#'   target = "I",
#'   init = list(S = 1e6, I = 10, R = 0),
#'   optim_method = "Nelder-Mead",
#'   method = "rk4",
#'   rtol = 1e-8,
#'   atol = 1e-10
#' )
#'
#' print(fit_I)
#' }
#'
#' @seealso
#' \code{\link{simulate_epi}}, \code{\link{predict.fit_epi_model}}
#'
#' @export
fit_epi_model <- function(x,
                          model = SIR_MODEL,
                          loss = c("logrmse", "rmse"),
                          init = NULL,
                          target = "incidence",
                          n_starts = 100,
                          control = list(maxit = 500),
                          seed = NULL,
                          eps = 1e-6,
                          verbose = TRUE,
                          optim_method = "L-BFGS-B",
                          ...) {

  stopifnot(inherits(model, "epi_model"))
  loss <- match.arg(loss)

  ## ---------------------------------------------------------------------------
  ## Validación de target
  ## ---------------------------------------------------------------------------
  if (!is.character(target) || length(target) != 1 || !target %in% model$outputs) {
    stop(
      "`target` must be one of: ",
      paste(model$outputs, collapse = ", "),
      call. = FALSE
    )
  }

  if (is.null(init)) {
    stop("Provide `init`.", call. = FALSE)
  }

  ## ---------------------------------------------------------------------------
  ## Estado inicial
  ## ---------------------------------------------------------------------------
  ini0 <- unlist(init)
  ini0 <- ini0[model$state_names]

  if (any(!is.finite(ini0)) || any(ini0 < 0)) {
    stop("`init` must contain finite, non-negative values.", call. = FALSE)
  }

  times <- seq_along(x) - 1

  ## ---------------------------------------------------------------------------
  ## Capturar control del integrador ODE
  ## ---------------------------------------------------------------------------
  ode_control <- list(...)

  ## ---------------------------------------------------------------------------
  ## Función objetivo
  ## ---------------------------------------------------------------------------
  fn <- objective_logrmse(
    model       = model,
    x           = x,
    ini0        = ini0,
    times       = times,
    target      = target,
    eps         = eps,
    ode_control = ode_control
  )

  ## ---------------------------------------------------------------------------
  ## Bounds en escala log
  ## ---------------------------------------------------------------------------
  lower <- log(model$lower[model$par_names])
  upper <- log(model$upper[model$par_names])

  uses_bounds <- optim_method %in% c("L-BFGS-B", "Brent")

  if (!is.null(seed)) set.seed(seed)

  ## ===========================================================================
  ## 1) MULTI-START
  ## ===========================================================================
  control_ms <- list(maxit = min(50L, control$maxit %||% 50L))
  best <- list(value = Inf, par = NULL)

  if (verbose) {
    message("Multi-start optimization\n------------------------")
  }

  for (i in seq_len(n_starts)) {

    ## Inicialización coherente con el optimizador
    if (uses_bounds) {
      par0 <- runif(length(lower), lower, upper)
    } else {
      par0 <- rnorm(length(lower), mean = 0, sd = 1)
    }

    ## Construir llamada a optim() correctamente
    opt_args <- list(
      par     = par0,
      fn      = fn,
      method  = optim_method,
      control = control_ms
    )

    if (uses_bounds) {
      opt_args$lower <- lower
      opt_args$upper <- upper
    }

    opt_i <- do.call(optim, opt_args)

    if (is.finite(opt_i$value) && opt_i$value < best$value) {
      best <- list(value = opt_i$value, par = opt_i$par)
      improved <- TRUE
    } else {
      improved <- FALSE
    }

    if (verbose) {
      message(sprintf(
        "  start %3d/%d | value = %.4g%s",
        i, n_starts, opt_i$value,
        if (improved) "  *best*" else ""
      ))
    }
  }

  if (!is.finite(best$value)) {
    stop(
      "Optimization failed: ODE integration failed for all starting points.",
      call. = FALSE
    )
  }

  ## ===========================================================================
  ## 2) OPTIMIZACIÓN FINAL
  ## ===========================================================================
  if (verbose) {
    message("\nFinal optimization\n------------------")
    message(sprintf("Starting value: %.4g", best$value))
  }

  opt_args <- list(
    par     = best$par,
    fn      = fn,
    method  = optim_method,
    control = control
  )

  if (uses_bounds) {
    opt_args$lower <- lower
    opt_args$upper <- upper
  }

  opt <- do.call(optim, opt_args)

  if (verbose) {
    message(sprintf("Final value: %.4g", opt$value))
    message("Convergence: ", opt$convergence)
    if (!is.null(opt$message) && nzchar(opt$message)) {
      message("Message: ", opt$message)
    }
  }

  par_hat <- setNames(exp(opt$par), model$par_names)

  ## ---------------------------------------------------------------------------
  ## Construir objeto resultado
  ## ---------------------------------------------------------------------------
  new_fit_epi_model(
    model       = model,
    par         = par_hat,
    optim       = opt,
    loss        = loss,
    target      = target,
    init        = init,
    ini0        = ini0,
    x           = x,
    ode_control = ode_control
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
  cat("  Model:  ", x$model$name, "\n", sep = "")
  cat("  Target: ", x$target, "\n", sep = "")
  cat("  Loss:   ", x$loss, "\n", sep = "")
  cat("  Value:  ", signif(x$value, 4), "\n", sep = "")
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
