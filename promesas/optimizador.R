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
    mu <- out[[target]][-1]

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
#' The fitting procedure is **model-agnostic**: all epidemiological structure,
#' parameter definitions, and (optional) parameter bounds are taken directly
#' from the supplied \code{\link{epi_model}} object.
#'
#' Parameter estimation is performed on the **log-scale**, ensuring positivity
#' of all model parameters. If the model defines lower and/or upper parameter
#' bounds, these are respected automatically during optimization. If no bounds
#' are defined, parameters are optimized over the entire positive real line.
#'
#' @details
#' ## Fitting strategy
#'
#' Model parameters are estimated by **trajectory matching**: the epidemic model
#' is solved as a system of ordinary differential equations (ODEs), and the
#' resulting model output specified by \code{target} is compared to the observed
#' data \code{x} using the selected loss function.
#'
#' The optimization is carried out in two stages:
#' \enumerate{
#'   \item A **multi-start phase**, in which the optimizer is run repeatedly from
#'   different initial parameter values to reduce sensitivity to local minima.
#'   \item A **final optimization phase**, initialized at the best solution found
#'   during the multi-start phase.
#' }
#'
#' By default, optimization is performed using \code{"L-BFGS-B"}, which allows
#' box constraints when parameter bounds are defined by the model. Models without
#' bounds are optimized over an unconstrained log-parameter space.
#'
#' ## Parameter bounds
#'
#' Parameter bounds are optional and must be defined as part of the
#' \code{\link{epi_model}} object via its \code{lower} and \code{upper} fields.
#'
#' \itemize{
#'   \item If bounds are provided, they are interpreted as **admissible
#'   epidemiological ranges** and enforced during optimization.
#'   \item If bounds are not provided, parameters are assumed to be strictly
#'   positive and no upper or lower limits are imposed beyond positivity.
#' }
#'
#' Bounds are **not** specified in \code{fit_epi_model()} itself; this ensures
#' full reproducibility and traceability of fitted models.
#'
#' ## Numerical integration
#'
#' The epidemic model is solved using \code{\link[deSolve]{ode}}. Additional
#' arguments passed via \code{...} are forwarded directly to the ODE solver,
#' allowing full control over numerical integration (e.g. \code{method},
#' \code{rtol}, \code{atol}).
#'
#' @param x Numeric vector of observed data used for fitting.
#'
#' @param model An object of class \code{"epi_model"} defining the epidemic model
#'   to be fitted.
#'
#' @param loss Character string specifying the loss function to minimize.
#'   One of \code{"logrmse"} (default) or \code{"rmse"}.
#'
#' @param init Named numeric vector or list specifying the initial state of the
#'   model. Names must match \code{model$state_names}.
#'
#' @param target Character string specifying which model output is matched to
#'   \code{x}. Must be one of \code{model$outputs} (e.g. \code{"incidence"},
#'   \code{"I"}, \code{"S"}).
#'
#' @param n_starts Integer. Number of random initializations used in the
#'   multi-start optimization phase.
#'
#' @param control List of control parameters passed to
#'   \code{\link{stats::optim}} for the final optimization step (e.g.
#'   \code{maxit}).
#'
#' @param seed Optional integer. If provided, sets the random seed for
#'   reproducible multi-start initialization.
#'
#' @param eps Small positive constant used for numerical stability in the
#'   log-RMSE loss.
#'
#' @param verbose Logical. If \code{TRUE}, progress information from the
#'   multi-start and final optimization phases is printed.
#'
#' @param optim_method Character string specifying the optimization algorithm
#'   passed to \code{\link{stats::optim}}. Defaults to \code{"L-BFGS-B"}.
#'
#' @param ... Additional arguments passed directly to
#'   \code{\link[deSolve]{ode}} to control numerical integration of the ODE system.
#'
#' @return
#' An object of class \code{"fit_epi_model"} containing the fitted model,
#' estimated parameters, optimization diagnostics, the fitting target, and
#' the numerical integration settings used during fitting.
#'
#' The returned object preserves full reproducibility and traceability between
#' the fitted parameters, the epidemic model definition, and subsequent
#' predictions obtained via \code{\link{predict.fit_epi_model}}.
#'
#' @examples
#' \dontrun{
#' ## ============================================================
#' ## Example 1: Fit a SIR model to simulated incidence data
#' ##            (model WITH parameter bounds)
#' ## ============================================================
#'
#' ## Simulate a SIR epidemic
#' sim <- simulate_epi(
#'   model = SIR_MODEL,
#'   times = 0:150,
#'   parms = c(beta = 0.50, gamma = 0.25),
#'   init  = c(S = 1e6 - 10, I = 10, R = 0),
#'   obs   = "poisson",
#'   seed  = 123
#' )
#' plot(sim)
#' plot(sim$incidence$inc)
#' ## Fit the model to observed incidence
#' fit_inc <- fit_epi_model(
#'   x      = sim$incidence$inc,
#'   model  = SIR_MODEL,
#'   target = "I",
#'   init   = c(S = 1e6 - 5, I = 10, R = 0),
#'   control = list(maxit = 5000),
#' )
#'
#' fit_inc
#'
#'
#' ## ============================================================
#' ## Example 2: Fit a state variable instead of incidence
#' ## ============================================================
#'
#' ## Fit to the infectious compartment I(t)
#' fit_I <- fit_epi_model(
#'   x      = sim$states$I,
#'   model  = SIR_MODEL,
#'   target = "I",
#'   init   = c(S = 1e6 - 5, I = 5, R = 0)
#' )
#'
#' fit_I
#'
#'
#' ## ============================================================
#' ## Example 3: Fit an SI model (model WITHOUT parameter bounds)
#' ## ============================================================
#'
#' ## Simulate an SI epidemic
#' sim_si <- simulate_epi(
#'   model = SI_MODEL,
#'   times = 0:100,
#'   parms = c(beta = 0.4),
#'   init  = SI_MODEL$init,
#'   obs   = "negbin",
#'   seed  = 42
#' )
#'
#' ## Fit the SI model to observed incidence
#' ## Note: beta is optimized over the entire positive real line
#' fit_si <- fit_epi_model(
#'   x      = sim_si$incidence$inc,
#'   model  = SI_MODEL,
#'   target = "incidence",
#'   init   = SI_MODEL$init
#' )
#'
#' fit_si
#'
#'
#' ## ============================================================
#' ## Example 4: Using a different optimizer and ODE solver
#' ## ============================================================
#'
#' fit_custom <- fit_epi_model(
#'   x      = sim$incidence$inc,
#'   model  = SIR_MODEL,
#'   target = "incidence",
#'   init   = c(S = 1e6 - 5, I = 5, R = 0),
#'   optim_method = "Nelder-Mead",
#'   method = "rk4",
#'   rtol   = 1e-8,
#'   atol   = 1e-10
#' )
#'
#' fit_custom
#' }

#'
#' @seealso
#' \code{\link{epi_model}},
#' \code{\link{new_epi_model}},
#' \code{\link{simulate_epi}},
#' \code{\link{predict.fit_epi_model}}
#'
#'
#'
#' @export
fit_epi_model <- function(x,
                          model = SIR_MODEL,
                          loss = c("logrmse", "rmse"),
                          init = NULL,
                          par_init = NULL,
                          target = "incidence",
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

  times <- seq_len(length(x) + 1) - 1

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
  ## Bounds en escala log (opcionales)
  ## ---------------------------------------------------------------------------
  n_par <- length(model$par_names)

  if (is.null(model$lower)) {
    lower <- rep(-Inf, n_par)
  } else {
    lower <- log(model$lower[model$par_names])
  }

  if (is.null(model$upper)) {
    upper <- rep( Inf, n_par)
  } else {
    upper <- log(model$upper[model$par_names])
  }

  uses_bounds <- optim_method %in% c("L-BFGS-B", "Brent")

  if (!is.null(seed)) set.seed(seed)

  ## ---------------------------------------------------------------------------
  ## Inicialización del optimizador
  ## par_init > model$defaults > fallback
  ## ---------------------------------------------------------------------------
  if (!is.null(par_init)) {

    stopifnot(
      is.numeric(par_init),
      !is.null(names(par_init))
    )

    missing <- setdiff(model$par_names, names(par_init))
    if (length(missing) > 0) {
      stop(
        "Missing initial parameters: ",
        paste(missing, collapse = ", "),
        call. = FALSE
      )
    }

    par0 <- log(par_init[model$par_names])

  } else if (!is.null(model$defaults)) {

    par0 <- log(model$defaults[model$par_names])

  } else {

    par0 <- rep(0, n_par)  # equivale a parámetros = 1
  }

  ## ---------------------------------------------------------------------------
  ## Optimización
  ## ---------------------------------------------------------------------------
  if (verbose) {
    message("Optimization\n------------")
  }

  opt_args <- list(
    par     = par0,
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

  ## ---------------------------------------------------------------------------
  ## Parámetros estimados (escala original)
  ## ---------------------------------------------------------------------------
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
