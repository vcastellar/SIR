
#' @keywords internal
#' @noRd
# objective using incidence_col directly (optimizing in log-parameter space)
#' @keywords internal
#' @noRd
# objective using incidence_col directly
# returns a function of log(theta)
objective_epi_incidence <- function(model, x, ini0,
                                    loss = c("mle", "rmse", "logrmse"),
                                    distr = c("poisson", "negbin"),
                                    eps = 1e-8,
                                    size = 10,
                                    times = NULL,
                                    rmse_eps = 0) {

  loss  <- match.arg(loss)
  distr <- match.arg(distr)

  if (is.null(times)) {
    # si x son conteos diarios, esto representa días 0,1,2,...,T-1
    times <- 0:(length(x) - 1)
  } else {
    stopifnot(length(times) == length(x))
  }

  inc_col <- (model$output$incidence_col %||% "incidence")

  # NOTA: la función devuelta espera log(theta)
  function(log_theta) {

    # Validación básica del punto en espacio log
    if (any(!is.finite(log_theta))) return(1e30)

    # Transformación a escala natural
    theta <- exp(log_theta)
    names(theta) <- model$par_names

    out <- try(
      deSolve::ode(
        y = ini0,
        times = times,
        func = model$rhs,
        parms = theta,
        method = "lsoda"
      ),
      silent = TRUE
    )
    if (inherits(out, "try-error")) return(1e30)

    out <- as.data.frame(out)
    if (!inc_col %in% names(out)) return(1e30)

    mu <- out[[inc_col]]

    # basic validity
    if (any(!is.finite(mu))) return(1e30)
    if (any(mu < 0)) return(1e30)

    # --- objective ---
    if (loss == "mle") {

      mu <- pmax(mu, eps)
      mu <- pmin(mu, 1e12)  # anti-overflow

      val <- switch(
        distr,
        poisson = -sum(stats::dpois(x, lambda = mu, log = TRUE)),
        negbin  = -sum(stats::dnbinom(x, mu = mu, size = size, log = TRUE))
      )

    } else if (loss == "rmse") { # loss == "rmse"

      # RMSE (opcionalmente con un pequeño eps si quieres estabilizar cuando mu~0)
      mu2 <- if (rmse_eps > 0) pmax(mu, rmse_eps) else mu
      # Nota: RMSE y SSE tienen el mismo óptimo; RMSE escala mejor para interpretar.
      val <- sqrt(mean((x - mu2)^2))

      # si quieres una versión más suave y típica para optimizadores:
      # val <- mean((x - mu2)^2)  # MSE (sin sqrt)
    } else if (loss == "logrmse") {
        val <- sqrt(mean((log1p(x) - log1p(mu))^2))
        }

    if (!is.finite(val)) 1e30 else val
  }
}




#' @keywords internal
#' @noRd
new_fit_epi_model <- function(model, par, optim, x, init, ini0, distr,
                              best_start = NULL, control = NULL) {
  structure(
    list(
      model = model,
      par = par,
      optim = optim,
      value = optim$value,
      convergence = optim$convergence,
      message = optim$message,
      distr = distr,
      init = init,
      ini0 = ini0,
      x = x,
      x_len = length(x),
      best_start = best_start,
      control = control
    ),
    class = "fit_epi_model"
  )
}


#' Fit an \code{epi_model} to incidence data by likelihood minimization
#' @name fit_epi_model
#' @description
#' Fits a deterministic compartmental epidemic model (an object of class
#' \code{"epi_model"}, such as \code{SIR_MODEL}, \code{SIRS_MODEL}, or \code{SEIR_MODEL})
#' to observed **incidence count data** by minimizing a **negative log-likelihood**.
#'
#' The function returns a fitted model object of class \code{"fit_epi_model"},
#' which contains the estimated parameters together with additional information
#' about the optimization process, the model definition, and the initial
#' conditions used for fitting.
#'
#' @details
#' ## Fitting approach
#' Internally, the function constructs an objective function (not exported) that:
#' \itemize{
#'   \item solves the model ODE system using \code{deSolve::ode()} and \code{model$rhs};
#'   \item extracts the model-implied incidence directly from the ODE output
#'     column specified by \code{model$output$incidence_col} (default: \code{"incidence"});
#'   \item evaluates a Poisson or Negative Binomial log-likelihood for the
#'     observed incidence vector \code{x}.
#' }
#'
#' To improve numerical stability and reduce sensitivity to initial values,
#' the optimization is initialized via an internal **multi-start** strategy
#' (random initial parameter vectors sampled within the bounds defined by the model),
#' and the best multi-start solution is used to initialize the final optimization.
#'
#' ## Data and time grid
#' The input vector \code{x} must represent incidence counts on an equally spaced
#' time grid (typically daily). The internal objective integrates the model on a
#' grid \code{times = 1:length(x)} so the extracted incidence time series aligns
#' with \code{x}.
#'
#' ## Initial conditions
#' Initial state values are constructed via \code{model$make_init()} using the
#' named list \code{init}. The function performs a minimal validation to ensure
#' that \code{init} includes all required (non-default) arguments of
#' \code{model$make_init()} and that the returned state vector contains the
#' required \code{model$state_names}.
#'
#' ## Model requirements
#' The supplied \code{model} must:
#' \itemize{
#'   \item be an object of class \code{"epi_model"};
#'   \item provide an ODE right-hand side function compatible with \code{deSolve::ode()};
#'   \item define parameter names via \code{model$par_names};
#'   \item define parameter bounds via \code{model$lower} and \code{model$upper};
#'   \item define \code{model$output$incidence_col} and return that column in the ODE output;
#'   \item provide \code{model$make_init()} returning a named numeric state vector
#'     containing at least \code{model$state_names}.
#' }
#'
#' @param x Numeric vector of observed incidence counts (e.g. daily cases).
#' @param model An \code{epi_model} object to be fitted (e.g. \code{SIR_MODEL},
#'   \code{SIRS_MODEL}, \code{SEIR_MODEL}). Defaults to \code{SIR_MODEL}.
#' @param distr Character string specifying the observation distribution used in
#'   the likelihood. One of \code{"poisson"} or \code{"negbin"}.
#' @param init Named list of arguments passed to \code{model$make_init()} to build
#'   the initial state vector (e.g. \code{list(N = 1e6, I0 = 10, R0 = 0)} or
#'   \code{list(N = 1e6, E0 = 0, I0 = 10, R0 = 0)} for SEIR).
#' @param size Numeric. Dispersion (size) parameter for the Negative Binomial
#'   likelihood (used when \code{distr = "negbin"}).
#' @param n_starts Integer. Number of random starting points used in the
#'   multi-start initialization.
#' @param control_ms List of control parameters passed to \code{optim()} during
#'   the multi-start phase.
#' @param control List of control parameters passed to \code{optim()} in the final
#'   optimization.
#' @param seed Optional integer seed for reproducible multi-start initialization.
#' @param eps Small positive constant used to avoid log-likelihood issues when
#'   the model-implied mean incidence is zero.
#' @param ... Reserved for future extensions.
#'
#' @return
#' An object of class \code{"fit_epi_model"} with components including:
#' \describe{
#'   \item{par}{Named numeric vector of fitted model parameters.}
#'   \item{model}{The \code{epi_model} object that was fitted.}
#'   \item{optim}{The full object returned by \code{stats::optim()}.}
#'   \item{value}{Final value of the negative log-likelihood.}
#'   \item{convergence}{Convergence code returned by \code{optim()}.}
#'   \item{message}{Optional convergence message from \code{optim()}.}
#'   \item{init}{List of initial condition arguments supplied by the user.}
#'   \item{ini0}{Named numeric vector of initial state values used internally.}
#'   \item{x}{Observed incidence vector used for fitting.}
#' }
#' @examples
#' \dontrun{
#' ## Simulate a SIR epidemic and fit the model to observed incidence
#' sim <- simulate_epi(
#'   model = SIRS_MODEL,
#'   n_days = 200,
#'   parms = SIRS_MODEL$default,
#'   init_args = list(N = 1e6, I0 = 20, R0 = 0),
#'   obs = "negbin",
#'   rho = 1,
#'   seed = 22
#' )
#' plot(sim)
#' x <- sim$incidence_obs$inc
#' plot(x)
#' fit <- fit_epi_model(
#'   x = x,
#'   model = SIRS_MODEL,
#'   loss = "logrmse",
#'   init = list(N = 1e6, I0 = 20, R0 = 0)
#' )
#' print(fit)
#'
#' ## Compare fitted incidence to observed incidence
#' pred <- predict(fit, n_days = length(x), init_args = list(N = 1e6, I0 = 20, R0 = 0))
#' plot(x, type = "l", xlab = "Day", ylab = "Incidence")
#' lines(pred$times_obs, pred$incidence, col = "red", lty = 2)
#' }
#'
#' fit2 <- fit_epi_model(
#'   x = x,
#'   model = SIRS_MODEL,
#'   loss = "mle",
#'   obs = "negbin",
#'   init = list(N = 1e6, I0 = 20, R0 = 0)
#' )
#' print(fit2)
#'
#' @seealso
#' \code{\link{simulate_epi}}, \code{\link{predict.fit_epi_model}},
#' \code{\link{SIR_MODEL}}, \code{\link{SIRS_MODEL}}, \code{\link{SEIR_MODEL}},
#' \code{\link[stats]{optim}}, \code{\link[deSolve]{ode}}
#'
#' @export
fit_epi_model <- function(x,
                          model = SIR_MODEL,
                          loss = c("mle", "rmse", "logrmse"),
                          distr = c("negbin", "poisson"),
                          init = list(N = 1e6, I0 = 10, R0 = 0),
                          size = 10,
                          control = list(maxit = 5000),
                          n_starts = 30,
                          control_ms = list(maxit = 50),
                          seed = 1,
                          eps = 1e-8,
                          ...) {

  stopifnot(inherits(model, "epi_model"))
  distr <- match.arg(distr)

  # 1) make_init required
  if (!is.function(model$make_init)) {
    stop("Model does not define `make_init()`. Please provide `model$make_init` or pass `ini0` externally.")
  }
  if (!is.list(init) || is.null(names(init))) {
    stop("`init` must be a named list with arguments for `model$make_init()`.")
  }

  # 2) minimal validation: required args of make_init() are present
  fml <- formals(model$make_init)
  req <- names(fml)[vapply(fml, function(z) identical(z, quote(expr = )), logical(1))]
  miss <- setdiff(req, names(init))
  if (length(miss) > 0) {
    stop("`init` is missing required argument(s) for `model$make_init()`: ",
         paste(miss, collapse = ", "))
  }

  # 3) build ini0 via make_init and validate state_names
  ini0 <- do.call(model$make_init, init)

  if (!is.numeric(ini0) || is.null(names(ini0))) {
    stop("`model$make_init()` must return a named numeric vector.")
  }
  miss_state <- setdiff(model$state_names, names(ini0))
  if (length(miss_state) > 0) {
    stop("`model$make_init()` output is missing state(s): ",
         paste(miss_state, collapse = ", "))
  }
  ini0 <- ini0[model$state_names]

  # 4) objective (expects log(theta))
  fn <- objective_epi_incidence(
    model = model,
    loss = loss,
    x = x,
    ini0 = ini0,
    distr = distr,
    eps = eps,
    size = size
  )

  # 5) bounds in natural space -> transform to log space
  lower <- unname(model$lower[model$par_names])
  upper <- unname(model$upper[model$par_names])

  if (any(!is.finite(lower)) || any(!is.finite(upper))) {
    stop("Model bounds (lower/upper) must be finite for log-optimization.")
  }
  if (any(lower <= 0)) {
    stop("All lower bounds must be > 0 to optimize in log-space.")
  }
  if (any(lower >= upper)) {
    stop("Invalid bounds: some lower >= upper.")
  }

  lower_log <- log(lower)
  upper_log <- log(upper)

  # 6) multi-start in log-space
  if (!is.null(seed)) set.seed(seed)

  best <- list(value = Inf, opt = NULL)

  for (i in seq_len(n_starts)) {

    # sample uniformly in natural space then log-transform (simple and ok)
    par0 <- stats::runif(length(model$par_names), min = lower, max = upper)
    par0_log <- log(par0)

    opt_i <- optim(
      par = par0_log,
      fn  = fn,
      method = "L-BFGS-B",
      lower = lower_log,
      upper = upper_log,
      control = control_ms
    )

    if (is.finite(opt_i$value) && opt_i$value < best$value) {
      best <- list(value = opt_i$value, opt = opt_i)
    }
  }

  if (is.null(best$opt)) {
    stop("Multi-start failed: no finite objective value found. Try different `init`, bounds, or `control_ms`.")
  }

  # 7) final optimization from best log-start
  opt <- optim(
    par = unname(best$opt$par),
    fn  = fn,
    method = "L-BFGS-B",
    lower = lower_log,
    upper = upper_log,
    control = control
  )

  # 8) back-transform to natural parameter scale
  par_hat <- setNames(exp(opt$par), model$par_names)

  # (opcional pero recomendable) guardar también el óptimo en log-escala
  opt_nat <- opt
  opt_nat$par <- par_hat

  new_fit_epi_model(
    model = model,
    par = par_hat,
    optim = opt_nat,
    x = x,
    init = init,
    ini0 = ini0,
    distr = distr,
    best_start = best,
    control = control
  )
}

#' Print a fitted epidemic model
#'
#' @name print.fit_epi_model
#'
#' @description
#' S3 print method for objects of class \code{"fit_epi_model"} as returned by
#' \code{\link{fit_epi_model}}. The function displays a concise, human-readable
#' summary of the fitted epidemic model, including the model name, likelihood
#' distribution, optimization diagnostics, and estimated parameters.
#'
#' This method is automatically invoked when a \code{"fit_epi_model"} object
#' is printed at the console.
#'
#' @param x An object of class \code{"fit_epi_model"}.
#' @param ... Further arguments passed to or from other methods (currently ignored).
#'
#' @return
#' Invisibly returns the input object \code{x}.
#'
#' @examples
#' \dontrun{
#' fit <- fit_epi_model(x, model = SIR_MODEL)
#' print(fit)
#' }
#'
#' @seealso
#' \code{\link{fit_epi_model}},
#' \code{\link{predict.fit_epi_model}},
#' \code{\link{summary.fit_epi_model}}
#'
#' @export
print.fit_epi_model <- function(x, ...) {
  cat("<fit_epi_model>\n")
  cat("  Model: ", x$model$name, "\n", sep = "")
  cat("  distr: ", x$distr, "\n", sep = "")
  cat("  value: ", x$value, "\n", sep = "")
  cat("  convergence: ", x$convergence, "\n", sep = "")
  if (!is.null(x$message) && nzchar(x$message)) {
    cat("  message: ", x$message, "\n", sep = "")
  }
  cat("  Parameters:\n")
  for (nm in names(x$par)) cat("    ", nm, ": ", x$par[[nm]], "\n", sep = "")
  invisible(x)
}


