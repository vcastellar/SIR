
#' @keywords internal
#' @noRd
# objective using incidence_col directly
objective_epi_incidence <- function(model, x, ini0,
                                    distr = c("poisson", "negbin"),
                                    eps = 1e-8,
                                    rho = 1,
                                    size = 10) {
  distr <- match.arg(distr)

  # time grid: ensure incidence vector aligns with x
  # We integrate on 1:length(x) so output length matches x
  times <- seq_len(length(x))

  inc_col <- (model$output$incidence_col %||% "incidence")

  function(theta) {
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

    # basic validity for counts
    if (any(!is.finite(mu))) return(1e30)
    if (any(mu < 0)) return(1e30)

    mu <- pmax(rho * mu, eps)

    val <- switch(
      distr,
      poisson = -sum(stats::dpois(x, lambda = mu, log = TRUE)),
      negbin  = -sum(stats::dnbinom(x, mu = mu, size = size, log = TRUE))
    )
    if (!is.finite(val)) 1e30 else val
  }
}

#' @keywords internal
#' @noRd
multi_start_optim <- function(model,
                              x, ini0,
                              n = 30,
                              distr = c("poisson", "negbin"),
                              control = list(maxit = 50),
                              seed = NULL) {

  if (!is.null(seed)) set.seed(seed)

  fn <- objective_epi(
    model = model,
    x = x,
    ini0 = ini0,
    distr = distr
  )

  lower <- unname(model$lower[model$par_names])
  upper <- unname(model$upper[model$par_names])

  best <- list(value = Inf, par = NULL, opt = NULL)

  for (i in seq_len(n)) {
    print(paste("paso: ", i))
    par0 <- stats::runif(length(model$par_names), min = lower, max = upper)

    opt <- optim(
      par = par0,
      fn  = fn,
      method = "L-BFGS-B",
      lower = lower,
      upper = upper,
      control = control
    )

    if (is.finite(opt$value) && opt$value < best$value) {
      best <- list(
        value = opt$value,
        par   = setNames(opt$par, model$par_names),
        opt   = opt
      )
    }
  }

  best
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
#' @param rho Numeric in \eqn{[0,1]}. Reporting fraction mapping true model incidence
#'   to expected observed incidence (\eqn{\mu(t) = \rho \, \text{incidence}(t)}).
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
#'
#' @examples
#' \dontrun{
#' ## Simulate a SIR epidemic and fit the model to observed incidence
#' sim <- simulate_epi(
#'   model = SIR_MODEL,
#'   n_days = 200,
#'   parms = c(beta = 0.30, gamma = 0.10),
#'   init_args = list(N = 1e6, I0 = 20, R0 = 0),
#'   obs = "poisson",
#'   rho = 1,
#'   seed = 22
#' )
#'
#' x <- sim$incidence_obs$inc
#' fit <- fit_epi_model(
#'   x = x,
#'   model = SIR_MODEL,
#'   distr = "poisson",
#'   init = list(N = 1e6, I0 = 20, R0 = 0)
#' )
#' print(fit)
#'
#' ## Compare fitted incidence to observed incidence
#' pred <- predict(fit, n_days = length(x), init_args = list(N = 1e6, I0 = 20, R0 = 0))
#' plot(x, type = "l", xlab = "Day", ylab = "Incidence")
#' lines(pred$incidence$time, pred$incidence$inc, col = "red", lty = 2)
#' }
#'
#' @seealso
#' \code{\link{simulate_epi}}, \code{\link{predict.fit_epi_model}},
#' \code{\link{SIR_MODEL}}, \code{\link{SIRS_MODEL}}, \code{\link{SEIR_MODEL}},
#' \code{\link[stats]{optim}}, \code{\link[deSolve]{ode}}
#'
#' @export
fit_epi_model <- function(x,
                          model = SIR_MODEL,
                          distr = c("negbin", "poisson"),
                          init = list(N = 1e6, I0 = 10, R0 = 0),
                          rho = 1,
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
    stop("`model$make_init()` output is missing state(s): ", paste(miss_state, collapse = ", "))
  }
  ini0 <- ini0[model$state_names]

  # 4) build objective
  fn <- objective_epi_incidence(
    model = model,
    x = x,
    ini0 = ini0,
    distr = distr,
    eps = eps,
    rho = rho,
    size = size
  )

  # 5) multi-start (assumes your existing multi_start_optim uses objective_epi-like signature)
  # Here we implement a minimal multi-start inline to avoid dependence on diff(C)
  if (!is.null(seed)) set.seed(seed)

  lower <- unname(model$lower[model$par_names])
  upper <- unname(model$upper[model$par_names])

  best <- list(value = Inf, opt = NULL)

  for (i in seq_len(n_starts)) {
    par0 <- stats::runif(length(model$par_names), min = lower, max = upper)

    opt_i <- optim(
      par = par0,
      fn  = fn,
      method = "L-BFGS-B",
      lower = lower,
      upper = upper,
      control = control_ms
    )

    if (is.finite(opt_i$value) && opt_i$value < best$value) {
      best <- list(value = opt_i$value, opt = opt_i)
    }
  }

  if (is.null(best$opt)) {
    stop("Multi-start failed: no finite objective value found. Try different `init`, bounds, or `control_ms`.")
  }

  # 6) final optimization from best start
  opt <- optim(
    par = unname(best$opt$par),
    fn  = fn,
    method = "L-BFGS-B",
    lower = lower,
    upper = upper,
    control = control
  )

  par_hat <- setNames(opt$par, model$par_names)

  new_fit_epi_model(
    model = model,
    par = par_hat,
    optim = opt,
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


