
#' @keywords internal
#' @noRd
objective_epi <- function(model, x, ini0,
                          distr = c("poisson", "negbin"),
                          eps = 1e-8) {
  distr <- match.arg(distr)
  times <- 0:length(x)

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
    C <- out[["C"]]
    if (any(!is.finite(C))) return(1e30)

    mu <- diff(C)

    # mu debe ser finita y >=0 para conteos
    if (any(!is.finite(mu))) return(1e30)
    if (any(mu < 0)) return(1e30)

    mu <- pmax(mu, eps)

    val <- switch(
      distr,
      negbin  = -sum(stats::dnbinom(x, mu = mu, size = 10, log = TRUE)),
      poisson = -sum(stats::dpois(x, lambda = mu, log = TRUE))
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
#' \code{"epi_model"}, such as \code{SIR_MODEL} or \code{SIRS_MODEL}) to observed
#' **incidence count data** by minimizing a **negative log-likelihood**.
#'
#' The function returns a fitted model object of class \code{"fit_epi_model"},
#' which contains the estimated parameters together with additional information
#' about the optimization process, the model definition, and the initial
#' conditions used for fitting.
#'
#' @details
#' ## Fitting approach
#' Internally, the function constructs an objective function using an internal
#' helper (see \code{objective_epi}, not exported) that:
#' \itemize{
#'   \item solves the model ODE system using \code{deSolve::ode()} and
#'     \code{model$rhs};
#'   \item computes expected incidence as increments of the cumulative state
#'     \eqn{C(t)} via \code{diff(C)};
#'   \item evaluates a Poisson or Negative Binomial log-likelihood for the
#'     observed incidence vector \code{x}.
#' }
#'
#' To improve numerical stability and reduce sensitivity to initial values,
#' the optimization is initialized via an internal **multi-start** strategy
#' (see \code{multi_start_optim}, not exported), which samples random initial
#' parameter vectors within the bounds defined by the model and retains the
#' best solution as the starting point for the final optimization.
#'
#' ## Data and time grid
#' The input vector \code{x} must represent incidence counts on an equally spaced
#' time grid (typically daily). The internal objective uses a time grid
#' \code{times = 0:length(x)} so that \code{diff(C)} has the same length as
#' \code{x}.
#'
#' ## Model requirements
#' The supplied \code{model} must:
#' \itemize{
#'   \item be an object of class \code{"epi_model"};
#'   \item provide an ODE right-hand side function compatible with
#'     \code{deSolve::ode()};
#'   \item define parameter names via \code{model$par_names};
#'   \item define parameter bounds via \code{model$lower} and \code{model$upper};
#'   \item include a cumulative state variable \code{C(t)} in the ODE output.
#' }
#'
#' @param x Numeric vector of observed incidence counts (e.g. daily cases).
#' @param model An \code{epi_model} object to be fitted (e.g. \code{SIR_MODEL} or
#'   \code{SIRS_MODEL}). Defaults to \code{SIR_MODEL}.
#' @param distr Character string specifying the observation distribution used in
#'   the likelihood. One of \code{"poisson"} or \code{"negbin"}.
#' @param init Named list defining the initial conditions through
#'   \code{init$N} (population size) and \code{init$I} (initial number of
#'   infectious individuals).
#' @param n_starts Integer. Number of random starting points used in the
#'   multi-start initialization.
#' @param control_ms List of control parameters passed to \code{optim()} during
#'   the multi-start phase.
#' @param control List of control parameters passed to \code{optim()} in the final
#'   optimization.
#' @param seed Optional integer seed for reproducible multi-start initialization.
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
#' }
#'
#' @examples
#' \dontrun{
#' ## Simulate a SIRS epidemic and fit the model to observed incidence
#' sim2 <- simulate_epi(
#'   model = SIR_MODEL,
#'   n_days = 200,
#'   parms = c(beta = 0.30, gamma = 0.10),
#'   init_args = list(N = 1e6, I0 = 20, R0 = 0),
#'   obs = "poisson",
#'   rho = 1,
#'   seed = 22
#' )
#' plot(sim2)
#'
#' x <- sim2$incidence_obs$inc
#' plot(x, type = "l", xlab = "Day", ylab = "Incidence")
#'
#' fit <- fit_epi_model(x, model = SIRS_MODEL, init = list(I = 10, N = 1e6))
#' fit
#' fit$par
#'
#' ## Compare fitted incidence to observed incidence
#' times <- 0:length(x)
#' ini0  <- c(S = 1e6 - 10, I = 10, R = 0, C = 10)
#'
#' out <- deSolve::ode(
#'   y = ini0,
#'   times = times,
#'   func = SIRS_MODEL$rhs,
#'   parms = fit$par,
#'   method = "lsoda"
#' )
#'
#' lines(as.data.frame(out)$incidence, col = "red")
#' }
#'
#' @seealso
#' \code{\link{simulate_epi}}, \code{\link{SIR_MODEL}}, \code{\link{SIRS_MODEL}},
#' \code{\link[stats]{optim}}, \code{\link[deSolve]{ode}}
#'
#' @export
fit_epi_model <- function(x,
                           model = SIR_MODEL,
                           distr = c("poisson", "negbin"),
                           init = list(I = 10, N = 1e6),
                           control = list(maxit = 5000),
                           n_starts = 30,
                           control_ms = list(maxit = 50),
                           seed = 1,
                           ...) {

  distr <- match.arg(distr)

  ini0 <- c(S = init$N - init$I, I = init$I, R = 0, C = init$I)

  best <- multi_start_optim(
    model = model,
    x = x,
    ini0 = ini0,
    n = n_starts,
    distr = distr,
    control = control_ms,
    seed = seed
  )

  fn <- objective_epi(
    model = model,
    x = x,
    ini0 = ini0,
    distr = distr
  )

  opt <- optim(
    par = unname(best$opt$par),
    fn  = fn,
    method = "L-BFGS-B",
    lower = unname(model$lower[model$par_names]),
    upper = unname(model$upper[model$par_names]),
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

#' @keywords internal
#' @noRd
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


