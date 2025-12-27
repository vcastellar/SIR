#' Fit an SIR model to incidence data via likelihood
#'
#' @name fit_sir_model
#'
#' @description
#' Estimates the transmission rate (\eqn{\beta}) and recovery/removal rate
#' (\eqn{\gamma}) of a deterministic SIR model by fitting the model-implied daily
#' incidence to observed incidence counts using maximum likelihood.
#'
#' The SIR dynamics are solved numerically with \code{deSolve::ode()} (method
#' \code{"lsoda"}). The model uses an auxiliary cumulative-infections state
#' \eqn{C(t)} whose daily increments \eqn{\Delta C(t)} are matched to the observed
#' daily incidence series \code{x}.
#'
#' @details
#' ## Model and mean incidence
#' The underlying ODE system is defined in \code{\link{sir}} and includes the
#' auxiliary variable \eqn{C(t)} with \eqn{dC/dt = \lambda(t)}. The expected
#' daily incidence is computed as:
#' \deqn{\mu_t = C(t) - C(t-1), \quad t = 2,\ldots,T}
#' In practice, \code{mu} is obtained with \code{diff(C)} on the ODE solution.
#'
#' A small positive floor (\code{eps = 1e-8}) is applied to \code{mu} to avoid
#' non-finite log-likelihood values when \eqn{\mu_t \le 0}.
#'
#' ## Observation model (likelihood)
#' The observation distribution is selected with \code{distr}:
#' \describe{
#'   \item{\code{"poisson"}}{Poisson likelihood:
#'     \eqn{X_t \sim \mathrm{Poisson}(\mu_t)}.}
#'   \item{\code{"negbin"}}{Negative binomial likelihood (overdispersed counts):
#'     \eqn{X_t \sim \mathrm{NegBin}(\mu_t, size)} using
#'     \code{stats::dnbinom(mu = mu, size = 10)} (fixed dispersion in this implementation).}
#' }
#'
#' ## Optimization
#' Parameters \eqn{(\beta,\gamma)} are estimated by minimizing the negative log-likelihood
#' using \code{optim(method = "L-BFGS-B")} with non-negativity bounds.
#'
#' To reduce sensitivity to local minima, the function performs a simple multi-start
#' scan: it draws \code{n} random initial points for \eqn{\beta} and \eqn{\gamma}
#' within user-defined ranges and keeps the best solution as the starting point for
#' the final optimization.
#'
#' @param x Numeric vector of observed daily incidence counts (non-negative).
#' @param distr Character. Observation model / likelihood. One of
#'   \code{"negbin"} or \code{"poisson"}.
#' @param init Named list with initial conditions and scale. Must include:
#'   \describe{
#'     \item{\code{I}}{Initial number of infectious individuals at day 0.}
#'     \item{\code{N}}{Total population size used to set \eqn{S(0) = N - I(0)}.}
#'   }
#'   The initial state vector is \code{c(S = N - I, I = I, R = 0, C = I)}.
#' @param ... Reserved for future extensions (currently unused).
#'
#' @return A named numeric vector with the maximum likelihood estimates:
#' \describe{
#'   \item{\code{beta}}{Estimated transmission rate (per day).}
#'   \item{\code{gamma}}{Estimated recovery/removal rate (per day).}
#' }
#'
#' @seealso
#' \code{\link{sir}}, \code{\link{simulate_epi}}
#'
#' @examples
#' # Simulate data and fit the model
#' sim <- simulate_epi(n_days = 200, model = "sir", N = 1e6,
#'                     beta = 0.30, gamma = 0.10, rho = 1,
#'                     obs = "negbin", size = 20, seed = 1)
#' x <- sim$incidence_obs$inc
#'
#' fit_pois <- fit_sir_model(x, distr = "poisson", init = list(I = 10, N = 1e6))
#' fit_nb   <- fit_sir_model(x, distr = "negbin",  init = list(I = 10, N = 1e6))
#'
#' fit_pois
#' fit_nb
#'
#' @export
fit_sir_model <- function(x, distr = c("negbin", "poisson"), init = list(I = 10, N = 1e6), ...) {

  distr <- match.arg(distr)

  times <- seq_along(x)
  ini0 <- c(S = init$N - init$I, I = init$I, R = 0, C = init$I)

  rss <- function(theta) {
    names(theta) <- c("beta", "gamma")

    out <- deSolve::ode(y = ini0,
                        times = times,
                        func = sir,
                        parms = theta,
                        method = "lsoda")
    mu <- diff(as.data.frame(out)[['C']])
    eps <- 1e-8
    mu <- pmax(mu, eps)

    verosimilitud <- switch(
      distr,
      negbin  = -sum(stats::dnbinom(x, mu = mu, size = 10, log = TRUE)),
      poisson = -sum(stats::dpois(x, lambda = mu, log = TRUE))
    )


    return(verosimilitud)

  }

  # escaneo Multi-start de parámetros beta y gamma iniciales
  multi_start_optim <- function(fn, n = 30,
                                beta_range  = c(0.01, 2),
                                gamma_range = c(0.02, 1)) {
    best <- list(value = Inf, par = NULL, opt = NULL)
    for (i in seq_len(n)) {
      par0 <- c(runif(1, beta_range[1], beta_range[2]),
                runif(1, gamma_range[1], gamma_range[2]))
      opt <- optim(par0, fn, method = "L-BFGS-B",
                   lower = c(1e-8, 1e-8),
                   upper = c(beta_range[2], gamma_range[2]))
      if (is.finite(opt$value) && opt$value < best$value) {
        best <- list(value=opt$value, par=opt$par, opt=opt)
      }
    return(best)
    }
  }

  best_init_parms <- multi_start_optim(fn = rss)

  opt <- optim(
    par = best_init_parms$opt$par,
    fn  = rss,
    method = "L-BFGS-B",
    lower = c(0, 0),
    control = list(maxit = 5000)
  )

  opt$message
  opt_par <- setNames(opt$par, c("beta", "gamma"))

  return(opt_par)
}
