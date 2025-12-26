
#' @examples
#' xsim <- simulate_sir(beta = 0.3, gamma = 0.1)
#' x <- xsim$incidence_obs$inc    # longitud n_days
#' plot(x)
#' fit <- fit_sir_model(x)
#' params <- fit
#' times <- seq_along(x)
#' init <- c(S = 1e6 -10,
#'           I = 10,
#'           R = 0,
#'           C = 10)
#' out <- deSolve::ode(y = init,
#'                     times = times,
#'                     func = sir_c,
#'                     parms = params,
#'                     method = "lsoda")
#' out <- as.data.frame(out)
#' plot(x)
#' lines(out$incidence, lty = 2, col = "red")
#' plot(cumsum(x))
#' lines(out$C, type = 'l', lty = 2, col = "red")
#' @export
fit_sir_model <- function(x, init = list(I = 10, N = 1e6)) {

  times <- seq_along(x)
  ini0 <- c(S = init$N - init$I,
            I = init$I,
            R = 0,
            C = init$I)

  rss <- function(theta) {
    names(theta) <- c("beta", "gamma")

    out <- deSolve::ode(y = ini0,
                        times = times,
                        func = sir_c,
                        parms = theta,
                        method = "lsoda")
    mu <- diff(as.data.frame(out)[['C']])
    eps <- 1e-8
    mu <- pmax(mu, eps)


    -sum(dpois(x, lambda = mu, log = TRUE))

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
