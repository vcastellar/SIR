
#' @examples
#' xsim <- simulate_sir(beta = 0.3, gamma = 0.15)
#' x <- round(diff(xsim$states$C ))    # longitud n_days
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
#' lines(out$incidence, lty = 2)
#' plot(cumsum(x))
#' lines(out$C, type = 'l', lty = 2)
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


  opt <- optim(
    par = c(0.32, 0.1),
    fn  = rss,
    method = "L-BFGS-B",
    lower = c(0, 0),
    control = list(maxit = 5000)
  )

  opt$message
  opt_par <- setNames(opt$par, c("beta", "gamma"))

  return(opt_par)
}
