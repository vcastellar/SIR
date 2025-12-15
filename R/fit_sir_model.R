
fit_sir_model <- function(infected, population, horizon_days) {

  inf_0 <- diff(infected)
  Day <- seq_along(inf_0)
  init <- c(S = population - inf_0[1], I = inf_0[1], R = 0)

  rss <- function(parameters) {

    names(parameters) <- c("beta", "gamma", "population")

    out <- deSolve::ode(y = init,
                        times = Day,
                        func = SIR,
                        parms = parameters)
    # plot(out)
    fit <- out[, 3]
    message("n obs: ", length(inf_0), " | n fit: ", length(fit))
    mean((inf_0 - fit) ^ 2)
  }

  # optimize with some sensible conditions
  opt <- optim(c(7, 3, 5e6),
               rss,
               method = "Nelder-Mead",
               control = list(maxit = 5000, factr = 1e4))

  opt <- optim(
    par = c(1, 1, 5e6),
    fn  = rss,
    method = "L-BFGS-B",
    lower = c(0, 0),
    upper = c(Inf, Inf),
    control = list(maxit = 5000)
  )
  opt$message
  opt_par <- setNames(opt$par, c("beta", "gamma", "population"))


  fit <- data.frame(ode(y = init,
                        times = Day,
                        func = SIR,
                        parms = opt_par))

  plot(inf_0)
  lines(fit$I, col = 'red')
  lines(fit$R)

  plot(inf_0)

  list(parameters = opt_par, fit = fit)
}
