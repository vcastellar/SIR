

objective_epi <- function(model, x, ini0,
                          distr = c("poisson", "negbin"),
                          eps = 1e-8) {
  distr <- match.arg(distr)

  # Para que diff(C) tenga la misma longitud que x:
  times <- 0:length(x)

  function(theta) {
    names(theta) <- model$par_names

    out <- deSolve::ode(
      y = ini0,
      times = times,
      func = model$rhs,
      parms = theta,
      method = "lsoda"
    )
    out <- as.data.frame(out)

    mu <- diff(out[["C"]])
    mu <- pmax(mu, eps)

    switch(
      distr,
      negbin  = -sum(stats::dnbinom(x, mu = mu, size = 10, log = TRUE)),
      poisson = -sum(stats::dpois(x, lambda = mu, log = TRUE))
    )
  }
}
multi_start_optim(model = SIR_MODEL, x = x, ini0 = ini0)
# escaneo Multi-start de parámetros beta y gamma iniciales
multi_start_optim <- function(model,
                              x, ini0,
                              n = 30,
                              distr = c("poisson", "negbin"),
                              control = list(maxit = 2000),
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

#-------------------------------------------------------------------------------

fit_sir_model <- function(x, model = NULL, distr = c("negbin", "poisson"),
                          init = list(I = 10, N = 1e6), ...) {

  distr <- match.arg(distr)

  times <- seq_along(x)
  ini0 <- c(S = init$N - init$I, I = init$I, R = 0, C = init$I)

  best <- multi_start_optim(model, x, ini0, n = 30,
                    distr = c("poisson", "negbin"),
                    control = list(maxit = 500),
                    seed = 1)

  fn <- objective_epi(
    model = model,
    x = x,
    ini0 = ini0,
    distr = distr
  )

  opt <- optim(
    par = best$opt$par,
    fn  = fn,
    method = "L-BFGS-B",
    lower = lower,
    upper = upper,
    control = control
  )

  opt$message
  opt_par <- setNames(opt$par, model$par_names)

  return(opt_par)
}

sim2 <- simulate_epi(
  model = SIRS_MODEL,
  n_days = 200,
  parms = c(beta = 0.30, gamma = 0.10, omega = 0.02),
  init_args = list(N = 1e6, I0 = 20, R0 = 0),
  obs = "poisson",
  rho = 1,
  seed = 22
)
plot(sim2)
sim2
x <- sim2$incidence_obs$inc
plot(x)
fit <- fit_sir_model(x, model = SIRS_MODEL, init = list(I = 10, N = 1e6))

out <- deSolve::ode(y = ini0,
                    times = times,
                    func = sir,
                    parms = fit,
                    method = "lsoda")
lines(as.data.frame(out)$incidence, col = "red")




