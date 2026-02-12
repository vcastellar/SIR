sird_rhs <- function(time, state, parms) {
  with(as.list(c(state, parms)), {
    N <- S + I + R
    lambda <- beta * S * I / N

    dS <- -lambda
    dI <-  lambda - gamma * I - mu * I
    dR <-  gamma * I
    dD <-  mu * I

    list(
      c(dS, dI, dR, dD),
      incidence = lambda
    )
  })
}

SIRD_MODEL <- epi_model(
  name      = "SIRD",
  rhs       = sird_rhs,
  states    = c("S", "I", "R", "D"),
  flows     = c("incidence"),
  par_names = c("beta", "gamma", "mu"),
  defaults  = c(beta = 0.3, gamma = 0.1, mu = 0.01),
  init      = c(S = 1e6, I = 10, R = 0, D = 0),
  roles     = list(susceptible = "S",
                   infectious  = "I",
                   recovered   = "R",
                   deceased    = "D"
                   )
)

register_epi_model(SIRD_MODEL, family = 'SIRD')
run_epi_app()

sim <- simulate_epi(model = SIRD_MODEL,
                    times = 0:200,
                    parms = SIRD_MODEL$defaults,
                    init =  SIRD_MODEL$init)
plot(sim)

plot(instantaneous_growth_rate(sim), type = 'l')
plot(doubling_time_ts(sim), type = 'l')
