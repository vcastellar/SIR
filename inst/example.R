sirv_rhs <- function(time, state, parms) {

  with(as.list(c(state, parms)), {

    dS <- -beta * S * I - nu * S
    dI <-  beta * S * I - gamma * I
    dR <-  gamma * I
    dV <-  nu * S

    list(c(dS, dI, dR, dV))
  })
}


SIRV_MODEL <- epi_model(
  name      = "SIRV",
  rhs       = sirv_rhs,
  par_names = c("beta", "gamma", "nu"),
  states    = c("S", "I", "R", "V")
)


register_epi_model(SIRV_MODEL, family = "SIR")
