## Demo: Simulación y visualización de un modelo SIR
##
## Este demo muestra cómo:
## 1) simular una epidemia SIR con simulate_sir()
## 2) integrar explícitamente el modelo SIR con sir_c() usando deSolve
## 3) visualizar infectados activos e incidencia observada

## Cargar paquetes necesarios
library(deSolve)

## Mensaje inicial
cat("Demo SIR: simulación e integración del modelo\n")
cat("=============================================\n\n")

## ------------------------------------------------------------
## 1) Simular datos con simulate_sir()
## ------------------------------------------------------------
cat("1) Simulando datos con simulate_sir()...\n")

sim <- simulate_sir(
  n_days = 300,
  N = 1e6,
  beta = 0.30,
  gamma = 0.10,
  I0 = 20,
  rho = 0.25,
  obs = "negbin",
  size = 15,
  seed = 1
)

cat("   Simulación completada.\n\n")

## Inspección rápida
cat("Primeras filas de la incidencia observada:\n")
print(head(sim$incidence_obs))

cat("\nPrimeras filas de los casos acumulados observados:\n")
print(head(sim$cumulative_obs))

## ------------------------------------------------------------
## 2) Integrar explícitamente el modelo SIR con sir_c()
## ------------------------------------------------------------
cat("\n2) Integrando el modelo SIR con sir_c()...\n")

times <- sim$states$time

init <- c(
  S = sim$params$N - sim$params$I0,
  I = sim$params$I0,
  R = 0,
  C = sim$params$I0
)

pars <- c(
  beta = sim$params$beta,
  gamma = sim$params$gamma
)

out <- ode(
  y = init,
  times = times,
  func = sir_c,
  parms = pars,
  method = "lsoda"
)

out <- as.data.frame(out)

cat("   Integración completada.\n\n")

## ------------------------------------------------------------
## 3) Gráficos
## ------------------------------------------------------------
cat("3) Mostrando gráficos...\n")

op <- par(no.readonly = TRUE)
par(mfrow = c(2, 1), mar = c(4, 4, 2, 1))

## Infectados activos
plot(
  out$time, out$I,
  type = "l",
  xlab = "Días",
  ylab = "I(t)",
  main = "SIR: Infectados activos"
)

## Incidencia observada
plot(
  sim$incidence_obs$time,
  sim$incidence_obs$inc,
  type = "h",
  xlab = "Días",
  ylab = "Incidencia observada",
  main = "Incidencia observada (conteos)"
)

par(op)

cat("\nDemo finalizada.\n")
cat("Puedes volver a ejecutarla con: demo('sir_demo')\n")
