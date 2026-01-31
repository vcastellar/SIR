# SIR

**SIR** is an R package for simulating and exploring deterministic
compartmental epidemiological models formulated as systems of ordinary
differential equations (ODEs).

The package focuses on **model definition and simulation**, providing
a small but flexible framework to study epidemic dynamics under different
assumptions and parameter values. An optional **Shiny application** is
included for interactive exploration.

---

## Features

- Deterministic simulation of epidemic models using ODEs
- Flexible definition of custom compartmental models
- Several classical models included:
  - SI
  - SIR
  - SIRS
  - SEIR
  - SEIRS
- Numerical integration based on `deSolve`
- Interactive Shiny app for visual exploration of model dynamics

> ⚠️ This package does **not** perform parameter estimation or statistical
> inference. Its scope is simulation and exploration only.

---

## Installation

Install the package from source:

```r
devtools::install("path/to/SIR")
library(SIR)

sim <- simulate_epi(
  model = SIR_MODEL,
  times = 0:200,
  parms = c(beta = 0.3, gamma = 0.1),
  init  = c(S = 1e6, I = 10, R = 0)
)

plot(sim)
plot(sim, what = "incidence")
```

## Defining a custom model

Users can define their own epidemic models by specifying the right-hand
side of the ODE system and creating an epi_model object.

```r
sir_rhs <- function(time, state, parms) {
  with(as.list(c(state, parms)), {
    N <- S + I + R
    lambda <- beta * S * I / N

    dS <- -lambda
    dI <-  lambda - gamma * I
    dR <-  gamma * I

    list(c(dS, dI, dR), incidence = lambda)
  })
}

my_sir <- new_epi_model(
  name        = "MySIR",
  rhs         = sir_rhs,
  state_names = c("S", "I", "R"),
  par_names   = c("beta", "gamma"),
  init        = c(S = 1e6, I = 10, R = 0),
  outputs     = c("S", "I", "R", "incidence")
)
```

This model can then be simulated using simulate_epi() like any built-in
model.

## Shiny application

The package includes a Shiny app for interactive exploration of epidemic
models.

```r
library(SIR)
run_epi_app()
```

The app allows users to:

- select a built-in model
- adjust parameters and initial conditions
- visualize state trajectories and incidence
- inspect the model equations

The Shiny app is optional; the package works without shiny installed.

## Scope and philosophy

SIR is designed as a simulation-oriented package:

- deterministic models
- explicit compartmental structure
- transparent dynamics
- minimal hidden state

It is intended for teaching, exploration, and rapid prototyping of
epidemiological models, rather than for statistical inference or
data-driven estimation.
