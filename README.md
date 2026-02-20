# SIR

**SIR** is an R package for simulating and exploring deterministic
compartmental epidemiological models formulated as systems of ordinary
differential equations (ODEs).

The package focuses on **model definition and simulation**, providing
a small but flexible framework to study epidemic dynamics under different
assumptions and parameter values. A **Shiny application** is
included for interactive exploration.

---

## Features

- Deterministic simulation of epidemic models using ODEs
- Flexible definition of custom compartmental models
- Several classical models included:
  - SI
  - SIR
  - SIR with vital dynamics
  - SIRS
  - SEIR
  - SEIRS
- Numerical integration based on `deSolve`
- Built-in post-simulation metrics (e.g., peak incidence, attack rate, growth rates)
- In-session model registry (`register_epi_model()`, `list_models()`, `get_model()`)
- Interactive Shiny app for visual exploration of model dynamics

> ⚠️ This package does **not** perform parameter estimation or statistical
> inference. Its scope is simulation and exploration only.

---

## Installation

Install the package from source:

```r
devtools::install_local("path/to/SIR")
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

You can also install from GitHub:

```r
remotes::install_github("vcastellar/SIR")
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

my_sir <- epi_model(
  name        = "MySIR",
  rhs         = sir_rhs,
  states      = c("S", "I", "R"),
  par_names   = c("beta", "gamma"),
  init        = c(S = 1e6, I = 10, R = 0)
)
```

This model can then be simulated using simulate_epi() like any built-in
model.

## Working with the model registry

Custom models can be registered in the package's in-memory registry so they
can be discovered and reused during the current R session.

```r
register_epi_model(my_sir)
list_models()
model <- get_model("MySIR")
```

To remove a model from the active session:

```r
unregister_epi_model("MySIR")
```

## Built-in summary metrics

After simulation, helper functions can be used to summarize epidemic
trajectories, for example:

```r
inc <- sim$derived$incidence
t   <- sim$states$time

peak_incidence(inc, t)
time_to_peak(inc, t)
attack_rate(inc)
```

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

The Shiny app is included in the package and uses the installed `shiny`
dependency.

## Scope and philosophy

SIR is designed as a simulation-oriented package:

- deterministic models
- explicit compartmental structure
- transparent dynamics
- minimal hidden state

It is intended for teaching, exploration, and rapid prototyping of
epidemiological models, rather than for statistical inference or
data-driven estimation.
