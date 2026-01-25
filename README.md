# **SIR** — Deterministic Compartmental Epidemic Models in R

**SIR** is an R package for defining, simulating, fitting, and predicting
**deterministic compartmental epidemic models** based on systems of ordinary
differential equations (ODEs).

The package is intentionally **model-agnostic**, **explicit**, and
**deterministic by design**.  
It focuses on transparent model definition, reproducible simulation, and
trajectory-matching parameter estimation, without imposing probabilistic
observation models or Bayesian inference frameworks.

---

## What makes this package different

Most epidemic-modeling packages tightly couple model structure, simulation,
and inference, often hiding assumptions behind complex APIs.

**SIR takes a different approach**:

- Models are **explicit objects** (`epi_model`) that fully declare:
  - states
  - parameters
  - equations
  - bounds
  - outputs
- Simulation, fitting, prediction, and plotting are **generic operations**
  defined *on* these models.
- Parameter estimation is performed by **trajectory matching**, not by
  likelihood-based inference.
- Predictions are **conditionally linked to fitted models**, preserving
  full traceability.

The result is a framework that is:
- easy to inspect,
- easy to extend,
- and well suited for teaching, prototyping, and deterministic analysis.

---

## Key Features

- **Explicit epidemic models** via the `epi_model` abstraction
- Deterministic simulation using `deSolve::ode()`
- Built-in SIR, SIRS, and SEIR models
- Generic model outputs (states, incidence, or user-defined observables)
- Trajectory-matching parameter estimation (RMSE / log-RMSE)
- Multi-start optimization with flexible optimizers
- Conditional prediction objects linking fits and forecasts
- Base R printing, plotting, and summary methods
- Minimal dependency footprint

---

## Installation

The package is under active development and not yet available on CRAN.

Install the development version from GitHub:

```r
# install.packages("remotes")
remotes::install_github("vcastellar/SIR")
```

## Dependencies

- **deSolve**
- **stats** (base R)

---

## Core Abstraction: `epi_model`

An `epi_model` object encapsulates all model-specific information:

- model name
- state variable names
- parameter names
- ODE right-hand side (`rhs`)
- parameter bounds and defaults
- declared model outputs (e.g. `"I"`, `"S"`, `"incidence"`)
- optional human-readable equations

Once defined, an `epi_model` can be passed unchanged to:

- `simulate_epi()` — deterministic simulation
- `fit_epi_model()` — trajectory-matching parameter estimation
- `predict()` — forward prediction from fitted models
- `plot()` / `summary()` — visualization and reporting

This cleanly separates **model definition** from **model use**.

---

## Built-in Models

The package currently provides:

### `SIR_MODEL`

Susceptible–Infectious–Recovered model.

### `SIRS_MODEL`

Extension of SIR with waning immunity (`R → S`).

### `SEIR_MODEL`

Susceptible–Exposed–Infectious–Recovered model with a latent period.

In all models, outputs are explicitly declared and can be used consistently
for fitting, plotting, and prediction.

---

## Simulation

```r
library(SIR)

sim <- simulate_epi(
  model = SIR_MODEL,
  times = 0:200,
  parms = c(beta = 0.3, gamma = 0.1),
  init  = list(S = 1e6, I = 20, R = 0)
)

plot(sim)
plot(sim, what = "I")
summary(sim)
```


Simulation returns a `sim_epi` object containing:

- state trajectories
- declared model outputs
- optional incidence if defined by the model

---

## Model Fitting (Trajectory Matching)

Model parameters can be estimated by minimizing a discrepancy between observed
data and simulated trajectories.

```r
fit <- fit_epi_model(
  x = sim$states$I,
  model = SIR_MODEL,
  target = "I",
  init = list(S = 1e6, I = 10, R = 0)
)

fit
```

## Fitting philosophy

- No observation likelihood is assumed
- No stochasticity is introduced
- The objective is **deterministic trajectory matching**
- Loss functions include RMSE and log-RMSE
- Multi-start optimization is used to mitigate local minima

This makes the fitting procedure:

- fast,
- transparent,
- and easy to interpret.

---

## Prediction from Fitted Models

Predictions are generated conditionally on a fitted model and returned as a
dedicated object linking the fit and the forecast.

```r
pred <- predict(
  fit,
  times = 201:400
)

pred
```

The returned object preserves:

- the fitted model
- estimated parameters
- initial conditions
- numerical integration settings
- predicted trajectories

Predictions are therefore fully reproducible and auditable.

---

## Plotting and Summaries

### Simulation objects (`sim_epi`)

- `plot(sim)` — all state trajectories
- `plot(sim, what = "I")` — a single output
- `plot(sim, what = "incidence")` — incidence (if defined)
- `summary(sim)` — epidemiological summaries

### Fitted models (`fit_epi_model`)

- `print(fit)` — parameter estimates and diagnostics
- `predict(fit, ...)` — conditional forward predictions

---

## Design Philosophy

- **Deterministic by default**
- Explicit model structure
- No hidden assumptions
- Minimal API surface
- Small dependency footprint
- Emphasis on clarity and reproducibility

This package is not intended to replace full probabilistic inference frameworks,
but to provide a clean, deterministic foundation for epidemic modelin
