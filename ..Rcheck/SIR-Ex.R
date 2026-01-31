pkgname <- "SIR"
source(file.path(R.home("share"), "R", "examples-header.R"))
options(warn = 1)
base::assign(".ExTimings", "SIR-Ex.timings", pos = 'CheckExEnv')
base::cat("name\tuser\tsystem\telapsed\n", file=base::get(".ExTimings", pos = 'CheckExEnv'))
base::assign(".format_ptime",
function(x) {
  if(!is.na(x[4L])) x[1L] <- x[1L] + x[4L]
  if(!is.na(x[5L])) x[2L] <- x[2L] + x[5L]
  options(OutDec = '.')
  format(x[1L:3L], digits = 7L)
},
pos = 'CheckExEnv')

### * </HEADER>
library('SIR')

base::assign(".oldSearch", base::search(), pos = 'CheckExEnv')
base::assign(".old_wd", base::getwd(), pos = 'CheckExEnv')
cleanEx()
nameEx("SEIRS_MODEL")
### * SEIRS_MODEL

flush(stderr()); flush(stdout())

base::assign(".ptime", proc.time(), pos = "CheckExEnv")
### Name: SEIRS_MODEL
### Title: SEIRS epidemic model with latent period and waning immunity
### Aliases: SEIRS_MODEL
### Keywords: datasets

### ** Examples

## Simulate a SEIRS epidemic
sim <- simulate_epi(
  model = SEIRS_MODEL,
  times = 0:300,
  parms = c(beta = 0.3, sigma = 0.2, gamma = 0.14, omega = 0.01),
  init  = c(S = 1e6, E = 0, I = 20, R = 0),
  obs   = "poisson"
)

plot(sim)

## Plot observed incidence
plot(sim, what = "incidence")

## Fit the model to observed incidence
fit_inc <- fit_epi_model(
  x = sim$incidence$inc,
  model = SEIRS_MODEL,
  init = SEIRS_MODEL$init,
  target = "incidence"
)

fit_inc




base::assign(".dptime", (proc.time() - get(".ptime", pos = "CheckExEnv")), pos = "CheckExEnv")
base::cat("SEIRS_MODEL", base::get(".format_ptime", pos = 'CheckExEnv')(get(".dptime", pos = "CheckExEnv")), "\n", file=base::get(".ExTimings", pos = 'CheckExEnv'), append=TRUE, sep="\t")
cleanEx()
nameEx("SEIR_MODEL")
### * SEIR_MODEL

flush(stderr()); flush(stdout())

base::assign(".ptime", proc.time(), pos = "CheckExEnv")
### Name: SEIR_MODEL
### Title: SEIR epidemic model with latent (exposed) period
### Aliases: SEIR_MODEL
### Keywords: datasets

### ** Examples

## Simulate a SEIR epidemic
sim <- simulate_epi(
  model = SEIR_MODEL,
  times = 0:200,
  parms = c(beta = 0.3, sigma = 0.2, gamma = 0.14),
  init  = c(S = 1e6, E = 5, I = 10, R = 0),
  obs   = "poisson"
)

plot(sim)

## Plot observed incidence
plot(sim, what = "incidence")

## Fit the model to observed incidence
fit_inc <- fit_epi_model(
  x = sim$incidence$inc,
  model = SEIR_MODEL,
  init = SEIR_MODEL$init,
  target = "incidence"
)

fit_inc




base::assign(".dptime", (proc.time() - get(".ptime", pos = "CheckExEnv")), pos = "CheckExEnv")
base::cat("SEIR_MODEL", base::get(".format_ptime", pos = 'CheckExEnv')(get(".dptime", pos = "CheckExEnv")), "\n", file=base::get(".ExTimings", pos = 'CheckExEnv'), append=TRUE, sep="\t")
cleanEx()
nameEx("SIRS_MODEL")
### * SIRS_MODEL

flush(stderr()); flush(stdout())

base::assign(".ptime", proc.time(), pos = "CheckExEnv")
### Name: SIRS_MODEL
### Title: SIRS epidemic model with waning immunity
### Aliases: SIRS_MODEL
### Keywords: datasets

### ** Examples

## Simulate a SIRS epidemic
sim <- simulate_epi(
  model = SIRS_MODEL,
  times = 0:200,
  parms = c(beta = 0.3, gamma = 0.1, omega = 0.02),
  init  = c(S = 1e6, I = 20, R = 0)
)

plot(sim)

## Plot observed incidence (if an observation model is used)
plot(sim, what = "incidence")

## Fit the model to observed incidence
fit_inc <- fit_epi_model(
  x = sim$incidence$inc,
  model = SIRS_MODEL,
  init = SIRS_MODEL$init,
  target = "incidence"
)

fit_inc




base::assign(".dptime", (proc.time() - get(".ptime", pos = "CheckExEnv")), pos = "CheckExEnv")
base::cat("SIRS_MODEL", base::get(".format_ptime", pos = 'CheckExEnv')(get(".dptime", pos = "CheckExEnv")), "\n", file=base::get(".ExTimings", pos = 'CheckExEnv'), append=TRUE, sep="\t")
cleanEx()
nameEx("SIR_MODEL")
### * SIR_MODEL

flush(stderr()); flush(stdout())

base::assign(".ptime", proc.time(), pos = "CheckExEnv")
### Name: SIR_MODEL
### Title: SIR epidemic model
### Aliases: SIR_MODEL
### Keywords: datasets

### ** Examples

## Simulate a SIR epidemic
sim <- simulate_epi(
  model = SIR_MODEL,
  times = 0:200,
  parms = c(beta = 0.3, gamma = 0.1),
  init  = c(S = 1e6, I = 10, R = 0)
)

plot(sim)

## Plot observed incidence (if an observation model is used)
plot(sim, what = "incidence")

## Fit the model to observed incidence
fit_inc <- fit_epi_model(
  x = sim$incidence$inc,
  model = SIR_MODEL,
  init = SIR_MODEL$init,
  target = "incidence"
)

fit_inc




base::assign(".dptime", (proc.time() - get(".ptime", pos = "CheckExEnv")), pos = "CheckExEnv")
base::cat("SIR_MODEL", base::get(".format_ptime", pos = 'CheckExEnv')(get(".dptime", pos = "CheckExEnv")), "\n", file=base::get(".ExTimings", pos = 'CheckExEnv'), append=TRUE, sep="\t")
cleanEx()
nameEx("SI_MODEL")
### * SI_MODEL

flush(stderr()); flush(stdout())

base::assign(".ptime", proc.time(), pos = "CheckExEnv")
### Name: SI_MODEL
### Title: SI epidemic model
### Aliases: SI_MODEL
### Keywords: datasets

### ** Examples

## Simulate an SI epidemic
sim <- simulate_epi(
  model = SI_MODEL,
  times = 0:100,
  parms = c(beta = 0.4),
  obs = "negbin",
  init = SI_MODEL$init
)

plot(sim)

## Plot observed incidence
plot(sim, what = "incidence")

## Fit the model to observed incidence
fit_inc <- fit_epi_model(
  x = sim$incidence$inc,
  model = SI_MODEL,
  init = SI_MODEL$init,
  target = "incidence"
)

fit_inc




base::assign(".dptime", (proc.time() - get(".ptime", pos = "CheckExEnv")), pos = "CheckExEnv")
base::cat("SI_MODEL", base::get(".format_ptime", pos = 'CheckExEnv')(get(".dptime", pos = "CheckExEnv")), "\n", file=base::get(".ExTimings", pos = 'CheckExEnv'), append=TRUE, sep="\t")
cleanEx()
nameEx("compute_growth_metrics")
### * compute_growth_metrics

flush(stderr()); flush(stdout())

base::assign(".ptime", proc.time(), pos = "CheckExEnv")
### Name: compute_growth_metrics
### Title: Compute short-horizon growth metrics from an infection time
###   series
### Aliases: compute_growth_metrics

### ** Examples

infected <- ts(c(100, 110, 121, 133, 146, 161))
out <- compute_growth_metrics(infected)
out$inc_pct
out$doubling_days



base::assign(".dptime", (proc.time() - get(".ptime", pos = "CheckExEnv")), pos = "CheckExEnv")
base::cat("compute_growth_metrics", base::get(".format_ptime", pos = 'CheckExEnv')(get(".dptime", pos = "CheckExEnv")), "\n", file=base::get(".ExTimings", pos = 'CheckExEnv'), append=TRUE, sep="\t")
cleanEx()
nameEx("fit_epi_model")
### * fit_epi_model

flush(stderr()); flush(stdout())

base::assign(".ptime", proc.time(), pos = "CheckExEnv")
### Name: fit_epi_model
### Title: Fit an epidemic model by trajectory matching
### Aliases: fit_epi_model

### ** Examples

## Not run: 
##D ## ------------------------------------------------------------
##D ## Example 1: Fit incidence using log-RMSE
##D ## ------------------------------------------------------------
##D sim <- simulate_epi(
##D   model = SIR_MODEL,
##D   times = 0:200,
##D   parms = SIR_MODEL$defaults,
##D   init  = list(S = 1e6, I = 20, R = 0),
##D   seed  = 22
##D )
##D 
##D x_inc <- sim$incidence$inc
##D 
##D fit_inc <- fit_epi_model(
##D   x = x_inc,
##D   model = SIR_MODEL,
##D   target = "incidence",
##D   loss = "logrmse",
##D   init = list(S = 1e6, I = 6, R = 0)
##D )
##D 
##D print(fit_inc)
##D 
##D ## ------------------------------------------------------------
##D ## Example 2: Fit a state variable using a different optimizer
##D ## ------------------------------------------------------------
##D fit_I <- fit_epi_model(
##D   x = sim$states$I,
##D   model = SIR_MODEL,
##D   target = "I",
##D   init = list(S = 1e6, I = 10, R = 0),
##D   optim_method = "Nelder-Mead",
##D   method = "rk4",
##D   rtol = 1e-8,
##D   atol = 1e-10
##D )
##D 
##D print(fit_I)
## End(Not run)




base::assign(".dptime", (proc.time() - get(".ptime", pos = "CheckExEnv")), pos = "CheckExEnv")
base::cat("fit_epi_model", base::get(".format_ptime", pos = 'CheckExEnv')(get(".dptime", pos = "CheckExEnv")), "\n", file=base::get(".ExTimings", pos = 'CheckExEnv'), append=TRUE, sep="\t")
cleanEx()
nameEx("new_epi_model")
### * new_epi_model

flush(stderr()); flush(stdout())

base::assign(".ptime", proc.time(), pos = "CheckExEnv")
### Name: new_epi_model
### Title: Create a new epidemic model object
### Aliases: new_epi_model

### ** Examples

sir_model <- new_epi_model(
  name = "SIR",
  rhs = sir_rhs,
  state_names = c("S", "I", "R"),
  par_names = c("beta", "gamma"),
  incidence = "incidence"
)




base::assign(".dptime", (proc.time() - get(".ptime", pos = "CheckExEnv")), pos = "CheckExEnv")
base::cat("new_epi_model", base::get(".format_ptime", pos = 'CheckExEnv')(get(".dptime", pos = "CheckExEnv")), "\n", file=base::get(".ExTimings", pos = 'CheckExEnv'), append=TRUE, sep="\t")
cleanEx()
nameEx("plot.sim_epi")
### * plot.sim_epi

flush(stderr()); flush(stdout())

base::assign(".ptime", proc.time(), pos = "CheckExEnv")
### Name: plot.sim_epi
### Title: Plot a simulated epidemic
### Aliases: plot.sim_epi

### ** Examples

sim <- simulate_epi(
  model = SIR_MODEL,
  times = 0:200,
  time_unit = "days",
  parms = c(beta = 0.3, gamma = 0.1),
  init  = c(S = 999990, I = 10, R = 0, C = 10),
  obs   = "poisson",
  seed  = 1
)

# Plot all model states
plot(sim)

# Plot observed incidence (requires an observation model)
plot(sim, what = "incidence")




base::assign(".dptime", (proc.time() - get(".ptime", pos = "CheckExEnv")), pos = "CheckExEnv")
base::cat("plot.sim_epi", base::get(".format_ptime", pos = 'CheckExEnv')(get(".dptime", pos = "CheckExEnv")), "\n", file=base::get(".ExTimings", pos = 'CheckExEnv'), append=TRUE, sep="\t")
cleanEx()
nameEx("predict.fit_epi_model")
### * predict.fit_epi_model

flush(stderr()); flush(stdout())

base::assign(".ptime", proc.time(), pos = "CheckExEnv")
### Name: predict.fit_epi_model
### Title: Predict epidemic dynamics from a fitted epidemic model
### Aliases: predict.fit_epi_model

### ** Examples

## Not run: 
##D sim <- simulate_epi(
##D   model = SIR_MODEL,
##D   n_days = 200,
##D   parms = c(beta = 0.30, gamma = 0.10),
##D   init = list(S = 1e6, I = 20, R = 0),
##D   obs = "poisson"
##D )
##D plot(sim)
##D inc_obs <- sim$incidence$inc
##D plot(seq_along(inc_obs) - 1, inc_obs,
##D      type = "l", xlab = "Day", ylab = "Incidence")
##D fit <- fit_epi_model(inc_obs,
##D                      loss = "logrmse",
##D                      model = SIRS_MODEL,
##D                      init = list(I0 = 6, N = 1e6))
##D init <- tail(sim$states, n = 1)[, -1]
##D pred <- predict(
##D   object = fit,
##D   n_days = 1000,
##D   init = init
##D )
##D plot(pred)
##D plot(sim)
##D plot(pred$states$time, pred$states$I)
##D plot(pred$states$I, col = "red", lty = 2)
##D lines(pred$states$S, col = "blue", lty = 2)
## End(Not run)




base::assign(".dptime", (proc.time() - get(".ptime", pos = "CheckExEnv")), pos = "CheckExEnv")
base::cat("predict.fit_epi_model", base::get(".format_ptime", pos = 'CheckExEnv')(get(".dptime", pos = "CheckExEnv")), "\n", file=base::get(".ExTimings", pos = 'CheckExEnv'), append=TRUE, sep="\t")
cleanEx()
nameEx("print.epi_model")
### * print.epi_model

flush(stderr()); flush(stdout())

base::assign(".ptime", proc.time(), pos = "CheckExEnv")
### Name: print.epi_model
### Title: Print an epidemic model object
### Aliases: print.epi_model

### ** Examples

SIR_MODEL




base::assign(".dptime", (proc.time() - get(".ptime", pos = "CheckExEnv")), pos = "CheckExEnv")
base::cat("print.epi_model", base::get(".format_ptime", pos = 'CheckExEnv')(get(".dptime", pos = "CheckExEnv")), "\n", file=base::get(".ExTimings", pos = 'CheckExEnv'), append=TRUE, sep="\t")
cleanEx()
nameEx("print.fit_epi_model")
### * print.fit_epi_model

flush(stderr()); flush(stdout())

base::assign(".ptime", proc.time(), pos = "CheckExEnv")
### Name: print.fit_epi_model
### Title: Print a fitted epidemic model
### Aliases: print.fit_epi_model

### ** Examples

## Not run: 
##D fit <- fit_epi_model(
##D   x = incidence_data,
##D   model = SIR_MODEL
##D )
##D 
##D fit
## End(Not run)




base::assign(".dptime", (proc.time() - get(".ptime", pos = "CheckExEnv")), pos = "CheckExEnv")
base::cat("print.fit_epi_model", base::get(".format_ptime", pos = 'CheckExEnv')(get(".dptime", pos = "CheckExEnv")), "\n", file=base::get(".ExTimings", pos = 'CheckExEnv'), append=TRUE, sep="\t")
cleanEx()
nameEx("print.sim_epi")
### * print.sim_epi

flush(stderr()); flush(stdout())

base::assign(".ptime", proc.time(), pos = "CheckExEnv")
### Name: print.sim_epi
### Title: Print a simulated epidemic
### Aliases: print.sim_epi

### ** Examples

sim <- simulate_epi(
  model = SIR_MODEL,
  times = 0:200,
  time_unit = "days",
  parms = c(beta = 0.30, gamma = 0.10),
  init  = c(S = 999990, I = 10, R = 0, C = 10),
  seed  = 1
)

sim




base::assign(".dptime", (proc.time() - get(".ptime", pos = "CheckExEnv")), pos = "CheckExEnv")
base::cat("print.sim_epi", base::get(".format_ptime", pos = 'CheckExEnv')(get(".dptime", pos = "CheckExEnv")), "\n", file=base::get(".ExTimings", pos = 'CheckExEnv'), append=TRUE, sep="\t")
cleanEx()
nameEx("simulate_epi")
### * simulate_epi

flush(stderr()); flush(stdout())

base::assign(".ptime", proc.time(), pos = "CheckExEnv")
### Name: simulate_epi
### Title: Simulate an epidemic model defined by an 'epi_model' object
### Aliases: simulate_epi

### ** Examples

## ------------------------------------------------------------------
## Example 1: SIR model with model-defined incidence
## ------------------------------------------------------------------

sim <- simulate_epi(
  model = SIR_MODEL,
  times = 0:200,
  time_unit = "week",
  parms = c(beta = 0.30, gamma = 0.10),
  init  = c(S = 1e6 - 10, I = 10, R = 0, C = 10),
  seed  = 1,
  method = 'lsoda'
)

# Plot state trajectories defined by the model
plot(sim)

# Plot observed incidence (requires the model to define incidence)
plot(sim, what = "incidence")


## ------------------------------------------------------------------
## Example 2: Using a different numerical integration method
## ------------------------------------------------------------------

sim_rk4 <- simulate_epi(
  model = SIR_MODEL,
  times = 0:200,
  parms = c(beta = 0.30, gamma = 0.10),
  init  = c(S = 1e6 - 10, I = 10, R = 0, C = 10),
  method = "rk4"
)

plot(sim_rk4)


## ------------------------------------------------------------------
## Example 3: Model without an incidence definition
## ------------------------------------------------------------------

sim <- simulate_epi(
  model = SI_MODEL,
  times = 30:100,
  parms = c(beta = 0.25),
  init  = c(S = 999, I = 1)
)

plot(sim)
plot(sim, what = "incidence")


## ------------------------------------------------------------------
## Example 4: Stochastic observation model
## ------------------------------------------------------------------

sim <- simulate_epi(
  model = SIR_MODEL,
  times = 0:150,
  parms = c(beta = 0.35, gamma = 0.12),
  init  = c(S = 1e5 - 5, I = 5, R = 0, C = 5),
  obs   = "negbin",
  size  = 20,
  seed  = 123
)

plot(sim, what = "incidence")




base::assign(".dptime", (proc.time() - get(".ptime", pos = "CheckExEnv")), pos = "CheckExEnv")
base::cat("simulate_epi", base::get(".format_ptime", pos = 'CheckExEnv')(get(".dptime", pos = "CheckExEnv")), "\n", file=base::get(".ExTimings", pos = 'CheckExEnv'), append=TRUE, sep="\t")
cleanEx()
nameEx("summary.sim_epi")
### * summary.sim_epi

flush(stderr()); flush(stdout())

base::assign(".ptime", proc.time(), pos = "CheckExEnv")
### Name: summary.sim_epi
### Title: Summarize a simulated epidemic
### Aliases: summary.sim_epi

### ** Examples

sim <- simulate_epi(n_days = 300, model = SIRS_MODEL, omega = 1/180, seed = 1)

summary(sim)




base::assign(".dptime", (proc.time() - get(".ptime", pos = "CheckExEnv")), pos = "CheckExEnv")
base::cat("summary.sim_epi", base::get(".format_ptime", pos = 'CheckExEnv')(get(".dptime", pos = "CheckExEnv")), "\n", file=base::get(".ExTimings", pos = 'CheckExEnv'), append=TRUE, sep="\t")
### * <FOOTER>
###
cleanEx()
options(digits = 7L)
base::cat("Time elapsed: ", proc.time() - base::get("ptime", pos = 'CheckExEnv'),"\n")
grDevices::dev.off()
###
### Local variables: ***
### mode: outline-minor ***
### outline-regexp: "\\(> \\)?### [*]+" ***
### End: ***
quit('no')
