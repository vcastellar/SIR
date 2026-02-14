library(testthat)
library(SIR)




# ------------------------------------------------------------------------------
test_that("build-in basic model", {
  seird_rhs <- function(time, state, parms) {
    with(as.list(c(state, parms)), {

      N <- S + E + I + R

      lambda <- beta * S * I / N

      dS <- -lambda
      dE <-  lambda - sigma * E
      dI <-  sigma * E - gamma * I - mu * I
      dR <-  gamma * I
      dD <-  mu * I

      list(
        c(dS, dE, dI, dR, dD),
        incidence = lambda,
        deaths    = mu * I
      )
    })
  }

  seird_model <- epi_model(
    name        = "SEIRD",
    rhs         = seird_rhs,
    states = c("S", "E", "I", "R", "D"),
    par_names   = c("beta", "sigma", "gamma", "mu"),
    flows       = c("incidence", "deaths"),
  )

  expect_s3_class(seird_model, "epi_model")
})
#-------------------------------------------------------------------------------



# ------------------------------------------------------------------------------
test_that("build-in basic model", {
  seird_rhs <- function(time, state, parms) {
    with(as.list(c(state, parms)), {

      N <- S + E + I + R

      lambda <- beta * S * I / N

      dS <- -lambda
      dE <-  lambda - sigma * E
      dI <-  sigma * E - gamma * I - mu * I
      dR <-  gamma * I
      dD <-  mu * I

      list(
        c(dS, dE, dI, dR, dD),
        incidence = lambda,
        deaths    = mu * I
      )
    })
  }

  seird_model <- epi_model(
    name        = "SEIRD",
    rhs         = seird_rhs,
    states = c("S", "E", "I", "R", "D"),
    par_names   = c("beta", "sigma", "gamma", "mu"),
    flows       = c("incidence", "deaths")
  )

  expect_s3_class(seird_model, "epi_model")
})
#-------------------------------------------------------------------------------


# ------------------------------------------------------------------------------
test_that("simulate_epi works for seird_model", {

  sim <- simulate_epi(
    model = seird_model,
    times = 0:200,
    parms = c(beta = 0.4, sigma = 0.2, gamma = 0.1, mu = 0.02),
    init  = c(S = 9999, E = 0, I = 1, R = 0, D = 0)
  )

  expect_s3_class(sim, "sim_epi")
  expect_equal(nrow(sim$states), 201)
})

# ------------------------------------------------------------------------------



# ------------------------------------------------------------------------------
test_that("print, plot and summary do not crash for seird_model", {

  expect_output(print(sim))

  expect_silent({
    grDevices::pdf(NULL)
    plot(sim)
    plot(sim, what = "I")
    grDevices::dev.off()
  })

  expect_true(is.list(summary(sim)))
})
# ------------------------------------------------------------------------------


# ------------------------------------------------------------------------------
test_that("metricks work for seird_model", {

  expect_error(print(peak_incidence(sim)))
  expect_error(print(time_to_peak(sim)))
  expect_error(print(peak_prevalence(sim)))


})
#-------------------------------------------------------------------------------



#-------------------------------------------------------------------------------
test_that("run_epi_app works for seird_model", {
  run_epi_app(models = list(SEIRD = seird_model))
})





# ------------------------------------------------------------------------------
test_that("built-in models exist", {

  expect_s3_class(SIR_MODEL,  "epi_model")
  expect_s3_class(SIRS_MODEL, "epi_model")
  expect_s3_class(SEIR_MODEL, "epi_model")
})

# ------------------------------------------------------------------------------
test_that("simulate_epi works for SIR", {

  sim <- simulate_epi(
    model = SIR_MODEL,
    times = 0:200,
    parms = c(beta = 0.3, gamma = 0.1),
    init = SIR_MODEL$init
  )

  expect_s3_class(sim, "sim_epi")
  expect_equal(nrow(sim$states), 201)
})

# ------------------------------------------------------------------------------
test_that("print, plot and summary do not crash", {

  expect_output(print(sim))

  expect_silent({
    grDevices::pdf(NULL)
    plot(sim)
    plot(sim, what = "I")
    grDevices::dev.off()
  })

  expect_output(summary(sim))
})

# ------------------------------------------------------------------------------
test_that("fit_epi_model works end-to-end", {

  init <- c(S = 99985, I = 15, R = 0, C = 15)

  sim <- simulate_epi(
    model = SIR_MODEL,
    init = init,
    parms = c(beta = 0.25, gamma = 0.1),
    times = 0:200,
    obs = "poisson",
    seed = 1
  )

  fit <- fit_epi_model(
    x = sim$incidence$inc,
    target = "incidence",
    model = SIR_MODEL,
    init = init,
    seed = 1
  )

  expect_s3_class(fit, "fit_epi_model")
  expect_true(all(is.finite(unlist(fit$par))))
})



# ------------------------------------------------------------------------------
test_that("invalid inputs raise errors", {

  expect_error(
    simulate_epi(
      model = SIR_MODEL,
      parms = c(beta = 0.3),
      init = c(S = 99990, I = 10, R = 0, C = 10),
      method = "rk4"
    )
  )

  expect_error(
    simulate_epi(
      model = SIR_MODEL,
      times = 0:200,
      time_unit = "week",
      parms = c(beta = 0.30, gamma = 0.10),
      init  = c(S = 1e6 - 10, I = 10, R = 0, C = 10),
      seed  = 1,
      method = 'bannana'
    )
  )
})

#-------------------------------------------------------------------------------
