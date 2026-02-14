library(testthat)
library(SIR)

make_seird_model <- function() {
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

  epi_model(
    name      = "SEIRD",
    rhs       = seird_rhs,
    states    = c("S", "E", "I", "R", "D"),
    par_names = c("beta", "sigma", "gamma", "mu"),
    flows     = c("incidence", "deaths")
  )
}

test_that("custom model constructor works", {
  seird_model <- make_seird_model()
  expect_s3_class(seird_model, "epi_model")
})

test_that("simulate_epi works for custom SEIRD model", {
  seird_model <- make_seird_model()

  sim <- simulate_epi(
    model = seird_model,
    times = 0:200,
    parms = c(beta = 0.4, sigma = 0.2, gamma = 0.1, mu = 0.02),
    init  = c(S = 9999, E = 0, I = 1, R = 0, D = 0)
  )

  expect_s3_class(sim, "sim_epi")
  expect_equal(nrow(sim$states), 201)
})

test_that("print, plot and summary work for custom SEIRD model", {
  seird_model <- make_seird_model()
  sim <- simulate_epi(
    model = seird_model,
    times = 0:200,
    parms = c(beta = 0.4, sigma = 0.2, gamma = 0.1, mu = 0.02),
    init  = c(S = 9999, E = 0, I = 1, R = 0, D = 0)
  )

  expect_output(print(sim))
  expect_silent({
    grDevices::pdf(NULL)
    plot(sim)
    plot(sim, what = "I")
    grDevices::dev.off()
  })

  expect_true(is.list(summary(sim)))
})

test_that("metrics work for custom SEIRD model", {
  seird_model <- make_seird_model()
  sim <- simulate_epi(
    model = seird_model,
    times = 0:200,
    parms = c(beta = 0.4, sigma = 0.2, gamma = 0.1, mu = 0.02),
    init  = c(S = 9999, E = 0, I = 1, R = 0, D = 0)
  )

  expect_true(is.numeric(peak_incidence(sim)))
  expect_true(is.numeric(time_to_peak(sim)))
  expect_true(is.numeric(peak_prevalence(sim)))
})

test_that("built-in models exist", {
  expect_s3_class(SIR_MODEL,  "epi_model")
  expect_s3_class(SIRS_MODEL, "epi_model")
  expect_s3_class(SEIR_MODEL, "epi_model")
})

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

test_that("print, plot and summary do not crash for SIR", {
  sim <- simulate_epi(
    model = SIR_MODEL,
    times = 0:200,
    parms = c(beta = 0.3, gamma = 0.1),
    init = SIR_MODEL$init
  )

  expect_output(print(sim))
  expect_silent({
    grDevices::pdf(NULL)
    plot(sim)
    plot(sim, what = "I")
    grDevices::dev.off()
  })

  expect_output(summary(sim))
})

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
      method = "bannana"
    )
  )
})
