library(testthat)
library(SIR)

# ------------------------------------------------------------------------------
test_that("built-in models exist", {

  expect_s3_class(SIR_MODEL,  "epi_model")
  expect_s3_class(SIRS_MODEL, "epi_model")
  expect_s3_class(SEIR_MODEL, "epi_model")
})

# ------------------------------------------------------------------------------
test_that("simulate_epi works for SIR", {

  init <- c(S = 1e5 - 10, I = 10, R = 0, C = 10)

  sim <- simulate_epi(
    model = SIR_MODEL,
    n_days = 20,
    parms = c(beta = 0.3, gamma = 0.1),
    init = init,
    obs = "none",
    method = "rk4"
  )

  expect_s3_class(sim, "sim_epi")
  expect_equal(nrow(sim$states), 21)
})

# ------------------------------------------------------------------------------
test_that("print, plot and summary do not crash", {

  init <- c(S = 99990, I = 10, R = 0, C = 10)

  sim <- simulate_epi(
    model = SIR_MODEL,
    n_days = 20,
    parms = c(beta = 0.3, gamma = 0.1),
    init = init,
    obs = "none",
    method = "rk4"
  )

  expect_output(print(sim))

  expect_silent({
    grDevices::pdf(NULL)
    plot(sim)
    plot(sim, what = "incidence")
    grDevices::dev.off()
  })

  expect_silent(invisible(summary(sim)))
})

# ------------------------------------------------------------------------------
test_that("fit_epi_model works end-to-end", {

  init <- c(S = 99985, I = 15, R = 0, C = 15)

  sim <- simulate_epi(
    model = SIR_MODEL,
    n_days = 40,
    parms = c(beta = 0.25, gamma = 0.1),
    init = init,
    obs = "poisson",
    seed = 1,
    method = "rk4"
  )

  fit <- fit_epi_model(
    x = sim$incidence_obs$inc,
    model = SIR_MODEL,
    init = list(I = 10, N = 1e5),
    n_starts = 3,
    seed = 1
  )

  expect_s3_class(fit, "fit_epi_model")
  expect_true(all(is.finite(unlist(fit$par))))
})

# ------------------------------------------------------------------------------
test_that("predict works on fitted model", {

  init <- c(S = 99990, I = 10, R = 0, C = 10)

  sim <- simulate_epi(
    model = SIR_MODEL,
    n_days = 30,
    parms = c(beta = 0.3, gamma = 0.1),
    init = init,
    obs = "poisson",
    seed = 1,
    method = "lsoda"
  )

  fit <- fit_epi_model(
    x = sim$incidence_obs$inc,
    model = SIR_MODEL,
    init = head(sim$states, n = 1)[,-1],
    n_starts = 2,
    seed = 2
  )

  pred <- predict(
    fit,
    n_days = 10,
    init = tail(sim$states, n = 1)[, -1]
  )

  expect_true(is.list(pred))
  expect_true("states" %in% names(pred))
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
})
