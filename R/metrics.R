get_role <- function(sim, role) {
  var <- sim$model$roles[[role]]
  if (is.null(var)) {
    stop("Model does not define role: ", role)
  }

  if (var %in% names(sim$states)) {
    sim$states[[var]]
  } else {
    sim$flows[[var]]
  }
}


## ============================================================================
## Epidemiological metrics based on roles
## ============================================================================

#' Peak incidence
#'
#' @description
#' Computes the maximum value of the incidence curve and the time at which
#' it occurs.
#'
#' This metric requires the epidemiological role \code{"incidence"}.
#'
#' @param sim An object of class \code{"sim_epi"}.
#'
#' @return
#' A named list with elements:
#' \itemize{
#'   \item \code{time}: Time at which peak incidence occurs.
#'   \item \code{value}: Maximum incidence value.
#' }
#'
#' @examples
#' sim <- simulate_epi(
#'   model = SI_MODEL,
#'   times = 0:100,
#'   parms = c(beta = 0.3),
#'   init  = c(S = 999, I = 1)
#' )
#'
#' peak_incidence(sim)
#'
#' @export
peak_incidence <- function(sim) {

  stopifnot(inherits(sim, "sim_epi"))

  inc <- get_role(sim, "incidence")
  t   <- sim$flows$time

  i <- which.max(inc)

  list(
    time  = t[i],
    value = inc[i]
  )
}


#' Time to peak incidence
#'
#' @description
#' Returns the time at which incidence reaches its maximum.
#'
#' This is a convenience wrapper around \code{peak_incidence()}.
#'
#' @param sim An object of class \code{"sim_epi"}.
#'
#' @return
#' A numeric scalar giving the time to peak incidence.
#'
#' @examples
#' sim <- simulate_epi(
#'   model = SI_MODEL,
#'   times = 0:100,
#'   parms = c(beta = 0.3),
#'   init  = c(S = 999, I = 1)
#' )
#'
#' time_to_peak(sim)
#'
#' @export
time_to_peak <- function(sim) {
  peak_incidence(sim)$time
}


#' Peak prevalence
#'
#' @description
#' Computes the maximum number of infectious individuals during the epidemic.
#'
#' This metric requires the epidemiological role \code{"infectious"}.
#'
#' @param sim An object of class \code{"sim_epi"}.
#'
#' @return
#' A numeric scalar giving the peak prevalence.
#'
#' @examples
#' sim <- simulate_epi(
#'   model = SIR_MODEL,
#'   times = 0:160,
#'   parms = c(beta = 0.3, gamma = 0.1),
#'   init  = c(S = 999, I = 1, R = 0)
#' )
#'
#' peak_prevalence(sim)
#'
#' @export
peak_prevalence <- function(sim) {

  stopifnot(inherits(sim, "sim_epi"))

  I <- get_role(sim, "infectious")

  max(I, na.rm = TRUE)
}


#' Attack rate
#'
#' @description
#' Computes the cumulative number of infections over the course of the epidemic,
#' defined as the integral (sum) of the incidence curve.
#'
#' This metric requires the epidemiological role \code{"incidence"}.
#'
#' @param sim An object of class \code{"sim_epi"}.
#'
#' @return
#' A numeric scalar giving the attack rate.
#'
#' @details
#' For models with permanent immunity (e.g. SIR, SEIR), this corresponds to the
#' final epidemic size. For models with waning immunity (e.g. SIRS, SEIRS), the
#' attack rate measures total infections and may exceed the population size.
#'
#' @examples
#' sim <- simulate_epi(
#'   model = SIR_MODEL,
#'   times = 0:160,
#'   parms = c(beta = 0.3, gamma = 0.1),
#'   init  = c(S = 999, I = 1, R = 0)
#' )
#'
#' attack_rate(sim)
#'
#' @export
attack_rate <- function(sim) {

  stopifnot(inherits(sim, "sim_epi"))

  inc <- get_role(sim, "incidence")

  sum(inc, na.rm = TRUE)
}


#' Initial exponential growth rate
#'
#' @description
#' Estimates the early exponential growth rate of the epidemic by fitting a
#' log-linear model to the initial segment of the incidence curve.
#'
#' This metric characterises the initial exponential phase of the epidemic.
#'
#' Requires the epidemiological role \code{"incidence"}.
#'
#' @param sim An object of class \code{"sim_epi"}.
#' @param n Integer. Number of initial time points to use for estimation.
#'
#' @return
#' A numeric scalar giving the initial exponential growth rate.
#'
#' @examples
#' sim <- simulate_epi(
#'   model = SI_MODEL,
#'   times = 0:50,
#'   parms = c(beta = 0.4),
#'   init  = c(S = 999, I = 1)
#' )
#'
#' initial_growth_rate(sim, n = 10)
#'
#' @export
initial_growth_rate <- function(sim, n = 7) {

  stopifnot(inherits(sim, "sim_epi"))
  stopifnot(n >= 2)

  inc <- get_role(sim, "incidence")[seq_len(n)]
  t   <- sim$flows$time[seq_len(n)]

  if (any(inc <= 0)) {
    stop("Incidence must be positive to estimate growth rate.")
  }

  fit <- stats::lm(log(inc) ~ t)

  coef(fit)[2]
}


#' Initial doubling time
#'
#' @description
#' Computes the epidemic doubling time during the initial exponential growth
#' phase.
#'
#' This is defined as \eqn{\log(2) / r}, where \eqn{r} is the initial exponential
#' growth rate.
#'
#' @param sim An object of class \code{"sim_epi"}.
#'
#' @return
#' A numeric scalar giving the initial doubling time.
#'
#' @examples
#' sim <- simulate_epi(
#'   model = SI_MODEL,
#'   times = 0:50,
#'   parms = c(beta = 0.4),
#'   init  = c(S = 999, I = 1)
#' )
#'
#' initial_doubling_time(sim)
#'
#' @export
initial_doubling_time <- function(sim) {

  r <- initial_growth_rate(sim)

  log(2) / r
}


#' Instantaneous growth rate
#'
#' @description
#' Computes the time-varying exponential growth rate of the epidemic based on
#' the incidence curve.
#'
#' This metric captures departures from exponential growth due to susceptible
#' depletion or other nonlinear effects.
#'
#' Requires the epidemiological role \code{"incidence"}.
#'
#' @param sim An object of class \code{"sim_epi"}.
#'
#' @return
#' A data.frame with columns:
#' \itemize{
#'   \item \code{time}: Time points.
#'   \item \code{r}: Instantaneous growth rate.
#' }
#'
#' @examples
#' sim <- simulate_epi(
#'   model = SI_MODEL,
#'   times = 0:100,
#'   parms = c(beta = 0.3),
#'   init  = c(S = 999, I = 1)
#' )
#'
#' head(instantaneous_growth_rate(sim))
#'
#' @export
instantaneous_growth_rate <- function(sim) {

  stopifnot(inherits(sim, "sim_epi"))

  inc <- get_role(sim, "incidence")
  t   <- sim$flows$time

  r <- diff(log(inc)) / diff(t)

  data.frame(
    time = t[-1],
    r = r
  )
}


#' Time-varying doubling time
#'
#' @description
#' Computes the time-varying epidemic doubling time based on the instantaneous
#' growth rate.
#'
#' As the epidemic approaches its peak, the doubling time increases and tends
#' to infinity when growth vanishes.
#'
#' Requires the epidemiological role \code{"incidence"}.
#'
#' @param sim An object of class \code{"sim_epi"}.
#'
#' @return
#' A data.frame with columns:
#' \itemize{
#'   \item \code{time}: Time points.
#'   \item \code{doubling_time}: Time-varying doubling time.
#' }
#'
#' @examples
#' sim <- simulate_epi(
#'   model = SI_MODEL,
#'   times = 0:100,
#'   parms = c(beta = 0.3),
#'   init  = c(S = 999, I = 1)
#' )
#'
#' head(doubling_time_ts(sim))
#'
#' @export
doubling_time_ts <- function(sim) {

  gr <- instantaneous_growth_rate(sim)

  Td <- log(2) / gr$r
  Td[gr$r <= 0] <- Inf

  data.frame(
    time = gr$time,
    doubling_time = Td
  )
}
