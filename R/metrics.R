## ============================================================================
## Helper: extract variable by epidemiological role
## ============================================================================

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
#' @details
#' Let \eqn{\lambda(t)} denote the incidence function.
#' The peak incidence is defined as
#' \deqn{
#' \lambda_{\max} = \max_t \lambda(t),
#' }
#' and the time to peak incidence is
#' \deqn{
#' t^* = \arg\max_t \lambda(t).
#' }
#'
#' This metric requires the epidemiological role \code{"incidence"}.
#'
#' @param sim An object of class \code{"sim_epi"}.
#'
#' @return
#' A named list with elements \code{time} and \code{value}.
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
#' @details
#' Let \eqn{\lambda(t)} be the incidence curve.
#' The time to peak incidence is defined as
#' \deqn{
#' t^* = \arg\max_t \lambda(t).
#' }
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
#'   times = 0:80,
#'   parms = c(beta = 0.25),
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
#' @details
#' Let \eqn{I(t)} denote the number of infectious individuals over time.
#' The peak prevalence is defined as
#' \deqn{
#' I_{\max} = \max_t I(t).
#' }
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
#' @details
#' Let \eqn{\lambda(t)} denote the incidence function.
#' The attack rate is defined as
#' \deqn{
#' AR = \int_0^T \lambda(t)\, dt.
#' }
#'
#' In discrete time simulations, this integral is approximated by
#' \deqn{
#' AR \approx \sum_{t=0}^T \lambda(t).
#' }
#'
#' For models with permanent immunity, this corresponds to the final epidemic
#' size.
#'
#' This metric requires the epidemiological role \code{"incidence"}.
#'
#' @param sim An object of class \code{"sim_epi"}.
#'
#' @return
#' A numeric scalar giving the attack rate.
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
#' @details
#' During the early epidemic phase, incidence is assumed to grow exponentially:
#' \deqn{
#' \lambda(t) \approx C e^{r t},
#' }
#' where \eqn{r} is the exponential growth rate.
#'
#' Taking logarithms yields
#' \deqn{
#' \log \lambda(t) = \log C + r t.
#' }
#'
#' Requires the epidemiological role \code{"incidence"}.
#'
#' @param sim An object of class \code{"sim_epi"}.
#' @param n Integer. Number of initial time points.
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
#' @details
#' If incidence grows exponentially, the doubling time is
#' \deqn{
#' T_d = \frac{\log 2}{r}.
#' }
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
#' Computes the time-varying exponential growth rate of the epidemic, with
#' optional smoothing via a moving average to reduce noise.
#'
#' @details
#' The instantaneous growth rate is calculated as the derivative of the
#' log-incidence. To handle noise, a moving average of size \code{window}
#' is applied to the incidence before the transformation:
#' \deqn{
#' \lambda_{smooth}(t) = \frac{1}{k} \sum_{j=-(k-1)/2}^{(k-1)/2} \lambda(t+j)
#' }
#' where \eqn{k} is the window size. The growth rate is then:
#' \deqn{
#' r(t_i) \approx \frac{\log \lambda_{smooth}(t_{i+1}) - \log \lambda_{smooth}(t_i)}{t_{i+1} - t_i}
#' }
#'
#' @param sim An object of class \code{"sim_epi"}.
#' @param window Integer. Size of the moving average window (must be odd for
#'   centered smoothing). Default is 1 (no smoothing).
#' @param offset Numeric. Small constant added to incidence to avoid \code{log(0)}.
#'
#' @return
#' A data.frame with columns \code{time} and \code{r}.
#'
#' @examples
#' sim <- simulate_epi(
#'   model = SI_MODEL,
#'   times = 0:50,
#'   parms = c(beta = 0.4),
#'   init  = c(S = 999, I = 1)
#' )
#' # Sin suavizado
#' plot(instantaneous_growth_rate(sim, window = 1))
#'
#' # Con suavizado (ventana de 7 días)
#' plot(instantaneous_growth_rate(sim, window = 7))
#'
#' @export
instantaneous_growth_rate <- function(sim, window = 1, offset = 0.5) {

  stopifnot(inherits(sim, "sim_epi"))

  inc <- get_role(sim, "incidence")
  t   <- sim$flows$time

  # 1. Aplicar suavizado (Moving Average)
  if (window > 1) {
    # Usamos stats::filter para una media móvil centrada
    inc <- as.numeric(stats::filter(inc, rep(1/window, window), sides = 2))

    # El filtrado produce NAs en los extremos, los eliminamos para el cálculo
    valid_idx <- !is.na(inc)
    inc <- inc[valid_idx]
    t   <- t[valid_idx]
  }

  # 2. Manejo de ceros mediante offset
  # Esto evita errores matemáticos sin detener la ejecución
  log_inc <- log(inc + offset)

  # 3. Cálculo de la tasa r
  r <- diff(log_inc) / diff(t)

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
#' @details
#' Given the instantaneous growth rate \eqn{r(t)}, the doubling time is
#' \deqn{
#' T_d(t) = \frac{\log 2}{r(t)}.
#' }
#'
#' When \eqn{r(t) \le 0}, the doubling time is set to infinity.
#'
#' Requires the epidemiological role \code{"incidence"}.
#'
#' @param sim An object of class \code{"sim_epi"}.
#'
#' @return
#' A data.frame with columns \code{time} and \code{doubling_time}.
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
