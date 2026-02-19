get_derived <- function(sim, variable) {
  derived_data <- sim$derived
  if (is.null(derived_data) || !"time" %in% names(derived_data)) {
    stop("Simulation does not define derived variables.")
  }

  if (!variable %in% names(derived_data)) {
    stop("Simulation does not define derived variable: ", variable)
  }

  derived_data[[variable]]
}




#' Peak incidence
#'
#' @description
#' Computes the maximum value of an incidence curve and the time at which
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
#' @param incidence Numeric vector giving the incidence curve.
#' @param time Numeric vector of the same length giving time points.
#'
#' @return
#' A named list with elements \code{time} and \code{value}.
#'
#' @examples
#' inc  <- c(1, 3, 7, 5, 2)
#' time <- 0:4
#' peak_incidence(inc, time)
#'
#' @export
peak_incidence <- function(incidence, time) {

  stopifnot(is.numeric(incidence),
            is.numeric(time),
            length(incidence) == length(time))

  i <- which.max(incidence)

  list(
    time  = time[i],
    value = incidence[i]
  )
}


#' Time to peak incidence
#'
#' @description
#' Returns the time at which incidence reaches its maximum.
#'
#' @param incidence Numeric vector giving the incidence curve.
#' @param time Numeric vector of time points.
#'
#' @return
#' A numeric scalar giving the time to peak incidence.
#'
#' @examples
#' inc  <- c(1, 4, 6, 3)
#' time <- 0:3
#' time_to_peak(inc, time)
#'
#' @export
time_to_peak <- function(incidence, time) {
  peak_incidence(incidence, time)$time
}


#' Peak prevalence
#'
#' @description
#' Computes the maximum value of a prevalence curve.
#'
#' @details
#' Let \eqn{I(t)} denote the prevalence over time.
#' The peak prevalence is defined as
#' \deqn{
#' I_{\max} = \max_t I(t).
#' }
#'
#' @param prevalence Numeric vector representing prevalence over time.
#'
#' @return
#' A numeric scalar giving the peak prevalence.
#'
#' @examples
#' I <- c(1, 5, 8, 4, 2)
#' peak_prevalence(I)
#'
#' @export
peak_prevalence <- function(prevalence) {

  stopifnot(is.numeric(prevalence))

  max(prevalence, na.rm = TRUE)
}



#' Attack rate
#'
#' @description
#' Computes the cumulative number of events as the sum of an incidence curve.
#'
#' @details
#' Let \eqn{\lambda(t)} denote the incidence function.
#' The attack rate is defined as
#' \deqn{
#' AR = \int_0^T \lambda(t)\, dt.
#' }
#'
#' In discrete time, this is approximated by:
#' \deqn{
#' AR \approx \sum_t \lambda(t).
#' }
#'
#' @param incidence Numeric vector giving the incidence curve.
#'
#' @return
#' A numeric scalar giving the attack rate.
#'
#' @examples
#' inc <- c(1, 2, 3, 4)
#' attack_rate(inc)
#'
#' @export
attack_rate <- function(incidence) {

  stopifnot(is.numeric(incidence))

  sum(incidence, na.rm = TRUE)
}


#' Initial exponential growth rate
#'
#' @description
#' Estimates the early exponential growth rate by fitting a log-linear model
#' to the initial segment of an incidence curve.
#'
#' @details
#' Assuming exponential growth:
#' \deqn{
#' \lambda(t) \approx C e^{r t}.
#' }
#'
#' Taking logarithms:
#' \deqn{
#' \log \lambda(t) = \log C + r t.
#' }
#'
#' @param incidence Numeric vector of incidence values.
#' @param time Numeric vector of time points.
#' @param n Integer. Number of initial time points used for estimation.
#'
#' @return
#' A numeric scalar giving the estimated growth rate.
#'
#' @export
initial_growth_rate <- function(incidence, time, n = 7) {

  stopifnot(is.numeric(incidence),
            is.numeric(time),
            length(incidence) == length(time),
            n >= 2)

  inc <- incidence[seq_len(n)]
  t   <- time[seq_len(n)]

  if (any(inc <= 0)) {
    stop("Incidence must be positive to estimate growth rate.")
  }

  fit <- stats::lm(log(inc) ~ t)

  coef(fit)[2]
}



#' Initial doubling time
#'
#' @description
#' Computes the epidemic doubling time during the initial exponential phase.
#'
#' @details
#' If incidence grows exponentially with rate \eqn{r},
#' the doubling time is:
#' \deqn{
#' T_d = \frac{\log 2}{r}.
#' }
#'
#' @param incidence Numeric vector of incidence values.
#' @param time Numeric vector of time points.
#' @param n Integer. Number of initial time points.
#'
#' @return
#' A numeric scalar giving the doubling time.
#'
#' @export
initial_doubling_time <- function(incidence, time, n = 7) {

  r <- initial_growth_rate(incidence, time, n)

  log(2) / r
}



#' Instantaneous growth rate
#'
#' @description
#' Computes the time-varying exponential growth rate from an incidence curve.
#'
#' @details
#' The growth rate is approximated by:
#' \deqn{
#' r(t_i) \approx \frac{\log \lambda(t_{i+1}) - \log \lambda(t_i)}
#' {t_{i+1} - t_i}.
#' }
#'
#' Optional smoothing via moving average can be applied.
#'
#' @param incidence Numeric vector of incidence values.
#' @param time Numeric vector of time points.
#' @param window Integer. Size of moving average window.
#' @param offset Numeric. Small constant added to avoid \code{log(0)}.
#'
#' @return
#' A data.frame with columns \code{time} and \code{r}.
#'
#' @export
instantaneous_growth_rate <- function(incidence,
                                      time,
                                      window = 1,
                                      offset = 0.5) {

  stopifnot(is.numeric(incidence),
            is.numeric(time),
            length(incidence) == length(time))

  inc <- incidence
  t   <- time

  if (window > 1) {
    inc <- as.numeric(stats::filter(inc, rep(1/window, window), sides = 2))
    valid <- !is.na(inc)
    inc <- inc[valid]
    t   <- t[valid]
  }

  log_inc <- log(inc + offset)

  r <- diff(log_inc) / diff(t)

  data.frame(
    time = t[-1],
    r = r
  )
}



#' Time-varying doubling time
#'
#' @description
#' Computes the time-varying doubling time from an incidence curve.
#'
#' @details
#' Given growth rate \eqn{r(t)}, the doubling time is:
#' \deqn{
#' T_d(t) = \frac{\log 2}{r(t)}.
#' }
#'
#' When \eqn{r(t) \le 0}, doubling time is set to \code{Inf}.
#'
#' @param incidence Numeric vector of incidence values.
#' @param time Numeric vector of time points.
#' @param window Integer. Moving average window.
#' @param offset Numeric. Small constant to avoid log(0).
#'
#' @return
#' A data.frame with columns \code{time} and \code{doubling_time}.
#'
#' @export
doubling_time_ts <- function(incidence,
                             time,
                             window = 1,
                             offset = 0.5) {

  gr <- instantaneous_growth_rate(incidence, time,
                                  window = window,
                                  offset = offset)

  Td <- log(2) / gr$r
  Td[gr$r <= 0] <- Inf

  data.frame(
    time = gr$time,
    doubling_time = Td
  )
}
