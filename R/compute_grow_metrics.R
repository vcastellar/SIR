#' Compute short-horizon growth metrics from an infection time series
#' @name compute_grow_metrics
#' @description
#' Computes two related growth metrics from a univariate time series of (typically)
#' **cumulative** counts:
#'
#' 1. **Daily percent increase** (`inc_pct`): for each time step \eqn{t},
#'    \deqn{inc\_pct_t = (x_{t+1} - x_t) / x_t}
#'    This is the one-step *relative* change, using the current value \eqn{x_t}
#'    as the denominator.
#'
#' 2. **Doubling time in days** (`doubling_days`): first computes a 4-point
#'    **left-aligned** moving average of `inc_pct` (call it \eqn{\alpha_t}),
#'    then converts that average growth rate into an implied doubling time
#'    under constant-rate exponential growth:
#'    \deqn{\alpha_t = mean(inc\_pct_t, inc\_pct_{t+1}, inc\_pct_{t+2}, inc\_pct_{t+3})}
#'    \deqn{doubling\_days_t = \log(2) / \log(1 + \alpha_t)}
#'
#' The doubling-time formula comes from solving \eqn{2 = (1+\alpha)^k} for \eqn{k}.
#' If \eqn{1+\alpha_t \le 0} or \eqn{\alpha_t = 0} then the logarithm or division
#' may yield non-finite values; infinite results are replaced with `NA`.
#'
#' @param infected A univariate time series (typically a `ts`) of counts ordered in
#'   time. In most epidemiological uses this is a **cumulative** count series, but
#'   the function will work for any numeric series with positive values.
#'
#' @return A named list with two numeric vectors:
#'   \describe{
#'     \item{inc_pct}{Vector of length `length(infected) - 1` with the one-step
#'       relative increases.}
#'     \item{doubling_days}{Vector of length `length(infected) - 4` with the implied
#'       doubling time (in days) computed from the 4-point moving average of `inc_pct`.}
#'   }
#'
#' @details
#' **Windowing / alignment**
#'
#' The percent increase is computed for adjacent pairs \eqn{(x_t, x_{t+1})}.
#' The 4-point moving average of `inc_pct` is **left-aligned**, meaning the value at
#' position \eqn{t} summarizes growth from \eqn{t} through \eqn{t+3}.
#'
#' **Practical notes**
#'
#' - If `infected` contains zeros, `inc_pct` will be `Inf` or `NaN` at those points.
#'   This will propagate into `alpha` and `doubling_days`.
#' - Negative or decreasing series can produce \eqn{\alpha < 0}, for which
#'   \eqn{\log(1+\alpha)} may be undefined when \eqn{1+\alpha \le 0}.
#'
#' @examples
#' infected <- ts(c(100, 110, 121, 133, 146, 161))
#' out <- compute_growth_metrics_base(infected)
#' out$inc_pct
#' out$doubling_days

compute_growth_metrics <- function(infected) {
  infected <- as.numeric(infected)

  n <- length(infected)
  if (n < 5) {
    stop("Se necesitan al menos 5 observaciones para calcular inc_pct (2) y alpha (4).")
  }

  # 1) Incremento porcentual diario: (x[t+1] - x[t]) / x[t], para t = 1..n-1
  inc_pct <- (infected[2:n] - infected[1:(n - 1)]) / infected[1:(n - 1)]

  # 2) Media móvil de 4 sobre inc_pct (equivalente a rollapply(width=4, mean, align='left'))
  #    Con embed: filas = ventanas, columnas = valores
  alpha <- rowMeans(embed(inc_pct, 4))

  # 3) Días de duplicación
  doubling_days <- log(2) / log(1 + alpha)
  doubling_days[is.infinite(doubling_days)] <- NA_real_

  list(
    inc_pct = inc_pct,
    doubling_days = doubling_days,
    alpha = alpha
  )
}
