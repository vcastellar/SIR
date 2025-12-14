

compute_growth_metrics_base <- function(infected) {
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
