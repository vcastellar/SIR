# Analisis del modelo SIR para COVID-19 en España
# El script descarga las series de casos confirmados desde el repositorio
# de la Universidad Johns Hopkins, calcula métricas de crecimiento y ajusta
# un modelo SIR y uno logístico simples.

suppressPackageStartupMessages({
  library(dplyr)
  library(lubridate)
  library(tidyr)
  library(reshape2)
  library(xts)
  library(zoo)
  library(deSolve)
})

DATA_URL <- "https://raw.githubusercontent.com/CSSEGISandData/COVID-19/master/csse_covid_19_data/csse_covid_19_time_series/time_series_covid19_confirmed_global.csv"

# Carga los datos de casos confirmados y elimina las columnas irrelevantes.
load_case_data <- function(url = DATA_URL) {
  read.csv(url, stringsAsFactors = FALSE) %>%
    select(-Lat, -Long)
}

# Convierte la tabla en formato largo y devuelve la serie del país solicitado.
prepare_country_series <- function(data, country = "Spain", min_cases = 100) {
  long <- melt(data, id.vars = c("Province.State", "Country.Region")) %>%
    mutate(
      variable = gsub("^X", "", variable),
      variable = gsub("\\.", "-", variable),
      variable = gsub("-20$", "-2020", variable),
      date = as.Date(as.character(variable), format = "%m-%d-%Y")
    )

  country_ts <- long %>%
    filter(Country.Region == country) %>%
    arrange(date) %>%
    transmute(date, cases = value) %>%
    filter(cases > min_cases)

  if (nrow(country_ts) == 0) {
    stop(paste0("No hay datos por encima de ", min_cases, " casos para ", country, "."))
  }

  list(
    series = country_ts,
    start_date = country_ts$date[1],
    infected = round(country_ts$cases)
  )
}

# Calcula métricas de crecimiento: porcentaje diario y días en que se duplica.
compute_growth_metrics <- function(infected) {
  inc_pct <- rollapply(infected, width = 2, FUN = function(x) (x[2] - x[1]) / x[1], align = "left")
  alpha <- rollapply(inc_pct, width = 4, FUN = mean)
  doubling_days <- log(2) / log(1 + alpha)
  doubling_days[is.infinite(doubling_days)] <- NA
  list(inc_pct = inc_pct, doubling_days = doubling_days)
}

# Función que define el sistema de ecuaciones SIR.
SIR <- function(time, state, parameters, N) {
  par <- as.list(c(state, parameters))
  with(par, {
    dS <- -beta / N * I * S
    dI <- beta / N * I * S - gamma * I
    dR <- gamma * I
    list(c(dS, dI, dR))
  })
}

# Ajusta el modelo SIR mediante mínimos cuadrados.
fit_sir_model <- function(infected, population, horizon_days) {
  Day <- seq_along(infected)
  init <- c(S = population - infected[1], I = infected[1], R = 0)

  rss <- function(parameters) {
    names(parameters) <- c("beta", "gamma")
    out <- ode(y = init, times = Day, func = function(time, state, parameters) {
      SIR(time, state, parameters, population)
    }, parms = parameters)
    fit <- out[, 3]
    mean((infected - fit) ^ 2)
  }

  opt <- optim(c(0.5, 0.5), rss, method = "Nelder-Mead", control = list(maxit = 500, factr = 1e4))
  opt_par <- setNames(opt$par, c("beta", "gamma"))

  t <- 1:(horizon_days + length(infected))
  fit <- data.frame(ode(y = init, times = t, func = function(time, state, parameters) {
    SIR(time, state, parameters, population)
  }, parms = opt_par))

  list(parameters = opt_par, fit = fit)
}

# Ajusta una curva logística simple por mínimos cuadrados.
fit_logistic_curve <- function(infected, dates) {
  Day <- seq_along(infected)
  logistic <- function(parms, x) {
    parms[3] / (1 + exp(-(x - parms[2]) / parms[1]))
  }

  rss <- function(parms) {
    out <- logistic(parms, x = Day)
    sum((infected - out) ^ 2)
  }

  opt <- optim(c(1, 1, max(infected) * 2), rss, method = "L-BFGS-B")
  fit <- logistic(opt$par, 1:length(dates))
  list(parameters = opt$par, fit = fit)
}

# Genera las gráficas del modelo SIR y de la curva logística.
plot_results <- function(fit_sir, fit_logistic, infected_dates, infected, population) {
  fechas_modelo <- seq.Date(from = infected_dates[1], length.out = nrow(fit_sir$fit), by = "day")
  fit_sir$fit$time <- fechas_modelo

  plot(fit_sir$fit$time, fit_sir$fit$S, type = "l", lwd = 2, ylim = c(0, population))
  lines(fit_sir$fit$time, fit_sir$fit$I, lwd = 2, col = "red")
  lines(fit_sir$fit$time, fit_sir$fit$R, lwd = 2, col = "green")
  points(infected_dates, infected, pch = 19)
  grid()
  legend("left", c("Susceptibles", "Infectados", "Recobrados"), lty = 1, lwd = 2, col = 1:3, inset = 0.05)

  plot(fit_sir$fit$time, fit_sir$fit$I, type = "l", lwd = 2, col = "red", pch = 19)
  points(infected_dates, infected, pch = 19)
  grid()

  plot(infected_dates, fit_logistic$fit, type = "l", lwd = 2, col = "blue",
       main = "Curva logística ajustada", ylab = "Casos", xlab = "Fecha")
  points(infected_dates, infected, pch = 19, col = "red")
  grid()
}

# Ejecuta el flujo completo: descarga datos, ajusta modelos y devuelve un resumen.
run_analysis <- function(country = "Spain", population = 47e6, min_cases = 100, horizon_days = 90) {
  data <- load_case_data()
  prepared <- prepare_country_series(data, country, min_cases)
  growth <- compute_growth_metrics(prepared$infected)

  sir <- fit_sir_model(prepared$infected, population, horizon_days)
  r0 <- sir$parameters["beta"] / sir$parameters["gamma"]
  logistic <- fit_logistic_curve(prepared$infected, prepared$series$date)

  plot_results(sir, logistic, prepared$series$date, prepared$infected, population)

  list(
    country = country,
    start_date = prepared$start_date,
    sir_parameters = sir$parameters,
    r0 = r0,
    logistic_parameters = logistic$parameters,
    growth = growth
  )
}

if (sys.nframe() == 0) {
  results <- run_analysis()
  print(results$sir_parameters)
  message(sprintf("R0 estimado: %.2f", results$r0))
}
