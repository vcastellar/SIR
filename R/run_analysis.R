#' Ejecutar el flujo completo de análisis SIR
#'
#' Descarga los datos de casos confirmados, prepara la serie para el país
#' seleccionado, calcula métricas de crecimiento, ajusta el modelo SIR y la
#' curva logística y genera las gráficas correspondientes.
#'
#' @param country País a analizar. Por defecto "Spain".
#' @param population Tamaño total de la población susceptible.
#' @param min_cases Número mínimo de casos a partir del cual se inicia la serie.
#' @param horizon_days Días de proyección adicionales para el modelo SIR.
#' @return Una lista con la información de país, fecha de inicio, parámetros
#'   del modelo SIR (`sir_parameters`), número reproductivo básico (`r0`),
#'   parámetros de la curva logística y métricas de crecimiento.
#' @examples
#' \donttest{
#' resultados <- run_analysis(country = "Spain", population = 47e6, min_cases = 100, horizon_days = 90)
#' resultados$r0
#' }
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
