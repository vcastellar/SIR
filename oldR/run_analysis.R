#
# run_analysis <- function(country = "Spain", population = 47e6, min_cases = 100, horizon_days = 90) {
#   data <- load_case_data()
#   prepared <- prepare_country_series(data, country, min_cases)
#   growth <- compute_growth_metrics(prepared$infected)
#
#   sir <- fit_sir_model(prepared$infected, population, horizon_days)
#   r0 <- sir$parameters["beta"] / sir$parameters["gamma"]
#   logistic <- fit_logistic_curve(prepared$infected, prepared$series$date)
#
#   plot_results(sir, logistic, prepared$series$date, prepared$infected, population)
#
#   list(
#     country = country,
#     start_date = prepared$start_date,
#     sir_parameters = sir$parameters,
#     r0 = r0,
#     logistic_parameters = logistic$parameters,
#     growth = growth
#   )
# }
#
# if (sys.nframe() == 0) {
#   results <- run_analysis()
#   print(results$sir_parameters)
#   message(sprintf("R0 estimado: %.2f", results$r0))
# }
