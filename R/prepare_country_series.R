#' Preparar la serie de un país
#'
#' Convierte la tabla de casos en formato largo y devuelve la serie temporal
#' para un país concreto filtrando por un número mínimo de casos.
#'
#' @param data Tabla de casos devuelta por [load_case_data()].
#' @param country Nombre del país a filtrar. Por defecto "Spain".
#' @param min_cases Umbral mínimo de casos para iniciar la serie.
#' @return Una lista con los elementos `series` (data.frame con fechas y casos),
#'   `start_date` (fecha de inicio) e `infected` (vector numérico).
#' @examples
#' data <- load_case_data()
#' serie <- prepare_country_series(data, country = "Spain", min_cases = 100)
#' head(serie$series)
#' @importFrom dplyr filter arrange transmute mutate
#' @importFrom reshape2 melt
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
