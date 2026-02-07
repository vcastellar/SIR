`%||%` <- function(x, y) if (is.null(x)) y else x

param_sliders_ui <- function(model) {

  lapply(model$par_names, function(p) {

    val   <- model$defaults[[p]] %||% 1
    lower <- model$lower[[p]]    %||% (val / 10)
    upper <- model$upper[[p]]    %||% (val * 10)

    shiny::sliderInput(
      inputId = paste0("par_", p),
      label   = p,
      min     = lower,
      max     = upper,
      value   = val
    )
  })
}

init_sliders_ui <- function(model) {

  lapply(model$states, function(s) {

    val <- model$init[[s]] %||% 0

    shiny::sliderInput(
      inputId = paste0("init_", s),
      label   = s,
      min     = 0,
      max     = max(1, val * 5),
      value   = val
    )
  })
}

rhs_text <- function(model) {
  paste(deparse(body(model$rhs)), collapse = "\n")
}
