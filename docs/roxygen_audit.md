# Auditoría de consistencia entre roxygen y funciones

Fecha: 2026-02-14

## Inconsistencias detectadas

1. **`remove_epi_model()` no documenta su argumento `name` en roxygen.**
   - La función recibe `name`, pero su bloque roxygen solo tiene título y `@export`.
   - Archivo: `R/registry_api.R`.

2. **`get_model()` no documenta su argumento `name` en roxygen.**
   - La función recibe `name`, pero su bloque roxygen solo tiene título y `@export`.
   - Archivo: `R/registry_api.R`.

3. **La documentación de `SIR_MODEL` omite un flujo declarado por la implementación (`recovery`).**
   - En la sección "Model variables" se documenta `incidence`, pero no `recovery`.
   - Sin embargo, `sir_rhs()` retorna ambos flujos (`incidence` y `recovery`) y `SIR_MODEL` declara ambos en `flows`.
   - Archivo: `R/MODEL_SIR.R`.

## Observación

- Los RHS internos (`*_rhs`) marcados con `@noRd` no tienen `@param`, pero al ser internos y sin Rd público, esto no bloquea la documentación pública del paquete.
