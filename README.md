# SIR
modelo epidemiológico SIR

## Requisitos
- R (>= 3.6)
- Paquetes: `dplyr`, `lubridate`, `tidyr`, `reshape2`, `xts`, `deSolve`

Puedes instalarlos con:

```r
install.packages(c("dplyr", "lubridate", "tidyr", "reshape2", "xts", "deSolve"))
```

## Cómo ejecutar el modelo
1. Abre R en la raíz del proyecto.
2. Ejecuta el script principal:

```r
source("R/modelo.R")
```

El script descarga los casos confirmados de COVID-19 y estima los parámetros del modelo SIR para España.

## Resultados esperados
- Gráficas de la serie de infectados y del ajuste del modelo SIR.
- Estimaciones de los parámetros `beta`, `gamma` y del número de reproducción básico (`R0`).
