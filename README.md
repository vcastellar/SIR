# SIR

Modelo epidemiológico SIR aplicado a la serie de casos confirmados de COVID-19.

## Requisitos
El script depende de los siguientes paquetes de R:
- dplyr
- lubridate
- tidyr
- reshape2
- xts
- zoo
- deSolve

Instala las dependencias desde una sesión de R si no están disponibles:

```r
install.packages(c("dplyr", "lubridate", "tidyr", "reshape2", "xts", "zoo", "deSolve"))
```

## Ejecución
El análisis principal se encuentra en `R/modelo.R` y descarga automáticamente la
serie de casos global desde el repositorio público de Johns Hopkins. Para
lanzarlo con los valores predeterminados (España, umbral de 100 casos y 90 días
de horizonte), ejecuta:

```bash
Rscript R/modelo.R
```

El script ajusta un modelo SIR y una curva logística básica, calcula el número
reproductivo básico estimado (`R0`) y genera gráficos de las curvas ajustadas.

## Parámetros
Puedes modificar el país, la población total, el umbral mínimo de casos o el
horizonte de proyección invocando la función `run_analysis()` desde R, por
 ejemplo:

```r
source("R/modelo.R")
resultados <- run_analysis(country = "Italy", population = 60e6, min_cases = 200)
```
