# SIR
modelo epidemiológico SIR

## ¿Qué es el modelo SIR?
El modelo **SIR** es un sistema de ecuaciones diferenciales:

$$
\begin{aligned}
\frac{dS}{dt} &= -\beta \cdot \frac{S \cdot I}{N}, \\
\frac{dI}{dt} &= \beta \cdot \frac{S \cdot I}{N} - \gamma \cdot I, \\
\frac{dR}{dt} &= \gamma \cdot I,
\end{aligned}
\qquad
N = S(t)+I(t)+R(t)
$$

que divide a la población en tres compartimentos: 

- susceptibles (S)
- infectados (I)
- recuperados (R).

Los susceptibles pueden infectarse con una tasa de transmisión `beta`, los infectados se recuperan con una tasa `gamma`, y la población total se mantiene constante. En este proyecto el script `R/modelo.R` resuelve el sistema con `deSolve`, estima `beta` y `gamma` mediante mínimos cuadrados y calcula el número de reproducción básico `R0 = beta / gamma` a partir de los casos confirmados de COVID-19 en España.

### Curva logística
Además del modelo SIR, el script ajusta una curva logística simple sobre la serie de infectados confirmados. La logística modela un crecimiento inicial exponencial que se desacelera al aproximarse a una capacidad máxima (el parámetro de saturación) y se usa como alternativa parsimoniosa para representar la trayectoria acumulada de casos. El ajuste se realiza por mínimos cuadrados (función `fit_logistic_curve`) y se grafica junto con los datos observados.

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
