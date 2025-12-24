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

- **susceptibles (S)**: individuos sin inmunidad al agente infeccioso, y que por tanto puede ser infectada si es expuesta al agente infeccioso
- **infectados (I)**: individuos que están infectados en un momento dado y pueden transmitir la infección a individuos de la población susceptible con la que entran en contacto
- **recuperados (R)**: individuos que son inmunes a la infección (o fallecidos), y consecuentemente no afectan a la transmisión cuando entran en contacto con otros individuos

- $\beta > 0$ representa el ratio de transmisión
- $\gamma > 0$ representa el ratio de recuperación. En otras palabras, $D = 1/\gamma$ representa la duración de la infeccion 

Los susceptibles pueden infectarse con una tasa de transmisión `beta`, los infectados se recuperan con una tasa `gamma`, y la población total se mantiene constante. En este proyecto el script `R/modelo.R` resuelve el sistema con `deSolve`, estima `beta` y `gamma` mediante mínimos cuadrados y calcula el número de reproducción básico `R0 = beta / gamma` a partir de los casos confirmados de COVID-19 en España.

## Cómo se ajustan los parámetros del modelo SIR
1. **Preparar la serie de infectados**: el script descarga los casos confirmados, los convierte a formato largo, filtra el país y aplica un umbral mínimo de casos; la serie resultante (`infected`) se redondea y sirve como dato de entrenamiento.
2. **Definir el sistema SIR**: se implementan las ecuaciones diferenciales con dos parámetros libres, `beta` (tasa de contagio) y `gamma` (tasa de recuperación), usando la población `N` para normalizar el contacto.
3. **Configurar condiciones iniciales**: se calcula el vector de estado inicial `init` con susceptibles `S = población − infectados_iniciales`, infectados iniciales `I` y recuperados `R = 0`, sobre una secuencia temporal igual a la longitud de la serie observada.
4. **Función objetivo (RSS)**: para un par candidato de parámetros, se resuelve el sistema SIR con `ode` y se extrae la serie simulada de infectados (`I`); luego se calcula el error cuadrático medio entre infectados observados y simulados. Esta función devuelve el valor a minimizar.
5. **Optimización numérica**: se llama a `optim` con método Nelder-Mead y punto de arranque `(0.5, 0.5)` para encontrar los valores de `beta` y `gamma` que minimizan el RSS. El resultado se nombra y se guarda como `opt_par`.
6. **Simulación con parámetros óptimos**: con los parámetros ajustados, se vuelve a integrar el modelo SIR en un horizonte extendido (`horizon_days + longitud_observada`) para generar trayectorias de `S`, `I` y `R` suavizadas y proyectadas.
7. **Cálculo de métricas derivadas**: el número básico de reproducción `R0` se calcula como `beta/gamma` a partir de los parámetros óptimos; estos resultados, junto con la curva ajustada, se devuelven y se usan en las gráficas.

### Curva logística
Además del modelo SIR, el script ajusta una curva logística simple sobre la serie de infectados confirmados. La logística modela un crecimiento inicial exponencial que se desacelera al aproximarse a una capacidad máxima (el parámetro de saturación) y se usa como alternativa parsimoniosa para representar la trayectoria acumulada de casos. El ajuste se realiza por mínimos cuadrados (función `fit_logistic_curve`) y se grafica junto con los datos observados.

## Requisitos
- R (>= 3.6)
- Paquetes: `deSolve`

Puedes instalarlos con:

```r
install.packages("deSolve")
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

## Referencias
- Wikipedia contributors. *Modelo SIR*. Wikipedia, The Free Encyclopedia.  
 [Modelo SIR](https://es.wikipedia.org/wiki/Modelo_SIR)
 [SIR Model Foundation](https://mat.uab.cat/matmat_antiga/PDFv2013/v2013n03.pdf)

