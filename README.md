# Procesamiento y análisis de espectros de emisión de rayos X

Este repositorio contiene una serie de scripts en Python desarrollados para el procesamiento y análisis de espectros experimentales y teóricos obtenidos mediante espectroscopía de emisión de rayos X. Los programas permiten combinar datos experimentales, calibrar espectros, generar espectros teóricos a partir de resultados obtenidos con el programa CTM4XAS, y comparar ambos conjuntos de datos.

## Descripción de los scripts

### suma.py

Este script genera la suma de los distintos archivos `.dat` obtenidos durante las mediciones experimentales de una misma muestra (en caso de que la medición haya sido realizada en partes).

A partir de múltiples archivos de datos, produce un único archivo final que contiene:

- los pasos del motor del cristal
- los pasos del motor del detector
- las cuentas medidas
- el tiempo de adiquisión
- las cuentas/segundo

Esto permite obtener un único espectro representativo para cada muestra a partir de varias mediciones experimentales.

---

### ajustefinal.py

Este script procesa los datos experimentales para generar un archivo calibrado en energía (el archivo de datos obtenido luego de una medición da los pasos de los motores del cristal y del detector, las cuentas obtenidas, el tiempo de adquisición y las cuentas por segundo. Es necesari).

El programa:

- procesa los datos experimentales originales
- genera un archivo `.dat` que contiene:
  - la energía de cada punto
  - el número de cuentas en cada punto
  - las cuentas normalizadas al máximo
  - el error en la energía
  - el error en la intensidad

El archivo generado se guarda con el nombre nombreoriginal_calib.dat

Además, el script permite procesar los datos teóricos al aplicar una función voigt a cada transición calculada y generar un gráfico donde se muestra:

- el espectro experimental
- el espectro teórico correspondiente

Este paso es necesario para obtener el archivo de datos experimentales calibrados, que luego será utilizado en los análisis posteriores.

---

### barras.py

Este script permite generar espectros teóricos a partir de las transiciones obtenidas con el programa CTM4XAS.

Entrada requerida:

Un archivo `.txt` con dos columnas:

1. energía de la transición  
2. intensidad de la transición  

Una vez ejecutado, el programa calcula y muestra el espectro teórico resultante. Cumple la misma función que la segunda parte del script ajustefinal.py

---

### comparar.py

Este script permite comparar espectros teóricos y experimentales.

Para su ejecución es necesario cargar:

- el archivo experimental calibrado (`*_calib.dat`)
- el archivo teórico que contiene la información de cada transición

El programa genera un gráfico donde se muestra la comparación entre:

- el espectro experimental
- el espectro teórico

---

### valores.py

Este script calcula dos magnitudes espectroscópicas importantes para cada muestra:

- IAD (Integrated Absolute Difference)
- M1 (primer momento energético)

El cálculo incluye propagación de errores a partir de los datos del espectro experimental.

---

### hs.py

Este script permite:

- generar el espectro de alto espín (HS) a partir de espectros experimentales
- calcular las poblaciones relativas de alto espín (HS) y bajo espín (LS) para las muestras de tulio

---

### graficar_hs.py

Este script genera un gráfico comparativo entre los distintos espectros experimentales, ademas del el espectro de alto espín (HS) calculado.

---

## Flujo de trabajo típico

Un flujo de trabajo común para el análisis de los datos es el siguiente:

1. **suma.py**  
   Combinar los archivos experimentales individuales en un único espectro por muestra.

2. **ajustefinal.py**  
   Calibrar el espectro experimental y generar el archivo `*_calib.dat`.

3. **barras.py**  
   Generar el espectro teórico a partir de las transiciones obtenidas con CTM4XAS.

4. **comparar.py**  
   Comparar el espectro teórico con el experimental calibrado.

5. **valores.py**  
   Calcular magnitudes espectroscópicas como el IAD y el primer momento energético.

6. **hs.py**  
   Obtener el espectro de alto espín y calcular las poblaciones de espín.

7. **graficar_hs.py**  
   Visualizar y comparar los espectros experimentales junto con el espectro HS calculado.
Esto permite visualizar directamente las diferencias entre las muestras y la contribución del estado de alto espín.

