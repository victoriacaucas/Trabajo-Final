# Generación de espectros continuos a partir de barras mediante perfiles Voigt.

import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from scipy.special import wofz

# ====================================================
# === DEFINICIÓN DE PERFILES ESPECTRALES ==============
# ====================================================

def voigt(x, x0, sigma, gamma0, gamma_slope, A, x_ref=None):
    """Perfil de Voigt normalizado, con ancho Lorentziano dependiente de la energía.

    Parámetros:
    x : array, eje de energía (eV)
    x0 : float, posición energética de la barra
    sigma : float, ancho Gaussiano (resolución instrumental)
    gamma0 : float, ancho Lorentziano base
    gamma_slope : float, variación lineal de gamma con la energía
    A : float, intensidad de la barra
    x_ref : float, energía de referencia para la variación de gamma

    Retorna el perfil Voigt evaluado en x"""

    # Si no se define referencia, se usa la posición de la barra
    if x_ref is None:
        x_ref = x0

    gamma_x = gamma0 + gamma_slope * (x - x_ref)    # ancho Lorentziano dependiente de energía
    gamma_x = np.maximum(gamma_x, 1e-3)             # el ancho no puede ser negativo

    z = ((x - x0) + 1j * gamma_x) / (sigma * np.sqrt(2))
    return A * np.real(wofz(z)) / (sigma * np.sqrt(2 * np.pi))

def barras_a_voigt(x, x0_list, A_list, sigma, gamma0, gamma_slope, x_ref=None):
    """Convierte un conjunto discreto de barras espectrales en un espectro continuo mediante convolución Voigt."""

    y = np.zeros_like(x)
    for x0, A in zip(x0_list, A_list):
        y += voigt(x, x0, sigma, gamma0, gamma_slope, A, x_ref)
    return y


# ====================================================
# === ARCHIVOS DE BARRAS ==============================
# ====================================================

archivos_barras = [
    # r"C:\Users\cauca\Documents\tesis\Datos\A2.txt",
    # r"C:\Users\cauca\Documents\tesis\Datos\B2.txt",
    r"C:\Users\cauca\Documents\tesis\Datos\C3.txt",
    # r"C:\Users\cauca\Documents\tesis\Datos\D2.txt",
]

# ====================================================
# === PARÁMETROS FÍSICOS DEL ENSANCHAMIENTO ===========
# ====================================================

FWHM_g = 3.0             # resolución instrumental
sigma = FWHM_g / 2.355   # ancho Gaussiano
gamma0 = 2.16            # ancho Lorentziano (efectos de vida/media)
gamma_slope = 0.000      # variación con energía

# Corrimiento global experimental, para coincidir con el espectro experimental [eV]
delta_exp = -2.13


# ====================================================
# === DEFINICIÓN DEL EJE DE ENERGÍA ===================
# ====================================================

# Se construye un eje común que cubra todas las barras
x_min, x_max = np.inf, -np.inf
for archivo in archivos_barras:
    data = np.loadtxt(archivo, skiprows=1)
    x_min = min(x_min, data[:, 0].min())
    x_max = max(x_max, data[:, 0].max())

# Margen amplio para evitar truncamientos del Voigt
x = np.linspace(x_min - 16, x_max + 16, 4000)


# ====================================================
# === CONSTRUCCIÓN DE LOS ESPECTROS ==================
# ====================================================

curvas = []
x0_lists, A_lists = [], []

# Colores solo para visualización
colores_barras = ["green", "deeppink"]
colores_curvas = ["darkgreen", "crimson"]

for i, archivo in enumerate(archivos_barras):

    data = np.loadtxt(archivo, skiprows=1)
    x0_list = data[:, 0]
    A_list = data[:, 1]

    x_ref = x0_list[-1]      # energía de referencia para la variación de gamma

    # Convolución de barras a espectro continuo
    y = barras_a_voigt(x, x0_list, A_list, sigma, gamma0, gamma_slope, x_ref)

    # Normalización al máximo
    y_max = np.nanmax(y)
    if y_max > 0:
        y *= A_list.max() / y_max

    curvas.append(y)
    x0_lists.append(x0_list)
    A_lists.append(A_list)


y_total = curvas[0]
y_total /= y_total.max()   


# ====================================================
# === CORRIMIENTO EXPERIMENTAL FINAL =================
# ====================================================

x_desplazada = x + delta_exp

# ====================================================
# === GUARDADO DEL ESPECTRO RESULTANTE ================
# ====================================================

# Se conserva la estructura esperada por otros scripts: columna 0 = energía, columna 2 = intensidad normalizada
ceros = np.zeros_like(y_total)
datos_suma = np.column_stack((x_desplazada, ceros, y_total))    # combina columnas: energía, ceros, intensidad

df_suma = pd.DataFrame(
    datos_suma,
    columns=["Energia_eV", "Col_dummy", "Intensidad_norm"]
)

nombres = [a.split("\\")[-1].replace(".txt", "") for a in archivos_barras]
nombre_salida = "sumaponderada_" + "_".join(nombres) + ".txt"
ruta_salida = rf"C:\Users\cauca\Documents\tesis\Datos\{nombre_salida}"

df_suma.to_csv(ruta_salida, sep="\t", index=False, float_format="%.8f")
print(f"Suma ponderada guardada en: {ruta_salida}")


# ====================================================
# === VISUALIZACIÓN ==================================
# ====================================================

fig, ax = plt.subplots(figsize=(8, 5))

# Barras + Voigts individuales
for x0_list, A_list, y, c_b, c_v, archivo in zip(
    x0_lists, A_lists, curvas, colores_barras, colores_curvas, archivos_barras
):
    label = archivo.split("\\")[-1].replace(".txt", "")

    A_max = np.nanmax(A_list)        
    if A_max > 0:
        A_norm = A_list / A_max    # intensidad = 0 -> no se normaliza
    else:
        A_norm = A_list

    markerline, stemlines, _ = ax.stem(x0_list, A_norm, basefmt=" ")
    plt.setp(markerline, color=c_b, marker="o", ms=3)
    plt.setp(stemlines, color=c_b, lw=0.8)

    ax.plot(x, y / y.max(), lw=1.5, color=c_v, label=f"Voigt {label}")

# Formato del gráfico
ax.set_xlabel("Energía (eV)")
ax.set_ylabel("Intensidad (norm.)")
ax.set_title("Convolución Voigt de espectros de barras")
ax.legend()
ax.grid(True)

plt.tight_layout()
plt.show()
