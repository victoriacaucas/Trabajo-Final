# Construcción del espectro HS y descomposición lineal de TmCoO3.

import numpy as np

# ====================================================
# === PARÁMETROS FÍSICOS DEL MODELO ==================
# ====================================================

# Coeficientes del modelo lineal:
# LaCoO3(E) = alpha*Cianuro(E) + gamma*HS(E)
alpha = 0.34
gamma = 0.74

# Incertidumbres asociadas a los coeficientes
sigma_alpha = 0.04
sigma_gamma = 0.04

# Corrimiento energético rígido aplicado al espectro de cianuro (alineación con el eje energético de LaCoO3) [eV]
delta_E = 0.42  


# ====================================================
# === RUTAS A LOS ARCHIVOS DE ENTRADA =================
# ====================================================

ruta_cian = r"C:\Users\cauca\Documents\tesis\Datos\cianuro_suma_calib.dat"
ruta_laco = r"C:\Users\cauca\Documents\tesis\Datos\LaCoO3_suma_calib.dat"
ruta_tm   = r"C:\Users\cauca\Documents\tesis\Datos\TmCoO3_suma_calib.dat"
ruta_hs   = r"C:\Users\cauca\Documents\tesis\Datos\HS.dat"


# ====================================================
# === CARGA DE DATOS EXPERIMENTALES ==================
# ====================================================
# Se cargan:
# - Energía (eV)
# - Intensidad normalizada
# - Error experimental de la intensidad

x_cian, y_cian, sigma_cian = np.loadtxt(
    ruta_cian, unpack=True, skiprows=1, usecols=(0, 2, 4)
)

x_laco, y_laco, sigma_laco = np.loadtxt(
    ruta_laco, unpack=True, skiprows=1, usecols=(0, 2, 4)
)

x_tm, y_tm, sigma_tm = np.loadtxt(
    ruta_tm, unpack=True, skiprows=1, usecols=(0, 2, 4)
)


# ====================================================
# === PREPROCESAMIENTO ENERGÉTICO ====================
# ====================================================

# Corrimiento rígido del eje energético del cianuro
x_cian = x_cian + delta_E

# Interpolación del espectro de cianuro (y su error) sobre la grilla energética del espectro de LaCoO3 (permite una combinación punto a punto de los espectros)
y_cian_i = np.interp(x_laco, x_cian, y_cian)           # Construye una función lineal C(x) a partir de los datos del cianuro y la evalúa en cada punto de x_laco
sigma_cian_i = np.interp(x_laco, x_cian, sigma_cian)

# En este punto:
# x_laco es la grilla común
# y_laco, y_cian_i están definidos punto a punto


# ====================================================
# === CONSTRUCCIÓN DEL ESPECTRO HS ===================
# ====================================================
# A partir del modelo LaCoO3 = alpha*Cianuro + gamma*HS se despeja:
# HS = (LaCoO3 − alpha*Cianuro) / gamma

y_hs = (y_laco - alpha * y_cian_i) / gamma


# ====================================================
# === PROPAGACIÓN DE ERRORES =========================
# ====================================================
# Se asume independencia entre variables y propagación gaussiana

sigma_hs = np.sqrt(
    (sigma_laco / gamma)**2 +
    (alpha * sigma_cian_i / gamma)**2 +
    (y_cian_i * sigma_alpha / gamma)**2 +
    (y_hs * sigma_gamma / gamma)**2
)


# ====================================================
# === GUARDADO DEL ARCHIVO HS ========================
# ====================================================
# Se conserva el formato de 5 columnas para compatibilidad con el resto del análisis

hs_out = np.zeros((len(x_laco), 5))
hs_out[:, 0] = x_laco               # Energía (eV)
hs_out[:, 1] = y_hs                 # Intensidad HS
hs_out[:, 2] = y_hs/y_hs.max()      # Intensidad HS normalizada
hs_out[:, 4] = sigma_hs/y_hs.max()  # Error de intensidad normalizado

np.savetxt(
    ruta_hs,
    hs_out,
    fmt="%.6e"
)


# ====================================================
# === PROYECCIÓN SOBRE LA GRILLA DE TmCoO3 ===========
# ====================================================
# En este caso se verifica que las grillas coincidan

if not np.allclose(x_tm, x_laco):
    raise ValueError(
        "Las grillas energéticas de TmCoO3 y LaCoO3 no coinciden."
    )

# Como coinciden, no es necesaria interpolación adicional
y_cian_tm = y_cian_i
sigma_cian_tm = sigma_cian_i

y_hs_tm = y_hs
sigma_hs_tm = sigma_hs


# ====================================================
# === AJUSTE LINEAL: TmCoO3 = A*Cianuro + B*HS =======
# ====================================================
# con la restricción física A + B = 1

# Se reescribe el problema como Y = A*X donde:
X = y_cian_tm - y_hs_tm
Y = y_tm - y_hs_tm

# Error efectivo en Y
sigma_Y = np.sqrt(sigma_tm**2 + sigma_hs_tm**2)

# Pesos para mínimos cuadrados ponderados
w = 1.0 / sigma_Y**2

# Estimador de A (mínimos cuadrados ponderados)
A = np.sum(w * X * Y) / np.sum(w * X**2)
B = 1.0 - A

# Incertidumbres
sigma_A = np.sqrt(1.0 / np.sum(w * X**2))
sigma_B = sigma_A


# ====================================================
# === RESULTADOS ====================================
# ====================================================

print("Resultados del ajuste:")
print("TmCoO3(E) = A*Cianuro(E) + B*HS(E)")
print(f"A = {A:.4f} ± {sigma_A:.4f}")
print(f"B = {B:.4f} ± {sigma_B:.4f}")
