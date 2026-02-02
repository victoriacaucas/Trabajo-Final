# Calibración energética, normalización de espectros y estimación de errores experimentales.

import numpy as np
import pandas as pd
import matplotlib.pyplot as plt

# ====================================================
# === CONFIGURACIÓN GENERAL ==========================
# ====================================================

# Parámetros de calibración
d_si = 0.8165779352      # Espaciado interplanar del cristal analizador (Ge(444)) [Å]
theta_0 = 82.96          # Ángulo de referencia del monocromador [°]
motor_a_ang = 1 / 200.0  # Conversión motor a ángulo
hc = 12398.42            # Constante hc [eV*Å]
DeltaE_const = 0.2       # Resolución instrumental en energía [eV]


# ====================================================
# === DEFINICIÓN DE TANDAS EXPERIMENTALES ============
# ====================================================
# Cada tanda incluye:
# - lista de archivos
# - valor de motor1_0: posición del máximo del Co metálico usado como referencia energética de esa tanda

tandas = [
    {
        "file_paths": [
            # r"C:\Users\cauca\Documents\tesis\Datos\Co_met4.dat",
            # r"C:\Users\cauca\Documents\tesis\Datos\Co3O4_2.dat",
        ],
        "motor1_0": -4.384,
    },
    {
        "file_paths": [
            # r"C:\Users\cauca\Documents\tesis\Datos\CoSO4_suma.dat",
        ],
        "motor1_0": -5.521,
    },
    {
        "file_paths": [
            r"C:\Users\cauca\Documents\tesis\Datos\TmCoO3_suma.dat",
            r"C:\Users\cauca\Documents\tesis\Datos\cianuro_suma.dat",
            r"C:\Users\cauca\Documents\tesis\Datos\LaCoO3_suma.dat",
            # r"C:\Users\cauca\Documents\tesis\Datos\TmFe04_suma.dat",
        ],
        "motor1_0": -6.941,
    },
    {
        "file_paths": [
            # r"C:\Users\cauca\Documents\tesis\Datos\TmFe06_suma.dat",
        ],
        "motor1_0": -1.452,
    },
]


# ====================================================
# === FUNCIÓN DE PROCESAMIENTO =======================
# ====================================================

def procesar_archivo(file_path, motor1_0, ax):
    """
    Procesa un archivo experimental:
    - Convierte motor a energía mediante ley de Bragg
    - Normaliza intensidades
    - Calcula errores experimentales (estadística de conteo)
    - Guarda archivo calibrado
    - Grafica el espectro con barras de error
    """

    # Lectura de datos
    data = pd.read_csv(file_path, skiprows=2, sep=r"\s+")
    motor1 = data.iloc[:, 0].to_numpy()
    cuentas = data.iloc[:, 2].to_numpy().astype(float)

    # Conversión a energía (Bragg)
    theta_bragg = theta_0 - (motor1 - motor1_0) * motor_a_ang
    lambda_m = 2 * d_si * np.sin(np.radians(theta_bragg))
    E = hc / lambda_m

    # Normalización de intensidad
    cuentas_max = cuentas.max()
    cuentas_norm = cuentas / cuentas_max

    # Cálculo de errores experimentales
    with np.errstate(divide="ignore", invalid="ignore"):
        DeltaI = cuentas_norm * np.sqrt((1 / cuentas) + (1 / cuentas_max))
    DeltaI[~np.isfinite(DeltaI)] = np.nan  # Los puntos con cuentas = 0 no tienen error bien definido, se asigna NaN para poder ignorarlos en análisis posteriores
    DeltaE = np.full_like(E, DeltaE_const)

    # Guardado del archivo calibrado
    df_calib = pd.DataFrame({
        "Energia_eV": E,
        "Cuentas": cuentas,
        "Cuentas_norm": cuentas_norm,
        "Error_Energia_eV": DeltaE,
        "Error_Intensidad": DeltaI,
    })

    output_calib = file_path.replace(".dat", "_calib.dat")
    df_calib.to_csv(output_calib, sep="\t", index=False, float_format="%.8f")
    print(f"Archivo calibrado guardado: {output_calib}")

    # Gráfica
    label = file_path.split("\\")[-1].replace(".dat", "")

    err = ax.errorbar(
        E, cuentas_norm,
        yerr=DeltaI,
        fmt="o", ms=2.5,
        elinewidth=0.8, capsize=2,
        alpha=0.8,
        label=label,
    )

    # Línea que une los puntos
    color = err[0].get_color()
    ax.plot(E, cuentas_norm, color=color, linewidth=1.2, alpha=0.9)

    # Energía del máximo experimental
    E_max = E[np.nanargmax(cuentas_norm)]
    print(f"Máximo experimental ({label}): {E_max:.2f} eV")


# ====================================================
# === PROCESAMIENTO DE DATOS EXPERIMENTALES ==========
# ====================================================

fig, ax_exp = plt.subplots(figsize=(8, 6))

for tanda in tandas:
    for file_path in tanda["file_paths"]:
        procesar_archivo(file_path, tanda["motor1_0"], ax_exp)


# ====================================================
# === DATOS HS (YA CALIBRADOS) =======================
# ====================================================

ruta_hs = r"C:\Users\cauca\Documents\tesis\Datos\HS.dat"
data_hs = pd.read_csv(ruta_hs, sep=r"\s+", header=None)

E_hs = data_hs.iloc[:, 0]
HS = data_hs.iloc[:, 1]
DeltaI_hs = data_hs.iloc[:, 4]

HS_norm = HS / HS.max()
DeltaI_hs_norm = DeltaI_hs / HS.max()

ax_exp.errorbar(
    E_hs, HS_norm,
    yerr=DeltaI_hs_norm,
    fmt="o", ms=2.5,
    elinewidth=0.8, capsize=2,
    alpha=0.8,
    label="HS",
)

ax_exp.plot(E_hs, HS_norm, linewidth=1.2, alpha=0.9)

print(f"Máximo HS: {E_hs.iloc[np.argmax(HS_norm)]:.2f} eV")


# ====================================================
# === FORMATO FINAL DE LA FIGURA ====================
# ====================================================

ax_exp.set_xlabel("Energía (eV)")
ax_exp.set_ylabel("Cuentas normalizadas")
ax_exp.set_title("Datos experimentales calibrados")
ax_exp.legend()
ax_exp.grid(True)

plt.tight_layout()
plt.show()
