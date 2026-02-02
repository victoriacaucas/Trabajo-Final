# Cálculo del IAD y del primer momento energético con propagación de errores.

import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from scipy.signal import savgol_filter


# ============================================================
#   FUNCIONES DE PROCESAMIENTO
# ============================================================

def normalizar_area(E, I, dI):
    """Normaliza el espectro de manera que el área total sea 1. Esta normalización es necesaria 
    para que el IAD compare solo diferencias de forma espectral, no de intensidad absoluta."""

    area = np.trapezoid(I, E)
    I_norm = I / area
    dI_norm = dI / area
    return I_norm, dI_norm


def centrar_M1(E, I):
    """Centra el espectro en su primer momento (M1). Esto elimina corrimientos 
    energéticos globales y permite una comparación directa entre espectros."""

    M1 = np.trapezoid(E * I, E) / np.trapezoid(I, E)
    return E - M1


def preparar_espectro(E, I, dI):
    """Aplica los pasos estándar previos al cálculo del IAD:
    1) Normalización por área
    2) Centrado energético por el primer momento M1"""

    I_norm, dI_norm = normalizar_area(E, I, dI)
    E_centered = centrar_M1(E, I_norm)
    return E_centered, I_norm, dI_norm


# ============================================================
#   CÁLCULO DEL IAD CON PROPAGACIÓN DE ERRORES
# ============================================================

def calcular_IAD_y_error(E1, I1, dE1, dI1,
                         E2, I2, dE2, dI2):

    """Calcula el IAD entre dos espectros previamente normalizados y centrados. El error se propaga considerando 
    únicamente las incertidumbres de intensidad, asumiendo un error energético efectivo constante."""

    # Interpolación del segundo espectro al eje energético del primero
    I2_interp = np.interp(E1, E2, I2)
    dI2_interp = np.interp(E1, E2, dI2)

    # Diferencia punto a punto
    Delta = I1 - I2_interp

    # Definición del IAD
    IAD = np.trapezoid(np.abs(Delta), E1)

    # Error energético efectivo (resolución experimental típica)
    dE_eff = 0.28  # eV

    # Propagación de errores
    dIAD = np.sqrt(np.sum((dE_eff ** 2) * (dI1**2 + dI2_interp**2)))

    return IAD, dIAD


# ============================================================
#   CÁLCULO DEL PRIMER MOMENTO M1 Y SU ERROR
# ============================================================

def calcular_M1_y_error(E, I, dE, dI):
    """Calcula el primer momento energético M1 restringido a la región del espectro por encima del 50% del máximo de intensidad."""

    Imax = np.max(I)
    umbral = 0.5 * Imax

    mask = I >= umbral
    E_sel = E[mask]
    I_sel = I[mask]
    dE_sel = dE[mask]
    dI_sel = dI[mask]

    if len(E_sel) == 0:
        raise ValueError("No hay puntos por encima del 50% del máximo")

    # Definición discreta de M1
    S_I = np.sum(I_sel)
    M1 = np.sum(E_sel * I_sel) / S_I

    # Propagación de errores
    termE = np.sum((I_sel / S_I) ** 2 * (dE_sel ** 2))
    termI = np.sum(((E_sel - M1) / S_I) ** 2 * (dI_sel ** 2))

    dM1 = np.sqrt(termE + termI)

    return M1, dM1, E_sel.min(), E_sel.max()


# ============================================================
#   LECTURA Y PROCESAMIENTO DE LOS ARCHIVOS EXPERIMENTALES
# ============================================================

file_paths = [
    r"C:\Users\cauca\Documents\tesis\Datos\Co_met4_calib.dat",
    r"C:\Users\cauca\Documents\tesis\Datos\CoSO4_suma_calib.dat",
    r"C:\Users\cauca\Documents\tesis\Datos\Co3O4_2_calib.dat",
    r"C:\Users\cauca\Documents\tesis\Datos\TmCoO3_suma_calib.dat",
    r"C:\Users\cauca\Documents\tesis\Datos\cianuro_suma_calib.dat",
    r"C:\Users\cauca\Documents\tesis\Datos\LaCoO3_suma_calib.dat",
    r"C:\Users\cauca\Documents\tesis\Datos\TmFe04_suma_calib.dat",
    r"C:\Users\cauca\Documents\tesis\Datos\TmFe06_suma_calib.dat",
]

espectros = {}
errores = {}

print("\n=== PROCESANDO ARCHIVOS ===")

for fp in file_paths:

    data = pd.read_csv(fp, skiprows=2, sep=r'\s+')

    E = data.iloc[:, 0].values
    I = data.iloc[:, 2].values
    dE = data.iloc[:, 3].values
    dI = data.iloc[:, 4].values

    # Normalización preliminar (solo para cálculo de M1)
    I_norm = I / I.max()

    espectros[fp] = {"E": E, "I": I_norm}
    errores[fp] = {"dE": dE, "dI": dI}

    # Cálculo de M1
    try:
        M1, dM1, Emin, Emax = calcular_M1_y_error(E, I_norm, dE, dI)
        print(f"{fp}:  M1 = {M1:.3f} ± {dM1:.3f} eV   rango [{Emin:.1f}, {Emax:.1f}]")
    except Exception as e:
        print(f"Error calculando M1 en {fp}: {e}")


# ============================================================
#   DEFINICIÓN DEL ESPECTRO DE REFERENCIA
# ============================================================

archivo_ref = r"C:\Users\cauca\Documents\tesis\Datos\cianuro_suma_calib.dat"

E_ref_raw = espectros[archivo_ref]["E"]
I_ref_raw = espectros[archivo_ref]["I"]
dE_ref = errores[archivo_ref]["dE"]
dI_ref_raw = errores[archivo_ref]["dI"]

E_ref, I_ref, dI_ref = preparar_espectro(E_ref_raw, I_ref_raw, dI_ref_raw)


# ============================================================
#   CÁLCULO DEL IAD RESPECTO AL REFERENTE
# ============================================================

print("\n=== CÁLCULO DE IAD ===")

for fp in file_paths:

    E2_raw = espectros[fp]["E"]
    I2_raw = espectros[fp]["I"]
    dE2 = errores[fp]["dE"]
    dI2_raw = errores[fp]["dI"]

    E2, I2, dI2 = preparar_espectro(E2_raw, I2_raw, dI2_raw)

    iad, diad = calcular_IAD_y_error(
        E_ref, I_ref, dE_ref, dI_ref,
        E2, I2, dE2, dI2
    )

    print(f"IAD({archivo_ref.split('\\')[-1]} vs {fp.split('\\')[-1]}) = {iad:.6f} ± {diad:.6f}")