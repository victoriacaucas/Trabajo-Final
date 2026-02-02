# Comparación visual de espectros experimentales y teóricos, y determinación de máximos energéticos.

import numpy as np
import pandas as pd
import matplotlib.pyplot as plt

# ====================================================
# === CONFIGURACIÓN DE ARCHIVOS ======================
# ====================================================

file_paths = [
    # Co2+ HS
    # r"C:\Users\cauca\Documents\tesis\Datos\CoSO4_suma_calib.dat",
    # r"C:\Users\cauca\Documents\tesis\Datos\sumaponderada_A2.txt",

    # # Co2+ LS
    # r"C:\Users\cauca\Documents\tesis\Datos\sumaponderada_D2.txt",

    # # Co3+ LS
    # r"C:\Users\cauca\Documents\tesis\Datos\cianuro_suma_calib.dat",
    # r"C:\Users\cauca\Documents\tesis\Datos\sumaponderada_B2.txt", 

    # # Co3+ HS
    r"C:\Users\cauca\Documents\tesis\Datos\HS.dat",
    r"C:\Users\cauca\Documents\tesis\Datos\sumaponderada_C3.txt",
]

# ====================================================
# === CREACIÓN DE LA FIGURA ==========================
# ====================================================

fig, ax = plt.subplots(figsize=(8, 6))

for file_path in file_paths:

    data = pd.read_csv(file_path, skiprows=2, sep=r'\s+')

    E = data.iloc[:, 0].to_numpy()
    cuentas = data.iloc[:, 2].to_numpy()

    # Normalización (solo para comparación visual)
    cuentas_norm = cuentas / np.nanmax(cuentas)

    # Identificación del máximo
    idx_max = np.nanargmax(cuentas_norm)
    E_max = E[idx_max]

    # Salida por pantalla
    nombre = file_path.split('\\')[-1]
    print(f"Máximo de {nombre}:")
    print(f"  Energía = {E_max:.2f} eV")

    ax.plot(E, cuentas_norm, linewidth=1.5, label=nombre)


# ====================================================
# === FORMATO FINAL DE LA FIGURA ====================
# ====================================================

ax.set_xlabel('Energía (eV)')
ax.set_ylabel('Intensidad normalizada')
ax.set_title('Comparación de espectros')
ax.legend()
ax.grid(True)

plt.tight_layout()
plt.show()
