import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from scipy.special import wofz

# ====================================================
# === CONFIGURACIÓN DE ARCHIVOS EXPERIMENTALES =======
# ====================================================
file_paths_exp1 = [
    # r"C:\Users\cauca\Documents\tesis\Datos\Co_met4.dat",
    # r"C:\Users\cauca\Documents\tesis\Datos\Co3O4_2.dat",
]

file_paths_exp2 = [
    r"C:\Users\cauca\Documents\tesis\Datos\CoSO4_suma.dat",
]

file_paths_exp3 = [
    # r"C:\Users\cauca\Documents\tesis\Datos\TmCoO3_suma.dat",
    # r"C:\Users\cauca\Documents\tesis\Datos\cianuro_suma.dat",
    # r"C:\Users\cauca\Documents\tesis\Datos\LaCoO3_suma.dat",
#     r"C:\Users\cauca\Documents\tesis\Datos\TmFe04_suma.dat",
]

# Parámetros de calibración
d_si = 0.8165779352  # Espaciado interplanar Ge(444)
motor1_0 = 0         # M1 de la medición de referencia (para cada tanda de mediciones)
theta_0 = 82.96      # Ángulo de Bragg

# ====================================================
# === FUNCIONES DE PERFIL DE VOIGT ===================
# ====================================================
def voigt(x, x0, sigma, gamma0, gamma_slope, A, x_ref=None):
    """Perfil de Voigt con gamma variable."""
    if x_ref is None:
        x_ref = x0
    gamma_x = gamma0 + gamma_slope * (x - x_ref)
    z = ((x - x0) + 1j * gamma_x) / (sigma * np.sqrt(2))
    return A * np.real(wofz(z)) / (sigma * np.sqrt(2*np.pi))

def barras_a_voigt(x, x0_list, A_list, sigma, gamma0, gamma_slope, x_ref=None):
    """Suma de perfiles Voigt a partir de barras discretas."""
    y = np.zeros_like(x)
    for x0, A in zip(x0_list, A_list):
        y += voigt(x, x0, sigma, gamma0, gamma_slope, A, x_ref)
    return y

# ====================================================
# === CONFIGURACIÓN DE ARCHIVOS DE BARRAS ============
# ====================================================
archivos_barras = [
    # r"C:\Users\cauca\Documents\tesis\Datos\A23.txt",
    # r"C:\Users\cauca\Documents\tesis\Datos\A1.txt",
    r"C:\Users\cauca\Documents\tesis\Datos\A2.txt",
    #r"C:\Users\cauca\Documents\tesis\Datos\A37.txt",
    #r"C:\Users\cauca\Documents\tesis\Datos\A63.txt",
    #r"C:\Users\cauca\Documents\tesis\Datos\experimento2.txt",

    # r"C:\Users\cauca\Documents\tesis\Datos\B1.txt",
    # r"C:\Users\cauca\Documents\tesis\Datos\B2.txt",
    # r"C:\Users\cauca\Documents\tesis\Datos\B4.txt",
    #r"C:\Users\cauca\Documents\tesis\Datos\C6.txt",
    # r"C:\Users\cauca\Documents\tesis\Datos\C14.txt",
    #r"C:\Users\cauca\Documents\tesis\Datos\D3.txt",
    # r"C:\Users\cauca\Documents\tesis\Datos\D10.txt",

    # r"C:\Users\cauca\Documents\tesis\Datos\C3.txt",

    # r"C:\Users\cauca\Documents\tesis\Datos\Delta4.txt"

    # r"C:\Users\cauca\Documents\tesis\Datos\D2.txt",
]

# Parámetros de ensanchamiento
FWHM_g = 3
sigma = FWHM_g / 2.355
gamma0 = 2.16
gamma_slope = -0.002       # Comportamiento lineal del ensanchamiento naatura

# Desplazamiento manual de la segunda curva teorica, para Co3O4 (en eV)
delta_E = -10.5  # <=== ajustar a mano, el programa genera curvas en distintos máximos

# ====================================================
# === PARTE 1: DATOS EXPERIMENTALES ==================
# ====================================================
fig, (ax_exp, ax_teo) = plt.subplots(2, 1, figsize=(8, 8), sharex=True)

for file_path in file_paths_exp1:

    motor1_01 = -4.384

    data = pd.read_csv(file_path, skiprows=2, sep=r'\s+')
    motor1 = data.iloc[:, 0]
    cuentas = data.iloc[:, 2]

    # Conversión a energía (Bragg)
    theta_bragg = theta_0 - (motor1 - motor1_01) / 200.0  
    lambda_m = 2 * d_si * np.sin(np.radians(theta_bragg))
    E = 12398.42 / lambda_m
    cuentas_norm = cuentas / cuentas.max()

    # === Cálculo del error experimental ===
    N_i = cuentas.astype(float)
    N_M = N_i.max()

    DeltaI = cuentas_norm * np.sqrt((1/N_i) + (1/N_M))
    DeltaE = 0.2

    # Reemplazar infinitos (cuando N_i = 0) por 0
    DeltaI = np.nan_to_num(DeltaI, nan=0.0, posinf=0.0, neginf=0.0)

    # Guardar archivo calibrado
    df_calib = pd.DataFrame({
        'Energía_eV': E,
        'Cuentas': cuentas,
        'Cuentas_norm': cuentas_norm,
        'Error energía': DeltaE,
        'Error intensidad': DeltaI
    })
    output_calib = file_path.replace('.dat', '_calib.dat')
    df_calib.to_csv(output_calib, sep='\t', index=False, float_format='%.8f')
    print(f"Archivo calibrado guardado: {output_calib}")

    # Graficar experimental
    label = file_path.split('\\')[-1].replace('.dat', '')
    # Graficar experimental con barras de error
    label = file_path.split('\\')[-1].replace('.dat', '')
    # ax_exp.errorbar(E, cuentas_norm, yerr=DeltaI, #xerr=DeltaE, 
    #                 fmt='o', ms=2.5, elinewidth=0.8, capsize=2, 
    #                 alpha=0.8, label=label)
    #ax_exp.plot(E, cuentas_norm, marker='o', ms=3, label=label, color='green')

    ax_exp.plot(E, cuentas_norm, lw=1.8, label="CoSO$_4$")
    # ax_exp.errorbar(E, cuentas_norm, yerr=DeltaI,
    #                 fmt='none',
    #                 elinewidth=0.7, capsize=1.5,
    #                 alpha=0.7)

for file_path in file_paths_exp2:

    motor1_02 = -5.521

    data = pd.read_csv(file_path, skiprows=2, sep=r'\s+')
    motor1 = data.iloc[:, 0]
    cuentas = data.iloc[:, 2]

    # Conversión a energía (Bragg)
    theta_bragg = theta_0 - (motor1 - motor1_02) / 200.0  
    lambda_m = 2 * d_si * np.sin(np.radians(theta_bragg))
    E = 12398.42 / lambda_m
    cuentas_norm = cuentas / cuentas.max()

    # === Cálculo del error experimental ===
    N_i = cuentas.astype(float)
    N_M = N_i.max()

    DeltaI = cuentas_norm * np.sqrt((1/N_i) + (1/N_M))
    DeltaE = 0.2

    # Reemplazar infinitos (cuando N_i = 0) por 0
    DeltaI = np.nan_to_num(DeltaI, nan=0.0, posinf=0.0, neginf=0.0)

    # Guardar archivo calibrado
    df_calib = pd.DataFrame({
        'Energía_eV': E,
        'Cuentas': cuentas,
        'Cuentas_norm': cuentas_norm,
        'Error energía': DeltaE,
        'Error intensidad': DeltaI
    })
    output_calib = file_path.replace('.dat', '_calib.dat')
    df_calib.to_csv(output_calib, sep='\t', index=False, float_format='%.8f')
    print(f"Archivo calibrado guardado: {output_calib}")

    # Graficar experimental
    label = file_path.split('\\')[-1].replace('.dat', '')
    # Graficar experimental con barras de error
    label = file_path.split('\\')[-1].replace('.dat', '')
    # ax_exp.errorbar(E, cuentas_norm, yerr=DeltaI, #xerr=DeltaE, 
    #                 fmt='o', ms=2.5, elinewidth=0.8, capsize=2, 
    #                 alpha=0.8, label=label)
    #ax_exp.plot(E, cuentas_norm, marker='o', ms=3, label=label, color='green')

    ax_exp.plot(E, cuentas_norm, lw=1.8, label="CoSO$_4$")
    # ax_exp.errorbar(E, cuentas_norm, yerr=DeltaI,
    #                 fmt='none',
    #                 elinewidth=0.7, capsize=1.5,
    #                 alpha=0.7)

for file_path in file_paths_exp3:

    motor1_03 = -6.941

    data = pd.read_csv(file_path, skiprows=2, sep=r'\s+')
    motor1 = data.iloc[:, 0]
    cuentas = data.iloc[:, 2]

    # Conversión a energía (Bragg)
    theta_bragg = theta_0 - (motor1 - motor1_03) / 200.0  
    lambda_m = 2 * d_si * np.sin(np.radians(theta_bragg))
    E = 12398.42 / lambda_m
    cuentas_norm = cuentas / cuentas.max()

    # === Cálculo del error experimental ===
    N_i = cuentas.astype(float)
    N_M = N_i.max()

    DeltaI = cuentas_norm * np.sqrt((1/N_i) + (1/N_M))
    DeltaE = 0.2

    # Reemplazar infinitos (cuando N_i = 0) por 0
    DeltaI = np.nan_to_num(DeltaI, nan=0.0, posinf=0.0, neginf=0.0)

    # Guardar archivo calibrado
    df_calib = pd.DataFrame({
        'Energía_eV': E,
        'Cuentas': cuentas,
        'Cuentas_norm': cuentas_norm,
        'Error energía': DeltaE,
        'Error intensidad': DeltaI
    })
    output_calib = file_path.replace('.dat', '_calib.dat')
    df_calib.to_csv(output_calib, sep='\t', index=False, float_format='%.8f')
    print(f"Archivo calibrado guardado: {output_calib}")

    # Graficar experimental
    label = file_path.split('\\')[-1].replace('.dat', '')
    # Graficar experimental con barras de error
    label = file_path.split('\\')[-1].replace('.dat', '')
    # ax_exp.errorbar(E, cuentas_norm, yerr=DeltaI, #xerr=DeltaE, 
    #                 fmt='o', ms=2.5, elinewidth=0.8, capsize=2, 
    #                 alpha=0.8, label=label)
    #ax_exp.plot(E, cuentas_norm, marker='o', ms=3, label=label, color='green')

    ax_exp.plot(E, cuentas_norm, lw=1.8, label="K$_3$[Co(CN)$_6$]")
    # ax_exp.errorbar(E, cuentas_norm, yerr=DeltaI,
    #                 fmt='none',
    #                 elinewidth=0.7, capsize=1.5,
    #                 alpha=0.7)

# Máximo experimental
E_max_exp = E[np.argmax(cuentas_norm)]
print(f"Máximo experimental: {E_max_exp:.2f} eV")

ax_exp.set_ylabel("Cuentas normalizadas")
# ax_exp.set_title("Datos experimentales calibrados")
ax_exp.legend()
ax_exp.grid(True)

# ====================================================
# === PARTE 2: CURVAS VOIGT Y SUMA PONDERADA =========
# ====================================================

# Determinar eje X global
x_min_global, x_max_global = np.inf, -np.inf
for archivo in archivos_barras:
    data = np.loadtxt(archivo, skiprows=1)
    x_min_global = min(x_min_global, data[:, 0].min())
    x_max_global = max(x_max_global, data[:, 0].max())
x = np.linspace(x_min_global - 12, x_max_global + 12, 4000)

barras_originales = []

# Generar curvas individuales
curvas = []
for i, archivo in enumerate(archivos_barras):
    data = np.loadtxt(archivo, skiprows=1)
    x0_list, A_list = data[:, 0], data[:, 1]

    barras_originales.append((x0_list.copy(), A_list.copy()))

    # Aplicar desplazamiento manual solo a la segunda curva
    if i == 1:
        x0_list = x0_list + delta_E

    x_ref = x0_list[-1]
    y = barras_a_voigt(x, x0_list, A_list, sigma, gamma0, gamma_slope, x_ref)
    y *= A_list.max() / y.max()
    curvas.append(y)

# Suma ponderada de las Voigts para Co3O4
if len(curvas) == 2:
    frac = [0.6666666, 0.3333333]
    areas = [np.trapezoid(y, x) for y in curvas]
    factores = [f / a for f, a in zip(frac, areas)]
    y_total = sum(f * y for f, y in zip(factores, curvas))
    y_total /= y_total.max()
elif len(curvas) == 1:
    y_total = curvas[0] / curvas[0].max()
else:
    raise ValueError("Configurado para 1 o 2 curvas de barras.")

# Máximo teórico
E_max_teo = x[np.argmax(y_total)]
delta_exp = E_max_exp - E_max_teo
print(f"Máximo teórico: {E_max_teo:.2f} eV")
print(f"Corrimiento experimental aplicado: {delta_exp:.2f} eV")

# Aplicar desplazamiento experimental a la curva total
alpha = 0.75
#x_shifted = (E_max_teo + alpha * (x - E_max_teo)) + delta_exp
x_shifted = x + delta_exp

for x0_list, A_list in barras_originales:
    ax_teo.vlines(x0_list + delta_exp, 0,
                  A_list / A_list.max(),
                  colors='purple', alpha=0.6, linewidth=1.2)


# ====================================================
# === GUARDAR SUMA PONDERADA =========================
# ====================================================
ceros = np.zeros_like(y_total)
datos_suma = np.column_stack((x_shifted, ceros, y_total))
df_suma = pd.DataFrame(datos_suma, columns=["Energia_eV", "Ceros", "Intensidad_norm"])

nombres_barras = [archivo.split("\\")[-1].replace(".txt", "") for archivo in archivos_barras]
nombre_salida = "sumaponderada_" + "_".join(nombres_barras) + ".txt"
ruta_salida = rf"C:\Users\cauca\Documents\tesis\Datos\{nombre_salida}"

df_suma.to_csv(ruta_salida, sep="\t", index=False, float_format="%.8f")
print(f"Suma ponderada guardada en: {ruta_salida}")

# ====================================================
# === GRÁFICO DE LA SUMA PONDERADA ===================
# ====================================================
ax_teo.plot(x_shifted, y_total, lw=2.5, color="indigo", label="Co$^{3+}$, HS")
ax_teo.set_xlabel("Energía (eV)")
ax_teo.set_ylabel("Intensidad [u.a]")
# ax_teo.set_title(f"Suma ponderada (ΔE = {delta_E} eV aplicado a la 2da curva)")
ax_teo.legend()
ax_teo.grid(True)

# Ajuste de ejes
xmin = min(ax_exp.get_xlim()[0], ax_teo.get_xlim()[0])
xmax = max(ax_exp.get_xlim()[1], ax_teo.get_xlim()[1])
ax_exp.set_xlim(xmin, xmax)
ax_teo.set_xlim(xmin, xmax)


plt.tight_layout()
plt.show()