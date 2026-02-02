# Suma de espectros experimentales

import pandas as pd
import glob
import matplotlib.pyplot as plt
import os
import re

# Lista para almacenar dataframes
dataframes = []

# Ruta de los archivos (ajusta según sea necesario)
#file_path_pattern = r"C:\Users\cauca\Documents\tesis\Datos\TmFe06_*.dat"

# Leer todos los archivos y almacenarlos en la lista
for file in glob.glob(file_path_pattern):
    # Leer el archivo, omitiendo la primera línea con la fecha
    df = pd.read_csv(file, skiprows=1, sep=' ', engine='python', names=['motor1', 'motor2', 'cuentas', 'tiempo', 'c/s'])
    # Convertir las columnas 'motor1', 'motor2', 'cuentas' y 'tiempo' a numéricas
    df['motor1'] = pd.to_numeric(df['motor1'], errors='coerce')
    df['motor2'] = pd.to_numeric(df['motor2'], errors='coerce')
    df['cuentas'] = pd.to_numeric(df['cuentas'], errors='coerce')
    df['tiempo'] = pd.to_numeric(df['tiempo'], errors='coerce')
    # Filtrar filas con valores no numéricos en 'motor1', 'motor2', 'cuentas' y 'tiempo'
    df = df.dropna(subset=['motor1', 'motor2', 'cuentas', 'tiempo'])
    df = df[df['tiempo'] != 0]  # Asegurar que no hay división por cero
    # Agregar una columna para identificar cada archivo
    df['archivo'] = file
    # Agregar el dataframe a la lista
    dataframes.append(df)

# Combinar todos los dataframes en uno solo
combined_df = pd.concat(dataframes, ignore_index=True)

# Agrupar por motor1 y motor2, y sumar las columnas 'cuentas' y 'tiempo'
grouped_df = combined_df.groupby(['motor1', 'motor2']).agg({'cuentas': 'sum', 'tiempo': 'sum'}).reset_index()

# Manejar la división por cero
grouped_df['c/s'] = grouped_df.apply(lambda row: row['cuentas'] / row['tiempo'] if row['tiempo'] != 0 else 0, axis=1)

# Ordenar por 'motor1' y 'motor2' después de agrupar
grouped_df = grouped_df.sort_values(by=['motor1', 'motor2']).reset_index(drop=True)

# Crear nombre del archivo de salida en la misma carpeta que los archivos originales

# Tomamos el primer archivo encontrado
first_file = glob.glob(file_path_pattern)[0]

# Ejemplo: C:\...\TmCoO3_001.dat → C:\...\TmCoO3.dat
base_file = first_file.split('_')[0] + '.dat'

# Ejemplo: C:\...\TmCoO3.dat → C:\...\TmCoO3_suma.dat
output_file = base_file.replace('.dat', '_suma.dat')

# Guardar el archivo combinado
grouped_df.to_csv(output_file, index=False, sep=' ')

print(f"Archivo guardado: {output_file}")

# Generar el gráfico de 'cuentas' en función de 'motor1'
plt.figure(figsize=(10, 6))
plt.plot(grouped_df['motor1'], grouped_df['cuentas'], marker='o')
plt.xlabel('motor1')
plt.ylabel('cuentas')
plt.title('Cuentas en función de motor1')
plt.grid(True)
plt.savefig('cuentas_vs_motor1_A.png')
plt.show()

# Generar el gráfico de 'motor1' vs 'cuentas' para todos los archivos juntos
plt.figure(figsize=(10, 6))
for file, group in combined_df.groupby('archivo'):
    plt.plot(group['motor1'], group['cuentas'], marker='o', label=file)
plt.xlabel('motor1')
plt.ylabel('cuentas')
plt.title('Cuentas en función de motor1 para todos los archivos')
plt.legend()
plt.grid(True)
plt.savefig('cuentas_vs_motor1_todos_archivos_A.png')
plt.show()