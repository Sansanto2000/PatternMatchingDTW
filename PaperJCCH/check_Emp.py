import os
import matplotlib.pyplot as plt
from utils import extract_lamp_info

# Especificar archivo teorico a emplear
act_dir = os.path.dirname(os.path.abspath(__file__))
file = os.path.join(act_dir, 'Lamparas_Meilan', 'w1466Blamp2.fits')

# Separar informacion
teo_x, teo_y, teo_head = extract_lamp_info(file, normalize=True)

# Determinar teorico real
try:
    teo_real_x = teo_x * teo_head['CD1_1'] + teo_head['CRVAL1']
except Exception as e:
    print(f"Error archivo {file} < Falta de headers")

print(teo_real_x)
print(len(teo_x))

plt.figure(figsize=(12, 4), dpi=600) # Ajuste de tamaño de la figura
plt.bar(teo_real_x, teo_y, width=1, label='Empirical', color='blue', align='edge', alpha=1) 
plt.savefig(os.path.join(act_dir, 'check.svg'))