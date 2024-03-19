import os
import sys
import matplotlib.pyplot as plt
sys.path.append(os.path.dirname(os.path.abspath(os.path.dirname(__file__))))
from utils import get_Data_FILE, get_Data_NIST, gaussianice, normalize_min_max, subconj_generator

# Directorio actual
act_dir = os.path.dirname(os.path.abspath(__file__))

# Definicion de constantes
FINDDIR = os.path.join(os.path.dirname(act_dir), 'WCOMPs')
FILENAME = "WCOMP01.fits"

# Obtencion de Teorico
filter:list=["He I", "Ar I", "Ar II"]
teo_x, teo_y = get_Data_NIST(dirpath=os.path.dirname(act_dir), filter=filter, normalize=True)

# Recorte de longitudes de onda de interes
teo_x, teo_y, _, _ = subconj_generator(teo_x, teo_y, 0, 22000)
teo_y = -teo_y

# Obtencion de empirico
obs_x, obs_y, obs_headers = get_Data_FILE(dirpath=FINDDIR, name=FILENAME, normalize=True)
obs_real_x = obs_x * obs_headers['CD1_1'] + obs_headers['CRVAL1']
obs_y /= 2.5

# grafico comparacion, empirico (observado) vs teorico (referencia NIST)
plt.figure(figsize=(12, 4), dpi=800)

plt.bar([0], [0], width=0, label='Referencia', color='blue', align='edge', alpha=1) # Teorico
for x, y in zip(teo_x, teo_y):
    plt.bar([x], [y], width=10, align='edge', color='blue', alpha=1)

# teo_x, teo_y = gaussianice(teo_x, teo_y, 80000, 20)
# teo_y, _, _ = normalize_min_max(teo_y)
# plt.plot(teo_x, teo_y, label='Datos de referencia', color='black', alpha=1) # Teorico
    
plt.bar([0], [0], width=0, label='Lampara Observada', color='black', align='edge', alpha=1) # Lamp Bruto
for x, y in zip(obs_real_x, obs_y):
    plt.bar([x], [y], width=3, align='edge', color='black', alpha=1)

plt.xlabel('Longitud de onda (Å)')
plt.ylabel('Intensity')

# Ocultar valores del eje Y
plt.yticks([]) 

# Ocultar lineas del recuadro
# ax = plt.gca()
# ax.spines['right'].set_visible(False)
# ax.spines['top'].set_visible(False)
# ax.spines['left'].set_visible(False)
# ax.spines['bottom'].set_visible(False)

plt.legend()

plt.savefig(os.path.join(act_dir, f"Comparacion_de_distancias.png"))
plt.close()