import os
import sys
import matplotlib.pyplot as plt
sys.path.append(os.path.dirname(
    os.path.dirname(os.path.abspath(os.path.dirname(__file__)))
    ))
from utils import get_Data_FILE, get_Data_NIST, gaussianice, normalize_min_max, subconj_generator

# Directorio actual
act_dir = os.path.dirname(os.path.abspath(__file__))

# Definicion de constantes
FINDDIR = r'C:\Users\santi\OneDrive\Documentos\Doctorado\PatternMatchingDTW\WCOMPs'
FILENAME = "WCOMP01.fits"
ROOTDIR = r'C:\Users\santi\OneDrive\Documentos\Doctorado\PatternMatchingDTW'

# Obtencion de Teorico
filter:list=["He I", "Ar I", "Ar II"]
teo_x, teo_y = get_Data_NIST(dirpath=ROOTDIR, filter=filter, normalize=True)

# Recorte de longitudes de onda de interes
teo_x, teo_y, _, _ = subconj_generator(teo_x, teo_y, 2000, 8000)
teo_y = -teo_y

# Obtencion de empirico
obs_x, obs_y, obs_headers = get_Data_FILE(dirpath=FINDDIR, name=FILENAME, normalize=True)
obs_real_x = obs_x * obs_headers['CD1_1'] + obs_headers['CRVAL1']
#obs_y /= 2.5

#teo_x, teo_y = gaussianice(teo_x, teo_y, 80000, 20)
#teo_y, _, _ = normalize_min_max(teo_y, min=0)
#obs_y, _, _ = normalize_min_max(obs_y, min=0)
#plt.plot(teo_x, teo_y, label='Datos de referencia', color='black', alpha=1) # Teorico

# grafico comparacion, empirico (observado) vs teorico (referencia NIST)
plt.figure(figsize=(24, 6), dpi=1200)

#plt.bar(teo_x, teo_y, width=10, label='Reference data', color='blue', align='edge', alpha=1) # Teorico
#plt.bar(obs_real_x, obs_y, width=3, label='Calibrated lamp', color='black', align='edge', alpha=1) # Lamp Bruto
plt.bar(obs_real_x+750, obs_y, width=3, label='Calculated', color='black', align='edge', alpha=1) # Calculated
plt.bar(obs_real_x, -obs_y, width=3, label='Real', color='green', align='edge', alpha=1) # Calculated

# Ajustar el espacio entre los ejes
plt.subplots_adjust(left=0.03, right=0.97, top=0.95, bottom=0.17)

plt.xticks(fontsize=18)
plt.yticks(fontsize=18)
plt.yticks([0]) # Ocultar valores del eje Y
plt.xlim(right=7000)

plt.xlabel('Wavelength (Å)', fontsize=20)
plt.ylabel('Intensity', fontsize=20)

# Ocultar lineas del recuadro
# ax = plt.gca()
# ax.spines['right'].set_visible(False)
# ax.spines['top'].set_visible(False)
# ax.spines['left'].set_visible(False)
# ax.spines['bottom'].set_visible(False)

plt.legend(fontsize=18)

plt.savefig(os.path.join(act_dir, f"Comparacion_de_distancias.svg"))
plt.close()