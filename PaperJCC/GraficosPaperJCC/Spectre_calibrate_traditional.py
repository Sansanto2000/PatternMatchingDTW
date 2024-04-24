import os
import sys
import numpy as npm
import matplotlib.pyplot as plt
sys.path.append(os.path.dirname(
    os.path.dirname(os.path.abspath(os.path.dirname(__file__)))
    ))
from utils import get_Data_FILE, normalize_min_max

# Directorio actual
act_dir = os.path.dirname(os.path.abspath(__file__))

# Definicion de constantes
LAMP_FINDDIR = r'C:\Users\santi\OneDrive\Documentos\Doctorado\PatternMatchingDTW\WCOMPs'
#LAMP_FINDDIR = r'/home/sponte/Documentos/Doctorado/PatternMatchingDTW/WCOMPs'
LAMP_FILENAME = "WCOMP01.fits"
SPEC_FINDDIR = r'C:\Users\santi\OneDrive\Documentos\Doctorado\PatternMatchingDTW\Ws'
#SPEC_FINDDIR = r'/home/sponte/Documentos/Doctorado/PatternMatchingDTW/Ws'
SPEC_FILENAME = "WOBJ01.fits"

# Obtencion de empirico
obs_x, obs_y, obs_headers = get_Data_FILE(dirpath=LAMP_FINDDIR, name=LAMP_FILENAME, normalize=False)
obs_real_x = obs_x * obs_headers['CD1_1'] + obs_headers['CRVAL1']

# Obtencion de espectro
spc_x, spc_y, spc_headers = get_Data_FILE(dirpath=SPEC_FINDDIR, name=SPEC_FILENAME, normalize=False)
spc_real_x = obs_x * spc_headers['CD1_1'] + spc_headers['CRVAL1']

# Postprocesado
#obs_y = normalize_min_max(obs_y, npm.min(spc_y))

# grafico de lamapara empirica en bruto (sin calibrar)
plt.figure(figsize=(24, 4), dpi=1200)

plt.bar(obs_real_x, spc_y, width=3, label='Spectre', color='green', align='edge', alpha=1) # Spectre

plt.bar(obs_real_x, obs_y, width=2, label='Lamp', color='black', align='edge', alpha=1) # Lamp Bruto

# Ajustar el espacio entre los ejes
plt.subplots_adjust(left=0.08, right=0.92, top=0.82, bottom=0.18)

# Ajustar el tamaño de la letra en las marcas de los ejes
# Añadir los tics del eje x en la parte superior
#plt.tick_params(axis='x', top=True, labeltop=True, bottom=True, labelbottom=True)
plt.xticks(fontsize=18)
plt.yticks(fontsize=18)

plt.xlabel('Nº pixel', fontsize=20)
plt.xlabel('Longitud de onda (Å)', fontsize=20)
plt.ylabel('Intensity', fontsize=20)

plt.legend(fontsize=18)

plt.savefig(os.path.join(act_dir, f"Spectre_calibrate_traditional.svg"))
plt.close()