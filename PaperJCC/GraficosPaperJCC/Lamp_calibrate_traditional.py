import os
import sys
import matplotlib.pyplot as plt
sys.path.append(os.path.dirname(
    os.path.dirname(os.path.abspath(os.path.dirname(__file__)))
    ))
from utils import get_Data_FILE

# Directorio actual
act_dir = os.path.dirname(os.path.abspath(__file__))

# Definicion de constantes
FINDDIR = r'C:\Users\santi\OneDrive\Documentos\Doctorado\PatternMatchingDTW\WCOMPs'
FILENAME = "WCOMP01.fits"

# Obtencion de empirico
obs_x, obs_y, obs_headers = get_Data_FILE(dirpath=FINDDIR, name=FILENAME, normalize=False)
obs_real_x = obs_x * obs_headers['CD1_1'] + obs_headers['CRVAL1']

# grafico de lamapara empirica en bruto (sin calibrar)
plt.figure(figsize=(24, 4), dpi=1200)

plt.bar(obs_real_x, obs_y, width=3, label='Calibrated lamp', color='black', align='edge', alpha=1) # Lamp Bruto

# Ajustar el espacio entre los ejes
plt.subplots_adjust(left=0.08, right=0.92, top=0.82, bottom=0.1)

# Ajustar el tamaño de la letra en las marcas de los ejes
plt.xticks(fontsize=18)
plt.yticks(fontsize=18)

plt.xlabel('Longitud de onda (Å)', fontsize=20)
plt.ylabel('Intensity', fontsize=20)

plt.legend(fontsize=18)

plt.savefig(os.path.join(act_dir, f"Lamp_calibrate_traditional.svg"))
plt.close()