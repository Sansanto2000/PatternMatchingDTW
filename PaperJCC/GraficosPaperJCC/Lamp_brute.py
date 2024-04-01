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
plt.figure(figsize=(12, 4), dpi=1200)

plt.bar(obs_x, obs_y, width=1.8, label='Extracted lamp', color='black', align='edge', alpha=1) # Lamp Bruto

plt.xlabel('Nº pixel')
plt.ylabel('Intensity')

plt.legend()

plt.savefig(os.path.join(act_dir, f"Lamp_brute.svg"))
plt.close()