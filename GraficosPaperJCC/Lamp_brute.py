import os
import sys
import matplotlib.pyplot as plt
sys.path.append(os.path.dirname(os.path.abspath(os.path.dirname(__file__))))
from utils import get_Data_FILE

# Directorio actual
act_dir = os.path.dirname(os.path.abspath(__file__))

# Definicion de constantes
FINDDIR = os.path.join(os.path.dirname(act_dir), 'WCOMPs')
FILENAME = "WCOMP01.fits"

# Obtencion de empirico
obs_x, obs_y, obs_headers = get_Data_FILE(dirpath=FINDDIR, name=FILENAME, normalize=False)
obs_real_x = obs_x * obs_headers['CD1_1'] + obs_headers['CRVAL1']

# grafico de lamapara empirica en bruto (sin calibrar)
plt.figure(figsize=(12, 4), dpi=800)

plt.bar([0], [0], width=0, label='Lampara Observada', color='black', align='edge', alpha=1) # Lamp Bruto
for x, y in zip(obs_x, obs_y):
    plt.bar([x], [y], width=1, align='edge', color='black', alpha=1)

plt.xlabel('Nº pixel')
plt.ylabel('Intensity')

plt.legend()

plt.savefig(os.path.join(act_dir, f"Lamp_brute.png"))
plt.close()