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
SPEC_FINDDIR = r'C:\Users\santi\OneDrive\Documentos\Doctorado\PatternMatchingDTW\Ws'
SPEC_FILENAME = "WOBJ01.fits"

# Obtencion de espectro
spc_x, spc_y, spc_headers = get_Data_FILE(dirpath=SPEC_FINDDIR, name=SPEC_FILENAME, normalize=False)

# grafico de lamapara empirica en bruto (sin calibrar)
plt.figure(figsize=(24, 4), dpi=1200)

plt.bar(spc_x, spc_y, width=1, label='Spectre', color='green', align='edge', alpha=1) # Spectre brute

# Ajustar el espacio entre los ejes
plt.subplots_adjust(left=0.08, right=0.92, top=0.82, bottom=0.18)

# Ajustar el tamaño de la letra en las marcas de los ejes
# Añadir los tics del eje x en la parte superior
plt.xticks(fontsize=18)
plt.yticks(fontsize=18)

plt.xlabel('Nº pixel', fontsize=20)
plt.ylabel('Intensity', fontsize=20)

plt.legend(fontsize=18)

plt.savefig(os.path.join(act_dir, f"Spectre_brute.svg"))
plt.close()