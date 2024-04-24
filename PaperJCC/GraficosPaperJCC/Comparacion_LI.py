import os
import sys
import numpy as np
import matplotlib.pyplot as plt
from scipy.signal import find_peaks
sys.path.append(os.path.dirname(
    os.path.dirname(os.path.abspath(os.path.dirname(__file__)))
    ))
from utils import get_Data_FILE, slice_with_range_step, zero_padding
from PaperJCC.Experimentos.Ejecucion import Config, find_best_calibration, get_Data_NIST, get_Data_LIBS

# Directorio actual
act_dir = os.path.dirname(os.path.abspath(__file__))

# Definicion de constantes
LAMP_FINDDIR = r'C:\Users\santi\OneDrive\Documentos\Doctorado\PatternMatchingDTW\WCOMPs'
LAMP_FILENAME = "WCOMP01.fits"
W_RANGE = 2000
W_STEP = 25


# Obtencion de empirico y empirico real
obs_x, obs_y, obs_headers = get_Data_FILE(dirpath=LAMP_FINDDIR, name=LAMP_FILENAME, normalize=False)
obs_real_x = obs_x * obs_headers['CD1_1'] + obs_headers['CRVAL1']

# Obtencion del teorico
# Obtencion de datos de referencia (Teorico) de donde corresponda
filter:list=["He I", "Ar I", "Ar II", "NeI"]
teo_x, teo_y = get_Data_NIST(dirpath=r'C:\Users\santi\OneDrive\Documentos\Doctorado\PatternMatchingDTW', 
                             name="Tabla(NIST)_Int_Long_Mat_Ref.csv", 
                             filter=filter)
# Aislado de picos teoricos(si corresponde)
if (True):
    indices, _ = find_peaks(teo_y, 
                            threshold=[0.0, np.inf],
                            height=[0.0, np.inf]
                            )
    teo_x = teo_x[indices]
    teo_y = teo_y[indices]

# Rellenado de ceros (si corresponde)
if (True):
    teo_x, teo_y = zero_padding(arr_x=teo_x, arr_y=teo_y, dist=10)

# Ventanas posibles
ranges, slices_x, slices_y = slice_with_range_step(arr_x=teo_x, arr_y=teo_y, LENGTH=W_RANGE, 
                                                   STEP=W_STEP, normalize=True)

#Filtrar aquellos arreglos que no tienen elementos
au_x = []
au_y = []
for i in range(len(slices_x)):
    if (len(slices_x[i]) > 0):
        au_x.append(slices_x[i])
        au_y.append(slices_y[i])
slices_x = np.array(au_x, dtype=object)
slices_y = np.array(au_y, dtype=object)

# Obtencion de calibrado con DTW
best_alignment, index = find_best_calibration(obs_y=obs_y, slices_y=slices_y, w_range=W_RANGE, w_step=W_STEP)
calibrado_x = np.full(len(best_alignment.index1), None)
for i in range(len(best_alignment.index1)): # Calibrado
    calibrado_x[best_alignment.index1[i]] = slices_x[index][best_alignment.index2[i]]
#print(best_alignment.index2)
    
# Obtencion de calibrado aplicando LI
c_inicio = slices_x[index][best_alignment.index2[0]] # inicio calibrado
c_fin = slices_x[index][best_alignment.index2[-1]]
#print(c_inicio, c_fin)
#print(obs_x.min(), obs_x.max())
calibrado_x_LI = np.interp(obs_x, (obs_x.min(), obs_x.max()), (c_inicio, c_fin))

# grafico de lamapara empirica en bruto (sin calibrar)
fig, axs = plt.subplots(1, 2, figsize=(24, 4), dpi=1200)

# Grafico izquierda
# axs[0].bar(obs_real_x, obs_y, width=3, label='Manual', color='black', align='edge', alpha=1)
axs[0].bar(calibrado_x, obs_y, width=3, label='Calibrated', color='red', align='edge', alpha=1)
axs[0].set_title('Best DTW', fontsize=24)
axs[0].set_xlabel('Longitud de onda (Å)', fontsize=20)
axs[0].set_ylabel('Intensity', fontsize=20)
axs[0].tick_params(axis='x', labelsize=18)
axs[0].tick_params(axis='y', labelsize=18)
axs[0].legend(fontsize=18)

# Grafico derecha
# axs[1].bar(obs_real_x, obs_y, width=3, label='Manual', color='black', align='edge', alpha=1)
axs[1].bar(calibrado_x_LI, obs_y, width=3, label='Calibrated', color='red', align='edge', alpha=1)
axs[1].set_title('Best DTW with LI', fontsize=24)
axs[1].set_xlabel('Longitud de onda (Å)', fontsize=20)
axs[1].set_ylabel('Intensity', fontsize=20)
axs[1].tick_params(axis='x', labelsize=18)
axs[1].tick_params(axis='y', labelsize=18)
axs[1].legend(fontsize=18)

# Ajustar el espacio entre los ejes
plt.subplots_adjust(left=0.06, right=0.94, top=0.84, bottom=0.2)

plt.savefig(os.path.join(act_dir, f"Comparacion_LI.svg"))
plt.close()