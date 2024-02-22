import os
import numpy as np
import matplotlib.pyplot as plt
from scipy.signal import find_peaks
from Calibration import Calibration
from NIST_Table_Interactor import NIST_Table_Interactor
from utils import normalize_min_max, getfileData, subconj_generator
from DTW import DTW
from IOU import IoU

# Datos y headers del observado
filename = "WCOMP01.fits"
script_dir = os.path.dirname(os.path.abspath(__file__))
filepath = os.path.join(script_dir, "WCOMPs")
filepath = os.path.join(filepath, filename)
obs_data, obs_headers = getfileData(filepath=filepath)

# Longitud de onda minima y maxima del observado calibrado
obs_long_min = obs_headers['CRVAL1'] + 0*obs_headers['CD1_1']
obs_long_max = obs_headers['CRVAL1'] + (len(obs_data)-1)*obs_headers['CD1_1']

# Separacion de datos observados para el eje X y el eje Y
obs_x = range(len(obs_data))
obs_y = obs_data

# Normalizado de los datos obserbados en el eje Y
obs_y, _, _ = normalize_min_max(obs_y)

# Busqueda de los picos empiricos
picos_x, _ = find_peaks(obs_y, height=0.025)
picos_y = obs_y[picos_x]

# Conversion necesaria para el graficado
obs_x = np.array(obs_x)
picos_x = np.array(picos_x, dtype=int)

# Graficar la señal y los picos
plt.plot(obs_x, obs_y, label='Señal')
plt.plot(picos_x, picos_y, 'ro', label='Picos', alpha=0.5)
plt.legend()
plt.savefig('peakFinder.png')
plt.show()
plt.clf()

# Datos de teoricos del NIST
script_dir = os.path.dirname(os.path.abspath(__file__))
csv_filename = os.path.join(script_dir, "Tabla(NIST)_Int_Long_Mat_Ref.csv")
nisttr = NIST_Table_Interactor(csv_filename=csv_filename)

# Especificacion del tipo de lampara a analizar
filter = ["He I", "Ar I", "Ar II"]

# Obtencion del dataframe
teorico_df = nisttr.get_dataframe(filter=filter)

# Separacion de datos teoricos para el eje X y el eje Y
teo_x = teorico_df['Wavelength(Ams)'].tolist()
teo_y = teorico_df['Intensity'].tolist()

# Recorte del segmento de longitudes de onda coincidente con el empirico
teo_x, teo_y, _, _ = subconj_generator(teo_x, teo_y, obs_long_min, obs_long_max)

# Normalizado de los datos teoricos en el eje Y
teo_y, _, _ = normalize_min_max(target=teo_y)

# Aplicación DTW del teorico respecto al recorte correcto del teorico
cal_X, cal_Y, NorAlgCos = DTW(picos_y, teo_y, teo_x)

# Determinación de la metrica IoU
Iou = IoU(teo_x[0], teo_x[-1], obs_long_min, obs_long_max)

# Agrupación de los datos relevantes de la calibración en un objeto y agregado a la lista calibrations
cal = Calibration(arr_X=cal_X, arr_Y=cal_Y, IoU=Iou, NaC=NorAlgCos)

print(f"IoU={cal.IoU}, NAC={cal.NaC}")

"""Hacer un iterador que aplique variaciones en pesos de penalizacion y borrado"""
#plt.bar(picos_x, picos_y, label='Empirico', color='yellow') 
#plt.bar(teo_x, teo_y, label='Teorico', color='green')
plt.bar(cal.arr_X, cal.arr_Y, label='calibrado', color='orange')
plt.legend()
plt.show()



