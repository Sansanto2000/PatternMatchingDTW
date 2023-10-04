from NIST_Table_Interactor import NIST_Table_Interactor

from astropy.io import fits

from fastdtw import fastdtw
from scipy.spatial.distance import euclidean
import matplotlib.pyplot as plt
import numpy as np

def getfileData(filepath):    
    
    hdul = fits.open(filepath) 
    
    if('WOBJ' in filepath):
        headers = hdul[0].header
        data = hdul[0].data[0][0]
    else:
        headers = hdul[0].header
        data = hdul[0].data
        
    return data

def gaussian(x, max_value, mean, std_deviation):
    return max_value * np.exp(-((x - mean) ** 2) / (2 * std_deviation ** 2))
    
filepath = "C:\\Users\\santi\\OneDrive\\Documentos\\Doctorado\\EFBTCOMP31.fits"
#filepath = "C:\\Users\\santi\\OneDrive\\Documentos\\Doctorado\\EFBTOBJ31.fits"

spectral_data = getfileData(filepath=filepath)

csv_filename = ".\\Tabla(NIST)_Int_Long_Mat_Ref.csv"
nisttr = NIST_Table_Interactor(csv_filename=csv_filename)

df = nisttr.get_dataframe(250)
x = df['Wavelength(Ams)'].tolist()
y = df['Intensity'].tolist()

plt.figure(figsize=(12, 6))
plt.plot(x, y)
plt.xlabel('Wavelength(Ams)')
plt.ylabel('Intensity')
plt.title('Mi Gráfico')
plt.show()

NIST_Reference = nisttr.get_gaussianized_data()
    # Extraer con pandas datos de Tabla(NIST)_Int_Long_Mat_Ref.csv
    # ANTES de extraer probablemente haya que sanitizar los datos.


# # OPCION 1: pip install fastdtw
# Supongamos que tienes dos listas de puntos o series temporales llamadas serie1 y serie2
#distancia, camino = fastdtw(spectral_data, serie2, dist=euclidean)

# distancia es la distancia DTW entre las series temporales
# camino contiene el alineamiento entre las series

# # OPCION 2: pip install tslearn
# from tslearn.metrics import dtw
# import numpy as np

# # Supongamos que tienes dos matrices numpy llamadas serie1 y serie2
# distancia, camino = dtw(serie1, serie2)

# # distancia es la distancia DTW entre las series temporales
# # camino contiene el alineamiento entre las series