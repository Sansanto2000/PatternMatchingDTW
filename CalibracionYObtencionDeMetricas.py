import IOU
import os
from DTW import DTW
from IOU import IoU
from tqdm import tqdm
from Calibration import Calibration
from utils import getfileData, normalize_min_max, gaussianice, slice_with_range_step, calibrate_with_observations
from NIST_Table_Interactor import NIST_Table_Interactor
import matplotlib.pyplot as plt
import numpy as np

def calibrate_and_obtain_metrics(filename:str, step:int, w_length:int, sigma:int) -> list:
    """Funcion que en base al nombre de un archivo FITS y un conjunto de parametros realiza iteradas calibraciones 
    DTW y determina las metricas obtenidas con cada una de estas

    Args:
        filename (str): Nombre del archivo FITS a analizar
        step (int): Cantidad de longitudes de onda entre cada inicio de ventana
        w_length (int): Cantidad de longitudes de onda que una ventana cubre a partir de su paso inicial
        sigma (int): Sigma a utilizar durante el proceso de gaussianizado

    Returns:
        Calibration: arreglo con las calibraciónes realizadas
    """
    # Datos y headers del observado
    script_dir = os.path.dirname(os.path.abspath(__file__))
    #filename = "EFBTCOMP31.fits"
    filepath = os.path.join(script_dir, "EFBTCOMPs")
    filepath = os.path.join(filepath, filename)
    obs_data, obs_headers = getfileData(filepath=filepath)

    # Longitud de onda minima y maxima del observado calibrado
    obs_long_min = obs_headers['CRVAL1'] + 0*obs_headers['CD1_1']
    obs_long_max = obs_headers['CRVAL1'] + (len(obs_data)-1)*obs_headers['CD1_1']

    # Separacion de datos observados para el eje X y el eje Y
    obs_x = range(len(obs_data))
    obs_y = obs_data

    # Normalizado de los datos en el eje Y
    obs_y, _, _ = normalize_min_max(obs_y)

    # Datos de teoricos del NIST
    csv_filename = os.path.join(script_dir, "Tabla(NIST)_Int_Long_Mat_Ref.csv")
    nisttr = NIST_Table_Interactor(csv_filename=csv_filename)

    # Especificacion del tipo de lampara a analizar
    filter = "He I"

    # Obtencion del dataframe
    teorico_df = nisttr.get_dataframe(filter=filter)

    # Separacion de datos teoricos para el eje X y el eje Y
    teo_x = teorico_df['Wavelength(Ams)'].tolist()
    teo_y = teorico_df['Intensity'].tolist()

    # Normalizado de los datos en el eje Y
    teo_y, _, _ = normalize_min_max(target=teo_y)

    # Definicion de constante a usar durante el gaussianizado
    SIGMA = sigma

    # Definicion de constantes para realizar la división del teorico en distintas ventanas
    STEP = step # Cantidad de longitudes de onda entre cada inicio de ventana
    W_RANGE = w_length # Rango de longitudes de onda que una ventana cubre a partir de su paso inicial

    # Recorte de los datos del teorico en conjuntos de ventanas
    ranges, slices_x, slices_y = slice_with_range_step(teo_x, teo_y, W_RANGE, STEP)

    Calibrations = []
    for k in tqdm(range(0, len(slices_x)), desc=filename):
        
        # Gaussianizado del recorte del teorico
        slice_x_gau, slice_y_gau = gaussianice(x=slices_x[k], y=slices_y[k], resolution=len(obs_x), sigma=SIGMA, rang=ranges[k])    
        
        # Normalizado del gaussinizado
        slice_y_gau, _, _ = normalize_min_max(slice_y_gau)
        
        # Aplicación DTW del observado respecto al gaussianizado
        cal_X, cal_Y, NorAlgCos = DTW(obs_y, slice_y_gau, slice_x_gau)
        
        # Determinación de la metrica IoU
        Iou = IoU(slice_x_gau[0], slice_x_gau[len(slice_x_gau)-1], obs_long_min, obs_long_max)
        
        # Agrupación de los datos relevantes de la calibración en un objeto y agregado a la lista calibrations
        Calibrations.append(Calibration(arr_X=cal_X, arr_Y=cal_Y, IoU=Iou, NaC=NorAlgCos))

    return Calibrations

# Archivos de datos observados para analizar
FILES = ["EFBTCOMP01.fits", "EFBTCOMP02.fits", "EFBTCOMP03.fits", "EFBTCOMP04.fits", "EFBTCOMP05.fits",
         "EFBTCOMP06.fits", "EFBTCOMP07.fits", "EFBTCOMP08.fits", "EFBTCOMP09.fits", "EFBTCOMP10.fits",
         "EFBTCOMP11.fits", "EFBTCOMP12.fits", "EFBTCOMP13.fits", "EFBTCOMP14.fits", "EFBTCOMP15.fits",
         "EFBTCOMP16.fits", "EFBTCOMP17.fits", "EFBTCOMP18.fits", "EFBTCOMP19.fits", "EFBTCOMP20.fits",
         "EFBTCOMP21.fits", "EFBTCOMP22.fits", "EFBTCOMP23.fits", "EFBTCOMP24.fits", "EFBTCOMP25.fits",
         "EFBTCOMP26.fits", "EFBTCOMP27.fits", "EFBTCOMP28.fits", "EFBTCOMP29.fits", "EFBTCOMP30.fits",
         "EFBTCOMP31.fits"]

# Definición de constantes
SIGMA = 50
STEP = 400 
W_RANGE = 1000

# inicializacion de variables para graficos posteriores
dist_entre_optima_y_determinada = {}

# Iteración para recorrer el arreglo de nombres de archivos
for filename in FILES:

    # Calibración del archivo
    calibrations = calibrate_and_obtain_metrics(filename=filename,step=STEP, w_length=W_RANGE, sigma=SIGMA)
    
    # Determinación de rango con mejor IoU y rango con mejor NAC (Normalize Aligment Cost)
    best_IoU_pos = 0
    best_NAC_pos = 0
    for i in range (0, len(calibrations)):
        
        if (calibrations[i].IoU > calibrations[best_IoU_pos].IoU):
            best_IoU_pos = i
            
        if (calibrations[i].NaC < calibrations[best_NAC_pos].NaC):
            best_NAC_pos = i
    
    # +1 en la posicion del arreglo de distancias que le corresponda a la diferencia entre best_IoU_pos y best_NAC_pos
    pos = abs(best_IoU_pos - best_NAC_pos)
    dist_entre_optima_y_determinada[pos] = dist_entre_optima_y_determinada.get(pos, 0) + 1

# Creación del arreglo a emplear para el grafico del histograma
arr_for_hist = []
dist_max = max(dist_entre_optima_y_determinada.keys()) + 1
for i in range(0, dist_max):
    arr_for_hist.append(dist_entre_optima_y_determinada.get(i, 0))

# Especifica los bins utilizando arange para asegurar de cubrir todos los enteros
bins = np.arange(max(arr_for_hist) + 2)

# Crear el diagrama de barras
plt.bar(range(len(arr_for_hist)), arr_for_hist, color='blue', align='center', alpha=0.7)

# Asegurar marcadores en todos los enteros del eje X
plt.xticks(range(len(arr_for_hist)))

# Añadir etiquetas y título
plt.xlabel('Cantidad de rangos de separación')
plt.ylabel('Cantidad')
plt.title('Dist. e/ mejor cal. segun IoU y NAC')

# Mostrar el histograma
plt.show()