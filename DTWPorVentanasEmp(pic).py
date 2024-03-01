"""
Aplicar proceso de deteccion por ventanas con Teo (suave), Emp (picos).
De cada archivo almacenar:
- ID de lámpara
- Cantidad de picos detectados
- Cantidad de ventanas probadas
- Parámetros de DTW
- NAC de la mejor ventana
- IoU mejor ventana
"""
import os
import numpy as np
from DTW import DTW
from IOU import IoU
import pandas as pd
from tqdm import tqdm
import matplotlib.pyplot as plt
from scipy.signal import find_peaks
from Calibration import Calibration
from utils import getfileData, normalize_min_max, slice_with_range_step, gaussianice_arr, subconj_generator, calibrate_with_observations
from NIST_Table_Interactor import NIST_Table_Interactor

def calibrate_for_windows(teo_slices_x:list, teo_slices_y:list, obs_y:list, obs_long_min:int, 
                          obs_long_max:int) -> list:
    """Funcion que en base a un conjunto de datos teorico y uno observado intenta realizar calibraciónes del observado
    en varias ventanas del teorico, de cada calibración realizada recopila varias metricas y las retorna como respuesta

    Args:
        teo_slices_x (list): Arreglo de segmentos de datos teorico para calibrar con el eje X
        teo_slices_y (list): Arreglo de arrglos de datos teorico para calibrar con el eje Y
        obs_y (list): Arreglo de datos observados para el eje Y
        obs_long_min (int): Valor minimo REAL que toman las longitudes de onda de los datos observados
        obs_long_max (int): Valor maximo REAL que toman las longitudes de onda de los datos observados

    Returns:
        list[Calibration]: Arreglo con las calibraciónes realizadas
    """
    
    Calibrations = []
    for i in range(0, len(teo_slices_x)):
        
        # Aplicación DTW del observado respecto al gaussianizado
        cal_X, cal_Y, NorAlgCos = DTW(obs_y, teo_slices_y[i], teo_slices_x[i])
        
        # Determinación de la metrica IoU
        Iou = IoU(teo_slices_x[i][0], teo_slices_x[i][-1], obs_long_min, obs_long_max)
        
        # Agrupación de los datos relevantes de la calibración en un objeto y agregado a la lista calibrations
        Calibrations.append(Calibration(arr_X=cal_X, arr_Y=cal_Y, IoU=Iou, NaC=NorAlgCos))

    return Calibrations

# Constantes
FILES = ["WCOMP01.fits", "WCOMP02.fits", "WCOMP03.fits", "WCOMP04.fits", "WCOMP05.fits",
         "WCOMP06.fits", "WCOMP07.fits", "WCOMP08.fits", "WCOMP09.fits", "WCOMP10.fits",
         "WCOMP11.fits", "WCOMP12.fits", "WCOMP13.fits", "WCOMP14.fits", "WCOMP15.fits",
         "WCOMP16.fits", "WCOMP17.fits", "WCOMP18.fits", "WCOMP19.fits", "WCOMP20_A.fits",
         "WCOMP20_A.fits", "WCOMP21.fits", "WCOMP22.fits", "WCOMP23.fits", "WCOMP24.fits", 
         "WCOMP25.fits", "WCOMP26.fits", "WCOMP27.fits", "WCOMP28.fits", "WCOMP29.fits", 
         "WCOMP30.fits", "WCOMP31.fits"] # Archivos a calibrar
HEIGHT = 0.025 # Altura minima a considerar para la busqueda de picos
W_STEP = 100 # Cantidad de longitudes de onda entre cada inicio de ventana
W_RANGE = 1900 # Rango de longitudes de onda que una ventana cubre
RESOLUTION = 300 # Resolucion con la que aplicar el suavizado del teorico 
SIGMA = 50 # Sigma a usar para el suavizado del teorico
SAVEPATH = os.path.join(os.path.dirname(os.path.abspath(__file__)), "DTWPorVentanasEmp(pic)") # Especificacion de carpeta para almacenar los graficos

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

# Achivado del teorico para reducir tiempos de testeo
teo_x, teo_y, _, _ = subconj_generator(teo_x, teo_y, 0, 25000)

# Normalizado de los datos en el eje Y
teo_y, _, _ = normalize_min_max(target=teo_y)

# Recorte de los datos del teorico en conjuntos de ventanas
ranges, slices_x, slices_y = slice_with_range_step(teo_x, teo_y, W_RANGE, W_STEP)

# Gaussianizado y normalizado de los recortes del teorico
slices_x, slices_y = gaussianice_arr(slices_x, slices_y, RESOLUTION, SIGMA, ranges, normalize=False)

# Declaración de diccionario donde se guardaran los datos a almacenar
metrics = {
    "ID": [],
    "Cant_Picos_Detectados": [],
    "Cant_Ventanas_Probadas": [],
    "SIGMA": [],
    "W_STEP": [],
    "W_RANGE": [],
    "RESOLUTION": [],
    "HEIGHT": [],
    "IoU_mejor_ventana": [],
    "NAC_mejor_ventana": []
}

for filename in tqdm(FILES, desc=f'Porcentaje de avance'):
    
    # Datos y headers del observado
    script_dir = os.path.dirname(os.path.abspath(__file__))
    filepath = os.path.join(script_dir, "WCOMPs")
    filepath = os.path.join(filepath, filename)
    obs_data, obs_headers = getfileData(filepath=filepath)

    # Longitud de onda minima y maxima del observado calibrado
    obs_long_min = obs_headers['CRVAL1'] + 0*obs_headers['CD1_1']
    obs_long_max = obs_headers['CRVAL1'] + (len(obs_data)-1)*obs_headers['CD1_1']
    obs_calib_real = calibrate_with_observations(obs_data, obs_headers['CRVAL1'], obs_headers['CD1_1'])
    
    # Separacion de datos observados para el eje X y el eje Y
    obs_x = range(len(obs_data))
    obs_y = obs_data

    # Normalizado de los datos obserbados en el eje Y
    obs_y, _, _ = normalize_min_max(obs_y)
    
    # Busqueda de los picos empiricos
    picos_x, _ = find_peaks(obs_y, height=HEIGHT)
    picos_y = obs_y[picos_x]
    
    # Calibrado por ventanas y obtencion de resultados
    calibrations = calibrate_for_windows(slices_x, slices_y, picos_y, obs_long_min, obs_long_max)
    
    # Busqueda de la mejor calibración
    best_calibration = None
    best_NAC = float('inf')
    for i in range(len(calibrations)):
        if (best_NAC > calibrations[i].NaC):
            best_NAC = calibrations[i].NaC
            best_calibration = calibrations[i]
    
    # Acomodado de los datos en el formato adecuado
    metrics["ID"].append(filename)
    metrics["Cant_Picos_Detectados"].append(len(picos_x))
    metrics["Cant_Ventanas_Probadas"].append(len(slices_x))
    metrics["SIGMA"].append(SIGMA)
    metrics["W_STEP"].append(W_STEP)
    metrics["W_RANGE"].append(W_RANGE)
    metrics["RESOLUTION"].append(RESOLUTION)
    metrics["HEIGHT"].append(HEIGHT)
    metrics["IoU_mejor_ventana"].append(best_calibration.IoU)
    metrics["NAC_mejor_ventana"].append(best_calibration.NaC)

    # # Conversion necesaria para el graficadoa
    # obs_x = np.array(obs_x)
    # picos_x = np.array(picos_x, dtype=int)
    
    # Ajusta el tamaño de la figura
    plt.figure(figsize=(10, 6))
    
    # Graficar la señal y los picos
    print([obs_long_min,obs_long_max])
    plt.bar(best_calibration.arr_X, best_calibration.arr_Y, label='Emp Calibrado', alpha=1, color='red', linewidth=2)
    plt.plot(obs_calib_real, obs_y, label='Emp Real', alpha=1, color='black', linewidth=2)
    # Falta la funcion
    #plt.plot(obs_x, obs_y, label='Empirico', alpha=1, color='black', linewidth=0.5, linestyle='--')

    plt.legend()
    save_location = os.path.join(SAVEPATH, f'{filename}_Graph.png')
    plt.savefig(save_location)
    plt.show()
    plt.close()

    
# Crear un DataFrame con los datos
df = pd.DataFrame(metrics)

# Nombre del archivo CSV
csv_name = "PruebaPorVentanasEmp(picos).csv"

# Escribir el DataFrame en el archivo CSV e informar por consola
df.to_csv(csv_name, index=False)
print(f"Datos de ejecución '{csv_name}' guardados exitosamente.")
    