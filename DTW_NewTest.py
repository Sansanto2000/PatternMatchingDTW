"""
Testeo para probar una nueva forma de ejecutar
"""
import os
import mplcursors
import numpy as np
from EM import EAM
from DTW import DTW
from IOU import IoU
import pandas as pd
from tqdm import tqdm
import matplotlib.pyplot as plt
from scipy.signal import find_peaks
from Calibration import Calibration
from utils import getfileData, normalize_min_max, slice_with_range_step, gaussianice_arr, subconj_generator, calibrate_with_observations
from NIST_Table_Interactor import NIST_Table_Interactor
from dtw import dtw, rabinerJuangStepPattern

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
#teo_x, teo_y, _, _ = subconj_generator(teo_x, teo_y, 0, 25000)

# Normalizado de los datos en el eje Y
teo_y, _, _ = normalize_min_max(target=teo_y)

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
    
    # Suavizado de las ventanas acorde a la cantidad de picos a emplear
    # slices_x_sua, slices_y_sua = gaussianice_arr(slices_x, slices_y, len(picos_y), SIGMA, ranges, normalize=False)
    
    # Calibrado por ventanas y obtencion de resultados
    #calibrations = calibrate_for_windows(slices_x_sua, slices_y_sua, picos_y, obs_long_min, obs_long_max)
    
    # Calibrado
    alignment = dtw(picos_y, teo_y, keep_internals=True)

    ## Display the warping curve, i.e. the alignment curve
    alignment.plot(type="threeway")
    
    ## Align and plot with the Rabiner-Juang type VI-c unsmoothed recursion
    dtw(picos_y, teo_y, keep_internals=True, 
        step_pattern=rabinerJuangStepPattern(6, "c"))\
        .plot(type="twoway",offset=-2)
    
    plt.show()
    
    # # Busqueda de la mejor calibración
    # best_calibration = None
    # best_EAM = float('inf')
    # for i in range(len(calibrations)):
    #     if (best_EAM > calibrations[i].EaM):
    #         best_EAM = calibrations[i].EaM
    #         best_calibration = calibrations[i]
    
    # # Ajusta el tamaño de figura para graficado
    # plt.figure(figsize=(10, 6))
    
    # # Grafica una seccion del teorico
    # min_teo_grap = best_calibration.arr_X[0] if best_calibration.arr_X[0] < obs_calib_real[0] else obs_calib_real[0]
    # min_teo_grap -= W_STEP
    # max_teo_grap = best_calibration.arr_X[-1] if best_calibration.arr_X[-1] > obs_calib_real[-1] else obs_calib_real[-1]
    # max_teo_grap += W_STEP
    # grap_teo_x, grap_teo_y, _, _ = subconj_generator(teo_x, teo_y, min_teo_grap, max_teo_grap)
    # # Iterar sobre los datos y dibujar barras
    # plt.bar([], [], width=0, label='Teorico', color='blue', align='edge', alpha=1)
    # for x, y in zip(grap_teo_x, grap_teo_y):
    #     # Dibujar la barra desde 0 hasta el punto (x, y)
    #     plt.bar([x], [y], width=8, align='edge', color='blue', alpha=0.7)
    # #plt.plot(grap_teo_x, grap_teo_y, label='Teorico', alpha=1, color='blue', linewidth=0.5, linestyle='--')
    
    # # Grafica el calibrado
    # plt.bar([0], [0], width=0, label='Emp Optimo', color='red', align='edge', alpha=1)
    # for x, y in zip(best_calibration.arr_X, best_calibration.arr_Y):
    #     # Dibujar la barra desde 0 hasta el punto (x, y)
    #     plt.bar([x], [y], width=2, align='edge', color='red', alpha=0.7)
    # #plt.plot(best_calibration.arr_X, best_calibration.arr_Y, label='Emp Calibrado', alpha=1, color='red', linewidth=0.5, linestyle='--')
    
    # # Grafica el empirico real
    # plt.bar([0], [0], width=0, label='Emp Real', color='black', align='edge', alpha=1)
    # for x, y in zip(obs_calib_real, obs_y):
    #     # Dibujar la barra desde 0 hasta el punto (x, y)
    #     plt.bar([x], [y], width=1, align='edge', color='black', alpha=0.7)
    # #plt.plot(obs_calib_real, obs_y, label='Emp Real', alpha=1, color='black', linewidth=0.5, linestyle='--')
        
    # mplcursors.cursor(hover=True) # Activar mplcursors
    # plt.legend()
    # save_location = os.path.join(SAVEPATH, f'{filename}_Graph.png')
    # plt.savefig(save_location)
    # #plt.show()
    # plt.close()
    
    # Corte temprano para etapas de prueba
    if (filename=="WCOMP02.fits"):
        break
    