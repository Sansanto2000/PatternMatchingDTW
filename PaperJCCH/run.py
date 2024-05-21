import os
import numpy as np
import pandas as pd
from config import Config
from utils import run_calibrations, concateAndExtractTeoricalFiles

# Obtencion de datos teoricos
EPSILON = 0.055
act_dir = os.path.dirname(os.path.abspath(__file__))
teo_files = [
    os.path.join(act_dir, 'Cartones_Info', 'W_Lamp01.fits'),
    os.path.join(act_dir, 'Cartones_Info', 'W_Lamp02.fits'),
    os.path.join(act_dir, 'Cartones_Info', 'W_Lamp03.fits'),
    os.path.join(act_dir, 'Cartones_Info', 'W_Lamp04.fits'),
    os.path.join(act_dir, 'Cartones_Info', 'W_Lamp05.fits'),
    os.path.join(act_dir, 'Cartones_Info', 'W_Lamp06.fits'),
    os.path.join(act_dir, 'Cartones_Info', 'W_Lamp07.fits'),
    os.path.join(act_dir, 'Cartones_Info', 'W_Lamp08.fits'),
    os.path.join(act_dir, 'Cartones_Info', 'W_Lamp09.fits'),
    os.path.join(act_dir, 'Cartones_Info', 'W_Lamp10.fits')
    ]
teo_x, teo_y = concateAndExtractTeoricalFiles(teo_files, EPSILON)
teo_x = teo_x[teo_y >= 0]
teo_y = teo_y[teo_y >= 0]

# Especificar lamparas a analizar
act_dir = os.path.dirname(os.path.abspath(__file__))
files = {   # Usamos las 3 lamparas Full 'Fe' calibradas por Meilan.
    os.path.join(act_dir, os.path.join("Lamparas_Meilan","w1466Blamp1.fits")),
    # os.path.join(act_dir, os.path.join("Lamparas_Meilan","w1466Blamp2.fits")),
    # os.path.join(act_dir, os.path.join("Lamparas_Meilan","w1467Blamp1.fits")),
    # os.path.join(act_dir, os.path.join("Lamparas_Meilan","w1467Blamp2.fits")),
    # os.path.join(act_dir, os.path.join("Lamparas_Meilan","w1576Blamp1.fits")),
    # os.path.join(act_dir, os.path.join("Lamparas_Meilan","w1576Blamp2.fits")),
    # os.path.join(act_dir, os.path.join("Lamparas_Meilan","w1581Blamp1.fits")),
    # os.path.join(act_dir, os.path.join("Lamparas_Meilan","w1581Blamp2.fits")),
    # os.path.join(act_dir, os.path.join("Lamparas_Meilan","w1600Blamp1.fits"))
}

# Definicion de direccion para el guardado de resultados
save_dir = os.path.join(act_dir, 'output')

# Especificacion de variables de configuracion
config = Config(FILES=files, SAVE_DIR=save_dir, WINDOW_STEP=75, 
                WINDOW_LENGTH=2000, GRAPH=False, OUTPUT_CSV_NAME="output.csv")

# Preparar CSV para persistencia de resultados
output_csv_path = os.path.join(config.SAVE_DIR, config.OUTPUT_CSV_NAME)

# Ejecutar calibraciones para todas las combinaciones de interes
# total_iteraciones = 2*2*2
# iteracion_actual = 1
# for detect_teorical_peaks in [True, False]:
#     for detect_empirical_peaks in [True, False]:
#         for normalize_windows in [True, False]:
#             #for zero_padding_bool in [True, False]:
#             print(f"{iteracion_actual}/{total_iteraciones}")
#             run_calibrations(
#                 teo_x=teo_x, 
#                 teo_y=teo_y, 
#                 files=config.FILES,
#                 window_length=config.WINDOW_LENGTH,
#                 window_step=config.WINDOW_STEP,
#                 detect_teorical_peaks=detect_teorical_peaks,
#                 detect_empirical_peaks=detect_empirical_peaks,
#                 zero_padding_bool=False,
#                 normalize_windows=normalize_windows,
#                 save_dir=save_dir,
#                 graph=config.GRAPH,
#                 output_csv_path=output_csv_path
#                 )
#             iteracion_actual += 1

run_calibrations(
                teo_x=teo_x, 
                teo_y=teo_y, 
                files=config.FILES,
                window_length=config.WINDOW_LENGTH,
                window_step=config.WINDOW_STEP,
                detect_teorical_peaks=False,
                detect_empirical_peaks=False,
                zero_padding_bool=True,
                normalize_windows=True,
                save_dir=save_dir,
                graph=config.GRAPH,
                output_csv_path=output_csv_path
                )