import os
import numpy as np
import pandas as pd
from config import Config
from utils import zero_padding, slice_with_range_step

# Especificar lamparas a analizar
act_dir = os.path.dirname(os.path.abspath(__file__))
FILES = {   # Usamos las 3 lamparas Full 'Fe' calibradas por Meilan.
    "materials": ["Fe I", "Fe II"],
    "files": [
        os.path.join(act_dir, os.path.join("LampData","Fe","...")),
        os.path.join(act_dir, os.path.join("LampData","Fe","...")),
        os.path.join(act_dir, os.path.join("LampData","Fe","..."))
    ]
}

teorical_path = os.path.join(act_dir, "TeoricalData", "...")
config = Config(FILES=FILES, TEORICAL_PATH=teorical_path, SAVE_DIR=act_dir,
                WINDOW_STEP=25, WINDOW_LENGTH=3000, NORMALIZE_WINDOWS=True,
                ZERO_PADDING=True, DETECT_TEORICAL_PEAKS=False, 
                DETECT_EMPIRICAL_PEAKS=False, GRAPH_BESTS=False,
                OUTPUT_CSV_NAME="output.csv")

# Preparar CSV para persistencia de resultados
output_csv_path = os.path.join(act_dir, Config.OUTPUT_CSV_NAME)
try: # Si existe lo lee
    df = pd.read_csv(output_csv_path)
except FileNotFoundError: # Si no existe lo crea
    df = pd.DataFrame(columns=[
        'peaks in teorical data',
        'peaks in empirical data', 
        'Normalized', 
        'Zero Padding',
        'W_STEP', 
        'W_LENGTH', 
        'Count of Windows', 
        '(AVG) Distance', 
        '(STD) Distance', 
        '(AVG) IoU', 
        '(STD) IoU', 
        '(AVG) Scroll Error',
        '(STD) Scroll Error'
        ])
    df.to_csv(output_csv_path, index=False)

# Como teorico usamos los datos del carton que mando Yael.
csvpath=os.path.join(act_dir, "TeoricalData", "...")
teo_x, teo_y = get_Data(csvpath=Config.TEORICAL_PATH, normalize=True)

# Rellenado de ceros (si corresponde)
if (Config.ZERO_PADDING):
    teo_x, teo_y = zero_padding(arr_x=teo_x, arr_y=teo_y, dist=10)
    
# Ventanado del Teorico
ranges, slices_x, slices_y = slice_with_range_step(teo_x, teo_y, Config.WINDOW_LENGTH, 
                                                   Config.WINDOW_STEP, Config.NORMALIZE_WINDOWS)
        
#Filtrar aquellos arreglos que no tienen elementos
au_x = []
au_y = []
for i in range(len(slices_x)):
    if (len(slices_x[i]) > 0):
        au_x.append(slices_x[i])
        au_y.append(slices_y[i])
slices_x = np.array(au_x, dtype=object)
slices_y = np.array(au_y, dtype=object)

# Procesado de los archivos del conjunto
Distances = np.array([])
IoUs = np.array([])
EAMs = np.array([])
for file in FILES["files"]:
        
    # Separar informacion
    emp_x, emp_y, emp_head = extract_lamp_info(file, normalize=True)

    # Determinar calibracion real
    try:
        emp_real_x = emp_x * emp_head['CD1_1'] + emp_head['CRVAL1']
    except Exception as e:
        print(f"Error archivo {file} < Falta de headers")
        continue

    # Determinar calibracion con DTW
    best_alignment, index = find_best_calibration(emp_y, slices_y, W_LENGTH, W_STEP)

    # ACOMODAR A PARTIR DE ACA
    # Determinación de la metrica IoU 
    c_inicio = slices_x[index][best_alignment.index2[0]] # inicio calibrado
    c_fin = slices_x[index][best_alignment.index2[-1]] # fin calibrado
    Iou = IoU(c_inicio, c_fin, emp_real_x[0], emp_real_x[-1]) # Segun mejor calibrado

    # Determinar subconjunto teorico de interes
    #grap_teo_x, grap_teo_y = subconj_generator(teo_x, teo_y, emp_real_x[0], emp_real_x[-1])


# Aplicar DTW sobre cada archivo

# Guardar resultados:
#   - ¿Hubo match?
#       - No:
#           - Considerar guardar mejor match sin restricciones
#       - Si:
#           - Guardar parametros y metricas