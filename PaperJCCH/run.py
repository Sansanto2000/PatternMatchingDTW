import os
import numpy as np
import pandas as pd
from config import Config
import matplotlib.pyplot as plt
from scipy.signal import find_peaks
from utils import zero_padding, slice_with_range_step, extract_lamp_info, find_best_calibration
from utils import IoU, EAM, subconj_generator

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
                DETECT_EMPIRICAL_PEAKS=False, GRAPH=True,
                OUTPUT_CSV_NAME="output.csv")

# Preparar CSV para persistencia de resultados
output_csv_path = os.path.join(act_dir, config.OUTPUT_CSV_NAME)
try: # Si existe lo lee
    df = pd.read_csv(output_csv_path)
except FileNotFoundError: # Si no existe lo crea
    df = pd.DataFrame(columns=[
        'Teorical peaks',
        'Empirical peaks', 
        'Normalized', 
        'Zero Padding',
        'W_STEP', 
        'W_LENGTH', 
        'Count of Windows', 
        'Distance', 
        'IoU', 
        'Scroll Error'
        ])
    df.to_csv(output_csv_path, index=False)

# Como teorico usamos los datos del carton que mando Yael.
csvpath=os.path.join(act_dir, "TeoricalData", "...")
teo_x, teo_y = get_Data(csvpath=config.TEORICAL_PATH, normalize=True)

# Aislado de picos (si corresponde)
if (config.DETECT_TEORICAL_PEAKS):
    indices, _ = find_peaks(teo_y, threshold=[0.0, np.inf], height= [0.0, np.inf])
    teo_x = teo_x[indices]
    teo_y = teo_y[indices]

# Rellenado de ceros (si corresponde)
if (config.ZERO_PADDING):
    teo_x, teo_y = zero_padding(arr_x=teo_x, arr_y=teo_y, dist=10)
    
# Ventanado del Teorico
ranges, slices_x, slices_y = slice_with_range_step(teo_x, teo_y, config.WINDOW_LENGTH, 
                                                   config.WINDOW_STEP, config.NORMALIZE_WINDOWS)
        
#Filtrar aquellos arreglos que no tienen elementos
au_x = []
au_y = []
for i in range(len(slices_x)):
    if (len(slices_x[i]) > 0):
        au_x.append(slices_x[i])
        au_y.append(slices_y[i])
slices_x = np.array(au_x, dtype=object)
slices_y = np.array(au_y, dtype=object)

# Procesado de los archivos Fe
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
    
    # Aislado de picos (si corresponde)
    if (config.DETECT_EMPIRICAL_PEAKS):
        indices, _ = find_peaks(obs_y, threshold=[0.0, np.inf], height= [0.0, np.inf])
        obs_x = obs_x[indices]
        obs_y = obs_y[indices]
        obs_real_x = obs_real_x[indices]

    # Determinar calibracion con DTW
    best_alignment, index = find_best_calibration(emp_y, slices_y, 
                                                  config.WINDOW_LENGTH, 
                                                  config.WINDOW_STEP)
    
    # Dispocición en vector de las longitudes de ondas calibradas
    calibrado_x = np.full(len(best_alignment.index1), None)
    for i in range(len(best_alignment.index1)): # Calibrado
        calibrado_x[best_alignment.index1[i]] = slices_x[index][best_alignment.index2[i]]

    # ACOMODAR A PARTIR DE ACA
    # Determinación de la metrica IoU 
    c_inicio = slices_x[index][best_alignment.index2[0]] # inicio calibrado
    c_fin = slices_x[index][best_alignment.index2[-1]] # fin calibrado
    Iou = IoU(c_inicio, c_fin, emp_real_x[0], emp_real_x[-1]) # Segun mejor calibrado
    
    # Agregado de metricas en arreglos de almacenamiento
    Distances = np.append(Distances, best_alignment.distance)
    IoUs = np.append(IoUs, Iou)
    EAMs = np.append(EAMs, EAM(calibrado_x, emp_real_x))

    if(config.GRAPH): 
        plt.figure(figsize=(12, 4), dpi=600) # Ajuste de tamaño de la figura

        min_teo_grap = calibrado_x[0] if calibrado_x[0] < emp_real_x[0] else emp_real_x[0] # Seccion del teorico
        min_teo_grap -= config.WINDOW_STEP
        max_teo_grap = calibrado_x[-1] if calibrado_x[-1] > emp_real_x[-1] else emp_real_x[-1]
        max_teo_grap += config.WINDOW_STEP
        grap_teo_x, grap_teo_y, _, _ = subconj_generator(teo_x, teo_y, min_teo_grap, max_teo_grap)
        plt.bar(grap_teo_x, -grap_teo_y, width=10, label='Teorical', color='blue', align='edge', alpha=0.7) 
        plt.bar(emp_real_x, emp_y, width=2, label='Emp Real', color='black', align='edge', alpha=0.7) # Real
        plt.bar(calibrado_x, emp_y, width=3, label='Emp Calibrated', color='red', align='edge', alpha=1) # Hallado
        plt.legend()
        fig_name = f"{os.path.splitext(os.path.basename(file))[0]}"
        fig_name += "_EP" if(config.DETECT_EMPIRICAL_PEAKS) else ""
        fig_name += "_TP" if(config.DETECT_TEORICAL_PEAKS) else ""
        fig_name += "_NOR" if(config.NORMALIZE_WINDOWS) else ""
        fig_name += "_ZP" if(config.ZERO_PADDING) else ""
        plt.savefig(os.path.join(config.SAVE_DIR, f"{fig_name}.svg"))
        plt.close()
        
# Guardar datos promedios de ejecucion en CSV
nueva_fila = { # Añadir la nueva fila al DataFrame
    'Teorical peaks':config.DETECT_TEORICAL_PEAKS, 
    'Empirical peaks':config.DETECT_EMPIRICAL_PEAKS, 
    'Normalized':config.NORMALIZE_WINDOWS, 
    'Zero Padding':config.ZERO_PADDING,
    'W_STEP':config.WINDOW_STEP, 
    'W_LENGTH':config.WINDOW_LENGTH, 
    'Count of Windows':len(slices_x),
    'Distance':f"{np.mean(Distances)} \u00B1({np.std(Distances)})", 
    'IoU':f"{np.mean(IoUs)} \u00B1({np.std(IoUs)})",
    'Scroll Error':f"{np.mean(EAMs)} \u00B1({np.std(EAMs)})"
}
df = df._append(nueva_fila, ignore_index=True)
df.to_csv(output_csv_path, index=False) # Guardar DataFrame actualizado


# Aplicar DTW sobre cada archivo

# Guardar resultados:
#   - ¿Hubo match?
#       - No:
#           - Considerar guardar mejor match sin restricciones
#       - Si:
#           - Guardar parametros y metricas