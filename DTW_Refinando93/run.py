import os
import numpy as np
import pandas as pd
from utils import inspect_files_comparison, extract_lamp_info, get_Data_NIST, subconj_generator
from utils import slice_with_range_step, zero_padding, find_best_calibration, IoU

# Especificar lamparas a analizar
act_dir = os.path.dirname(os.path.abspath(__file__))
FILES = {
    "Ar": {
        "materials": ["Ar I", "Ar II"],
        "files": [
            os.path.join(act_dir, os.path.join("LampData","Ar","whn1.fits")),
            os.path.join(act_dir, os.path.join("LampData","Ar","whn7.fits"))
        ]
    },
    "CuNeAr": {
        "materials": ["Cu I", "Cu II", "Ne I", "Ne II", "Ar I", "Ar II"],
        "files": [
            os.path.join(act_dir, os.path.join("LampData","CuNeAr","cdcu-16.ec.fits")),
            os.path.join(act_dir, os.path.join("LampData","CuNeAr","cdcu-21.ec.fits")),
            os.path.join(act_dir, os.path.join("LampData","CuNeAr","cdcu-26.ec.fits")),
            os.path.join(act_dir, os.path.join("LampData","CuNeAr","cdcu-30.ec.fits")),
            os.path.join(act_dir, os.path.join("LampData","CuNeAr","cdcu-35.ec.fits")),
            os.path.join(act_dir, os.path.join("LampData","CuNeAr","cdcu-135.ec.fits")),
            os.path.join(act_dir, os.path.join("LampData","CuNeAr","cdcu-185.ec.fits")),
            os.path.join(act_dir, os.path.join("LampData","CuNeAr","cdcu-235.ec.fits")),
            os.path.join(act_dir, os.path.join("LampData","CuNeAr","cdcu-285.ec.fits")),
            os.path.join(act_dir, os.path.join("LampData","CuNeAr","cdcu-325.ec.fits")),
            os.path.join(act_dir, os.path.join("LampData","CuNeAr","cdcu16.fits")),
            os.path.join(act_dir, os.path.join("LampData","CuNeAr","cdcu21.fits")),
            os.path.join(act_dir, os.path.join("LampData","CuNeAr","cdcu26.fits")),
            os.path.join(act_dir, os.path.join("LampData","CuNeAr","cdcu30.fits")),
            os.path.join(act_dir, os.path.join("LampData","CuNeAr","cdcu35.fits")),
            os.path.join(act_dir, os.path.join("LampData","CuNeAr","cdcu135.fits")),
            os.path.join(act_dir, os.path.join("LampData","CuNeAr","cdcu185.ecc.fits")),
            os.path.join(act_dir, os.path.join("LampData","CuNeAr","cdcu185.fits")),
            os.path.join(act_dir, os.path.join("LampData","CuNeAr","cdcu235.fits")),
            os.path.join(act_dir, os.path.join("LampData","CuNeAr","cdcu275.fits")),
            os.path.join(act_dir, os.path.join("LampData","CuNeAr","cdcu325.fits")),
            os.path.join(act_dir, os.path.join("LampData","CuNeAr","cunear.fits")),
            os.path.join(act_dir, os.path.join("LampData","CuNeAr","cunear.R500.fits"))
        ]
    },
    "HeArNe": {
        "materials": ["Hg I", "Hg II", "He I", "He II", "Cd I", "Cd II"],
        "files": [
            os.path.join(act_dir, os.path.join("LampData","HeArNe","cdla.19.ec.fits")),
            os.path.join(act_dir, os.path.join("LampData","HeArNe","cdla.23.ec.fits")),
            os.path.join(act_dir, os.path.join("LampData","HeArNe","cdla.43.ec.fits")),
            os.path.join(act_dir, os.path.join("LampData","HeArNe","cdlam.35.ec.fits"))
        ]
    },
    "HgHeCd": {
        "materials": ["Hg I", "Hg II", "He I", "He II", "Cd I", "Cd II"],
        "files": [
            os.path.join(act_dir, os.path.join("LampData","HgHeCd","lamp2_FORS.fits"))
        ]
    },
    "NeAr": {
        "materials": ["Ne I", "Ne II", "Ar I", "Ar II"],
        "files": [
            os.path.join(act_dir, os.path.join("LampData","NeAr","wc01.ec.fits")),
            os.path.join(act_dir, os.path.join("LampData","NeAr","whn16.fits"))
        ]
    },
    "ThAr": {
        "materials": ["Th I", "Th II", "Ar I", "Ar II"],
        "files": [
            os.path.join(act_dir, os.path.join("LampData","ThAr","cdto16.fits")),
            os.path.join(act_dir, os.path.join("LampData","ThAr","cdto21.fits")),
            os.path.join(act_dir, os.path.join("LampData","ThAr","cdto26.fits")),
            os.path.join(act_dir, os.path.join("LampData","ThAr","cdto30.fits")),
            os.path.join(act_dir, os.path.join("LampData","ThAr","cdto30a.fits")),
            os.path.join(act_dir, os.path.join("LampData","ThAr","cdto30b.fits")),
            os.path.join(act_dir, os.path.join("LampData","ThAr","cdto35a.1.fits")),
            os.path.join(act_dir, os.path.join("LampData","ThAr","cdto135.fits")),
            os.path.join(act_dir, os.path.join("LampData","ThAr","cdto185.fits")),
            os.path.join(act_dir, os.path.join("LampData","ThAr","cdto235.fits")),
            os.path.join(act_dir, os.path.join("LampData","ThAr","cdto275.fits")),
            os.path.join(act_dir, os.path.join("LampData","ThAr","cdto325.fits")),
            os.path.join(act_dir, os.path.join("LampData","ThAr","cdto325a.fits")),
            os.path.join(act_dir, os.path.join("LampData","ThAr","thar_sophia.fits")),
            os.path.join(act_dir, os.path.join("LampData","ThAr","toar-DC-R1000.fits")),
            os.path.join(act_dir, os.path.join("LampData","ThAr","toar-DC-R2000.fits")),
            os.path.join(act_dir, os.path.join("LampData","ThAr","toar-DC.fits")),
            os.path.join(act_dir, os.path.join("LampData","ThAr","twrcomp03.fits")),
            os.path.join(act_dir, os.path.join("LampData","ThAr","wfbc03.ec.fits"))
        ]
    }
}

# inspect_files_comparison(act_dir=act_dir, files=FILES)
# # Calcular resolucion de forma dinamica, usando min max longitud de onda

# file = FILES["CuNeAr"]["files"][0]
# #file = FILES["Ar"]["files"][0]
# emp_x, emp_y, emp_head = extract_lamp_info(file, normalize=True)
# emp_real_x = emp_x * emp_head['CD1_1'] + emp_head['CRVAL1']
# print(emp_x)
# print(emp_real_x)
# print(emp_head['CD1_1'])
# print(emp_head['CRVAL1'])

# ---------------------------------------------------------------------
# Preparar CSV para persistencia de los datos
CSV_NAME = 'fullDataOutput.csv'
csv_path = os.path.join(act_dir, CSV_NAME)
try: # Si existe lo lee
    df = pd.read_csv(csv_path)
except FileNotFoundError: # Si no existe lo crea
    df = pd.DataFrame(columns=[
        'Source of reference', 
        'peaks in teorical data', 
        'Treshold teorical data', 
        'Height teorical data',
        'peaks in empirical data', 
        'Treshold empirical data', 
        'Height empirical data',
        'Normalized', 
        'Zero Padding',
        'Linear Interpolated',
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
    df.to_csv(csv_path, index=False)

# Recorre los conjuntos de lamparas de cada material
for material_set in FILES.keys():
        
    # Determinar teorico que le corresponde al conjunto
    csvpath=os.path.join(act_dir, "NIST", "Tabla(NIST)_Int_Long_Mat_Ref.csv")
    teo_x, teo_y = get_Data_NIST(csvpath=csvpath, filter=FILES[material_set]["materials"], 
                                 normalize=True)
    
    # Rellenado de ceros
    teo_x, teo_y = zero_padding(arr_x=teo_x, arr_y=teo_y, dist=10)
    
    # Ventanado del Teorico
    W_LENGTH = 2000
    W_STEP = 25
    NORMALIZE_WINDOW = True
    ranges, slices_x, slices_y = slice_with_range_step(teo_x, teo_y, W_LENGTH, W_STEP, NORMALIZE_WINDOW)
        
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
    for file in FILES[material_set]["files"]:
            
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