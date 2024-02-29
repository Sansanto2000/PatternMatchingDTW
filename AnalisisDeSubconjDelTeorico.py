import os
import random
from tqdm import tqdm
from Calibration import Calibration
from NIST_Table_Interactor import NIST_Table_Interactor
from utils import normalize_min_max, slice_with_range_step, gaussianice, subconj_generator
from DTW import DTW
from IOU import IoU
import pandas as pd

def calibrate_and_obtain_metrics(teo_x:list, teo_y:list, obs_x:list, obs_y:list, obs_long_min:int, 
                                 obs_long_max:int,step:int,  w_length:int, sigma:int) -> list:
    """Funcion que en base a un conjunto de datos teorico y uno observado intenta realizar calibraciónes del observado
    en varias ventanas del teorico, de cada calibración realizada recopila varias metricas y las retorna como respuesta

    Args:
        teo_x (list): Arreglo de datos teorico para el eje X
        teo_y (list): Arreglo de datos teorico para el eje Y
        obs_x (list): Arreglo de datos observados para el eje X
        obs_y (list): Arreglo de datos observados para el eje Y
        obs_long_min (int): Valor minimo REAL que toman las longitudes de onda de los datos observados
        obs_long_max (int): Valor maximo REAL que toman las longitudes de onda de los datos observados
        step (int): Cantidad de longitudes de onda entre cada inicio de ventana
        w_length (int): Cantidad de longitudes de onda que una ventana cubre a partir de su paso inicial
        sigma (int): Sigma a utilizar durante el proceso de gaussianizado

    Returns:
        list[Calibration]: Arreglo con las calibraciónes realizadas
    """
    # Recorte de los datos del teorico en conjuntos de ventanas
    ranges, slices_x, slices_y = slice_with_range_step(teo_x, teo_y, w_length, step)

    Calibrations = []
    for k in range(0, len(slices_x)):
        
        # Gaussianizado del recorte del teorico
        slice_x_gau, slice_y_gau = gaussianice(x=slices_x[k], y=slices_y[k], resolution=len(obs_x), sigma=sigma, rang=ranges[k])    
        
        # Normalizado del gaussinizado
        slice_y_gau, _, _ = normalize_min_max(slice_y_gau)
        
        # Aplicación DTW del observado respecto al gaussianizado
        cal_X, cal_Y, NorAlgCos = DTW(obs_y, slice_y_gau, slice_x_gau)
        
        # Determinación de la metrica IoU
        Iou = IoU(slice_x_gau[0], slice_x_gau[-1], obs_long_min, obs_long_max)
        
        # Agrupación de los datos relevantes de la calibración en un objeto y agregado a la lista calibrations
        Calibrations.append(Calibration(arr_X=cal_X, arr_Y=cal_Y, IoU=Iou, NaC=NorAlgCos))

    return Calibrations
    
# Datos de teoricos del NIST
script_dir = os.path.dirname(os.path.abspath(__file__))
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

# Definición de constantes
SIGMA = 50 # Sigma a considerar durante el gaussianizado
STEP = 300 # Medida de pasos entre inicios de ventanas
W_RANGE = 2500 # Tamaño de ventanas
ITERATIONS = 10 # Cantidad de iteraciones a considerar para la evaluación de subconjuntos

# Declaración de diccionario donde se guardaran las metricas
metrics = {
    "Pos de +IoU": [],
    "+IoU": [],
    "NAC": [],
    "Pos de +NAC": [],
    "IoU": [],
    "+NAC": []
}

for i in tqdm(range(ITERATIONS), desc=f'Porcentaje de avance'):
    
    # Determinacion de longitudes de onda minima y maxima del subconjunto a considerar
    sub_min = random.randint(0, 25000)
    sub_max = sub_min + random.randint(W_RANGE - STEP, W_RANGE + STEP)

    # Determinación de subconjunto del teorico a usar como observado
    sub_x, sub_y, _, _ = subconj_generator(teo_x, teo_y, sub_min, sub_max)

    # Tratamiento de los datos del eje X para que parescan observados
    sub_x, sub_y = gaussianice(x=sub_x, y=sub_y, resolution=300, sigma=SIGMA, rang=(sub_min, sub_max))    
    sub_x = range(len(sub_x))

    # Normalizado de los datos del subconjunto
    sub_y, _, _ = normalize_min_max(sub_y)
    
    # Calibración y obtencion de metricas de los datos del subconjunto respecto a los datos del teorico
    calibrations = calibrate_and_obtain_metrics(teo_x, teo_y, sub_x, sub_y, sub_min, sub_max, STEP, W_RANGE, SIGMA)
    
    # Declaración de metricas a resguardar en el dataframe
    best_IoU_range:int
    best_IoU:float = float('-inf')
    NAC_of_best_IoU:float
    best_NAC_range:int
    best_NAC:float = float('inf')
    IoU_of_best_NAC:float
    
    # Recorrer arreglo de calibrados para recuperar las metricas de interes
    for j in range(len(calibrations)):
        
        if (best_IoU < calibrations[j].IoU):
            best_IoU = calibrations[j].IoU
            best_IoU_range = j
            NAC_of_best_IoU = calibrations[j].NaC
            
        if (best_NAC > calibrations[j].NaC):
            best_NAC = calibrations[j].NaC
            best_NAC_range = j
            IoU_of_best_NAC = calibrations[j].IoU
    
    # Acomodado de los datos en el formato adecuado
    metrics["Pos de +IoU"].append(best_IoU_range)
    metrics["+IoU"].append(best_IoU)
    metrics["NAC"].append(NAC_of_best_IoU)
    metrics["Pos de +NAC"].append(best_NAC_range)
    metrics["IoU"].append(IoU_of_best_NAC)
    metrics["+NAC"].append(best_NAC)
    
# Crear un DataFrame con los datos
df = pd.DataFrame(metrics)

# Nombre del archivo CSV
csv_name = "analicisDeSubconjuntosTeoricos.csv"

# Escribir el DataFrame en el archivo CSV e informar por consola
df.to_csv(csv_name, index=False)
print(f"Datos de metricas '{csv_name}' guardados exitosamente.")