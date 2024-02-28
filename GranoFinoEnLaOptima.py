import os
import numpy as np
from tqdm import tqdm
import matplotlib.pyplot as plt
from scipy.signal import find_peaks
from Calibration import Calibration
from NIST_Table_Interactor import NIST_Table_Interactor
from utils import normalize_min_max, getfileData, subconj_generator, gaussianice
from DTW import DTW
from IOU import IoU
import pandas as pd

def iteration_matrix_values(iteration:int, total:int):
    """Funcion para recuperar matriz de penalizado a usar respecto una iteracion 
    sobre un total

    Args:
        iteration (int): Numero de iteracion
        total (int): Cantidad total de iteraciones

    Returns:
        numpy.array: Arreglo de 3 elementos correspodientes a [penalizacion de coincidencia, 
        penalizacion de insertado, penalizacion de borrado]
    """
    
    # Declaracion de valores por defecto
    pty_match = 1
    pty_insert = 1
    pty_deletion = 1
    
    # Ajuste de los valores segun el valor de iteracion actual
    max_per_case = total / 7
    
    # Calculo de valor a remplazar para esta iteracion
    value = 1 + 0.5 * (iteration % max_per_case)
    
    if (iteration < max_per_case):
        pty_match = value
        
    elif (iteration < max_per_case*2):
        pty_insert = value
        
    elif (iteration < max_per_case*3):
        pty_deletion = value
        
    elif (iteration < max_per_case*4):
        pty_match = value
        pty_insert = value
        
    elif (iteration < max_per_case*5):
        pty_match = value
        pty_deletion = value
        
    elif (iteration < max_per_case*6):
        pty_insert = value
        pty_deletion = value
        
    elif (iteration < max_per_case*7):
        pty_match = value
        pty_insert = value
        pty_deletion = value
    
    return np.array([pty_match, pty_insert, pty_deletion])

# Constantes para flujo del algoritmo
PICOS_EMPIRICO = True
SUAVIZADO_TEORICO = False

# Definicion de constante a usar durante el suavisado del teorico en caso de que corresponda
SIGMA = 50

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

# Busqueda de los picos empiricos en caso de que sea necesario
if (PICOS_EMPIRICO):
    obs_x, _ = find_peaks(obs_y, height=0.025)
    obs_y = obs_y[obs_x]

    # # Conversion necesaria para el graficadoa
    # obs_x = np.array(obs_x)
    # picos_x = np.array(picos_x, dtype=int)

# # Graficar la señal y los picos
# plt.plot(obs_x, obs_y, label='Señal')
# plt.plot(picos_x, picos_y, 'ro', label='Picos', alpha=0.5)
# plt.legend()
# plt.savefig('peakFinder.png')
# #plt.show()
# plt.clf()

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

# Suavizado del teorico en caso de que corresponda
if(SUAVIZADO_TEORICO):
    teo_x, teo_y = gaussianice(x=teo_x, y=teo_y, resolution=len(obs_x), sigma=SIGMA)    

# Declaración de diccionario donde se guardaran las metricas
metrics = {
    "Pty_Match": [],
    "Pty_Insert": [],
    "Pty_Deletion": [],
    "IoU": [],
    "NAC": [],
}

# Definicion de constantes
ITERACIONES = 100*7

for i in tqdm(range(ITERACIONES), desc=f'Porcentaje de avance'):

    # Calculo de matriz de penalizado a usar
    penalty_matrix = iteration_matrix_values(i, ITERACIONES)

    # Aplicación DTW del teorico respecto al recorte correcto del teorico
    cal_X, cal_Y, NorAlgCos = DTW(obs_y, teo_y, teo_x, penalty_matrix)

    # Determinación de la metrica IoU
    Iou = IoU(teo_x[0], teo_x[-1], obs_long_min, obs_long_max)

    # Agrupación de los datos relevantes de la calibración en un objeto y agregado a la lista calibrations
    cal = Calibration(arr_X=cal_X, arr_Y=cal_Y, IoU=Iou, NaC=NorAlgCos) 

    metrics["Pty_Match"].append(penalty_matrix[0])
    metrics["Pty_Insert"].append(penalty_matrix[1])
    metrics["Pty_Deletion"].append(penalty_matrix[2])
    metrics["IoU"].append(cal.IoU)
    metrics["NAC"].append(cal.NaC)

# Crear un DataFrame con los datos
df = pd.DataFrame(metrics)

#plt.bar(picos_x, picos_y, label='Empirico', color='yellow') 
#plt.bar(teo_x, teo_y, label='Teorico', color='green')
# plt.bar(cal.arr_X, cal.arr_Y, label='calibrado', color='orange')
# plt.legend()
# plt.show()

# Constante para daefinir que cantidad de datos de los mejores resultados se guardaran
TOP_CANT = 5

# Arreglos donde se almacenaran los tops de mejores resultados
bests_IoUs = [{ "Pty_Match": None, "Pty_Insert": None, 
               "Pty_Deletion": None, "IoU": -float('inf'), "NAC": None }] * TOP_CANT
bests_NACs = [{ "Pty_Match": None, "Pty_Insert": None, 
               "Pty_Deletion": None, "IoU": None, "NAC": float('inf') }] * TOP_CANT

# Recorre los datos almacenados para determinar los mejores IoU y los de mejores NAC
for indice, fila in df.iterrows():
    
    # Variables auxiliares para almacenar los valores intermedios durante los desplazamientos en el top
    ant_IoU = None
    ant_NAC = None
    
    # Iteracion segun la cantidad de elementos que se busca almacenar del top
    for i in range(TOP_CANT):
        
        # Busqueda y guardado de los mejores valores de IoU
        if (ant_IoU is None) and (bests_IoUs[i]['IoU'] < fila['IoU']):
            # Revisa que no tenga los mismos valores
            if (bests_IoUs[i]['Pty_Match'] != fila['Pty_Match']) and (bests_IoUs[i]['Pty_Insert'] != fila['Pty_Insert']) \
                and (bests_IoUs[i]['Pty_Deletion'] != fila['Pty_Deletion']):
                ant_IoU = bests_IoUs[i]
                bests_IoUs[i] = fila.to_dict()
        elif not (ant_IoU is None):
            aux = bests_IoUs[i]
            bests_IoUs[i] = ant_IoU
            ant_IoU = aux
            
        # Busqueda y guardado de los mejores valores de NAC
        if (ant_NAC is None) and (bests_NACs[i]['NAC'] >= fila['NAC']):
            # Revisa que no tenga los mismos valores
            if (bests_NACs[i]['Pty_Match'] != fila['Pty_Match']) and (bests_NACs[i]['Pty_Insert'] != fila['Pty_Insert']) \
                and (bests_NACs[i]['Pty_Deletion'] != fila['Pty_Deletion']):
                ant_NAC = bests_NACs[i]
                bests_NACs[i] = fila.to_dict()
        elif not (ant_NAC is None):
            aux = bests_NACs[i]
            bests_NACs[i] = ant_NAC
            ant_NAC = aux
        

# Imprimir en consola top IoU
print(f"---------TOP {TOP_CANT} IoU-------------")
for i in range(TOP_CANT):
    print(f"{i}. | IoU={bests_IoUs[i]['IoU']} | NAC={bests_IoUs[i]['NAC']} | Match={bests_IoUs[i]['Pty_Match']} | \
        Insert={bests_IoUs[i]['Pty_Insert']} | Del={bests_IoUs[i]['Pty_Deletion']}")
    
print(f"---------TOP {TOP_CANT} NAC-------------")
for i in range(TOP_CANT):
    print(f"{i}. | IoU={bests_NACs[i]['IoU']} | NAC={bests_NACs[i]['NAC']} | Match={bests_NACs[i]['Pty_Match']} | \
        Insert={bests_NACs[i]['Pty_Insert']} | Del={bests_NACs[i]['Pty_Deletion']}")
    
# Nombre del archivo CSV donde se almacenaran los datos de ejecución
csv_name = "PenalizacionesGranoFino.csv"

# Escribir el DataFrame en el archivo CSV e informar por consola
df.to_csv(csv_name, index=False)
print(f"Datos de metricas '{csv_name}' guardados exitosamente.")