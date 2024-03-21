import os
from dtw import *
import numpy as np
import pandas as pd
sys.path.append(os.path.dirname(os.path.abspath(os.path.dirname(__file__))))
from EM import EAM
from IOU import IoU
from tqdm import tqdm
import matplotlib.pyplot as plt
from scipy.signal import find_peaks
from NIST_Table_Interactor import NIST_Table_Interactor
from utils import getfileData, normalize_min_max, slice_with_range_step, subconj_generator, get_Data_LIBS, get_Data_FILE, get_Data_NIST

class Config:
    
    def __init__(self, finddir:str = os.path.dirname(os.path.abspath(__file__)), files:list=['WCOMP01.fits'], 
                 treshold_emp:float=0.0, treshold_teo:float=0.0, w_step:int=100, w_range:int=1900, 
                 picos_empirico:bool = False, picos_teorico:bool = False, 
                 savepath:str=os.path.dirname(os.path.abspath(__file__)), 
                 csv_name:str="DTW_AlternativoXVentanas.csv", use_nist_data:bool=True,
                 normalize_window:bool=False, teo_name:str="Tabla(NIST)_Int_Long_Mat_Ref.csv"):
        """Funcion de inicializacion de la clase Config

        Args:
            finddir (str, optional): Especificación de carpeta donde buscar los archivos. Defaults to 
            os.path.dirname(os.path.abspath(__file__)).
            files (list, optional): Arreglo con los nombres de los archivos a calibrar. Defaults to ['WCOMP01.fits'].
            treshold_emp (float, optional): Altura minima para la busqueda de picos en el empirico. Defaults to 0.0.
            treshold_teo (float, optional): Altura minima para la busqueda de picos en el teorico. Defaults to 0.0.
            picos_empirico (bool, optional): Condicion boleana par saber si se recortara el empirico segun los picos
            disponibles o no. Defaults to False.
            picos_teorico (bool, optional): Condicion boleana par saber si se recortara el empirico segun los picos
            disponibles o no. Defaults to False.
            w_step (int, optional): Cantidad de longitudes de onda entre inicios de ventanas. Defaults to 100.
            w_range (int, optional): Rango de longitudes de onda que una ventana cubre. Defaults to 1900.
            savepath (str, optional): Especificación de carpeta para almacenar los graficos. Defaults to 
            os.path.dirname(os.path.abspath(__file__)).
            csv_name (str, optional): Nombre del archivo CSV. Defaults to "DTW_AlternativoXVentanas.csv".
            use_nist_data (bool, optional): Condicion boleana par saber si se usaran los datos del NIST o de otra 
            fuente alternativa.
            normalize_window (bool, optional): Condicion boleana par saber si se normalizaran las ventanas teoricas 
            posterior a su segmentado.
            teo_name (str, optional): Nombre del archivo de donde se extraeran los datos teoricos.
        """
        self.FINDDIR = finddir
        self.FILES = files
        self.TRESHOLD_EMP = treshold_emp
        self.TRESHOLD_TEO = treshold_teo
        self.W_STEP = w_step
        self.W_RANGE = w_range
        self.PICO_EMPIRICO = picos_empirico
        self.PICO_TEORICO = picos_teorico
        self.SAVEPATH = savepath
        self.CSV_NAME = csv_name
        self.USE_NIST_DATA = use_nist_data
        self.NORMALIZE_WINDOW = normalize_window
        self.TEO_NAME = teo_name
    
    def __str__(self) -> str:
        """Metodo de especificacion de cadena (str) representativa de instancia de la clase

        Returns:
            str: Cadena representativa de la instancia de la clase
        """
        return f"FILE={self.FILE}, \
            TRESHOLD={self.TRESHOLD}, \
            W_STEP={self.W_STEP}, \
            W_RANGE={self.W_RANGE}, \
            W_SAVEPATH={self.W_SAVEPATH}, \
            W_CSV_NAME={self.W_CSV_NAME}]"

def find_best_calibration(obs_y:np.ndarray, slices_y:np.ndarray, w_range:int, w_step:int):
    """Funcion para hallar la mejor calibracion de un conjunto de datos observados respecto a un conjunto
    de posibles alternativas

    Args:
        obs_y (np.ndarray): Arreglo de datos a calibrar
        slices_y (np.ndarray): Conjunto de arreglos entre los que se debe encontrar la mejor coincidencia
        w_range (int, optional): Rango de longitudes de onda que una ventana cubre. Defaults to 1900
        w_step (int, optional): Cantidad de longitudes de onda entre inicios de ventanas. Defaults to 100

    Returns:
        DTW: Objeto DTW con detalles de la mejor coincidencia encontrada entre la funcion a comparar y las 
        posibles funciones objetivo
        int: Indice de la ventana que se corresponde con la mejor calibracion
        float: IoU correspondiente a la mejor calibración
    """

    # Plantillas para acumulado de resultados
    alignments = np.array([])
    IoUs = np.array([])
        
    for i in range(0, len(slices_y)):
        
        # Aplicación DTW del observado respecto al gaussianizado
        alignment = dtw(obs_y, slices_y[i], keep_internals=True, step_pattern=asymmetric, open_begin=True, open_end=True)
        
        # Determinación de la metrica IoU
        w_inicio = w_step*i
        w_fin = w_inicio + w_range
        # Iou = IoU(w_inicio, w_fin, obs_real_x[0], obs_real_x[-1]) # Segun mejor ventana
        c_inicio = slices_x[i][alignment.index2[0]]
        c_fin = slices_x[i][alignment.index2[-1]]
        Iou = IoU(c_inicio, c_fin, obs_real_x[0], obs_real_x[-1]) # Segun mejor calibrado
        
        # Agrega datos a arreglo
        alignments = np.append(alignments, alignment)
        IoUs = np.append(IoUs, Iou)
    
    # Busca el alineado con mejor metrica de distancia
    index = np.argmin([alig.distance for alig in alignments])
    best_alignment = alignments[index]
    
    return best_alignment, index, IoUs[index]


# Directorio actual
act_dir = os.path.dirname(os.path.abspath(__file__))

# Variables de configuracion
findDir = os.path.join(os.path.dirname(act_dir), 'WCOMPs')
files = ["WCOMP01.fits", "WCOMP02.fits", "WCOMP03.fits", "WCOMP04.fits", "WCOMP05.fits",
         "WCOMP06.fits", "WCOMP07.fits", "WCOMP08.fits", "WCOMP09.fits", "WCOMP10.fits",
         "WCOMP11.fits", "WCOMP12.fits", "WCOMP13.fits", "WCOMP14.fits", "WCOMP15.fits",
         "WCOMP16.fits", "WCOMP17.fits", "WCOMP18.fits", "WCOMP19.fits", "WCOMP20_A.fits",
         "WCOMP20_A.fits", "WCOMP21.fits", "WCOMP22.fits", "WCOMP23.fits", "WCOMP24.fits", 
         "WCOMP25.fits", "WCOMP26.fits", "WCOMP27.fits", "WCOMP28.fits", "WCOMP29.fits", 
         "WCOMP30.fits", "WCOMP31.fits"]
CONFIG = Config(files=files, 
                finddir=findDir,
                #teo_name="Tabla(NIST)_Int_Long_Mat_Ref.csv",
                #teo_name="LIBS_He_Ar_Ne_Resolution=100.csv",
                teo_name="LIBS_He_Ar_Ne_Resolution=1000.csv",
                use_nist_data=False,
                w_step=25, 
                picos_empirico=False, 
                picos_teorico=True, 
                treshold_emp=0.025,
                treshold_teo=0.000001,
                normalize_window=True)

# Preparar CSV para persistencia de los datos
df = pd.DataFrame(columns=['File', 'Cant_Ventanas_Probadas', 'W_STEP', 'W_RANGE', 
                                 'Distance_mejor_ventana', 'IoU_mejor_ventana', 'Error_de_desplazamiento'])
csv_path = os.path.join(CONFIG.SAVEPATH, CONFIG.CSV_NAME)
df.to_csv(csv_path, index=False)

# Obtencion de Teorico (de donde corresponda)
filter:list=["He I", "Ar I", "Ar II"]
if (CONFIG.USE_NIST_DATA):
    teo_x, teo_y = get_Data_NIST(dirpath=os.path.dirname(act_dir), name=CONFIG.TEO_NAME, 
                filter=filter)
else:
    teo_x, teo_y = get_Data_LIBS(dirpath=os.path.dirname(act_dir), name=CONFIG.TEO_NAME)

# Aislado de picos del Teorico. Solo si corresponde
if (CONFIG.PICO_TEORICO):
    #indices, _ = find_peaks(teo_y, height=CONFIG.TRESHOLD_TEO)
    indices, _ = find_peaks(teo_y, threshold=[CONFIG.TRESHOLD_TEO, np.inf])
    teo_x = teo_x[indices]
    teo_y = teo_y[indices]

# Ventanado del Teorico
ranges, slices_x, slices_y = slice_with_range_step(teo_x, teo_y, CONFIG.W_RANGE, CONFIG.W_STEP, CONFIG.NORMALIZE_WINDOW)
    
#Filtrar aquellos arreglos que no tienen elementos
au_x = []
au_y = []
for i in range(len(slices_x)):
    if (len(slices_x[i]) > 0):
        au_x.append(slices_x[i])
        au_y.append(slices_y[i])
slices_x = np.array(au_x, dtype=object)
slices_y = np.array(au_y, dtype=object)

# --------PROCESADO DE ARCHIVO-----------
for file in tqdm(CONFIG.FILES, desc=f'Porcentaje de avance'):

    # Obtencion de empirico
    obs_x, obs_y, obs_headers = get_Data_FILE(dirpath=CONFIG.FINDDIR, name = file)
    obs_real_x = obs_x * obs_headers['CD1_1'] + obs_headers['CRVAL1']
    
    # Aislado de picos del Empirico. Solo si corresponde
    if (CONFIG.PICO_EMPIRICO):
        indices, _ = find_peaks(obs_y, height=CONFIG.TRESHOLD_EMP)
        obs_x = obs_x[indices]
        obs_y = obs_y[indices]
        obs_real_x = obs_real_x[indices]

    # Busqueda de la mejor calibracion
    best_alignment, index, Iou = find_best_calibration(obs_y, slices_y, CONFIG.W_RANGE, CONFIG.W_STEP)

    # Dispocición en vector de las longitudes de ondas calibradas para obs
    calibrado_x = np.full(len(best_alignment.index1), None)
    for i in range(len(best_alignment.index1)): # Calibrado
        calibrado_x[best_alignment.index1[i]] = slices_x[index][best_alignment.index2[i]]

    # Display the warping curve, i.e. the alignment curve
    best_alignment.plot(type="threeway")
    plt.savefig(os.path.join(CONFIG.SAVEPATH, f"{file}_threeway.png"))
    plt.close()

    # grafico de correlacion de puntos 
    best_alignment.plot(type="twoway",offset=-2)
    plt.savefig(os.path.join(CONFIG.SAVEPATH, f"{file}_twoway.png"))
    plt.close()

    # grafico del espectro calibrado en contraste con su lampara de compaonracion

    plt.figure(figsize=(10, 6), dpi=600) # Ajuste de tamaño de la figura

    min_teo_grap = calibrado_x[0] if calibrado_x[0] < obs_real_x[0] else obs_real_x[0] # Seccion del teorico
    min_teo_grap -= CONFIG.W_STEP
    max_teo_grap = calibrado_x[-1] if calibrado_x[-1] > obs_real_x[-1] else obs_real_x[-1]
    max_teo_grap += CONFIG.W_STEP
    grap_teo_x, grap_teo_y, _, _ = subconj_generator(teo_x, teo_y, min_teo_grap, max_teo_grap)
    plt.bar([], [], width=0, label='Teorico', color='blue', align='edge', alpha=1) 
    for x, y in zip(grap_teo_x, grap_teo_y):
        plt.bar([x], [-y], width=10, align='edge', color='blue', alpha=0.7)

    plt.bar([0], [0], width=0, label='Emp Real', color='black', align='edge', alpha=1) # Real
    for x, y in zip(obs_real_x, obs_y):
        plt.bar([x], [-y], width=2, align='edge', color='black', alpha=0.7)
    
    plt.bar([0], [0], width=0, label='Emp Optimo', color='red', align='edge', alpha=1) # Optimo
    for x, y in zip(calibrado_x, obs_y):
        plt.bar([x], [y], width=3, align='edge', color='red', alpha=1)
        
    plt.legend()

    plt.savefig(os.path.join(CONFIG.SAVEPATH, f"{file}_calibrado.png"))
    plt.close()

    # Guardar datos de ejecucion en CSV
    nueva_fila = { # Añadir la nueva fila al DataFrame
        'File': file, 
        'Cant_Ventanas_Probadas': len(slices_x), 
        'W_STEP': CONFIG.W_STEP,
        'W_RANGE': CONFIG.W_RANGE,
        'Distance_mejor_ventana': best_alignment.distance,
        'IoU_mejor_ventana': Iou,
        'Error_de_desplazamiento': EAM(calibrado_x, obs_real_x)
        }
    df = df._append(nueva_fila, ignore_index=True)
    # df = pd.concat([df, pd.DataFrame(list(nueva_fila))], ignore_index=True)

    df.to_csv(csv_path, index=False) # Guardar DataFrame actualizado