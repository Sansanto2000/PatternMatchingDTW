"""Script para generar experimentos para Paper a presentar en la JCC 2024 """
import os
import sys
import dtw
import math
import numpy as np
import pandas as pd
from tqdm import tqdm
import matplotlib.pyplot as plt
from scipy.signal import find_peaks
sys.path.append(
    os.path.dirname( # Root (folder)
        os.path.dirname( # PaperJCC (folder)
            os.path.abspath( # Experimentos (folder | absolute path)
                os.path.dirname(__file__) # Experimentos (folder | relative path)
                ) 
            )
        )
    )
from utils import get_Data_LIBS, get_Data_NIST, get_Data_FILE, zero_padding, slice_with_range_step, subconj_generator
from IOU import IoU
from EM import EAM

class Config:
    
    def __init__(self, FILES_DIR:str=os.path.dirname(os.path.abspath(__file__)), FILES:list=['WCOMP01.fits'],
                 USE_NIST_DATA:bool=True, TEORICAL_DATA_FOLDER:str=os.path.dirname(os.path.abspath(__file__)),
                 TEORICAL_DATA_CSV_NAME:str="Tabla(NIST)_Int_Long_Mat_Ref.csv",
                 SAVE_DIR:str=os.path.dirname(os.path.abspath(__file__)), OUTPUT_CSV_NAME:str="output.csv", 
                 WINDOW_STEP:int=100, WINDOW_LENGTH:int=1900, NORMALIZE_WINDOWS:bool=False, 
                 ZERO_PADDING:bool=False, DETECT_TEORICAL_PEAKS:bool=False, 
                 TRESHOLD_TEORICAL_PEAKS:float=0.0, HEIGHT_TEORICAL_PEAKS:bool=False,
                 DETECT_EMPIRICAL_PEAKS:bool=False, TRESHOLD_EMPIRICAL_PEAKS:float=0.0, 
                 HEIGHT_EMPIRICAL_PEAKS:bool=False, GRAPH_BESTS:bool=False):
        """Funcion de inicializacion de la clase Config

        Args:
            FILES_DIR (str, optional): Especificación de carpeta donde buscar los archivos. Defaults to 
            os.path.dirname(os.path.abspath(__file__)).
            FILES (list, optional): Arreglo con los nombres de los archivos a calibrar. Defaults to ['WCOMP01.fits'].
            SAVE_DIR (str, optional): Especificación de carpeta para almacenar los graficos. Defaults to 
            os.path.dirname(os.path.abspath(__file__)).
            OUTPUT_CSV_NAME (str, optional): Nombre del archivo CSV. Defaults to "DTW_AlternativoXVentanas.csv".
            TEORICAL_DATA_CSV_NAME (str, optional): Nombre del archivo de donde se extraeran los datos teoricos.
            USE_NIST_DATA (bool, optional): Condicion boleana par saber si se usaran los datos del NIST o de otra 
            fuente alternativa.
            WINDOW_STEP (int, optional): Cantidad de longitudes de onda entre inicios de ventanas. Defaults to 100.
            WINDOW_LENGTH (int, optional): Rango de longitudes de onda que una ventana cubre. Defaults to 1900.
            NORMALIZE_WINDOWS (bool, optional): Condicion boleana par saber si se normalizaran las ventanas teoricas 
            posterior a su segmentado. Defaults to False.
            DETECT_TEORICAL_PEAKS (bool, optional): Condicion boleana par saber si se recortara el teorico segun los picos
            disponibles o no. Defaults to False.
            TRESHOLD_TEORICAL_PEAKS (float, optional): Diferencia minima con picos vecinos para la busqueda de picos en el 
            teorico. Defaults to 0.0.
            HEIGHT_TEORICAL_PEAKS (bool, optional): Altura minima para la busqueda de picos en el teorico. Defaults to 0.0.
            DETECT_EMPIRICAL_PEAKS (bool, optional): Condicion boleana par saber si se recortara el empirico segun los picos
            disponibles o no. Defaults to False.
            TRESHOLD_EMPIRICAL_PEAKS (float, optional): Diferencia minima con picos vecinos para la busqueda de picos en el 
            empirico. Defaults to 0.0.
            HEIGHT_EMPIRICAL_PEAKS (bool, optional): Altura minima para la busqueda de picos en el empirico. Defaults to 0.0.
            GRAPH_BESTS (bool, optional): Condicion boleana par saber si se generaran y guardaran graficos de referencia de 
            las mejores calibraciones o no. Defaults to False.
            """
        self.FILES_DIR = FILES_DIR
        self.FILES = FILES
        self.WINDOW_STEP = WINDOW_STEP
        self.OUTPUT_CSV_NAME = OUTPUT_CSV_NAME
        self.USE_NIST_DATA = USE_NIST_DATA
        self.TEORICAL_DATA_CSV_NAME = TEORICAL_DATA_CSV_NAME
        self.TEORICAL_DATA_FOLDER = TEORICAL_DATA_FOLDER
        self.WINDOW_LENGTH = WINDOW_LENGTH
        self.SAVE_DIR = SAVE_DIR
        self.NORMALIZE_WINDOWS = NORMALIZE_WINDOWS
        self.ZERO_PADDING = ZERO_PADDING
        self.DETECT_TEORICAL_PEAKS = DETECT_TEORICAL_PEAKS
        self.TRESHOLD_TEORICAL_PEAKS = TRESHOLD_TEORICAL_PEAKS
        self.HEIGHT_TEORICAL_PEAKS = HEIGHT_TEORICAL_PEAKS
        self.DETECT_EMPIRICAL_PEAKS = DETECT_EMPIRICAL_PEAKS
        self.TRESHOLD_EMPIRICAL_PEAKS = TRESHOLD_EMPIRICAL_PEAKS
        self.HEIGHT_EMPIRICAL_PEAKS = HEIGHT_EMPIRICAL_PEAKS
        self.GRAPH_BESTS = GRAPH_BESTS

def find_best_calibration(obs_y:np.ndarray, slices_y:np.ndarray, w_range:int, w_step:int):
    """Funcion para hallar las calibraciones correspondientes a todas los segmentos de una ventana.
    Devuelve la mejor calibracion encontrada.

    Args:
        obs_y (np.ndarray): Arreglo de datos a calibrar
        slices_y (np.ndarray): Conjunto de arreglos entre los que se debe encontrar la mejor coincidencia
        w_range (int, optional): Rango de longitudes de onda que una ventana cubre. Defaults to 1900
        w_step (int, optional): Cantidad de longitudes de onda entre inicios de ventanas. Defaults to 100

    Returns:
        DTW: Objeto DTW con detalles de la mejor coincidencia encontrada entre la funcion a comparar y las 
        posibles funciones objetivo
        int: Indice de la ventana que se corresponde con la mejor calibracion
    """

    # Plantillas para acumulado de resultados
    alignments = np.array([])
        
    for i in range(0, len(slices_y)):
        
        # Aplicación DTW del observado respecto al gaussianizado
        alignment = dtw.dtw(obs_y, 
                        slices_y[i], 
                        keep_internals=True, 
                        step_pattern=dtw.asymmetric, 
                        #window_type="sakoechiba", 
                        #window_args={'window_size':1000},
                        open_begin=True, 
                        open_end=True
                        )
        
        # Agrega datos a arreglo
        alignments = np.append(alignments, alignment)
    
    # Busca el alineado con mejor metrica de distancia
    best_index = np.argmin([alig.distance for alig in alignments])
    
    return alignments[best_index], best_index

def run_calibrations(CONFIG:Config):
    """Funcion que realiza un conjunto de calibraciones completo sobre un conjunto de archivos.
    Almacena datos estadisticos varios de la ejecucion en un archivo CSV. Tambien guardar puede un 
    grafico de muestra para comprobacion del usuario.

    Args:
        CONFIG (Config): Archivo con datos de configuracion que detalla como se deben hacer las
        calibraciones del grupo de archivos a procesar.
    """
    # Verifica existencia de la carpeta de guardado de resultados, si no existe la crea
    if not os.path.exists(CONFIG.SAVE_DIR):
        os.makedirs(CONFIG.SAVE_DIR)

    # Preparar CSV para persistencia de los datos
    csv_path = os.path.join(CONFIG.SAVE_DIR, CONFIG.OUTPUT_CSV_NAME)
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
            'W_STEP', 
            'W_LENGTH', 
            'Count of Windows', 
            '(AVG) Distance', 
            '(AVG) IoU', 
            '(AVG) Scroll Error'
            ])
        df.to_csv(csv_path, index=False)


    # Obtencion de datos de referencia (Teorico) de donde corresponda
    if (CONFIG.USE_NIST_DATA):
        filter:list=["He I", "Ar I", "Ar II"]
        teo_x, teo_y = get_Data_NIST(dirpath=CONFIG.TEORICAL_DATA_FOLDER, 
                                    name=CONFIG.TEORICAL_DATA_CSV_NAME, 
                                    filter=filter)
    else:
        teo_x, teo_y = get_Data_LIBS(dirpath=CONFIG.TEORICAL_DATA_FOLDER, 
                                    name=CONFIG.TEORICAL_DATA_CSV_NAME)

    # Aislado de picos (si corresponde)
    if (CONFIG.DETECT_TEORICAL_PEAKS):
        indices, _ = find_peaks(teo_y, 
                                threshold=[CONFIG.TRESHOLD_TEORICAL_PEAKS, np.inf],
                                height= [CONFIG.HEIGHT_TEORICAL_PEAKS, np.inf]
                                )
        teo_x = teo_x[indices]
        teo_y = teo_y[indices]

    # Rellenado de ceros (si corresponde)
    if (CONFIG.ZERO_PADDING):
        teo_x, teo_y = zero_padding(arr_x=teo_x, arr_y=teo_y, dist=10)

    # Ventanado de datos de referencia(teorico)
    ranges, slices_x, slices_y = slice_with_range_step(teo_x, teo_y, 
                                                    CONFIG.WINDOW_LENGTH, 
                                                    CONFIG.WINDOW_STEP, 
                                                    CONFIG.NORMALIZE_WINDOWS)

    #Filtrar aquellos arreglos que no tienen elementos
    au_x = []
    au_y = []
    for i in range(len(slices_x)):
        if (len(slices_x[i]) > 0):
            au_x.append(slices_x[i])
            au_y.append(slices_y[i])
    slices_x = np.array(au_x, dtype=object)
    slices_y = np.array(au_y, dtype=object)

    # --------PROCESADO DE ARCHIVOS-----------
    Distances = np.array([])
    IoUs = np.array([])
    EAMs = np.array([])
    for file in tqdm(CONFIG.FILES, desc=f'Porcentaje de avance'):

        # Obtencion de empirico
        obs_x, obs_y, obs_headers = get_Data_FILE(dirpath=CONFIG.FILES_DIR, name=file)
        obs_real_x = obs_x * obs_headers['CD1_1'] + obs_headers['CRVAL1']
        
        # Aislado de picos (si corresponde)
        if (CONFIG.DETECT_EMPIRICAL_PEAKS):
            indices, _ = find_peaks(obs_y, 
                                    threshold=[CONFIG.TRESHOLD_EMPIRICAL_PEAKS, np.inf],
                                    height= [CONFIG.HEIGHT_EMPIRICAL_PEAKS, np.inf]
                                    )
            obs_x = obs_x[indices]
            obs_y = obs_y[indices]
        
        # Busqueda de la mejor calibracion y obtencion de metricas
        best_alignment, index = find_best_calibration(obs_y, slices_y, 
                                                           CONFIG.WINDOW_LENGTH, CONFIG.WINDOW_STEP)
        
        # Determinación de la metrica IoU
        c_inicio = slices_x[index][best_alignment.index2[0]] # inicio calibrado
        c_fin = slices_x[index][best_alignment.index2[-1]] # fin calibrado
        Iou = IoU(c_inicio, c_fin, obs_real_x[0], obs_real_x[-1]) # Segun mejor calibrado
        
        # Dispocición en vector de las longitudes de ondas calibradas para obs
        calibrado_x = np.full(len(best_alignment.index1), None)
        for i in range(len(best_alignment.index1)): # Calibrado
            calibrado_x[best_alignment.index1[i]] = slices_x[index][best_alignment.index2[i]]
            
        # Agregado de metricas en arreglos de almacenamiento
        Distances = np.append(Distances, best_alignment.distance)
        IoUs = np.append(IoUs, Iou)
        EAMs = np.append(EAMs, EAM(calibrado_x, obs_real_x))

    # Genera y almacena grafico de ultimo archivo en caso de que corresponda
    if (CONFIG.GRAPH_BESTS):
        plt.figure(figsize=(10, 6), dpi=800) # Ajuste de tamaño de la figura

        min_teo_grap = calibrado_x[0] if calibrado_x[0] < obs_real_x[0] else obs_real_x[0] # Seccion del teorico
        min_teo_grap -= CONFIG.WINDOW_STEP
        max_teo_grap = calibrado_x[-1] if calibrado_x[-1] > obs_real_x[-1] else obs_real_x[-1]
        max_teo_grap += CONFIG.WINDOW_STEP
        grap_teo_x, grap_teo_y, _, _ = subconj_generator(teo_x, teo_y, min_teo_grap, max_teo_grap)
        plt.bar(grap_teo_x, -grap_teo_y, width=10, label='Teorico', color='blue', align='edge', alpha=0.7) 

        plt.bar(obs_real_x, obs_y, width=2, label='Emp Real', color='black', align='edge', alpha=0.7) # Real
        
        plt.bar(calibrado_x, obs_y, width=3, label='Emp Optimo', color='red', align='edge', alpha=1) # Optimo
            
        plt.legend()
        fig_name = "NIST" if(CONFIG.USE_NIST_DATA) else "LIBS"
        if (fig_name=="LIBS"):
            fig_name += "100" if("100" in CONFIG.TEORICAL_DATA_CSV_NAME) else "260"
        fig_name += "_EP" if(CONFIG.DETECT_EMPIRICAL_PEAKS) else ""
        fig_name += "_TP" if(CONFIG.DETECT_TEORICAL_PEAKS) else ""
        fig_name += "_NOR" if(CONFIG.NORMALIZE_WINDOWS) else ""
        fig_name += "_ZP" if(CONFIG.ZERO_PADDING) else ""
        plt.savefig(os.path.join(CONFIG.SAVE_DIR, f"{fig_name}.svg"))
        plt.close()
            
    # Guardar datos promedios de ejecucion en CSV
    nueva_fila = { # Añadir la nueva fila al DataFrame
        'Source of reference':'NIST' if (CONFIG.USE_NIST_DATA) else 'LIBS', 
        'peaks in teorical data':CONFIG.DETECT_TEORICAL_PEAKS, 
        'Treshold teorical data':CONFIG.TRESHOLD_TEORICAL_PEAKS if (CONFIG.DETECT_TEORICAL_PEAKS) else '-', 
        'Height teorical data':CONFIG.HEIGHT_TEORICAL_PEAKS if (CONFIG.DETECT_TEORICAL_PEAKS) else '-',
        'peaks in empirical data':CONFIG.DETECT_EMPIRICAL_PEAKS, 
        'Treshold empirical data':CONFIG.TRESHOLD_EMPIRICAL_PEAKS if (CONFIG.DETECT_EMPIRICAL_PEAKS) else '-', 
        'Height empirical data':CONFIG.HEIGHT_EMPIRICAL_PEAKS if (CONFIG.DETECT_EMPIRICAL_PEAKS) else '-',
        'Normalized':CONFIG.NORMALIZE_WINDOWS, 
        'Zero Padding':CONFIG.ZERO_PADDING,
        'W_STEP':CONFIG.WINDOW_STEP, 
        'W_LENGTH':CONFIG.WINDOW_LENGTH, 
        'Count of Windows':len(slices_x), 
        '(AVG) Distance':np.mean(Distances), 
        '(AVG) IoU':np.mean(IoUs), 
        '(AVG) Scroll Error':np.mean(EAMs)
        }
    df = df._append(nueva_fila, ignore_index=True)

    df.to_csv(csv_path, index=False) # Guardar DataFrame actualizado

act_dir = os.path.dirname(os.path.abspath(__file__)) # Directorio actual
files = ["WCOMP01.fits", "WCOMP02.fits", "WCOMP03.fits", "WCOMP04.fits", "WCOMP05.fits", # Archivos a procesar
         "WCOMP06.fits", "WCOMP07.fits", "WCOMP08.fits", "WCOMP09.fits", "WCOMP10.fits",
         "WCOMP11.fits", "WCOMP12.fits", "WCOMP13.fits", "WCOMP14.fits", "WCOMP15.fits",
         "WCOMP16.fits", "WCOMP17.fits", "WCOMP18.fits", "WCOMP19.fits", "WCOMP20_A.fits",
         "WCOMP20_A.fits", "WCOMP21.fits", "WCOMP22.fits", "WCOMP23.fits", "WCOMP24.fits", 
         "WCOMP25.fits", "WCOMP26.fits", "WCOMP27.fits", "WCOMP28.fits", "WCOMP29.fits", 
         "WCOMP30.fits", "WCOMP31.fits"]

total_iteraciones = 3*2*2*2*2
iteracion_actual = 1

# Ejecutar calibraciones para todas las combinaciones de interes
for teorical_data_csv_name in ["Tabla(NIST)_Int_Long_Mat_Ref.csv", 
                      "LIBS_He_Ar_Ne_Resolution=100.csv", 
                      "LIBS_He_Ar_Ne_Resolution=260.csv"]:
    for detect_empirical_peaks in [True, False]:
        for detect_teorical_peaks in [True, False]:
            for normalize_windows in [True, False]:
                for zero_padding_bool in [True, False]:
                    print(f"{iteracion_actual}/{total_iteraciones}")
                    config = Config(
                        FILES_DIR=os.path.join(os.path.dirname(os.path.dirname(act_dir)), 'WCOMPs'),
                        FILES=files,
                        SAVE_DIR=os.path.join(act_dir, 'output'),
                        OUTPUT_CSV_NAME="output.csv",
                        
                        TEORICAL_DATA_FOLDER=os.path.dirname(os.path.dirname(act_dir)),
                        TEORICAL_DATA_CSV_NAME=teorical_data_csv_name,
                        USE_NIST_DATA= True if("NIST" in teorical_data_csv_name) else False,
                        
                        WINDOW_STEP=25,
                        WINDOW_LENGTH=2000,
                        NORMALIZE_WINDOWS=normalize_windows,
                        ZERO_PADDING=zero_padding_bool,
                        DETECT_TEORICAL_PEAKS=detect_teorical_peaks,
                        TRESHOLD_TEORICAL_PEAKS=0.0,
                        HEIGHT_TEORICAL_PEAKS=0.0,
                        DETECT_EMPIRICAL_PEAKS=detect_empirical_peaks,
                        TRESHOLD_EMPIRICAL_PEAKS=0.0,
                        HEIGHT_EMPIRICAL_PEAKS=0.0,
                        
                        GRAPH_BESTS=True
                    )
                    run_calibrations(CONFIG=config)
                    iteracion_actual += 1