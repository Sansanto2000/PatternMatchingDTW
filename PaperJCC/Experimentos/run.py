import os
import sys
import numpy as np
import pandas as pd
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
from utils import get_Data_LIBS, get_Data_NIST

class Config:
    
    def __init__(self, FILES_DIR:str=os.path.dirname(os.path.abspath(__file__)), FILES:list=['WCOMP01.fits'],
                 USE_NIST_DATA:bool=True, TEORICAL_DATA_FOLDER:str=os.path.dirname(os.path.abspath(__file__)),
                 TEORICAL_DATA_CSV_NAME:str="Tabla(NIST)_Int_Long_Mat_Ref.csv",
                 SAVE_DIR:str=os.path.dirname(os.path.abspath(__file__)), OUTPUT_CSV_NAME:str="output.csv", 
                 WINDOW_STEP:int=100, WINDOW_LENGTH:int=1900, NORMALIZE_WINDOWS:bool=False, 
                 DETECT_EMPIRICAL_PEAKS:bool=False, TRESHOLD_EMPIRICAL_PEAKS:float=0.0, HEIGHT_EMPIRICAL_PEAKS:bool=False):
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
            posterior a su segmentado. Defaults to False
            DETECT_EMPIRICAL_PEAKS (bool, optional): Condicion boleana par saber si se recortara el empirico segun los picos
            disponibles o no. Defaults to False.
            TRESHOLD_EMPIRICAL_PEAKS (float, optional): Diferencia minima con picos vecinos para la busqueda de picos en el 
            empirico. Defaults to 0.0.
            HEIGHT_EMPIRICAL_PEAKS (bool, optional): Altura minima para la busqueda de picos en el empirico. Defaults to 0.0.
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
        self.DETECT_EMPIRICAL_PEAKS = DETECT_EMPIRICAL_PEAKS
        self.TRESHOLD_EMPIRICAL_PEAKS = TRESHOLD_EMPIRICAL_PEAKS
        self.HEIGHT_EMPIRICAL_PEAKS = HEIGHT_EMPIRICAL_PEAKS

act_dir = os.path.dirname(os.path.abspath(__file__)) # Directorio actual
files = ["WCOMP01.fits", "WCOMP02.fits", "WCOMP03.fits", "WCOMP04.fits", "WCOMP05.fits", # Archivos a procesar
         "WCOMP06.fits", "WCOMP07.fits", "WCOMP08.fits", "WCOMP09.fits", "WCOMP10.fits",
         "WCOMP11.fits", "WCOMP12.fits", "WCOMP13.fits", "WCOMP14.fits", "WCOMP15.fits",
         "WCOMP16.fits", "WCOMP17.fits", "WCOMP18.fits", "WCOMP19.fits", "WCOMP20_A.fits",
         "WCOMP20_A.fits", "WCOMP21.fits", "WCOMP22.fits", "WCOMP23.fits", "WCOMP24.fits", 
         "WCOMP25.fits", "WCOMP26.fits", "WCOMP27.fits", "WCOMP28.fits", "WCOMP29.fits", 
         "WCOMP30.fits", "WCOMP31.fits"]

CONFIG = Config(    # Constantes de configuracion
    FILES_DIR=os.path.join(os.path.dirname(os.path.dirname(act_dir)), 'WCOMPs'),
    FILES=files,
    SAVE_DIR=os.path.join(act_dir, 'output'),
    OUTPUT_CSV_NAME="output.csv",
    
    TEORICAL_DATA_FOLDER=os.path.dirname(os.path.dirname(act_dir)),
    TEORICAL_DATA_CSV_NAME="Tabla(NIST)_Int_Long_Mat_Ref.csv",
    #TEORICAL_DATA_CSV_NAME="LIBS_He_Ar_Ne_Resolution=100.csv",
    #TEORICAL_DATA_CSV_NAME="LIBS_He_Ar_Ne_Resolution=260.csv",
    USE_NIST_DATA=True,
    
    WINDOW_STEP=25,
    WINDOW_LENGTH=2000,
    NORMALIZE_WINDOWS=False,
    DETECT_EMPIRICAL_PEAKS = False,
    TRESHOLD_EMPIRICAL_PEAKS=0.0,
    HEIGHT_EMPIRICAL_PEAKS=0.0
)

# Verifica existencia de la carpeta de guardado de resultados, si no existe la crea
if not os.path.exists(CONFIG.SAVE_DIR):
    os.makedirs(CONFIG.SAVE_DIR)

# Preparar CSV para persistencia de los datos
df = pd.DataFrame(columns=['Source of reference', 
                           'peaks in data', 'Treshold', 'Height'
                           'Normalized', 
                           'W_STEP', 'W_RANGE', 'Count of Windows', 
                           '(AVG) Distance', 
                           '(AVG) IoU_mejor_ventana', 
                           '(AVG) Error_de_desplazamiento'])
csv_path = os.path.join(CONFIG.SAVE_DIR, CONFIG.OUTPUT_CSV_NAME)
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
if (CONFIG.DETECT_EMPIRICAL_PEAKS):
    indices, _ = find_peaks(teo_y, 
                            threshold=[CONFIG.TRESHOLD_EMPIRICAL_PEAKS, np.inf],
                            height= [CONFIG.HEIGHT_EMPIRICAL_PEAKS, np.inf]
                            )
    teo_x = teo_x[indices]
    teo_y = teo_y[indices]

print ('----------------')
print (f'LEN(teo_x)={len(teo_x)}\
    [{teo_x[0]}, ...,{teo_x[-1]}]')
print (f'LEN(teo_y)={len(teo_y)} \
    [{teo_y[0]}, ...,{teo_y[-1]}]')
print ('----------------')