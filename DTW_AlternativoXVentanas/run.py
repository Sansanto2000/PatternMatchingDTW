import os
from dtw import *
import numpy as np
sys.path.append(os.path.dirname(os.path.abspath(os.path.dirname(__file__))))
from IOU import IoU
import matplotlib.pyplot as plt
from NIST_Table_Interactor import NIST_Table_Interactor
from utils import getfileData, normalize_min_max, slice_with_range_step

class Config:
    
    def __init__(self, file:str='WCOMPs/WCOMP01.fits', treshold:float=0.0, w_step:int=100, w_range:int=1900, 
                 savepath:str=os.path.dirname(os.path.abspath(__file__)), csv_name:str="DTW_AlternativoXVentanas.csv" ):
        """Funcion de inicializacion de la clase Config

        Args:
            file (str, optional): Nombre del archivo a calibrar. Defaults to 'WCOMP01.fits'.
            treshold (float, optional): Altura minima para la busqueda de picos. Defaults to 0.0.
            w_step (int, optional): Cantidad de longitudes de onda entre inicios de ventanas. Defaults to 100.
            w_range (int, optional): Rango de longitudes de onda que una ventana cubre. Defaults to 1900.
            savepath (str, optional): Especificación de carpeta para almacenar los graficos. Defaults to 
            os.path.dirname(os.path.abspath(__file__)).
            csv_name (str, optional): Nombre del archivo CSV. Defaults to "DTW_AlternativoXVentanas.csv".
        """
        self.FILE = file
        self.TRESHOLD = treshold
        self.W_STEP = w_step
        self.W_RANGE = w_range
        self.SAVEPATH = savepath
        self.CSV_NAME = csv_name
    
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

def get_Data_FILE(dirpath:str=os.path.dirname(os.path.abspath(__file__)), name:str='WCOMP01.fits', normalize:bool=True):
    """Funcion para obtener los datos de un archivo correspondiente a una lampara de comparación

    Args:
        dirpath (str, optional): Direccion de la carpeta contenedora del archivo. Defaults to 
        os.path.dirname(os.path.abspath(__file__)).
        name (str, optional): Nombre del archivo. Defaults to 'WCOMP01.fits'.
        normalize (bool, optional): Booleano para saber si los datos de respuesta deben estar normalizados o no. Defaults to True.

    Returns:
        numpy.ndarray: Datos de la lampara correspondientes al eje X
        numpy.ndarray: Datos de la lampara correspondientes al eje Y
        list: Headers adjuntos al archivo
    """
    # Datos y headers del observado
    filepath = os.path.join(dirpath, name)
    obs_data, obs_headers = getfileData(filepath=filepath)
    
    # Separacion de datos observados para el eje X y el eje Y
    obs_x = np.array(range(len(obs_data)))
    obs_y = obs_data
    
    # Normalizado de los datos obserbados en el eje Y
    if (normalize):
        obs_y, _, _ = normalize_min_max(obs_y)
    
    return obs_x, obs_y, obs_headers

def get_Data_NIST(dirpath:str=os.path.dirname(os.path.abspath(__file__)), name:str='Tabla(NIST)_Int_Long_Mat_Ref.csv', 
                  filter:list=["He I", "Ar I", "Ar II"], normalize:bool=True):
    """Funcion para obtener los datos teoricos a analizar

    Args:
        dirpath (str, optional): Direccion de la carpeta contenedora del archivo. Defaults to os.path.dirname(os.path.abspath(__file__)).
        name (str, optional): Nombre del archivo. Defaults to 'Tabla(NIST)_Int_Long_Mat_Ref.csv'.
        filter (list, optional): Elementos quimicos de los que se quieren los picos. Defaults to ["He I", "Ar I", "Ar II"].
        normalize (bool, optional): Booleano para saber si los datos de respuesta deben estar normalizados o no. Defaults to True.

    Returns:
        numpy.ndarray: Datos del teorico correspondientes al eje X
        numpy.ndarray: Datos del teorico correspondientes al eje Y
    """
    
    # Datos de teoricos del NIST
    filepath = os.path.join(dirpath, name)
    nisttr = NIST_Table_Interactor(csv_filename=filepath)

    # Obtencion del dataframe
    teorico_df = nisttr.get_dataframe(filter=filter)

    # Separacion de datos teoricos para el eje X y el eje Y
    teo_x = np.array(teorico_df['Wavelength(Ams)'])
    teo_y = np.array(teorico_df['Intensity'])
    
    # Normalizado de los datos en el eje Y
    if (normalize):
        teo_y, _, _ = normalize_min_max(target=teo_y)
    
    return teo_x, teo_y


CONFIG = Config(w_step=25)

# Directorio actual
act_dir = os.path.dirname(os.path.abspath(__file__))

# Obtencion de empirico
obs_x, obs_y, obs_headers = get_Data_FILE(dirpath=os.path.dirname(act_dir), name = CONFIG.FILE)
obs_real_x = obs_x * obs_headers['CD1_1'] + obs_headers['CRVAL1']

# Obtencion de Teorico
teo_x, teo_y = get_Data_NIST(dirpath=os.path.dirname(act_dir))

# Ventanado
ranges, slices_x, slices_y = slice_with_range_step(teo_x, teo_y, CONFIG.W_RANGE, CONFIG.W_STEP)
    
#Filtrar aquellos arreglos que tienen al menos un elemento
au_x = []
au_y = []
for i in range(len(slices_x)):
    if (len(slices_x[i]) > 0):
        au_x.append(slices_x[i])
        au_y.append(slices_y[i])
slices_x = np.array(au_x, dtype=object)
slices_y = np.array(au_y, dtype=object)

# Plantillas para acumulado de resultados
alignments = np.array([])
IoUs = np.array([])

for i in range(0, len(slices_x)):
        
    # Aplicación DTW del observado respecto al gaussianizado
    alignment = dtw(obs_y, slices_y[i], keep_internals=True, step_pattern=asymmetric, open_begin=True, open_end=True)
    #print(alignment.distance)
    
    # Determinación de la metrica IoU
    w_inicio = CONFIG.W_STEP*i
    w_fin = w_inicio + CONFIG.W_RANGE
    Iou = IoU(w_inicio, w_fin, obs_real_x[0], obs_real_x[-1])
    
    # Agrega datos a arreglo
    alignments = np.append(alignments, alignment)
    IoUs = np.append(IoUs, Iou)

# Busca el alineado con mejor metrica de distancia
index = np.argmin([alig.distance for alig in alignments])
best_alignment = alignments[index]

# Informa el IoU del mencionado alineamiento
print(IoUs[index])

# Display the warping curve, i.e. the alignment curve
best_alignment.plot(type="threeway")

# Correlacion
best_alignment.plot(type="twoway",offset=-2)

plt.show()
