import os
import math
from enum import Enum
from astropy.io import fits
from scipy.ndimage import gaussian_filter1d
from scipy.signal import savgol_filter
import numpy as np
from scipy.stats import norm
from NIST_Table_Interactor import NIST_Table_Interactor
import pandas as pd

def zero_padding(arr_x:np.ndarray, arr_y:np.ndarray, dist:int):
    """Funcion para rellenar de ceros un arreglo de datos Y conforme falten valores
    intermedios en un arrglo X ordenado. Retorna el arreglo de datos X con los nuevos
    valores intermedios y el arreglo de datos Y con las incerciones de 0 correspondientes

    Args:
        arr_x (np.ndarray): Arreglo de datos X ordenado, de este se revizara si hay 
        intervalos de valores que falten agregar.
        arr_y (np.ndarray): Arreglo d edatos Y. Acorde a los datos faltantes de X se 
        agregaran valores en el mismo.
        dist (int): Cantidad de unidades cada cual se deben agregar inserciones.

    Returns:
        numpy.ndarray: Arreglo de datos X actualizado.
        numpy.ndarray: Arreglo de datos Y actualizado.
    """

    # Encontrar las diferencias entre valores consecutivos
    diffs = np.diff(arr_x)

    # Encontrar las posiciones donde se deben insertar ceros
    zeros_pos = np.where(diffs > dist)[0]

    # Insertar ceros en las posiciones encontradas
    acumulate_count = 0
    for i, pos in enumerate(zeros_pos):
        # Calcular cuántos ceros se deben insertar en esta posición
        zero_count = math.floor((diffs[pos] - 1) / dist)

        # Calcula arreglos de datos a insertar en cada arreglo recibido
        news_x = np.arange(
            arr_x[pos+acumulate_count]+dist, 
            arr_x[pos+acumulate_count+1],
            dist, dtype=arr_x.dtype
            )
        
        news_y = np.zeros(len(news_x), dtype=arr_y.dtype)
        
        # Insertar nuevas longitudes de onda
        arr_x = np.insert(arr_x, pos+acumulate_count+1, news_x)

        # Insertar los ceros correspondientes
        arr_y = np.insert(arr_y, pos+acumulate_count+1, news_y)
        
        # Actualiza la cantidad acumulada para poder insertar conociendo los 
        # desplazamientos previos
        acumulate_count += len(news_x)
    
    return arr_x, arr_y

def get_Data_LIBS(dirpath:str=os.path.dirname(os.path.abspath(__file__)), name:str='LIBS_He_Ar_Ne_Resolution=1000.csv', 
                  normalize:bool=True):
    """Funcion para obtener los datos teoricos a analizar

    Args:
        dirpath (str, optional): Direccion de la carpeta contenedora del archivo. Defaults to os.path.dirname(os.path.abspath(__file__)).
        name (str, optional): Nombre del archivo. Defaults to 'Tabla(NIST)_Int_Long_Mat_Ref.csv'.
        normalize (bool, optional): Booleano para saber si los datos de respuesta deben estar normalizados o no. Defaults to True.

    Returns:
        numpy.ndarray: Datos del teorico correspondientes al eje X
        numpy.ndarray: Datos del teorico correspondientes al eje Y
    """
    
    # Obtencion del dataframe
    filepath = os.path.join(dirpath, name)
    libs_df = pd.read_csv(filepath)

    # Separacion de datos teoricos para el eje X y el eje Y
    teo_x = np.array(libs_df['Wavelength (Å)'])
    #libs_df['Sum'] = libs_df['Sum'].str.replace(',', '.')
    teo_y = np.array(pd.to_numeric(libs_df['Sum'], errors='coerce'))
    teo_y = np.nan_to_num(teo_y, nan=0)
    
    # Normalizado de los datos en el eje Y
    if (normalize):
        teo_y, _, _ = normalize_min_max(target=teo_y)
    return teo_x, teo_y

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
                  filter:list=["He I", "Ar I", "Ar II"]
                  , normalize:bool=True):
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

def getfileData(filepath:str): 
    """Funcion para obtener los datos de un archivo

    Args:
        filepath (str): Path del archivo

    Returns:
        numpy.ndarray: Datos del archivo
    """
    hdul = fits.open(filepath) 
    if('WOBJ' in filepath):
        headers = hdul[0].header
        data = hdul[0].data[0][0]
    else:
        headers = hdul[0].header
        data = hdul[0].data
    return data, headers

def normalize_min_max(target, min:float=None, max:float=None):
    """Dada un arreglo de datos objetivo se normalizan sus valores entre cero y uno

    Args:
        target (numpy.ndarray): Arreglo de datos a normalizar
        min (float, optional): Valor minimo a considerar como referencia para valor 0 post normalizado. Defaults to None.
        max (float, optional): Valor maximo a considerar como referencia para valor 1 post normalizado. Defaults to None.

    Returns:
        numpy.ndarray: Arreglo de datos normalizados entre 0 y 1
        float, optional: valor minimo empleado para la normalización
        float, optional: valor maximo empleado para la normalización
    """
    if(not min):
        min = np.min(target)
        
    if(not max):
        max = np.max(target)
        
    if(max==min):
        return target, min, max
    
    nor_target = (target - min) / (max - min)
    return nor_target, min, max

def gaussianice(x:list, y:list, resolution:int, sigma:float, rang=None): 
    """Funcion que gaussianiza conjuntos X, Y de datos

    Args:
        x (list): Arreglo de datos correspondientes al eje X
        y (list): Arreglo de datos correspondientes al eje Y
        resolution (int): Cantidad de elementos que tiene que demostrar la funcion gaussianizada resultante
        sigma (float): Sigma empleado
        rang (_type_, optional): rango de valores entre los que se quiere que se informe el gaussianizado. Defaults to None.

    Returns:
        numpy.ndarray: Arreglos con datos gaussianizados para el eje X
        numpy.ndarray: Arreglos con datos gaussianizados para el eje Y
    """
    if(rang != None):
        x_au = np.linspace(rang[0], rang[1], resolution)
    else:
        x_au = np.linspace(np.min(x), np.max(x), resolution)
    y_au = np.zeros_like(x_au)
    for i in range(len(y)):
        mu = x[i]
        pdf = norm.pdf(x_au, mu, sigma)
        pdf, _, _ = normalize_min_max(pdf)
        y_au += (pdf * y[i])        
    return x_au, y_au

def gaussianice_arr(xs:list, ys:list, resolution:int, sigma:float, ranges:list, normalize:bool=False): 
    """Funcion que gaussianiza arreglos de conjuntos conjuntos X, Y de datos

    Args:
        x (list): Arreglo de arreglos con datos correspondientes al eje X
        y (list): Arreglo de arreglos con datos correspondientes al eje Y
        resolution (int): Cantidad de elementos que tiene que demostrar la funcion gaussianizada resultante
        sigma (float): Sigma empleado
        rang (list): Arreglo de rango de valores entre los que se quiere que se informe el 
        gaussianizado.
        normalize (bool): Condicion booleana para saber si se deben normalizar los datos gauccianizados. Por
        defecto es False
        

    Returns:
        numpy.ndarray: Arreglo de arreglos con datos gaussianizados para el eje X
        numpy.ndarray: Arreglo de arreglos con datos gaussianizados para el eje Y
    """
    gaussianised_xs = []
    gaussianised_ys = []
    
    for i in range(len(xs)):
        
        aux_xs, aux_ys = gaussianice(xs[i], ys[i], resolution, sigma, ranges[i])
        
        if(normalize):
            aux_ys, _, _ = normalize_min_max(target=aux_ys)
            
        gaussianised_xs.append(aux_xs)
        gaussianised_ys.append(aux_ys)
        
    return gaussianised_xs, gaussianised_ys

def slice_with_range_step(arr_x, arr_y, W_RANGE, STEP, normalize:bool=False):
    """Divide en varios subarreglos los datos correspondientes al eje X y el eje Y de una funcion

    Args:
        arr_x (numpy.ndarray): Arreglo con los datos de la funcion correspondientes al eje X
        arr_y (numpy.ndarray): Arreglo con los datos de la funcion correspondientes al eje Y
        W_RANGE (_type_): Rango de valores que una ventana cubre a partir de su inicio
        STEP (_type_): Cantidad de valores a considerar entre cada inicio de recorte
        normalize (bool): Condicion boleana para indicar si las ventanas deben ser normalizadas o no

    Returns:
        numpy.ndarray: Arreglo de tuplas en el que cada tupla contiene los datos correspondientes al inicio y fin en el eje
        X de cada recorte
        numpy.ndarray: Arreglo con los datos de la funcion resultante correspondientes al eje X
        numpy.ndarray: Arreglo con los datos de la funcion resultante correspondientes al eje Y
    """
    ranges = []
    sub_arrs_x = []
    sub_arrs_y = []
    inicio = 0 #arr_x[0]
    
    if len(arr_x)==0:
        raise ValueError("No se puede rebanar un arreglo vacio")
    
    while inicio < (arr_x[len(arr_x)-1]):
        fin = inicio + W_RANGE
        arr_aux=[]
        arr_auy=[]
        i = 0
        
        while i<len(arr_x) and arr_x[i]<fin:
            if (arr_x[i]>=inicio):
                arr_aux.append(arr_x[i])
                arr_auy.append(arr_y[i])
            i+=1
        
        ranges.append((inicio, fin))
        inicio += STEP
        sub_arrs_x.append(arr_aux)
        if(normalize and len(arr_aux)>0): # En caso de que corresponda normaliza el la ventana calculada
            arr_auy, _, _ = normalize_min_max(np.array(arr_auy))
        sub_arrs_y.append(arr_auy)
        
    return ranges, sub_arrs_x, sub_arrs_y

def calibrate_with_observations(obs_data:np.ndarray, crval1:float, crpix1:float) -> list:
    """Función que calibra un conjunto de obserbaciones en longitud de onda segun los parametros
    de observación disponibles

    Args:
        obs_data (np.ndarray): Arreglo de datos con las intensidades*pixel obserbadas
        crval1 (float): Valor de longitud de onda inicial a considerar para la calibración
        crpix1 (float): Valor de longitud de onda por paso a considerar para la calibración

    Returns:
        list: Arreglo de datos con las intensidades*longitudDeOnda calibrados
    """
    
    # Calibracion de cada dato del arreglo observado respecto a los parametros recividos
    wav_arr = []
    for i in range(len(obs_data)):
        wav_arr.append(crval1+i*crpix1)
        
    return wav_arr

def subconj_generator(conj_x:np.ndarray, conj_y:np.ndarray, value_min:int, value_max:int):
    """Funcion que en base un subconjunto de datos correspondientes a una funcion genera
    un subconjunto de los mismos teniendo en cuenta determinados valores max y min que 
    puede tomar el eje X del subconjunto

    Args:
        conj_x (numpy.ndarray): Arreglo de datos del eje X
        conj_y (numpy.ndarray): arreglo de datos del eje Y
        value_min (int): Valor minimo que puede tener el subconjunto en el eje X
        value_max (int): Valor maximo que puede tener el subconjunto en el eje X

    Returns:
        numpy.ndarray: Arreglo de datos del subconjunto para el eje X
        numpy.ndarray: Arreglo de datos del subconjunto para el eje Y
        int: Valor minimo real que tiene el subconjunto en el eje X
        int: Valor maximo real que tiene el subconjunto en el eje X
    """
    # Determinación de subconjunto del teorico a usar como observado
    sub_x = []
    sub_y = []
    for i in range(len(conj_x)):
        if (value_min <= conj_x[i] and conj_x[i] <= value_max):
            sub_x.append(conj_x[i])
            sub_y.append(conj_y[i])
        elif (conj_x[i] > value_max):
            break
        
    # Redeterminacion de min y max del subconjunto con valores reales
    sub_min = conj_x[0]
    sub_max = conj_x[-1]
    
    return np.array(sub_x), np.array(sub_y), sub_min, sub_max
    

class Processor:
    """(Legacy Class)
    """
    class FuntionType(Enum):
        LINEAL = "LINEAL"
        GAUSSIAN = "GAUSSIAN"
        SG = "SG"
    def equalized_and_smooth(self, reference, target, function = FuntionType.LINEAL, sigma:float = None, window_length=10, polyorder=4):
        target_equalized_lineal = np.interp(np.linspace(0, 1, len(reference)), 
                                            np.linspace(0, 1, len(target)), target)
        if(function is self.FuntionType.LINEAL):
            return target_equalized_lineal
        elif(function is self.FuntionType.GAUSSIAN):
            target_equalized_gaussian = gaussian_filter1d(target_equalized_lineal, sigma=sigma)
            return target_equalized_gaussian
        elif(function is self.FuntionType.SG):
            target_equalized_SG = savgol_filter(target_equalized_lineal, 
                                                window_length=window_length, 
                                                polyorder=polyorder)
            return target_equalized_SG
        else:
            return None