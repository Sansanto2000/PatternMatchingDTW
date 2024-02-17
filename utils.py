import os
from enum import Enum
from astropy.io import fits
from scipy.ndimage import gaussian_filter1d
from scipy.signal import savgol_filter
import numpy as np
from scipy.stats import norm

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
    return data

def normalize_min_max(target, min:int=None, max:int=None):
    """Dada un arreglo de datos objetivo se normalizan sus valores entre cero y uno

    Args:
        target (numpy.ndarray): Arreglo de datos a normalizar
        min (int, optional): Valor minimo a considerar como referencia para valor 0 post normalizado. Defaults to None.
        max (int, optional): Valor maximo a considerar como referencia para valor 1 post normalizado. Defaults to None.

    Returns:
        numpy.ndarray: Arreglo de datos normalizados entre 0 y 1
        int, optional: valor minimo empleado para la normalización
        int, optional: valor maximo empleado para la normalización
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
    if(rang):
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

def slice_with_range_step(arr_x, arr_y, W_RANGE, STEP):
    """Divide en varios subarreglos los datos correspondientes al eje X y el eje Y de una funcion

    Args:
        arr_x (numpy.ndarray): Arreglo con los datos de la funcion correspondientes al eje X
        arr_y (numpy.ndarray): Arreglo con los datos de la funcion correspondientes al eje Y
        W_RANGE (_type_): Rango de valores que una ventana cubre a partir de su inicio
        STEP (_type_): Cantidad de valores a considerar entre cada inicio de recorte

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
        inicio=inicio+STEP
        sub_arrs_x.append(arr_aux)
        sub_arrs_y.append(arr_auy)
        
    return ranges, sub_arrs_x, sub_arrs_y

class Processor:
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