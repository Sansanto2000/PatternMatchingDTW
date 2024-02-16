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
        print(x_au)
    else:
        x_au = np.linspace(np.min(x), np.max(x), resolution)
    y_au = np.zeros_like(x_au)
    for i in range(len(y)):
        mu = x[i]
        pdf = norm.pdf(x_au, mu, sigma)
        pdf, _, _ = normalize_min_max(pdf)
        y_au += (pdf * y[i])        
    return x_au, y_au

def dp(dist_mat):
    """
    Find minimum-cost path through matrix `dist_mat` using dynamic programming.

    The cost of a path is defined as the sum of the matrix entries on that
    path. See the following for details of the algorithm:

    - http://en.wikipedia.org/wiki/Dynamic_time_warping
    - https://www.ee.columbia.edu/~dpwe/resources/matlab/dtw/dp.m

    The notation in the first reference was followed, while Dan Ellis's code
    (second reference) was used to check for correctness. Returns a list of
    path indices and the cost matrix.
    """
    N, M = dist_mat.shape
    # Initialize the cost matrix
    cost_mat = np.zeros((N + 1, M + 1))
    for i in range(1, N + 1):
        cost_mat[i, 0] = np.inf
    for i in range(1, M + 1):
        cost_mat[0, i] = np.inf
    # Fill the cost matrix while keeping traceback information
    traceback_mat = np.zeros((N, M))
    for i in range(N):
        for j in range(M):
            penalty = [
                cost_mat[i, j],      # match (0)
                cost_mat[i, j + 1],  # insertion (1) * PENALTY
                cost_mat[i + 1, j]]  # deletion (2) * PENALTY /////////////////
            i_penalty = np.argmin(penalty)
            cost_mat[i + 1, j + 1] = dist_mat[i, j] + penalty[i_penalty]
            traceback_mat[i, j] = i_penalty
    # Traceback from bottom right
    i = N - 1
    j = M - 1
    path = [(i, j)]
    while i > 0 or j > 0:
        tb_type = traceback_mat[i, j]
        if tb_type == 0:
            # Match
            i = i - 1
            j = j - 1
        elif tb_type == 1:
            # Insertion
            i = i - 1
        elif tb_type == 2:
            # Deletion
            j = j - 1
        path.append((i, j))
    # Strip infinity edges from cost_mat before returning
    cost_mat = cost_mat[1:, 1:]
    return (path[::-1], cost_mat)

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