import numpy as np
from astropy.io import fits

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

def extract_lamp_info(filepath:str, normalize:bool=False):
    """Funcion para obtener los datos de un archivo correspondiente a una lampara de comparación

    Args:
        filepath (str, optional): Direccion del archivo.
        normalize (bool, optional): Booleano para saber si los datos de respuesta deben estar 
        normalizados o no. Defaults to False.

    Returns:
        numpy.ndarray: Datos de la lampara correspondientes al eje X
        numpy.ndarray: Datos de la lampara correspondientes al eje Y
        list: Headers adjuntos al archivo
    """
    # Extraer datos y headers del archivo
    hdul = fits.open(filepath) 
    if('WOBJ' in filepath): # Espectro calibrado
        headers = hdul[0].header
        data = hdul[0].data[0][0]
    else:
        headers = hdul[0].header
        data = hdul[0].data
    
    # Separ datos en X e Y
    obs_x = np.array(range(len(data)))
    obs_y = data
    
    # Normalizado de los datos obserbados en el eje Y
    if (normalize):
        obs_y, _, _ = normalize_min_max(obs_y)
    
    return obs_x, obs_y, headers