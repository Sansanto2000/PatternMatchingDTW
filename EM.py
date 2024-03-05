import math
import numpy as np
    
def ELM(arr1:np.ndarray, arr2:np.ndarray) -> float:
    """Funcion que en base a los datos de 2 funciones calcula el error lineal medio entre
    las 2 respecto al eje Y
    
    Args:
        arr1 (numpy.ndarray): Arreglo de datos correspondientes a la 1ra funcion a comparar (Eje Y)
        arr2 (numpy.ndarray): Arreglo de datos correspondientes a la 2da funcion a comparar (Eje Y)
    
    Returns:
        float: Metrica 'Error Lineal Medio'. Un ELM cercano a cero indica una mejor similitud 
        entre las series temporales, tambien sugiere que las dos secuencias están mejor alineadas. Un 
        valor lejano a 0 indica que las series temporales son distintas y que no se alinean de forma optima
    """
    
    # Variable para acumulaciona de diferencias
    sum:float = 0
    
    # Recorrer los arreglos y acumulado de diferencias lineales
    for y1, y2 in zip(arr1, arr2):
        sum += abs(y1 - y2)
    
    # Obtencion de la diferencia promedio entre los valores
    elm = sum / len(arr1)
    
    return elm

def ERM(arr1:np.ndarray, arr2:np.ndarray) -> float:
    """Funcion que en base a los datos de 2 funciones calcula el error racional medio entre
    las 2 respecto al eje Y. Esta funcion penaliza un poco menos diferencias fuertes de valores. 
    Pensada para aminorar el impacto de valores desviados.
    
    Args:
        arr1 (numpy.ndarray): Arreglo de datos correspondientes a la 1ra funcion a comparar (Eje Y)
        arr2 (numpy.ndarray): Arreglo de datos correspondientes a la 2da funcion a comparar (Eje Y)
    
    Returns:
        float: Metrica 'Error Racional Medio'. Un ERM cercano a cero indica una mejor similitud 
        entre las series temporales, tambien sugiere que las dos secuencias están mejor alineadas. Un 
        valor lejano a 0 indica que las series temporales son distintas y que no se alinean de forma optima
    """
    
    # Variable para acumulaciona de diferencias
    sum:float = 0
    
    # Recorrer los arreglos y acumulado de diferencias racionalizadas
    for y1, y2 in zip(arr1, arr2):
        sum += math.sqrt(abs(y1 - y2))
    
    # Obtencion de la diferencia promedio entre los valores
    erm = sum / len(arr1)
    
    return erm