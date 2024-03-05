import math
import numpy as np
    
def EAM(arr1:np.ndarray, arr2:np.ndarray) -> float:
    """Funcion que en base a los datos de 2 funciones calcula el error lineal medio entre
    las 2 respecto al eje Y
    
    Args:
        arr1 (numpy.ndarray): Arreglo de datos correspondientes a la 1ra funcion a comparar (Eje Y)
        arr2 (numpy.ndarray): Arreglo de datos correspondientes a la 2da funcion a comparar (Eje Y)
    
    Returns:
        float: Metrica 'Error Absoluto Medio'. Un EAM cercano a cero indica una mejor similitud 
        entre las series temporales, tambien sugiere que las dos secuencias están mejor alineadas. Un 
        valor lejano a 0 indica que las series temporales son distintas y que no se alinean de forma optima
    """
    
    # Calculo de diferencias entre las pocisiones de los arreglos y promediado
    EaM = np.mean(np.abs(arr1 - arr2))
    
    return EaM

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
    
    # Calculo de diferencias racionalizadas entre las pocisiones de los arreglos y promediado
    ErM = np.mean(np.sqrt(np.abs(arr1 - arr2)))
    
    return ErM