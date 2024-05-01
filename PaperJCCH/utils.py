import math
import numpy as np

def zero_padding(arr_x:np.ndarray, arr_y:np.ndarray, dist:int):
    """Funcion para rellenar de ceros un arreglo de datos Y conforme falten valores
    intermedios en un arrglo X ordenado. Retorna el arreglo de datos X con los nuevos
    valores intermedios y el arreglo de datos Y con las incerciones de 0 correspondientes

    Args:
        arr_x (np.ndarray): Arreglo de datos X ordenado, de este se revizara si hay 
        intervalos de valores que falten agregar.
        arr_y (np.ndarray): Arreglo d edatos Y. Acorde a los datos faltantes de X se 
        agregaran valores en el mismo.
        dist (int): Cantidad de unidades cada cual se deben insertar datos.

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

def slice_with_range_step(arr_x, arr_y, LENGTH, STEP, normalize:bool=False):
    """Divide en varios subarreglos los datos correspondientes al eje X y el eje Y de una funcion

    Args:
        arr_x (numpy.ndarray): Arreglo con los datos de la funcion correspondientes al eje X
        arr_y (numpy.ndarray): Arreglo con los datos de la funcion correspondientes al eje Y
        LENGTH (_type_): Rango de valores que una ventana cubre a partir de su inicio
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
        fin = inicio + LENGTH
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

