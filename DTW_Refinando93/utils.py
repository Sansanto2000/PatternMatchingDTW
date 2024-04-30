import os
import re
import dtw
import math
import numpy as np
import pandas as pd
from astropy.io import fits
import matplotlib.pyplot as plt

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

class NIST_Table_Interactor:
    """Clase que centraliza la logica necesaria para procesar los datos csv del NIST
    """
    
    df = None
    
    def __init__(self, csv_filename):
        self.df = pd.read_csv(csv_filename)
        self._sanitize()
        
    def _sanitize(self):
        """Funcion para sanitizar los campos especificos del dataset y convertirlos
        en datos del tipo requerido.
        """
        self.df['Intensity'] = self.df['Intensity'].str.replace(r'[^0-9]+', '', regex=True)
        self.df['Intensity'] = pd.to_numeric(self.df['Intensity'])
        
        self.df['Wavelength(Ams)'] = self.df['Wavelength(Ams)'].apply(lambda x: re.sub(r'[^\d.]', '', x))
        self.df['Wavelength(Ams)'] = pd.to_numeric(self.df['Wavelength(Ams)'])
    
    def get_dataframe(self, cant:int = None, filter=None) -> pd.DataFrame:
        """Funcion para recuperar datos del conjunto de datos analizado

        Args:
            cant (int, optional): Cantidad de filas que se quieren recuperar. Defaults to None.
            filter (str, list, optional): Filtro de que tipo de materiales se quieren recuperar. Defaults to None.

        Returns:
            pd.DataFrame: DataFrame correspondiente.
        """
        df_aux = None
            
        if (cant):
            df_aux = self.df.head(cant) 
        else:
            df_aux = self.df
            
        if (filter):
            if (type(filter) == str):
                df_aux = df_aux[df_aux['Spectrum'] == filter]
            else:
                df_aux = df_aux[df_aux['Spectrum'].isin(filter)]
            
        return df_aux

def get_Data_NIST(csvpath:str, filter:list, normalize:bool=False):
    """Funcion para obtener las lineas de intensidad de ciertos materiales.

    Args:
        csvpath (str, optional): Direccion del csv del que se extraeran los datos.
        filter (list, optional): Elementos quimicos de los que se quieren los picos. Defaults to ["He I", "Ar I", "Ar II"].
        normalize (bool, optional): Booleano para saber si los datos de respuesta deben estar normalizados o no. Defaults to Falses.

    Returns:
        numpy.ndarray: Longitudes de onda para al eje X
        numpy.ndarray: Intensidades para el eje Y
    """
    
    # Datos de teoricos del NIST
    nisttr = NIST_Table_Interactor(csv_filename=csvpath)

    # Obtencion del dataframe
    lines_df = nisttr.get_dataframe(filter=filter)

    # Separacion de las lineas de intensidad para el eje X y el eje Y
    teo_x = np.array(lines_df['Wavelength(Ams)'])
    teo_y = np.array(lines_df['Intensity'])
    
    # Normalizado de los datos en el eje Y
    if (normalize):
        teo_y, _, _ = normalize_min_max(target=teo_y)
    
    return teo_x, teo_y

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
    
    return np.array(sub_x), np.array(sub_y)

def inspect_files_comparison(act_dir:str, files:dict):
    """Inspecciona secuencialmente distintos archivos de lampara y genera un grafico de 
    contraste con sus correspondientes lineas de intensidad teoricas.

    Args:
        act_dir (str): directorio para buscar datos y almacenar graficos.
        files (dict): conjunto de archivos a comparar.
    """
    
    for material_set in files.keys():

        # Determinar teorico que le corresponde
        csvpath=os.path.join(act_dir, "NIST", "Tabla(NIST)_Int_Long_Mat_Ref.csv")
        teo_x, teo_y = get_Data_NIST(csvpath=csvpath, filter=files[material_set]["materials"], 
                                     normalize=True)

        for file in files[material_set]["files"]:

            # Separar informacion
            emp_x, emp_y, emp_head = extract_lamp_info(file, normalize=True)

            # Determinar calibracion real
            try:
                emp_real_x = emp_x * emp_head['CD1_1'] + emp_head['CRVAL1']
            except Exception as e:
                print(f"Error archivo {file} < Falta de headers")
                continue

            # Determinar subconjunto teorico de interes
            grap_teo_x, grap_teo_y = subconj_generator(teo_x, teo_y, emp_real_x[0], emp_real_x[-1])

            # Graficar
            plt.figure(figsize=(24, 4), dpi=600)
            plt.bar(emp_real_x, emp_y, width=2, label='Real', color='red', align='edge', alpha=1)
            plt.bar(grap_teo_x, grap_teo_y, width=4, label='Teorical', color='blue', align='edge', alpha=0.9)
            plt.subplots_adjust(left=0.05, right=0.95, top=0.98, bottom=0.16)
            plt.xlabel('Wavelength (Å)', fontsize=20)
            plt.ylabel('Intensity', fontsize=20)
            plt.savefig(os.path.join(act_dir,"LampGraphs", material_set, f"{os.path.basename(file)}.svg"))
            plt.close()

def find_best_calibration(obs_y:np.ndarray, slices_y:np.ndarray, w_range:int, w_step:int):
    """Funcion para hallar las calibraciones correspondientes a todas los segmentos de una ventana.
    Devuelve la mejor calibracion encontrada.

    Args:
        obs_y (np.ndarray): Arreglo de datos a calibrar.
        slices_y (np.ndarray): Conjunto de arreglos entre los que se debe encontrar la mejor coincidencia.
        w_range (int, optional): Rango de longitudes de onda que una ventana cubre. Defaults to 1900.
        w_step (int, optional): Cantidad de longitudes de onda entre inicios de ventanas. Defaults to 100.

    Returns:
        DTW: Objeto DTW con detalles de la mejor coincidencia encontrada entre la funcion a comparar y las 
        posibles funciones objetivo.
        int: Indice de la ventana que se corresponde con la mejor calibracion.
    """

    # Plantillas para acumulado de resultados
    alignments = np.array([])
        
    for i in range(0, len(slices_y)):
        
        # Aplicación DTW del observado respecto al gaussianizado
        alignment = dtw.dtw(obs_y, 
                        slices_y[i], 
                        keep_internals=True, 
                        step_pattern=dtw.asymmetric, 
                        #distance_only=only_distance, # No hizo nada
                        #window_type="sakoechiba", # Resamplear
                        #window_args={'window_size':1000},
                        open_begin=True, 
                        open_end=True
                        )
        
        # Agrega datos a arreglo
        alignments = np.append(alignments, alignment)
    
    # Busca el alineado con mejor metrica de distancia
    best_index = np.argmin([alig.distance for alig in alignments])
    
    return alignments[best_index], best_index

def IoU(obs_min:float, obs_max:float, real_min:float, real_max:float) -> float:
    """Funcion para realizar el calculo de [Interseccion / Union] dados dos rangos de valores

    Args:
        obs_min (float): Valor minimo del rango 1
        obs_max (float): Valor maximo del rango 1
        real_min (float): Valor minimo del rango 2
        real_max (float): Valor maximo del rango 2

    Returns:
        float: Interseccion/Union
    """

    if (obs_min > real_max or real_min > obs_max):
        interseccion = 0
        # Obs:        |----| 
        # Rea: |----|
        #      or
        # Obs: |----| 
        # Rea:        |----|
        
    elif (real_min <= obs_min and real_max >= obs_max):
        interseccion = obs_max - obs_min 
        # Obs:   |----| 
        # Rea: |--------|
        
    elif (real_min >= obs_min and real_max <= obs_max):
        interseccion = real_max - real_min
        # Obs: |--------| 
        # Rea:   |----|
        
    elif (obs_min <= real_min and real_min <= obs_max and obs_max <= real_max):
        interseccion = obs_max - real_min
        # Obs: |----| 
        # Rea:    |----|
        
    elif (real_min <= obs_min and obs_min <= real_max and real_max <= obs_max):
        interseccion = real_max - obs_min
        # Obs:    |----| 
        # Rea: |----|

    union = real_max - real_min + obs_max - obs_min - interseccion
    iou =  interseccion / union
    return iou



