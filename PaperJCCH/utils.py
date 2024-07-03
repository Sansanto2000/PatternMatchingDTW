import os
import dtw
import math
import time
import numpy as np
import pandas as pd
from config import Config
from astropy.io import fits
import matplotlib.pyplot as plt
from scipy.signal import find_peaks

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

def get_Data_NIST(dirpath:str, name:str, filter:list=None, normalize:bool=True):
    """Funcion para obtener los datos teoricos a analizar desde los datos del NIST.

    Args:
        dirpath (str, optional): direccion de la carpeta contenedora del archivo.
        name (str, optional): nombre del archivo.
        filter (list, optional): elementos quimicos de los que se quieren los picos. Defaults to None.
        normalize (bool, optional): booleano para saber si los datos de respuesta deben estar normalizados o no. Defaults to True.

    Returns:
        numpy.ndarray: Datos del teorico correspondientes al eje X
        numpy.ndarray: Datos del teorico correspondientes al eje Y
    """
    
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
        
        def get_dataframe(self, cant:int = None, filter:str=None) -> pd.DataFrame:
            """Funcion para recuperar datos del conjunto de datos analizado

            Args:
                cant (int, optional): Cantidad de filas que se quieren recuperar. Defaults to None.
                filter (str, optional): Filtro de que tipo de materiales se quieren recuperar. Defaults to None.

            Returns:
                pd.DataFrame: DataFrame correspondiente.
            """
            df = None
                
            if (cant):
                df = self.df.head(cant) 
            else:
                df = self.df
                
            if (filter):
                if (type(filter) == str):
                    df = df[df['Spectrum'] == filter]
                else:
                    df = df[df['Spectrum'].isin(filter)]
                
            return df
        
    # Datos de teoricos del NIST
    filepath = os.path.join(dirpath, name)
    nisttr = NIST_Table_Interactor(csv_filename=filepath)

    # Obtencion del dataframe
    if(filter is None):
        teorico_df = nisttr.get_dataframe()
    else:
        teorico_df = nisttr.get_dataframe(filter=filter)

    # Separacion de datos teoricos para el eje X y el eje Y
    teo_x = np.array(teorico_df['Wavelength(Ams)'])
    teo_y = np.array(teorico_df['Intensity'])
    
    # Normalizado de los datos en el eje Y
    if (normalize):
        teo_y, _, _ = normalize_min_max(target=teo_y)
    
    return teo_x, teo_y

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
    best_distance = np.inf
    best_aligment = None
    best_index = None
        
    for i in range(0, len(slices_y)):
        
        # Aplicación DTW del observado respecto al gaussianizado
        alignment = dtw.dtw(obs_y, 
                        slices_y[i], 
                        keep_internals=False, 
                        step_pattern=dtw.asymmetric, 
                        #distance_only=only_distance, # No hizo nada
                        #window_type="sakoechiba", # Resamplear
                        #window_args={'window_size':1000},
                        open_begin=True, 
                        open_end=True
                        )
        
        if (alignment.distance < best_distance):
            best_distance = alignment.distance
            best_aligment = alignment
            best_index = i
    
    return best_aligment, best_index

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

def EAM(arr1:np.ndarray, arr2:np.ndarray) -> float:
    """Funcion que en base a los datos de 2 arreglos calcula el error lineal medio entre
    las 2
    
    Args:
        arr1 (numpy.ndarray): Arreglo de datos correspondientes a la 1ra funcion a comparar
        arr2 (numpy.ndarray): Arreglo de datos correspondientes a la 2da funcion a comparar
    
    Returns:
        float: Metrica 'Error Absoluto Medio'. Un EAM cercano a cero indica una mejor similitud 
        entre las series temporales, tambien sugiere que las dos secuencias están mejor alineadas. Un 
        valor lejano a 0 indica que las series temporales son distintas y que no se alinean de forma optima
    """
    
    # Calculo de diferencias entre las pocisiones de los arreglos y promediado
    EaM = np.mean(np.abs(arr1 - arr2))
    
    return EaM

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

def run_calibrations(teo_x:np.ndarray, teo_y:np.ndarray, files:np.ndarray, window_length:int, 
                     window_step:int, detect_teorical_peaks:bool, detect_empirical_peaks:bool, 
                     zero_padding_bool:bool, normalize_windows:bool, save_dir:str, graph:bool, 
                     output_csv_path:str, teorical_weed_out:bool, empirical_weed_out:bool, 
                     minimal_data_for_weed_out:int):
    """Funcion que realiza un conjunto de calibraciones completo sobre un conjunto de archivos.
    Almacena datos estadisticos varios de la ejecucion en un archivo CSV. Tambien guarda un  grafico 
    de muestra para comprobacion del usuario.

    Args:
        teo_x (np.ndarray): Datos teoricos para el eje X.
        teo_y (np.ndarray): Datos teoricos para el eje Y.
        files (np.ndarray): Conjunto de archivos de lampara a calibrar.
        window_length (int): Rango de longitudes de onda que una ventana cubre.
        window_step (int): Cantidad de longitudes de onda entre inicios de ventanas.
        detect_teorical_peaks (bool): Condicion booleana para aislar picos de los datos teoricos.
        detect_empirical_peaks (bool): Condicion booleana para aislar picos de los datos empiricoss.
        zero_padding_bool (bool): Condicion booleana para rellenar con ceros los espacios sin registros 
        de intensidades en los datos teoricos a usar.
        normalize_windows (bool): Condicion booleana para normalizar individualmente las ventanas 
        extraidas de los datos teoricos.
        save_dir (str): Path a la carpeta donde guardar los resultados.
        graph (bool): Condicion booleana para saber si se deben generar graficos de las calibraciones 
        generadas o no.
        output_csv_path (str): Path al CSV donde guardar los resultados estadisticos.
        teorical_weed_out (bool): condicion booleana para saber si se debe realizar un desmalezado de 
        datos teoricos o no.
        empirical_weed_out (bool): condicion booleana para saber si se debe realizar un desmalezado de 
        datos empiricos o no.
        minimal_data_for_weed_out (int): valor minimo a considerar para detener el desmalezado de un 
        conjunto de datos.
    """

    # Preparar CSV para persistencia de resultados
    output_csv_path = os.path.join(save_dir, os.path.basename(output_csv_path))
    try: # Si existe lo lee
        df = pd.read_csv(output_csv_path)
    except FileNotFoundError: # Si no existe lo crea
        df = pd.DataFrame(columns=[
            'Teorical peaks',
            'Empirical peaks', 
            'Normalized', 
            'Zero Padding',
            'Teorical WeedOut',
            'Empirical WeedOut',
            'WeedOut Min Value',
            'W_STEP', 
            'W_LENGTH', 
            'Count of Windows', 
            'Distance', 
            'IoU', 
            'Scroll Error',
            'Time'
            ])
        df.to_csv(output_csv_path, index=False)
    
    # Aislado de picos (si corresponde)
    if (detect_teorical_peaks):
        teo_y, _, _ = normalize_min_max(teo_y)
        indices, _ = find_peaks(teo_y, threshold=[0.0, np.inf], height= [0.01, np.inf])
        #print(len(indices))
        #print('-----------')
        teo_x = teo_x[indices]
        teo_y = teo_y[indices]

    # Rellenado de ceros (si corresponde)
    if (zero_padding_bool):
        teo_x, teo_y = zero_padding(arr_x=teo_x, arr_y=teo_y, dist=1.5)
    
    if (teorical_weed_out):
        teo_x, teo_y, _ = weedOut(teo_x, teo_y, minimal_data_for_weed_out, 8)
        
    # Ventanado del Teorico
    ranges, slices_x, slices_y = slice_with_range_step(teo_x, teo_y, window_length, 
                                                    window_step, normalize_windows)
            
    #Filtrar aquellos arreglos que no tienen elementos
    au_x = []
    au_y = []
    for i in range(len(slices_x)):
        if (len(slices_x[i]) > 0):
            au_x.append(slices_x[i])
            au_y.append(slices_y[i])
    slices_x = np.array(au_x, dtype=object)
    slices_y = np.array(au_y, dtype=object)

    # Procesado de los archivos Fe
    Distances = np.array([])
    IoUs = np.array([])
    EAMs = np.array([])
    Times = np.array([])
    for file in files:
            
        # Separar informacion
        emp_x, emp_y, emp_head = extract_lamp_info(file, normalize=True)

        # Determinar calibracion real
        try:
            emp_real_x = emp_x * emp_head['CD1_1'] + emp_head['CRVAL1']
        except Exception as e:
            print(f"Error archivo {file} < Falta de headers")
            continue
        
        # Aislado de picos (si corresponde)
        if (detect_empirical_peaks):
            indices, _ = find_peaks(emp_y, threshold=[0.0, np.inf], height= [0.2, np.inf])
            #print(len(indices))
            emp_x = emp_x[indices]
            emp_y = emp_y[indices]
            emp_real_x = emp_real_x[indices]
            
        if (empirical_weed_out):
            emp_real_x, emp_y, _ = weedOut(emp_real_x, emp_y, minimal_data_for_weed_out, 8)
        
        # print('OK1')
        
        # Registro del tiempo de inicio
        start = time.time()
        
        # Determinar calibracion con DTW
        best_alignment, index = find_best_calibration(emp_y, slices_y, window_length, window_step)
        
        # Registro del tiempo transcurrido
        end = time.time()        
        transcurred_time = end - start
        
        # print('OK2')

        # print('-------------------------------')
        # print(f"best_alignment.index1=[{best_alignment.index1[0]}, "+
        #       f"{best_alignment.index1[1]}, ...{best_alignment.index1[-1]}]")
        # print(f"len={len(best_alignment.index1)}")
        # print('-------------------------------')
        # print(f"best_alignment.index2=[{best_alignment.index2[0]}, "+
        #       f"{best_alignment.index2[1]}, ...{best_alignment.index2[-1]}]")
        # print(f"len={len(best_alignment.index2)}")
        # print('-------------------------------')
        
        # Dispocición en vector de las longitudes de ondas calibradas
        calibrado_x = np.full(len(best_alignment.index1), None)
        for i in range(len(best_alignment.index1)): # Calibrado
            calibrado_x[best_alignment.index1[i]] = slices_x[index][best_alignment.index2[i]]

        # ACOMODAR A PARTIR DE ACA
        # Determinación de la metrica IoU 
        c_inicio = slices_x[index][best_alignment.index2[0]] # inicio calibrado
        c_fin = slices_x[index][best_alignment.index2[-1]] # fin calibrado
        Iou = IoU(c_inicio, c_fin, emp_real_x[0], emp_real_x[-1]) # Segun mejor calibrado
        
        # Agregado de metricas en arreglos de almacenamiento
        Distances = np.append(Distances, best_alignment.distance)
        IoUs = np.append(IoUs, Iou)
        EAMs = np.append(EAMs, EAM(calibrado_x, emp_real_x))
        Times = np.append(Times, transcurred_time)

        if(graph): 
            plt.figure(figsize=(12, 4), dpi=600) # Ajuste de tamaño de la figura

            min_teo_grap = calibrado_x[0] if calibrado_x[0] < emp_real_x[0] else emp_real_x[0] # Seccion del teorico
            min_teo_grap -= window_step
            max_teo_grap = calibrado_x[-1] if calibrado_x[-1] > emp_real_x[-1] else emp_real_x[-1]
            max_teo_grap += window_step
            grap_teo_x, grap_teo_y, _, _ = subconj_generator(teo_x, teo_y, min_teo_grap, max_teo_grap)
            plt.bar(grap_teo_x, -grap_teo_y, width=10, label='Teorical', color='blue', align='edge', alpha=0.7) 
            plt.bar(emp_real_x, emp_y, width=2, label='Emp Real', color='black', align='edge', alpha=0.7) # Real
            plt.bar(calibrado_x, emp_y, width=3, label='Emp Calibrated', color='red', align='edge', alpha=1) # Hallado
            plt.legend()
            fig_name = f"{os.path.splitext(os.path.basename(file))[0]}"
            fig_name += "_EP" if(detect_empirical_peaks) else ""
            fig_name += "_TP" if(detect_teorical_peaks) else ""
            fig_name += "_NOR" if(normalize_windows) else ""
            fig_name += "_ZP" if(zero_padding_bool) else ""
            plt.savefig(os.path.join(save_dir, f"{fig_name}.svg"))
            plt.close()
            
    # Guardar datos promedios de ejecucion en CSV
    nueva_fila = { # Añadir la nueva fila al DataFrame
        'Teorical peaks':detect_teorical_peaks, 
        'Empirical peaks':detect_empirical_peaks, 
        'Normalized':normalize_windows, 
        'Zero Padding':zero_padding_bool,
        'Teorical WeedOut': teorical_weed_out,
        'Empirical WeedOut': empirical_weed_out,
        'WeedOut Min Value': minimal_data_for_weed_out,
        'W_STEP':window_step, 
        'W_LENGTH':window_length, 
        'Count of Windows':len(slices_x),
        'Distance':f"{np.mean(Distances)} \u00B1({np.std(Distances)})", 
        'IoU':f"{np.mean(IoUs)} \u00B1({np.std(IoUs)})",
        'Scroll Error':f"{np.mean(EAMs)} \u00B1({np.std(EAMs)})",
        'Time':f"{np.mean(Times)} \u00B1({np.std(Times)})"
        
    }
    df = df._append(nueva_fila, ignore_index=True)
    df.to_csv(output_csv_path, index=False) # Guardar DataFrame actualizado
    
def concateAndExtractTeoricalFiles(files:list, epsilon):
    """Funcion para extraer y concaternar los datos de varios archivos
    teoricos. Si multiples archivos se superponen entre si entonces 
    resuelve sus conflictos eligiendo los valores maximos.

    Args:
        files (_type_): Conjunto de archivos teoricos.
        epsilon (_type_): Distancia minima que tiene que haber entre
        2 datos teoricos. Si dos datos estan a menos distancia entonces
        se deja el que tenga un valor adjunto mayor.

    Returns:
        numpy.ndarray: Valores para X de la concatenacion realizada.
        numpy.ndarray: Valores para Y de la concatenacion realizada.
    """

    total_x = np.array([])
    total_y = np.array([])
    for i in range(len(files)):
        
        # Separar informacion
        teo_x, teo_y, teo_head = extract_lamp_info(files[i], normalize=False)

        # Determinar teorico real
        try:
            teo_real_x = teo_x * teo_head['CD1_1'] + teo_head['CRVAL1']
        except Exception as e:
            print(f"Error archivo {files[i]} < Falta de headers")
        
        if (i==0):
            total_x = np.array(teo_real_x)
            total_y = np.array(teo_y)
        else:
            inter_total_index = np.where(total_x >= teo_real_x[0])
            x1 = total_x[inter_total_index]
            y1 = total_y[inter_total_index]
            inter_real_index = np.where(teo_real_x <= total_x[-1])
            x2 = teo_real_x[inter_real_index]
            y2 = teo_y[inter_real_index]
            
            # Combinar arreglos y ordenarlos
            x_combinado = np.concatenate((x1, x2))
            y_combinado = np.concatenate((y1, y2))
            indices_orden = np.argsort(x_combinado)
            x_combinado = x_combinado[indices_orden]
            y_combinado = y_combinado[indices_orden]
            
            # Filtrar segun epsilo
            if(len(x_combinado)==0):
                x_filtrado = x_combinado
                y_filtrado = y_combinado
            else:
                indices_seleccionados = []
                for j in range(len(x_combinado)-1):
                    
                    diff = x_combinado[j+1] - x_combinado[j]
                    
                    if diff <= epsilon:
                        # Obtener el índice del elemento con el mayor valor en y
                        indice_max_y = np.argmax(y_combinado[j:j+2]) + j
                        indices_seleccionados.append(indice_max_y)
                    else:
                        indices_seleccionados.append(j)
                indices_seleccionados.append(len(x_combinado)-1)
                
                # Crear los arreglos filtrados
                x_filtrado = x_combinado[indices_seleccionados]
                y_filtrado = y_combinado[indices_seleccionados]
            
            # Concatenamos los arreglos finales
            indices = np.arange(len(total_x))
            masc = np.logical_not(np.isin(indices, inter_total_index))
            x1 = total_x[masc]
            y1 = total_y[masc]
            indices = np.arange(len(teo_real_x))
            masc = np.logical_not(np.isin(indices, inter_real_index))
            x2 = teo_real_x[masc]
            y2 = teo_y[masc]
            total_x = np.concatenate([x1, x_filtrado, x2])
            total_y = np.concatenate([y1, y_filtrado, y2])

    return total_x, total_y

def weedOut(arr_x:np.ndarray , arr_y:np.ndarray, min_data:int, mowing_function_grade:int):
    """Funcion para 'desmalezar' conjuntos de pares de datos.

    Args:
        arr_x (numpy.ndarray): arreglo de datos en el eje X.
        arr_y (numpy.ndarray): arreglo de datos en el eje Y.
        min_data (int): cantidad minima de valores a considerar.
        mowing_function_grade (int): grado de la funcion de desmalezado a entrenar.

    Returns:
        numpy.ndarray: arreglo X desmalezado.
        numpy.ndarray: arreglo Y desmalezado.
        int: cantidad de iteraciones realizadas hasta recortar suficientes datos.
    """
    
    iter = 0
    while (len(arr_x) > min_data):
        iter += 1
        coeficientes = np.polyfit(arr_x, arr_y, mowing_function_grade)
        polinomio = np.poly1d(coeficientes)
        emp_fit_y = polinomio(arr_x)
        indices = arr_y > emp_fit_y
        arr_x = arr_x[indices]
        arr_y = arr_y[indices]
    
    return arr_x, arr_y, iter

def process_batch(teo_x:np.ndarray, teo_y:np.ndarray, files:dict):
    """Funcion que realiza el procesado de un conjunto de archivos variando los valores
    de sus posibles parametros de ejecución.

    Args:
        teo_x (numpy.ndarray): datos teoricos a usar. Correspondientes al eje X.
        teo_y (numpy.ndarray): datos teoricos a usar. Correspondientes al eje Y.
        files (dict): conjunto de paths de los archivos a analizar.
    """
    
    # Definicion de direccion para el guardado de resultados
    act_dir = os.path.dirname(os.path.abspath(__file__))
    save_dir = os.path.join(act_dir, 'output')

    # Especificacion de variables de configuracion
    config = Config(FILES=files, SAVE_DIR=save_dir, WINDOW_STEP=75, 
                    WINDOW_LENGTH=3000, GRAPH=False, OUTPUT_CSV_NAME="output.csv")

    # Preparar CSV para persistencia de resultados
    output_csv_path = os.path.join(config.SAVE_DIR, config.OUTPUT_CSV_NAME)

    # Ejecutar calibraciones para todas las combinaciones de interes
    total_iteraciones = 2*2*2*2*2*2*3
    iteracion_actual = 1
    for detect_teorical_peaks in [True, False]:
        for detect_empirical_peaks in [True, False]:
            for normalize_windows in [True, False]:
                for zero_padding_bool in [True, False]:
                    for teorical_weed_out in [True, False]:
                        for empirical_weed_out in [True, False]:
                            for minimal_data_for_weed_out in [20, 100, 1000]:
                                print(f"{iteracion_actual}/{total_iteraciones}")
                                run_calibrations(
                                    teo_x=teo_x, 
                                    teo_y=teo_y, 
                                    files=config.FILES,
                                    window_length=config.WINDOW_LENGTH,
                                    window_step=config.WINDOW_STEP,
                                    detect_teorical_peaks=detect_teorical_peaks,
                                    detect_empirical_peaks=detect_empirical_peaks,
                                    zero_padding_bool=zero_padding_bool,
                                    normalize_windows=normalize_windows,
                                    save_dir=save_dir,
                                    graph=config.GRAPH,
                                    output_csv_path=output_csv_path,
                                    teorical_weed_out=teorical_weed_out, 
                                    empirical_weed_out=empirical_weed_out, 
                                    minimal_data_for_weed_out=minimal_data_for_weed_out
                                    )
                                iteracion_actual += 1