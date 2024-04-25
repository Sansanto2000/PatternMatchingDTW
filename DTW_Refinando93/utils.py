import re
import numpy as np
import pandas as pd
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