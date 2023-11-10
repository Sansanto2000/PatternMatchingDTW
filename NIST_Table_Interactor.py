import pandas as pd

import re

class NIST_Table_Interactor:
    
    df = None
    
    def __init__(self, csv_filename):
        self.df = pd.read_csv(csv_filename)
        self._sanitize()
        
    def _sanitize(self):
        self.df['Intensity'] = self.df['Intensity'].str.replace(r'[^0-9]+', '', regex=True)
        self.df['Intensity'] = pd.to_numeric(self.df['Intensity'])
        
        self.df['Wavelength(Ams)'] = self.df['Wavelength(Ams)'].apply(lambda x: re.sub(r'[^\d.]', '', x))
        self.df['Wavelength(Ams)'] = pd.to_numeric(self.df['Wavelength(Ams)'])
    
    def get_dataframe(self, cant:int = None, filter:str=None) -> pd.DataFrame:
        df = None
        if (cant):
            df = self.df.head(cant)
        else:
            df = self.df
        if(filter):
            df = df[df['Spectrum'] == filter]
        return df
    
    def _subdivide(self, begin, end, resolution:int) -> list:
        # begin: Longitud de inicio inicial
        # end: Longitud de inicio final
        # resolution: Cantidad de elementos por unidad de longitud de onda
        diff = end - begin
        cant = diff * resolution
        step_size = diff / cant
        elems = []
        for i in range (cant):
            elems.append(begin + i*step_size)
        elems.append(end)
        return elems
    
    def _find_middle_values(self, begin, end) -> list:
        # Funcion para determinar el valor medio entre 2 valores
        diff = end - begin
        middle = begin + diff/2
        return middle
    
    def get_gaussianized_data(self) -> list:
        # Cada elementos del datafreame subdividir la cantida de elementos
        # hallar el punto medio entre los mismos y aplicar a cada lado 
        # la funcion gaussiana que corresponda
        return []
    