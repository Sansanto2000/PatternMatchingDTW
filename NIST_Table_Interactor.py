import pandas as pd

import re

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