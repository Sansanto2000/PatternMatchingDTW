"""
Testeo para probar una nueva forma de ejecutar DTW
"""
import os
from dtw import *
import numpy as np
from NIST_Table_Interactor import NIST_Table_Interactor
from utils import getfileData, normalize_min_max

def get_Data_FILE():
    # Datos y headers del observado
    script_dir = os.path.dirname(os.path.abspath(__file__))
    filepath = os.path.join(script_dir, "WCOMPs")
    filepath = os.path.join(filepath, 'WCOMP01.fits')
    obs_data, obs_headers = getfileData(filepath=filepath)
    
    # Separacion de datos observados para el eje X y el eje Y
    obs_x = np.array(range(len(obs_data)))
    obs_y = obs_data
    
    # Normalizado de los datos obserbados en el eje Y
    #obs_y, _, _ = normalize_min_max(obs_y)
    
    return obs_x, obs_y


def get_Data_NIST():
    # Datos de teoricos del NIST
    script_dir = os.path.dirname(os.path.abspath(__file__))
    csv_filename = os.path.join(script_dir, "Tabla(NIST)_Int_Long_Mat_Ref.csv")
    nisttr = NIST_Table_Interactor(csv_filename=csv_filename)
    
    # Especificacion del tipo de lampara a analizar
    filter = ["He I", "Ar I", "Ar II"]

    # Obtencion del dataframe
    teorico_df = nisttr.get_dataframe(filter=filter)

    # Separacion de datos teoricos para el eje X y el eje Y
    teo_x = teorico_df['Wavelength(Ams)'].tolist()
    teo_y = teorico_df['Intensity'].tolist()
    
    # Normalizado de los datos en el eje Y
    #teo_y, _, _ = normalize_min_max(target=teo_y)
    
    return teo_x, teo_y

# Datos de empirico para usar como query
obs_x, obs_y = get_Data_FILE()

print(type(obs_x))

## A noisy sine wave as query
idx = np.linspace(0,6.28,num=100)
query = np.sin(idx) 

## A cosine is for template; sin and cos are offset by 25 samples
template = np.cos(idx)

## Find the best match with the canonical recursion formula
#alignment = dtw(obs_y, template, keep_internals=True)

## Display the warping curve, i.e. the alignment curve
#alignment.plot(type="threeway")

dtw(obs_y, template, keep_internals=True, 
    step_pattern=rabinerJuangStepPattern(6, "c"))\
    .plot(type="twoway",offset=-2)

import matplotlib.pyplot as plt

plt.show()