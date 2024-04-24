import os
import numpy as np
import matplotlib.pyplot as plt
from scipy.signal import find_peaks
from utils import getfileData, normalize_min_max

# Archivos de los que se quieren determinar los picos
FILES = ["WCOMP01.fits", "WCOMP02.fits", "WCOMP03.fits", "WCOMP04.fits", "WCOMP05.fits",
         "WCOMP06.fits", "WCOMP07.fits", "WCOMP08.fits", "WCOMP09.fits", "WCOMP10.fits",
         "WCOMP11.fits", "WCOMP12.fits", "WCOMP13.fits", "WCOMP14.fits", "WCOMP15.fits",
         "WCOMP16.fits", "WCOMP17.fits", "WCOMP18.fits", "WCOMP19.fits", "WCOMP20_A.fits",
         "WCOMP20_A.fits", "WCOMP21.fits", "WCOMP22.fits", "WCOMP23.fits", "WCOMP24.fits", 
         "WCOMP25.fits", "WCOMP26.fits", "WCOMP27.fits", "WCOMP28.fits", "WCOMP29.fits", 
         "WCOMP30.fits", "WCOMP31.fits"]

# Constantes
HEIGHT = 0.025 # Altura minima a considerar
SAVEPATH = os.path.join(os.path.dirname(os.path.abspath(__file__)), "Graficos") # Especificacion de carpeta para almacenar los graficos

for filename in FILES:
    
    # Datos y headers del observado
    script_dir = os.path.dirname(os.path.abspath(__file__))
    filepath = os.path.join(script_dir, "WCOMPs")
    filepath = os.path.join(filepath, filename)
    obs_data, obs_headers = getfileData(filepath=filepath)

    # Longitud de onda minima y maxima del observado calibrado
    obs_long_min = obs_headers['CRVAL1'] + 0*obs_headers['CD1_1']
    obs_long_max = obs_headers['CRVAL1'] + (len(obs_data)-1)*obs_headers['CD1_1']
    
    # Separacion de datos observados para el eje X y el eje Y
    obs_x = range(len(obs_data))
    obs_y = obs_data

    # Normalizado de los datos obserbados en el eje Y
    obs_y, _, _ = normalize_min_max(obs_y)
    
    # Busqueda de los picos empiricos
    picos_x, _ = find_peaks(obs_y, height=HEIGHT)
    picos_y = obs_y[picos_x]

    # Conversion necesaria para el graficadoa
    obs_x = np.array(obs_x)
    picos_x = np.array(picos_x, dtype=int)
    
    # Ajusta el tamaño de la figura
    plt.figure(figsize=(10, 6))
    
    # Graficar la señal y los picos
    plt.bar(picos_x, picos_y, label='Picos', alpha=1, color='red', linewidth=2)
    plt.plot(obs_x, obs_y, label='Empirico', alpha=1, color='black', linewidth=0.5, linestyle='--')

    plt.legend()
    save_location = os.path.join(SAVEPATH, f'{filename}_PeakFinder.png')
    plt.savefig(save_location)
    #plt.show()
    plt.close()