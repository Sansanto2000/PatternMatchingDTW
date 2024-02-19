import IOU
import os
from DTW import DTW
from IOU import IoU
from utils import getfileData, normalize_min_max, gaussianice, slice_with_range_step, calibrate_with_observations
from NIST_Table_Interactor import NIST_Table_Interactor

# Datos y headers del observado
script_dir = os.path.dirname(os.path.abspath(__file__))
filepath = os.path.join(script_dir, "EFBTCOMPs")
filepath = os.path.join(filepath, "EFBTCOMP31.fits")
obs_data, obs_headers = getfileData(filepath=filepath)

# Longitud de onda minima y maxima del observado calibrado
obs_long_min = obs_headers['CRVAL1'] + 0*obs_headers['CD1_1']
obs_long_max = obs_headers['CRVAL1'] + (len(obs_data)-1)*obs_headers['CD1_1']

# Separacion de datos observados para el eje X y el eje Y
obs_x = range(len(obs_data))
obs_y = obs_data

# Normalizado de los datos en el eje Y
obs_y, _, _ = normalize_min_max(obs_y)

# Datos de teoricos del NIST
csv_filename = os.path.join(script_dir, "Tabla(NIST)_Int_Long_Mat_Ref.csv")
nisttr = NIST_Table_Interactor(csv_filename=csv_filename)

# Especificacion del tipo de lampara a analizar
filter = "He I"

# Obtencion del dataframe
teorico_df = nisttr.get_dataframe(filter=filter)

# Separacion de datos teoricos para el eje X y el eje Y
teo_x = teorico_df['Wavelength(Ams)'].tolist()
teo_y = teorico_df['Intensity'].tolist()

# Normalizado de los datos en el eje Y
teo_y, _, _ = normalize_min_max(target=teo_y)

# Gaussianizado del teorico (innecesario aca)
SIGMA = 50
teo_x_gau, teo_y_gau = gaussianice(x=teo_x, y=teo_y, resolution=len(obs_x), sigma=SIGMA)

# Normalizado del gaussianizado (innecesario aca)
teo_y_gau, _, _ = normalize_min_max(teo_y_gau)

# Aplicación de DTW del observado respecto al teorico (innecesario aca)
DTW(obs_y, teo_y_gau, teo_x_gau)

# Definicion de constantes para realizar la división del teorico en distintas ventanas
STEP = 300 # Cantidad de longitudes de onda entre cada inicio de ventana
W_RANGE = 2000 # Rango de longitudes de onda que una ventana cubre a partir de su paso inicial

# Recorte de los datos del teorico en conjuntos de ventanas
ranges, slices_x, slices_y = slice_with_range_step(teo_x, teo_y, W_RANGE, STEP)

for k in range(0, len(slices_x)):
    
    # Gaussianizado del recorte del teorico
    slice_x_gau, slice_y_gau = gaussianice(x=slices_x[k], y=slices_y[k], 
                                                    resolution=len(obs_x), sigma=SIGMA,
                                                    rang=ranges[k])    
    
    # Normalizado del gaussinizado
    slice_y_gau, _, _ = normalize_min_max(slice_y_gau)
    
    # Aplicación DTW del observado respecto al gaussianizado
    _, _, NorAlgCos = DTW(obs_y, slice_y_gau, slice_x_gau)
    
    # Informe por consola de metricas resultantes de cada calibración tentativa
    print("-----------------")
    print(f"IoU={IoU(slice_x_gau[0], slice_x_gau[len(slice_x_gau)-1], obs_long_min, obs_long_max)}") # IoU 
    print("Normalized alignment cost={:.4f}".format(NorAlgCos)) # Normalized alignment cost
    print("-----------------")