import IOU
import os
from utils import getfileData, normalize_min_max, gaussianice
from NIST_Table_Interactor import NIST_Table_Interactor

# Datos del observado
script_dir = os.path.dirname(os.path.abspath(__file__))
filepath = os.path.join(script_dir, "EFBTCOMP31.fits")
obs_data = getfileData(filepath=filepath)

# Separacion de datos observados para el eje X y el eje Y
obs_x = range(len(obs_data))
obs_y = obs_data

# Normalizado de los datos en el eje Y
lamp_y, _, _ = normalize_min_max(obs_y)

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
teorico_y, _, _ = normalize_min_max(target=teo_y)

# Gaussianizado del teorico
SIGMA = 50
teo_x_gau, teo_y_gau = gaussianice(x=teo_x, y=teo_y, resolution=len(obs_x), sigma=SIGMA)    
print(type(teo_y))

# Normalizado del gaussianizado
teo_y_gau, _, _ = normalize_min_max(teo_y_gau)