import os
from NIST_Table_Interactor import NIST_Table_Interactor
from utils import normalize_min_max, getfileData, subconj_generator

# Datos y headers del observado
filename = "WCOMP01.fits"
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
#Scipy peak_find


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

# Recorte del segmento de longitudes de onda coincidente con el empirico
teo_x, teo_y, _, _ = subconj_generator(teo_x, teo_y, obs_long_min, obs_long_max)

# Normalizado de los datos teoricos en el eje Y
teo_y, _, _ = normalize_min_max(target=teo_y)

print(len(teo_x))

