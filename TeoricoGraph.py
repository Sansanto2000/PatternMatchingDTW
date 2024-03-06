import os
import numpy as np
from NIST_Table_Interactor import NIST_Table_Interactor
from utils import normalize_min_max
import matplotlib.pyplot as plt
import plotly.graph_objects as go

# Datos de teoricos del NIST
script_dir = os.path.dirname(os.path.abspath(__file__))
csv_filename = os.path.join(script_dir, "Tabla(NIST)_Int_Long_Mat_Ref.csv")
nisttr = NIST_Table_Interactor(csv_filename=csv_filename)

# Especificacion del tipo de lampara a analizar
filter = ["He I"]
#filter = ["He I", "Ar I", "Ar II"]

# Obtencion del dataframe
teorico_df = nisttr.get_dataframe(filter=filter)

# Separacion de datos teoricos para el eje X y el eje Y
teo_x = np.array(teorico_df['Wavelength(Ams)'])
teo_y = np.array(teorico_df['Intensity'])

# Normalizado de los datos en el eje Y
#teo_y, _, _ = normalize_min_max(target=teo_y)

# Rango de valores a graficar
min = 2000
max = 20000
grap_teo_x = teo_x[(teo_x >= min) & (teo_x <= max)]
grap_teo_y = teo_y[(teo_x >= min) & (teo_x <= max)]

# Crear figura
fig = go.Figure()

# Añadir un gráfico de barras
fig.add_trace(go.Scatter(x=grap_teo_x, y=grap_teo_y, mode='markers+text', name='Puntos'))
fig.add_trace(go.Bar(x=grap_teo_x, y=grap_teo_y, marker=dict(line=dict(width=200)), name='Barras'))


# Diseño del diseño del gráfico
fig.update_layout(title='Diagrama de Barras',
                  xaxis_title='Longitud de onda',
                  yaxis_title='Intensidad')

# Mostrar el gráfico
fig.show()