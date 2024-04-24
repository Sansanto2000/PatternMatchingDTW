import os
from utils import extract_lamp_info


# Especificar lamparas a analizar
act_dir = act_dir = os.path.dirname(os.path.abspath(__file__))
FILES = [
    os.path.join(act_dir, r"/LampData/Ar/whn1.fits"),
    os.path.join(act_dir, r"/LampData/Ar/whn7.fits"),
]

# Recuperar informacion de funcion



# Aplicar DTW sobre cada archivo

# Guardar resultados:
#   - ¿Hubo match?
#       - No:
#           - Considerar guardar mejor match sin restricciones
#       - Si:
#           - Guardar parametros y metricas