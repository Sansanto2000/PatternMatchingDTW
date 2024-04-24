import os
from utils import extract_lamp_info

# Especificar lamparas a analizar
act_dir = os.path.dirname(os.path.abspath(__file__))
print(act_dir)
FILES = [
    os.path.join(act_dir, "LampData\\Ar\\whn1.fits"),
    os.path.join(act_dir, "LampData\\Ar\\whn7.fits"),
]

emp_x, emp_y, emp_head = extract_lamp_info(FILES[0])

print(len(emp_head))
print(emp_head)

# Recuperar informacion de funcion



# Aplicar DTW sobre cada archivo

# Guardar resultados:
#   - ¿Hubo match?
#       - No:
#           - Considerar guardar mejor match sin restricciones
#       - Si:
#           - Guardar parametros y metricas