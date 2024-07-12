import os 
import pandas as pd
import matplotlib.pyplot as plt
from utils import extract_lamp_info, get_Data_NIST, subconj_generator, get_Data_IRAF

act_dir = os.path.dirname(os.path.abspath(__file__))
file = os.path.join(act_dir, os.path.join("Lamparas_Meilan","w1466Blamp1.fits"))

# Separar informacion
emp_x, emp_y, emp_head = extract_lamp_info(file, normalize=False)

# Determinar empirico real
try:
    emp_real_x = emp_x * emp_head['CD1_1'] + emp_head['CRVAL1']
except Exception as e:
    print(f"Error archivo {file} < Falta de headers")
    
print(f"Rango empirico: [{emp_real_x[0]}, {emp_real_x[-1]}]")

plt.figure(figsize=(20, 2), dpi=600) # Ajuste de tamaño de la figura
plt.subplots_adjust(left=0, right=1, top=1, bottom=0.2)
plt.bar(emp_real_x, emp_y, width=1, label='Empirical', color='green', align='edge', alpha=1) 
plt.savefig(os.path.join(act_dir, 'Fe_Lamp.png'))

# # Obtencion de datos teoricos (NIST)
# act_dir = os.path.dirname(os.path.abspath(__file__))
# teo_x, teo_y = get_Data_NIST(
#     dirpath=act_dir, 
#     name='Tabla(NIST)_Int_Long_Mat_Ref.csv', 
#     filter=["Fe I", "Fe II"], 
#     normalize=True)
# teo_x, teo_y = subconj_generator(teo_x, teo_y, 0, 10000)

# Obtencion de datos teoricos (IRAF)
act_dir = os.path.dirname(os.path.abspath(__file__))
teo_x, teo_y = get_Data_IRAF(
    dirpath=os.path.join(act_dir, 'IRAF_Lamp'), 
    name='fear.dat', 
    normalize=True)
#teo_x, teo_y = subconj_generator(teo_x, teo_y, 0, 10000)

print(f"Rango datos de referencia: [{teo_x[0]}, {teo_x[-1]}]")

plt.figure(figsize=(40, 2), dpi=600) # Ajuste de tamaño de la figura
plt.subplots_adjust(left=0, right=1, top=1, bottom=0.2)
plt.bar(teo_x, teo_y, width=3, label='Teorical', color='blue', align='edge', alpha=1) 
plt.savefig(os.path.join(act_dir, 'Fear_IRAF.png'))

