import os
import matplotlib.pyplot as plt
from utils import extract_lamp_info

# Especificar archivo teorico a emplear
act_dir = os.path.dirname(os.path.abspath(__file__))
files = [
    os.path.join(act_dir, 'Cartones_Info', 'W_Lamp01.fits'),
    os.path.join(act_dir, 'Cartones_Info', 'W_Lamp02.fits'),
    #os.path.join(act_dir, 'Cartones_Info', 'W_Lamp03.fits')
    ]
colors = ['red', 'blue', 'green', 'yellow']
EPSILON = 0.06

plt.figure(figsize=(12, 4), dpi=600) # Ajuste de tamaño de la figura

total_teo_x = []
total_teo_y = []
for i in range(len(files)):



    # Separar informacion
    teo_x, teo_y, teo_head = extract_lamp_info(files[i], normalize=False)

    # Determinar teorico real
    try:
        teo_real_x = teo_x * teo_head['CD1_1'] + teo_head['CRVAL1']
    except Exception as e:
        print(f"Error archivo {files[i]} < Falta de headers")

    total_teo_x.extend(teo_real_x)
    total_teo_y.extend(teo_y)

    print(teo_head['CD1_1'])
    print(len(teo_x))
    print(teo_real_x)

    plt.bar(teo_real_x, teo_y, width=1, label='Teorical', color=colors[i], align='edge', alpha=1) 

plt.savefig(os.path.join(act_dir, 'check_Teo.svg'))