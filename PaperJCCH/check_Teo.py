import os
import numpy as np
import matplotlib.pyplot as plt
from scipy.signal import find_peaks
from utils import extract_lamp_info, concateAndExtractTeoricalFiles, normalize_min_max

# Especificar archivo teorico a emplear
act_dir = os.path.dirname(os.path.abspath(__file__))
files = [
    # os.path.join(act_dir, 'Cartones_Info', 'W_Lamp01.fits'),
    # os.path.join(act_dir, 'Cartones_Info', 'W_Lamp02.fits'),
    # os.path.join(act_dir, 'Cartones_Info', 'W_Lamp03.fits'),
    # os.path.join(act_dir, 'Cartones_Info', 'W_Lamp04.fits'),
    os.path.join(act_dir, 'Cartones_Info', 'W_Lamp05.fits'),
    os.path.join(act_dir, 'Cartones_Info', 'W_Lamp06.fits'),
    os.path.join(act_dir, 'Cartones_Info', 'W_Lamp07.fits'),
    os.path.join(act_dir, 'Cartones_Info', 'W_Lamp08.fits'),
    os.path.join(act_dir, 'Cartones_Info', 'W_Lamp09.fits'),
    os.path.join(act_dir, 'Cartones_Info', 'W_Lamp10.fits')
    ]

EPSILON = 0.055

total_x, total_y = concateAndExtractTeoricalFiles(files, EPSILON)

print(len(total_x))
print(total_x[-1]-total_x[0])

# Tratamiento
indices = np.where(total_y >= 0)
total_x = total_x[indices]
total_y = total_y[indices]
total_y, _, _ = normalize_min_max(total_y)
indices, _ = find_peaks(total_y, threshold=[0.0, np.inf], height= [0.08, np.inf])
total_x = total_x[indices]
total_y = total_y[indices]

print(len(total_x))

# Graficado
plt.figure(figsize=(12, 4), dpi=600) # Ajuste de tamaño de la figura
plt.bar(total_x, total_y, width=1, label='Teorical', color='blue', align='edge', alpha=1) 
plt.savefig(os.path.join(act_dir, 'check_Teo.svg'))