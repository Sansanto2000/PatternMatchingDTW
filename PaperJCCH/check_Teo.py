import os
import numpy as np
import matplotlib.pyplot as plt
from utils import extract_lamp_info, concateAndExtractTeoricalFiles

# Especificar archivo teorico a emplear
act_dir = os.path.dirname(os.path.abspath(__file__))
files = [
    os.path.join(act_dir, 'Cartones_Info', 'W_Lamp01.fits'),
    os.path.join(act_dir, 'Cartones_Info', 'W_Lamp02.fits'),
    os.path.join(act_dir, 'Cartones_Info', 'W_Lamp03.fits'),
    os.path.join(act_dir, 'Cartones_Info', 'W_Lamp04.fits'),
    os.path.join(act_dir, 'Cartones_Info', 'W_Lamp05.fits'),
    os.path.join(act_dir, 'Cartones_Info', 'W_Lamp06.fits'),
    os.path.join(act_dir, 'Cartones_Info', 'W_Lamp07.fits'),
    os.path.join(act_dir, 'Cartones_Info', 'W_Lamp08.fits'),
    os.path.join(act_dir, 'Cartones_Info', 'W_Lamp09.fits'),
    os.path.join(act_dir, 'Cartones_Info', 'W_Lamp10.fits')
    ]

EPSILON = 0.055

plt.figure(figsize=(12, 4), dpi=600) # Ajuste de tamaño de la figura

total_x = np.array([])
total_y = np.array([])
for i in range(len(files)):
    
    # Separar informacion
    teo_x, teo_y, teo_head = extract_lamp_info(files[i], normalize=False)

    # Determinar teorico real
    try:
        teo_real_x = teo_x * teo_head['CD1_1'] + teo_head['CRVAL1']
    except Exception as e:
        print(f"Error archivo {files[i]} < Falta de headers")
    
    if (i==0):
        total_x = np.array(teo_real_x)
        total_y = np.array(teo_y)
    else:
        inter_total_index = np.where(total_x >= teo_real_x[0])
        x1 = total_x[inter_total_index]
        y1 = total_y[inter_total_index]
        inter_real_index = np.where(teo_real_x <= total_x[-1])
        x2 = teo_real_x[inter_real_index]
        y2 = teo_y[inter_real_index]
        
        # Combinar arreglos y ordenarlos
        x_combinado = np.concatenate((x1, x2))
        y_combinado = np.concatenate((y1, y2))
        indices_orden = np.argsort(x_combinado)
        x_combinado = x_combinado[indices_orden]
        y_combinado = y_combinado[indices_orden]
        
        # Filtrar segun EPSILON
        if(len(x_combinado)==0):
            x_filtrado = x_combinado
            y_filtrado = y_combinado
        else:
            indices_seleccionados = []
            for j in range(len(x_combinado)-1):
                
                diff = x_combinado[j+1] - x_combinado[j]
                
                if diff <= EPSILON:
                    # Obtener el índice del elemento con el mayor valor en y
                    indice_max_y = np.argmax(y_combinado[j:j+2]) + j
                    indices_seleccionados.append(indice_max_y)
                else:
                    indices_seleccionados.append(j)
            indices_seleccionados.append(len(x_combinado)-1)
            
            # Crear los arreglos filtrados
            x_filtrado = x_combinado[indices_seleccionados]
            y_filtrado = y_combinado[indices_seleccionados]
        
        # Concatenamos los arreglos finales
        indices = np.arange(len(total_x))
        masc = np.logical_not(np.isin(indices, inter_total_index))
        x1 = total_x[masc]
        y1 = total_y[masc]
        indices = np.arange(len(teo_real_x))
        masc = np.logical_not(np.isin(indices, inter_real_index))
        x2 = teo_real_x[masc]
        y2 = teo_y[masc]
        total_x = np.concatenate([x1, x_filtrado, x2])
        total_y = np.concatenate([y1, y_filtrado, y2])

plt.bar(total_x, total_y, width=1, label='Teorical', color='blue', align='edge', alpha=1) 

plt.savefig(os.path.join(act_dir, 'check_Teo.svg'))