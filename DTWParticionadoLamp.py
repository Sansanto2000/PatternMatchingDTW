from NIST_Table_Interactor import NIST_Table_Interactor

import os

import matplotlib.pyplot as plt
import numpy as np

from utils import dp, getfileData, Processor

processor = Processor()

# Datos del espectro a utilizar
script_dir = os.path.dirname(os.path.abspath(__file__))
filepath = os.path.join(script_dir, "EFBTCOMP31.fits")
spectral_data = getfileData(filepath=filepath)
x1=range(len(spectral_data))
y1=spectral_data
min_y1 = np.min(y1)
max_y1 = np.max(y1)
nor_y1 = (y1 - min_y1) / (max_y1 - min_y1)

# Datos de picos del NIST
csv_filename = os.path.join(script_dir, "Tabla(NIST)_Int_Long_Mat_Ref.csv")
nisttr = NIST_Table_Interactor(csv_filename=csv_filename)
filter = "He I"
# filter = "Ar II"
# filter = "Ne II"
# filter = "Rb II"
df = nisttr.get_dataframe(filter=filter)

#///////////////////////////////////////
#cant_subarr = 6
# subarr_size = 100 #len(x1) // cant_subarr
# x1_arr = [x1[i:i + subarr_size] for i in range(0, len(x1), subarr_size)]
# y1_arr = [y1[rang] for rang in x1_arr]
# nor_y1_arr = [nor_y1[rang] for rang in x1_arr]
# nor_y1_min = np.min(nor_y1)
# nor_y1_max = np.max(nor_y1)

#///////////////////////////////////////
batch = 25
x2 = df['Wavelength(Ams)'].tolist()
x2_arr = [x2[i:i + batch] for i in range(0, len(x2), batch)]
for i in range(0, len(x2_arr)):
    x2_arr[i] = processor.equalized_and_smooth(reference=x1, target=x2_arr[i], 
                                    function=Processor.FuntionType.LINEAL)
# x2_equalized = processor.equalized_and_smooth(reference=x1, target=x2, 
#                                     function=Processor.FuntionType.LINEAL)
# x2_arr = [x2_equalized[i:i + batch] for i in range(0, len(x2_equalized), batch)]
y2 = df['Intensity'].tolist()
#y2_equalized = processor.equalized_and_smooth(reference=y1, target=y2, 
#                                    function=Processor.FuntionType.SG, 
#                                    window_length=10, polyorder=4)
#nor_y2, _, _ = processor.normalize_min_max(y2_equalized)
nor_y2_arr = []
for i in range(0, len(x2_arr)):
    x = x2_arr[i]
    y = y2[(i*batch):(i*batch+batch)]
    y_equalized = processor.equalized_and_smooth(reference=y1, target=y, 
                                    function=Processor.FuntionType.SG, 
                                    window_length=10, polyorder=4)
    nor_y, _, _ = processor.normalize_min_max(y_equalized)
    #nor_y, _, _ = processor.normalize_min_max(y_equalized, max=np.max(y2), min=np.min(y2))
    nor_y2_arr.append(nor_y)
    
num_rows = 2
num_cols = 2
fig, axs = plt.subplots(num_rows, num_cols, figsize=(8, 5), 
                        subplot_kw={'yticks': []})
fig.suptitle(f"{filter} Ref")
for i in range(num_rows):
    for j in range(num_cols):
        axs[i, j].plot(x2_arr[i*num_cols+j], nor_y2_arr[i*num_cols+j], 
                       color="orange", alpha=1, label=f"Part {i}")    
        axs[i, j].set_ylim(0, 1)
        if j==0:
            axs[i, j].set_yticks(np.linspace(0, 1, 3))
        if(i==2):
            axs[i, j].set_xlabel('\u03BB')
plt.tight_layout()
plt.savefig("partionsOfLamp_graph.png")
plt.show()

fig, axs = plt.subplots(num_rows, num_cols, figsize=(8, 5), 
                        subplot_kw={'yticks': []})
fig.suptitle(f"{filter} Ref")
for row in range(num_rows):
    for col in range(num_cols):
        x = x2_arr[row*num_cols+col]
        y = nor_y2_arr[row*num_cols+col]
        
        N = nor_y1.shape[0]
        M = len(y)
        dist_mat = np.zeros((N, M))
        for i in range(N):
            for j in range(M):
                dist_mat[i, j] = abs(nor_y1[i] - y[j])

        # DTW
        path, cost_mat = dp(dist_mat)
        
        new_x1 = [x[tupla[1]] for tupla in path]
        new_nor_y1 = [nor_y1[tupla[0]] for tupla in path]
        
        axs[row, col].plot(new_x1, new_nor_y1, color="blue", 
                                     alpha=1, label=f"Espectro")
        axs[row, col].plot(x, y, color="orange", 
                                     alpha=0.6, label=f"{filter} calibrado")    
        axs[row, col].set_ylim(0, 1)
        
        if j==0:
            axs[i, j].set_yticks(np.linspace(0, 1, 3))
        if(i==num_rows):
            axs[i, j].set_xlabel('\u03BB')

        axs[row, col].set_title("Normalized alignment cost: {:.4f}".format(cost_mat[N - 1, M - 1]/(N + M)))

fig.suptitle(f"Calibrado")
plt.tight_layout()
plt.savefig("partionsOfLampCalib_graph.png")
plt.show()