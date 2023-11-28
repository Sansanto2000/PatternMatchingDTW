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

batch = 150
x1_arr = [x1[i:i + batch] for i in range(0, len(x1), batch)]
#y1_arr = [y1[rang] for rang in x1_arr]
nor_y1_arr = [nor_y1[rang] for rang in x1_arr]

# Datos de picos del NIST
csv_filename = os.path.join(script_dir, "Tabla(NIST)_Int_Long_Mat_Ref.csv")
nisttr = NIST_Table_Interactor(csv_filename=csv_filename)
filter = "He I"
# filter = "Ar II"
# filter = "Ne II"
# filter = "Rb II"
df = nisttr.get_dataframe(filter=filter)
x2 = df['Wavelength(Ams)'].tolist()
y2 = df['Intensity'].tolist()
x2_arr = []
y2_arr = []
for i in range(0, len(x1_arr)):
    x = processor.equalized_and_smooth(reference=x1_arr[i], target=x2, 
                                       function=Processor.FuntionType.LINEAL)
    y = processor.equalized_and_smooth(reference=nor_y1_arr[i], target=y2, 
                                       function=Processor.FuntionType.SG, 
                                       window_length=10, polyorder=4)
    y, _, _ = processor.normalize_min_max(y)
    x2_arr.append(x)
    y2_arr.append(y)

#////////////////////////////////////////
num_rows = 2
num_cols = 2
fig, axs = plt.subplots(num_rows, num_cols, figsize=(8, 5), 
                        subplot_kw={'yticks': []})
fig.suptitle(f"Particionado del esptectro")
for i in range(num_rows):
    for j in range(num_cols):
        axs[i, j].plot(x1_arr[i*num_cols+j], nor_y1_arr[i*num_cols+j], 
                       color="orange", alpha=1, label=f"Part {i}")    
        axs[i, j].set_ylim(0, 1)
        if j==0:
            axs[i, j].set_yticks(np.linspace(0, 1, 3))
        if(i==2):
            axs[i, j].set_xlabel('\u03BB')
plt.tight_layout()
plt.savefig("partionsOfSpec_graph.png")
plt.show()

# fig, axs = plt.subplots(num_rows, num_cols, figsize=(8, 5), 
#                         subplot_kw={'yticks': []})
# fig.suptitle(f"{filter} Ref")
# for row in range(num_rows):
#     for col in range(num_cols):
#         x = x2_arr[row*num_cols+col]
#         y = nor_y2_arr[row*num_cols+col]
        
#         N = nor_y1.shape[0]
#         M = len(y)
#         dist_mat = np.zeros((N, M))
#         for i in range(N):
#             for j in range(M):
#                 dist_mat[i, j] = abs(nor_y1[i] - y[j])

#         # DTW
#         path, cost_mat = dp(dist_mat)
        
#         new_x1 = [x[tupla[1]] for tupla in path]
#         new_nor_y1 = [nor_y1[tupla[0]] for tupla in path]
        
#         axs[row, col].plot(new_x1, new_nor_y1, color="blue", 
#                                      alpha=1, label=f"Espectro")
#         axs[row, col].plot(x, y, color="orange", 
#                                      alpha=0.6, label=f"{filter} calibrado")    
#         axs[row, col].set_ylim(0, 1)
        
#         if j==0:
#             axs[i, j].set_yticks(np.linspace(0, 1, 3))
#         if(i==num_rows):
#             axs[i, j].set_xlabel('\u03BB')

#         axs[row, col].set_title("Normalized alignment cost: {:.4f}".format(cost_mat[N - 1, M - 1]/(N + M)))

# fig.suptitle(f"Calibrado")
# plt.tight_layout()
# plt.savefig("partionsOfLampCalib_graph.png")
# plt.show()