from NIST_Table_Interactor import NIST_Table_Interactor
from DTW import dp
import os
import matplotlib.pyplot as plt
import numpy as np
from utils import dp, getfileData, Processor, normalize_min_max, slice_with_range_step, gaussianice
processor = Processor()

# Datos del espectro a utilizar
script_dir = os.path.dirname(os.path.abspath(__file__))
filepath = os.path.join(script_dir, "EFBTCOMP31.fits")
lamp_data = getfileData(filepath=filepath)
lamp_x = range(len(lamp_data))
lamp_y = lamp_data
# min_lamp_y = np.min(lamp_y)
# max_lamp_y = np.max(lamp_y)
# nor_lamp_y = (lamp_y - min_lamp_y) / (max_lamp_y - min_lamp_y)
lamp_y, _, _ = normalize_min_max(lamp_y)

# plt.figure(figsize=(8, 5))
# plt.plot(lamp_x, lamp_y, color="green", alpha=1, label=f"Lampara")
# plt.ylim(0, 1)
# plt.xlabel('\u03BB')
# plt.ylabel('Intensity')
# plt.title("Lampara")
# plt.legend()
# plt.tight_layout()
# plt.savefig(f"Lampara.png")
# plt.show()

# Datos de picos del NIST
csv_filename = os.path.join(script_dir, "Tabla(NIST)_Int_Long_Mat_Ref.csv")
nisttr = NIST_Table_Interactor(csv_filename=csv_filename)
filter = "He I"
# filter = "Ar II"
# filter = "Ne II"
# filter = "Rb II"
teorico_df = nisttr.get_dataframe(filter=filter)
teorico_x = teorico_df['Wavelength(Ams)'].tolist()
teorico_y = teorico_df['Intensity'].tolist()
teorico_y, _, _ = normalize_min_max(target=teorico_y)

sigma= 50
teorico_x_au, teorico_y_au = processor.gaussianice(x=teorico_x, y=teorico_y, 
                                                   resolution=len(lamp_x), sigma=sigma)    
teorico_y_au, _, _ = processor.normalize_min_max(teorico_y_au)

# plt.figure(figsize=(8, 5))
# plt.vlines(teorico_x, 0, teorico_y, colors='orange', linewidth=2)
# plt.plot(teorico_x_au, teorico_y_au, color="blue", linestyle='-', alpha=0.6)
# plt.xlabel('\u03BB')
# plt.ylabel('Intensity')
# plt.title(f'Desviación Estándar={sigma}')
# plt.savefig(f"TeoricoHeIGaussS={sigma}.png")
# plt.show()

#/////////////////////////////////////// 
# DTW
N = lamp_y.shape[0]
M = len(teorico_y_au)
dist_mat = np.zeros((N, M))
for i in range(N):
    for j in range(M):
        dist_mat[i, j] = abs(lamp_y[i] - teorico_y_au[j])
path, cost_mat = dp(dist_mat)

new_lamp_x = [teorico_x_au[tupla[1]] for tupla in path]
new_lamp_y = [lamp_y[tupla[0]] for tupla in path]

# plt.figure(figsize=(8, 5))
# plt.plot(new_lamp_x, new_lamp_y, color="green", alpha=1, label=f"Lampara")
# plt.plot(teorico_x_au, teorico_y_au, color="orange", alpha=0.6, label=f"{filter} calibrado")    
# plt.ylim(0, 1)
# plt.xlabel('\u03BB')
# plt.ylabel('Intensity')
# plt.title("Normalized alignment cost: {:.4f}, SIGMA={}".format(cost_mat[N - 1, M - 1]/(N + M), sigma))
# plt.legend()
# plt.tight_layout()
# plt.savefig(f"CalibradoConHeIGaussS={sigma}.png")
# plt.show()

#/////////////////////////////////////// 
# DTW por partes

# min_interior = rango_lamp[0] if rango_teorico[0] > rango_lamp[0] else rango_teorico[0]
# max_interior = rango_lamp[1] if rango_lamp[1] < rango_teorico[1] else rango_teorico[1]
# interseccion = max_interior - min_interior
# union = None
# IoU = interseccion / union

W_RANGE = 2000
STEP = 300
SIGMA = 50

ranges, slices_x, slices_y = slice_with_range_step(teorico_x, teorico_y, W_RANGE, STEP)
slices_x_au = []
slices_y_au = []
print("lamp_X=",len(lamp_x))
for k in range(0, len(slices_x)): # <== Meter graficos de las rebanadas
    slice_x_au, slice_y_au = gaussianice(x=slices_x[k], y=slices_y[k], 
                                                    resolution=len(lamp_x), sigma=sigma,
                                                    rang=ranges[k])    
    slice_y_au, _, _ = normalize_min_max(slice_y_au)
    slices_x_au.append(slice_x_au)
    slices_y_au.append(slice_y_au)

for k in range(0, len(slices_x)):
    N = lamp_y.shape[0]
    M = len(slices_y[k])
    dist_mat = np.zeros((N, M))
    for i in range(N):
        for j in range(M):
            dist_mat[i, j] = abs(lamp_y[i] - slices_y[k][j])
    if dist_mat.size == 0:
        continue
    path, cost_mat = dp(dist_mat)

    new_lamp_x = [slices_x[k][tupla[1]] for tupla in path]
    new_lamp_y = [lamp_y[tupla[0]] for tupla in path]
    
    plt.figure(figsize=(8, 5))
    plt.plot(new_lamp_x, new_lamp_y, color="green", alpha=1, label=f"Lampara")
    plt.plot(slices_x[k], slices_y[k], color="orange", alpha=0.6, label=f"{filter} calibrado")    
    print(len(new_lamp_y), len(slices_y[k]))
    plt.ylim(0, 1)
    plt.xlabel('\u03BB')
    plt.ylabel('Intensity')
    plt.title(f"Range = {ranges[k]}")
    plt.legend()
    plt.tight_layout()
    plt.savefig(f"CalibradoConHeIGaussS={sigma}.png")
    plt.show()

