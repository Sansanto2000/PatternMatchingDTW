from NIST_Table_Interactor import NIST_Table_Interactor

import os

import matplotlib.pyplot as plt
import numpy as np

from utils import dp, getfileData, Processor
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
lamp_y, _, _ = processor.normalize_min_max(lamp_y)

plt.figure(figsize=(8, 5))
plt.plot(lamp_x, lamp_y, color="green", alpha=1, label=f"Lampara")
plt.ylim(0, 1)
plt.xlabel('\u03BB')
plt.ylabel('Intensity')
plt.title("Lampara")
plt.legend()
plt.tight_layout()
plt.savefig(f"Lampara.png")
plt.show()

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
teorico_y, _, _ = processor.normalize_min_max(target=teorico_y)

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

plt.figure(figsize=(8, 5))
plt.plot(new_lamp_x, new_lamp_y, color="green", alpha=1, label=f"Lampara")
plt.plot(teorico_x_au, teorico_y_au, color="orange", alpha=0.6, label=f"{filter} calibrado")    
plt.ylim(0, 1)
plt.xlabel('\u03BB')
plt.ylabel('Intensity')
plt.title("Normalized alignment cost: {:.4f}, SIGMA={}".format(cost_mat[N - 1, M - 1]/(N + M), sigma))
plt.legend()
plt.tight_layout()
plt.savefig(f"CalibradoConHeIGaussS={sigma}.png")
plt.show()

