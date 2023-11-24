from NIST_Table_Interactor import NIST_Table_Interactor

from astropy.io import fits
import os
from enum import Enum

#from fastdtw import fastdtw
#from scipy.spatial.distance import euclidean
import matplotlib.pyplot as plt
import numpy as np
from scipy.ndimage import gaussian_filter1d
from scipy.signal import savgol_filter

def getfileData(filepath):    
    hdul = fits.open(filepath) 
    if('WOBJ' in filepath):
        headers = hdul[0].header
        data = hdul[0].data[0][0]
    else:
        headers = hdul[0].header
        data = hdul[0].data
    return data

def dp(dist_mat):
    """
    Find minimum-cost path through matrix `dist_mat` using dynamic programming.

    The cost of a path is defined as the sum of the matrix entries on that
    path. See the following for details of the algorithm:

    - http://en.wikipedia.org/wiki/Dynamic_time_warping
    - https://www.ee.columbia.edu/~dpwe/resources/matlab/dtw/dp.m

    The notation in the first reference was followed, while Dan Ellis's code
    (second reference) was used to check for correctness. Returns a list of
    path indices and the cost matrix.
    """
    N, M = dist_mat.shape
    # Initialize the cost matrix
    cost_mat = np.zeros((N + 1, M + 1))
    for i in range(1, N + 1):
        cost_mat[i, 0] = np.inf
    for i in range(1, M + 1):
        cost_mat[0, i] = np.inf
    # Fill the cost matrix while keeping traceback information
    traceback_mat = np.zeros((N, M))
    for i in range(N):
        for j in range(M):
            penalty = [
                cost_mat[i, j],      # match (0)
                cost_mat[i, j + 1],  # insertion (1)
                cost_mat[i + 1, j]]  # deletion (2)
            i_penalty = np.argmin(penalty)
            cost_mat[i + 1, j + 1] = dist_mat[i, j] + penalty[i_penalty]
            traceback_mat[i, j] = i_penalty
    # Traceback from bottom right
    i = N - 1
    j = M - 1
    path = [(i, j)]
    while i > 0 or j > 0:
        tb_type = traceback_mat[i, j]
        if tb_type == 0:
            # Match
            i = i - 1
            j = j - 1
        elif tb_type == 1:
            # Insertion
            i = i - 1
        elif tb_type == 2:
            # Deletion
            j = j - 1
        path.append((i, j))
    # Strip infinity edges from cost_mat before returning
    cost_mat = cost_mat[1:, 1:]
    return (path[::-1], cost_mat)
    
# Datos del espectro a utilizar
script_dir = os.path.dirname(os.path.abspath(__file__))
filepath = os.path.join(script_dir, "EFBTCOMP31.fits")
#filepath = "C:\\Users\\santi\\OneDrive\\Documentos\\Doctorado\\EFBTCOMP31.fits"
#filepath = "C:\\Users\\santi\\OneDrive\\Documentos\\Doctorado\\EFBTOBJ31.fits"
spectral_data = getfileData(filepath=filepath)
#print(spectral_data)
x1=range(len(spectral_data))
y1=spectral_data
min_y1 = np.min(y1)
max_y1 = np.max(y1)
nor_y1 = (y1 - min_y1) / (max_y1 - min_y1)

# Datos de picos del NIST
#csv_filename = ".\\Tabla(NIST)_Int_Long_Mat_Ref.csv"
csv_filename = os.path.join(script_dir, "Tabla(NIST)_Int_Long_Mat_Ref.csv")
nisttr = NIST_Table_Interactor(csv_filename=csv_filename)
filter = "He I"
# filter = "Ar II"
# filter = "Ne II"
# filter = "Rb II"
df = nisttr.get_dataframe(filter=filter)
#df = nisttr.get_dataframe()

class Processor:
    class FuntionType(Enum):
        LINEAL = "LINEAL"
        GAUSSIAN = "GAUSSIAN"
        SG = "SG"
    def get_procesed(reference, target, tipe = [FuntionType.LINEAL]):
        x2_equalized_lineal = np.interp(np.linspace(0, 1, len(x1)), np.linspace(0, 1, len(x2)), x2)


SIGMA = 7
x2 = df['Wavelength(Ams)'].tolist()
x2_equalized_lineal = np.interp(np.linspace(0, 1, len(x1)), np.linspace(0, 1, len(x2)), x2)
y2 = df['Intensity'].tolist()
y2_equalized_lineal = np.interp(np.linspace(0, 1, len(y1)), np.linspace(0, 1, len(y2)), y2)
y2_equalized_gaussian = gaussian_filter1d(y2_equalized_lineal, sigma=SIGMA)
y2_equalized_SG = savgol_filter(y2_equalized_lineal, window_length=10, polyorder=4)

min_y2 = np.min(y2)
max_y2 = np.max(y2)
nor_y2 = (y2 - min_y2) / (max_y2 - min_y2)
nor_y2_equalized_lineal = (y2_equalized_lineal - min_y2) / (max_y2 - min_y2)
nor_y2_equalized_gaussian = (y2_equalized_gaussian - min_y2) / (max_y2 - min_y2)
nor_y2_equalized_SG = (y2_equalized_SG - min_y2) / (max_y2 - min_y2)

plt.figure(figsize=(8, 5))
plt.plot(x2, nor_y2, color="blue", label=f"{filter}")
#plt.plot(x2_equalized_lineal, nor_y2_equalized_lineal, color="orange", alpha=0.8, label=f"{filter}++")
#plt.plot(x2_equalized_lineal, nor_y2_equalized_gaussian, color="orange", alpha=0.8, label=f"{filter}++")
plt.plot(x2_equalized_lineal, nor_y2_equalized_SG, color="orange", alpha=0.8, label=f"{filter}++")
plt.title(f'{filter} Ref')
plt.legend()
# plt.savefig(f"{filter}_interpoladoLineal++.png")
plt.savefig(f"{filter}_interpoladoGaussiano_SIGMA={SIGMA}++.png")
plt.show()
x2 = x2_equalized_lineal # Comentar si no se quiere suavisar
nor_y2 = nor_y2_equalized_SG # Comentar si no se quiere suavisar

# plt.figure(figsize=(8, 5))
# plt.plot(x1, nor_y1, color="blue", label="Espectro")
# plt.title('Espectro')
# plt.savefig("espectro.png")
# plt.show()

N = nor_y1.shape[0]
M = len(nor_y2)
dist_mat = np.zeros((N, M))
for i in range(N):
    for j in range(M):
        dist_mat[i, j] = abs(nor_y1[i] - nor_y2[j])

# DTW
path, cost_mat = dp(dist_mat)

#\\\\\\\\\\\\\\\\\\\\\\\\\\
# Informacion general
# print(f"N: {N}")
# print(f"M: {M}")
# print(f"dist_mat.shape: {dist_mat.shape}")
# print("Alignment cost: {:.4f}".format(cost_mat[N - 1, M - 1]))
# print("Normalized alignment cost: {:.4f}".format(cost_mat[N - 1, M - 1]/(N + M)))
# print(len(path))
# print(len(x1))
# print(len(x2))

#\\\\\\\\\\\\\\\\\\\\\\\\\\
# Grafico de Correspondencia de Entradas.
# plt.figure(figsize=(5, 5))
# x_path = [i[0] for i in path]
# y_path = [i[1] for i in path]
# plt.plot(y_path, x_path, label='cost_path')
# plt.ylabel('Espectro')
# plt.xlabel(f'{filter} Ref')
# plt.title(f'Correspondencia de Entradas ({filter})')
# plt.savefig("corr_img.png")
# plt.show()

#\\\\\\\\\\\\\\\\\\\\\\\\\\
# Subgráfico 1: Esquina superior izquierda
plt.figure(figsize=(8, 5))
plt.subplot(2, 2, 1)
plt.plot(x1, nor_y1, color="blue", label="Espectro")
plt.title('Espectro')

# Subgráfico 2: Esquina superior derecha
plt.subplot(2, 2, 2)
plt.plot(x2, nor_y2, color="orange", alpha=0.6, label=f"{filter} calibrado")
plt.title(f"{filter} Ref")

# Subgráfico 3: Ocupa los dos recuadros restantes
new_x1 = [x2[tupla[1]] for tupla in path]
new_nor_y1 = [nor_y1[tupla[0]] for tupla in path]
plt.subplot(2, 2, (3, 4))
plt.plot(new_x1, new_nor_y1, color="blue", label="Espectro")
plt.plot(x2, nor_y2, color="orange", alpha=0.6, label=f"{filter} calibrado")

# new_x2 = [x1[tupla[0]] for tupla in path]
# new_nor_y2 = [nor_y2[tupla[1]] for tupla in path]
# plt.subplot(2, 2, (3, 4))
# plt.plot(x1, nor_y1, color="blue", label="Espectro")
# plt.plot(new_x2, new_nor_y2, color="orange", alpha=0.6, label=f"{filter} calibrado")
#plt.legend()
plt.xlabel('\u03BB')
plt.ylabel('Pixel_Val')
plt.title("Normalized alignment cost: {:.4f}".format(cost_mat[N - 1, M - 1]/(N + M)))

plt.tight_layout()
plt.savefig("3part_graph.png")
plt.show()


#NIST_Reference = nisttr.get_gaussianized_data()
    # Extraer con pandas datos de Tabla(NIST)_Int_Long_Mat_Ref.csv
    # ANTES de extraer probablemente haya que sanitizar los datos.


# # OPCION 1: pip install fastdtw
# Supongamos que tienes dos listas de puntos o series temporales llamadas serie1 y serie2
#distancia, camino = fastdtw(spectral_data, serie2, dist=euclidean)

# distancia es la distancia DTW entre las series temporales
# camino contiene el alineamiento entre las series

# # OPCION 2: pip install tslearn
# from tslearn.metrics import dtw
# import numpy as np

# # Supongamos que tienes dos matrices numpy llamadas serie1 y serie2
# distancia, camino = dtw(serie1, serie2)

# # distancia es la distancia DTW entre las series temporales
# # camino contiene el alineamiento entre las series