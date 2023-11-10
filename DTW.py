from NIST_Table_Interactor import NIST_Table_Interactor

from astropy.io import fits
import os

#from fastdtw import fastdtw
#from scipy.spatial.distance import euclidean
import matplotlib.pyplot as plt
import numpy as np

def getfileData(filepath):    
    
    hdul = fits.open(filepath) 
    
    if('WOBJ' in filepath):
        headers = hdul[0].header
        data = hdul[0].data[0][0]
    else:
        headers = hdul[0].header
        data = hdul[0].data
        
    return data

def gaussian(x, max_value, mean, std_deviation):
    return max_value * np.exp(-((x - mean) ** 2) / (2 * std_deviation ** 2))

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
# plt.figure(figsize=(8, 4))
# #plt.plot(x1, y1)
# plt.plot(x1, nor_y1)
# plt.xlabel('Index')
# plt.ylabel('Pixel_Val')
# plt.title('Mi Gráfico')
# plt.show()

# Datos de picos del NIST
#csv_filename = ".\\Tabla(NIST)_Int_Long_Mat_Ref.csv"
csv_filename = os.path.join(script_dir, "Tabla(NIST)_Int_Long_Mat_Ref.csv")
nisttr = NIST_Table_Interactor(csv_filename=csv_filename)
#df = nisttr.get_dataframe(filter="He I")
#df = nisttr.get_dataframe(filter="Ar II")
#df = nisttr.get_dataframe(filter="Ne II")
df = nisttr.get_dataframe(filter="Rb II")
#df = nisttr.get_dataframe()
x2 = df['Wavelength(Ams)'].tolist()
y2 = df['Intensity'].tolist()
min_y2 = np.min(y2)
max_y2 = np.max(y2)
nor_y2 = (y2 - min_y2) / (max_y2 - min_y2)
# plt.figure(figsize=(8, 4))
# plt.plot(x2, y2)
# #plt.plot(x2, nor_y2)
# plt.xlabel('Wavelength(Ams)')
# plt.ylabel('Intensity')
# plt.title('Mi Gráfico')
# #plt.show()

# Distance matrix
# N = y1.shape[0]
# M = len(y2)
# dist_mat = np.zeros((N, M))
# for i in range(N):
#     for j in range(M):
#         dist_mat[i, j] = abs(y1[i] - y2[j])
N = nor_y1.shape[0]
M = len(nor_y2)
dist_mat = np.zeros((N, M))
for i in range(N):
    for j in range(M):
        dist_mat[i, j] = abs(nor_y1[i] - nor_y2[j])


# DTW
path, cost_mat = dp(dist_mat)
print("Alignment cost: {:.4f}".format(cost_mat[N - 1, M - 1]))
print("Normalized alignment cost: {:.4f}".format(cost_mat[N - 1, M - 1]/(N + M)))

plt.figure(figsize=(6, 4))
plt.subplot(121)
plt.title("Distance matrix")
plt.imshow(dist_mat, cmap=plt.cm.binary, interpolation="nearest", origin="lower")
plt.subplot(122)
plt.title("Cost matrix")
plt.imshow(cost_mat, cmap=plt.cm.binary, interpolation="nearest", origin="lower")
x_path, y_path = zip(*path)
plt.plot(y_path, x_path)



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