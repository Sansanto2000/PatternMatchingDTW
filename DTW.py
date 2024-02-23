import numpy as np
    
def dp(dist_mat, penalty_matrix = np.array([1,1,1])):
    """
    Find minimum-cost path through matrix `dist_mat` using dynamic programming.

    The cost of a path is defined as the sum of the matrix entries on that
    path. See the following for details of the algorithm:

    - http://en.wikipedia.org/wiki/Dynamic_time_warping
    - https://www.ee.columbia.edu/~dpwe/resources/matlab/dtw/dp.m

    The notation in the first reference was followed, while Dan Ellis's code
    (second reference) was used to check for correctness. Returns a list of
    path indices and the cost matrix.
    
    Args:
        arr1 (_type_): Matriz de distancia
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
                cost_mat[i, j + 1],  # insertion (1) * PENALTY
                cost_mat[i + 1, j]]  # deletion (2) * PENALTY /////////////////
            penalty*=penalty_matrix
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

def DTW(arr1:np.ndarray, arr2:np.ndarray, arrRef:np.ndarray, penalty_matrix = np.array([1,1,1])):
    """Funcion que en base a los datos de 2 funciones y una refencia en eje X aplica DTW para
    obtener los datos X e Y de la funcion councidente
    
    Args:
        arr1 (numpy.ndarray): Arreglo de datos correspondientes a la 1ra funcion a comparar 
        (Usualmente se emplean los arreglos de datos correspondientes al eje Y)
        arr2 (numpy.ndarray): Arreglo de datos correspondientes a la 2da funcion a comparar 
        (Usualmente se emplean los arreglos de datos correspondientes al eje Y)
        arrRef (numpy.ndarray): Arreglo de datos para usar como referencia para distribuir los
        resultados del analisis en el eje X (Usualmente se emplean los datos de arr2 
        correspondientes al eje X)
        penalty_matrix (numpy.array): Arreglo de pesos a considerar para multiplicar las penalizaciones
        de [coincidencia, insercion, borrado]
    
    Returns:
        np.ndarray: Arreglo resultante para el eje X
        np.ndarray: Arreglo resultante para el eje Y
        numpy.float64: Metrica 'Normalized Alignment Cost'. Un NAC cercano a cero indica una mejor similitud 
        entre las series temporales, tambien sugiere que las dos secuencias están alineadas de manera óptima 
        y son más similares. Un valor lejano a 0 indica que las series temporales son distintas y que no se 
        alinean de forma optima
    """
    if (type(arr1)==list):
        N = len(arr1)
    else:
        N = arr1.shape[0]
    M = len(arr2)
    
    dist_mat = np.zeros((N, M))
    for i in range(N):
        for j in range(M):
            dist_mat[i, j] = abs(arr1[i] - arr2[j])
            
    path, cost_mat = dp(dist_mat, penalty_matrix)

    new_arr_x = [arrRef[tupla[1]] for tupla in path]
    new_arr_y = [arr1[tupla[0]] for tupla in path]
    NorAlgCos = cost_mat[N - 1, M - 1]/(N + M)
    
    return new_arr_x, new_arr_y, NorAlgCos