from enum import Enum
from astropy.io import fits
from scipy.ndimage import gaussian_filter1d
from scipy.signal import savgol_filter
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
                cost_mat[i, j + 1],  # insertion (1) * PENALTY
                cost_mat[i + 1, j]]  # deletion (2) * PENALTY /////////////////
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

class Processor:
    class FuntionType(Enum):
        LINEAL = "LINEAL"
        GAUSSIAN = "GAUSSIAN"
        SG = "SG"
    def equalized_and_smooth(self, reference, target, function = FuntionType.LINEAL, sigma:float = None, window_length=10, polyorder=4):
        target_equalized_lineal = np.interp(np.linspace(0, 1, len(reference)), 
                                            np.linspace(0, 1, len(target)), target)
        if(function is self.FuntionType.LINEAL):
            return target_equalized_lineal
        elif(function is self.FuntionType.GAUSSIAN):
            target_equalized_gaussian = gaussian_filter1d(target_equalized_lineal, sigma=sigma)
            return target_equalized_gaussian
        elif(function is self.FuntionType.SG):
            target_equalized_SG = savgol_filter(target_equalized_lineal, 
                                                window_length=window_length, 
                                                polyorder=polyorder)
            return target_equalized_SG
        else:
            return None
    def normalize_min_max(self, target, min=None, max=None):
        if(not min):
            min = np.min(target)
        if(not max):
            max = np.max(target)
        nor_target = (target - min) / (max - min)
        return nor_target, min, max