import matplotlib.pyplot as plt
import numpy as np
from scipy.signal import find_peaks
from utils import get_Data_LIBS, get_Data_NIST

libs_x, libs_y = get_Data_LIBS(name="LIBS_He_Ar_Ne_Resolution=1000.csv")
if (True): # Aislado de picos
    indices, _ = find_peaks(libs_y, 
                            height=[0.025,np.inf], 
                            #threshold=[0.0000001,np.inf]
                            )
    libs_x = libs_x[indices]
    libs_y = libs_y[indices]
    
nist_x, nist_y = get_Data_NIST()


plt.figure(figsize=(10, 6), dpi=800)

#plt.bar(nist_x, -nist_y, width=6, label='Emp Real', color='purple', align='edge', alpha=1) # NIST

plt.bar(libs_x, libs_y, width=6, label='Emp Real', color='black', align='edge', alpha=1) # LIBS

plt.savefig("Contraste_NIBS_LIBS.png")

plt.close()