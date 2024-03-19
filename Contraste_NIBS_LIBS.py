import matplotlib.pyplot as plt
from utils import get_Data_LIBS, get_Data_NIST

libs_x, libs_y = get_Data_LIBS()

nist_x, nist_y = get_Data_NIST()

plt.figure(figsize=(10, 6), dpi=800)

plt.bar(nist_x, -nist_y, width=6, label='Emp Real', color='purple', align='edge', alpha=1) # NIST

plt.bar(libs_x, libs_y, width=6, label='Emp Real', color='black', align='edge', alpha=1) # LIBS

plt.savefig("prueba.png")

plt.close()