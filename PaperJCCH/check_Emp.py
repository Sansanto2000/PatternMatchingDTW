import os
import numpy as np
import matplotlib.pyplot as plt
from scipy.signal import find_peaks
from utils import extract_lamp_info

# Especificar archivo teorico a emplear
act_dir = os.path.dirname(os.path.abspath(__file__))
file = os.path.join(act_dir, 'Lamparas_Meilan', 'w1466Blamp2.fits')

# Separar informacion
emp_x, emp_y, emp_head = extract_lamp_info(file, normalize=True)

# Determinar teorico real
try:
    emp_real_x = emp_x * emp_head['CD1_1'] + emp_head['CRVAL1']
except Exception as e:
    print(f"Error archivo {file} < Falta de headers")

print(len(emp_x))
print (emp_real_x[-1]-emp_real_x[0])

plt.figure(figsize=(12, 4), dpi=600) # Ajuste de tamaño de la figura
plt.bar(emp_real_x, emp_y, width=1, label='Empirical', color='blue', align='edge', alpha=1) 
plt.savefig(os.path.join(act_dir, 'check_Emp1.svg'))
plt.close()
    
# Aislado de picos
indices, _ = find_peaks(emp_y, threshold=[0.0, np.inf], height= [0.1, np.inf])
emp_x = emp_x[indices]
emp_y = emp_y[indices]
emp_real_x = emp_real_x[indices]

print(len(emp_x))

plt.figure(figsize=(12, 4), dpi=600) # Ajuste de tamaño de la figura
plt.bar(emp_real_x, emp_y, width=1, label='Empirical', color='blue', align='edge', alpha=1) 
plt.savefig(os.path.join(act_dir, 'check_Emp2.svg'))

# # Aislado de picos (OTRA VEZ)
# indices, _ = find_peaks(emp_y, threshold=[0.0, np.inf], height= [0.0, np.inf])
# emp_x = emp_x[indices]
# emp_y = emp_y[indices]
# emp_real_x = emp_real_x[indices]

# print(len(emp_x))

# plt.figure(figsize=(12, 4), dpi=600) # Ajuste de tamaño de la figura
# plt.bar(emp_real_x, emp_y, width=1, label='Empirical', color='blue', align='edge', alpha=1) 
# plt.savefig(os.path.join(act_dir, 'check_Emp3.svg'))