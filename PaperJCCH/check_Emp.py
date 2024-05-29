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

print(f"Cantidad de datos original: {len(emp_x)}")
print(f"Rango original: {emp_real_x[-1]-emp_real_x[0]}")

#-----------------------------------------------------------
# # Grafico TEorico real
# plt.figure(figsize=(12, 4), dpi=600) # Ajuste de tamaño de la figura
# plt.bar(emp_real_x, emp_y, width=1, label='Empirical', color='blue', align='edge', alpha=1) 
# plt.savefig(os.path.join(act_dir, 'check_Emp1.svg'))
# plt.legend()
# plt.close()

#-----------------------------------------------------------    
# Aislado de picos
# indices, _ = find_peaks(emp_y, threshold=[0.0, np.inf], height= [0.0, np.inf])
# emp_x = emp_x[indices]
# emp_y = emp_y[indices]
# emp_real_x = emp_real_x[indices]

# print(len(emp_x))

# plt.figure(figsize=(12, 4), dpi=600) # Ajuste de tamaño de la figura
# plt.bar(emp_real_x, emp_y, width=1, label='Empirical', color='blue', align='edge', alpha=1) 
# plt.savefig(os.path.join(act_dir, 'check_Emp2.svg'))

#-----------------------------------------------------------
# # Polinomio Fit
# coeficientes = np.polyfit(emp_real_x, emp_y, 8)
# polinomio = np.poly1d(coeficientes)
# emp_fit_y = polinomio(emp_real_x)

# plt.figure(figsize=(12, 4), dpi=600) # Ajuste de tamaño de la figura
# plt.bar(emp_real_x, emp_y, width=1, label='Empirical', color='blue', align='edge', alpha=1) 
# plt.plot(emp_real_x, emp_fit_y, label='Reaper Line', color='red', alpha=1)
# plt.savefig(os.path.join(act_dir, 'check_Emp3.svg'))
# plt.legend()
# plt.close()

#-----------------------------------------------------------
# # Recorte de maleza (UNICO)
# coeficientes = np.polyfit(emp_real_x, emp_y, 8)
# polinomio = np.poly1d(coeficientes)
# emp_fit_y = polinomio(emp_real_x)
# indices = emp_y > emp_fit_y
# emp_x = emp_x[indices]
# emp_real_x_old = emp_real_x
# emp_real_x = emp_real_x[indices]
# emp_y = emp_y[indices]

# plt.figure(figsize=(12, 4), dpi=600) # Ajuste de tamaño de la figura
# plt.bar(emp_real_x, emp_y, width=1, label='Empirical', color='blue', align='edge', alpha=1) 
# plt.plot(emp_real_x_old, emp_fit_y, label='Reaper Line', color='red', alpha=1)
# plt.savefig(os.path.join(act_dir, 'check_Emp4.svg'))
# plt.legend()
# plt.close()

#-----------------------------------------------------------
# Recorte de maleza (MULTIPLE)
MIN_D = 20 
iteracion = 0
while (len(emp_real_x) > MIN_D):
    iteracion += 1
    coeficientes = np.polyfit(emp_real_x, emp_y, 8)
    polinomio = np.poly1d(coeficientes)
    emp_fit_y = polinomio(emp_real_x)
    indices = emp_y > emp_fit_y
    emp_x = emp_x[indices]
    emp_real_x_old = emp_real_x
    emp_real_x = emp_real_x[indices]
    emp_y = emp_y[indices]

print(f"Cantidad de datos post desmalezado: {len(emp_x)}")
print(f"Cantidad de iteraciones = {iteracion}")

plt.figure(figsize=(12, 4), dpi=600)
plt.bar(emp_real_x, emp_y, width=1, label='Empirical', color='blue', align='edge', alpha=1) 
plt.plot(emp_real_x_old, emp_fit_y, label=f'Reaper Line ({iteracion})', color='red', alpha=1)
plt.savefig(os.path.join(act_dir, 'check_Emp5.svg'))