from NIST_Table_Interactor import NIST_Table_Interactor

import os

import matplotlib.pyplot as plt

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
nor_lamp_y, _, _ = processor.normalize_min_max(lamp_y)

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
plt.figure(figsize=(8, 5))
plt.vlines(teorico_x, 0, teorico_y, colors='orange', linewidth=2)
plt.plot(teorico_x_au, teorico_y_au, color="blue", linestyle='-', alpha=0.6)
plt.xlabel('\u03BB')
plt.ylabel('Intensity')
plt.title(f'Desviación Estándar={sigma}')
plt.savefig(f"TeoricoHeIGaussS={sigma}.png")
plt.show()
