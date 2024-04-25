import os
import matplotlib.pyplot as plt
from utils import extract_lamp_info, get_Data_NIST, subconj_generator

# Especificar lamparas a analizar
act_dir = os.path.dirname(os.path.abspath(__file__))
print(act_dir)
FILES = {
    "Ar": {
        "materials": ["Ar I", "Ar II"],
        "files": [
            os.path.join(act_dir, "LampData\\Ar\\whn1.fits"),
            os.path.join(act_dir, "LampData\\Ar\\whn7.fits")
        ]
    },
    "CuNeAr": {
        "materials": ["Cu I", "Cu II", "Ne I", "Ne II", "Ar I", "Ar II"],
        "files": [
            os.path.join(act_dir, "LampData\\CuNeAr\\cdcu-16.ec.fits"),
            os.path.join(act_dir, "LampData\\CuNeAr\\cdcu-21.ec.fits"),
            os.path.join(act_dir, "LampData\\CuNeAr\\cdcu-26.ec.fits"),
            os.path.join(act_dir, "LampData\\CuNeAr\\cdcu-30.ec.fits"),
            os.path.join(act_dir, "LampData\\CuNeAr\\cdcu-35.ec.fits"),
            os.path.join(act_dir, "LampData\\CuNeAr\\cdcu-135.ec.fits"),
            os.path.join(act_dir, "LampData\\CuNeAr\\cdcu-185.ec.fits"),
            os.path.join(act_dir, "LampData\\CuNeAr\\cdcu-235.ec.fits"),
            os.path.join(act_dir, "LampData\\CuNeAr\\cdcu-285.ec.fits"),
            os.path.join(act_dir, "LampData\\CuNeAr\\cdcu-325.ec.fits"),
            os.path.join(act_dir, "LampData\\CuNeAr\\cdcu16.fits"),
            os.path.join(act_dir, "LampData\\CuNeAr\\cdcu21.fits"),
            os.path.join(act_dir, "LampData\\CuNeAr\\cdcu26.fits"),
            os.path.join(act_dir, "LampData\\CuNeAr\\cdcu30.fits"),
            os.path.join(act_dir, "LampData\\CuNeAr\\cdcu35.fits"),
            os.path.join(act_dir, "LampData\\CuNeAr\\cdcu135.fits"),
            os.path.join(act_dir, "LampData\\CuNeAr\\cdcu185.ecc.fits"),
            os.path.join(act_dir, "LampData\\CuNeAr\\cdcu185.fits"),
            os.path.join(act_dir, "LampData\\CuNeAr\\cdcu235.fits"),
            os.path.join(act_dir, "LampData\\CuNeAr\\cdcu275.fits"),
            os.path.join(act_dir, "LampData\\CuNeAr\\cdcu325.fits"),
            os.path.join(act_dir, "LampData\\CuNeAr\\cunear.fits"),
            os.path.join(act_dir, "LampData\\CuNeAr\\cunear.R500.fits")
        ]
    },
    "NeAr": {
        "materials": ["Ne I", "Ne II", "Ar I", "Ar II"],
        "files": [
            os.path.join(act_dir, "LampData\\NeAr\\wc01.ec.fits"),
            os.path.join(act_dir, "LampData\\NeAr\\whn16.fits")
        ]
    }
}

for material_set in FILES.keys():

    for file in FILES[material_set]["files"]:

        # Separar informacion
        emp_x, emp_y, emp_head = extract_lamp_info(file, normalize=True)

        # Determinar calibracion real
        try:
            emp_real_x = emp_x * emp_head['CD1_1'] + emp_head['CRVAL1']
        except Exception as e:
            print(f"Error archivo {file} < Falta de headers")
            continue

        # Determinar teorico que le corresponde
        filter:list=["He I", "Ar I", "Ar II", "NeI"]
        csvpath=os.path.join(act_dir, "NIST\\Tabla(NIST)_Int_Long_Mat_Ref.csv")
        teo_x, teo_y = get_Data_NIST(csvpath=csvpath, filter=filter, normalize=True)
        grap_teo_x, grap_teo_y = subconj_generator(teo_x, teo_y, emp_real_x[0], emp_real_x[-1])

        # Graficar
        plt.figure(figsize=(24, 4), dpi=400)
        plt.bar(emp_real_x, emp_y, width=2, label='Real', color='red', align='edge', alpha=1)
        plt.bar(grap_teo_x, grap_teo_y, width=4, label='Teorical', color='blue', align='edge', alpha=0.9)
        plt.subplots_adjust(left=0.05, right=0.95, top=0.98, bottom=0.16)
        plt.xlabel('Wavelength (Å)', fontsize=20)
        plt.ylabel('Intensity', fontsize=20)
        plt.savefig(os.path.join(act_dir,"LampGraphs", material_set, f"{os.path.basename(file)}.svg"))
        plt.close()

# Aplicar DTW sobre cada archivo

# Guardar resultados:
#   - ¿Hubo match?
#       - No:
#           - Considerar guardar mejor match sin restricciones
#       - Si:
#           - Guardar parametros y metricas