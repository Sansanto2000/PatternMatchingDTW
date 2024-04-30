import os
import pandas as pd
from utils import inspect_files_comparison, extract_lamp_info

# Especificar lamparas a analizar
act_dir = os.path.dirname(os.path.abspath(__file__))
FILES = {
    "Ar": {
        "materials": ["Ar I", "Ar II"],
        "files": [
            os.path.join(act_dir, os.path.join("LampData","Ar","whn1.fits")),
            os.path.join(act_dir, os.path.join("LampData","Ar","whn7.fits"))
        ]
    },
    "CuNeAr": {
        "materials": ["Cu I", "Cu II", "Ne I", "Ne II", "Ar I", "Ar II"],
        "files": [
            os.path.join(act_dir, os.path.join("LampData","CuNeAr","cdcu-16.ec.fits")),
            os.path.join(act_dir, os.path.join("LampData","CuNeAr","cdcu-21.ec.fits")),
            os.path.join(act_dir, os.path.join("LampData","CuNeAr","cdcu-26.ec.fits")),
            os.path.join(act_dir, os.path.join("LampData","CuNeAr","cdcu-30.ec.fits")),
            os.path.join(act_dir, os.path.join("LampData","CuNeAr","cdcu-35.ec.fits")),
            os.path.join(act_dir, os.path.join("LampData","CuNeAr","cdcu-135.ec.fits")),
            os.path.join(act_dir, os.path.join("LampData","CuNeAr","cdcu-185.ec.fits")),
            os.path.join(act_dir, os.path.join("LampData","CuNeAr","cdcu-235.ec.fits")),
            os.path.join(act_dir, os.path.join("LampData","CuNeAr","cdcu-285.ec.fits")),
            os.path.join(act_dir, os.path.join("LampData","CuNeAr","cdcu-325.ec.fits")),
            os.path.join(act_dir, os.path.join("LampData","CuNeAr","cdcu16.fits")),
            os.path.join(act_dir, os.path.join("LampData","CuNeAr","cdcu21.fits")),
            os.path.join(act_dir, os.path.join("LampData","CuNeAr","cdcu26.fits")),
            os.path.join(act_dir, os.path.join("LampData","CuNeAr","cdcu30.fits")),
            os.path.join(act_dir, os.path.join("LampData","CuNeAr","cdcu35.fits")),
            os.path.join(act_dir, os.path.join("LampData","CuNeAr","cdcu135.fits")),
            os.path.join(act_dir, os.path.join("LampData","CuNeAr","cdcu185.ecc.fits")),
            os.path.join(act_dir, os.path.join("LampData","CuNeAr","cdcu185.fits")),
            os.path.join(act_dir, os.path.join("LampData","CuNeAr","cdcu235.fits")),
            os.path.join(act_dir, os.path.join("LampData","CuNeAr","cdcu275.fits")),
            os.path.join(act_dir, os.path.join("LampData","CuNeAr","cdcu325.fits")),
            os.path.join(act_dir, os.path.join("LampData","CuNeAr","cunear.fits")),
            os.path.join(act_dir, os.path.join("LampData","CuNeAr","cunear.R500.fits"))
        ]
    },
    "HeArNe": {
        "materials": ["Hg I", "Hg II", "He I", "He II", "Cd I", "Cd II"],
        "files": [
            os.path.join(act_dir, os.path.join("LampData","HeArNe","cdla.19.ec.fits")),
            os.path.join(act_dir, os.path.join("LampData","HeArNe","cdla.23.ec.fits")),
            os.path.join(act_dir, os.path.join("LampData","HeArNe","cdla.43.ec.fits")),
            os.path.join(act_dir, os.path.join("LampData","HeArNe","cdlam.35.ec.fits"))
        ]
    },
    "HgHeCd": {
        "materials": ["Hg I", "Hg II", "He I", "He II", "Cd I", "Cd II"],
        "files": [
            os.path.join(act_dir, os.path.join("LampData","HgHeCd","lamp2_FORS.fits"))
        ]
    },
    "NeAr": {
        "materials": ["Ne I", "Ne II", "Ar I", "Ar II"],
        "files": [
            os.path.join(act_dir, os.path.join("LampData","NeAr","wc01.ec.fits")),
            os.path.join(act_dir, os.path.join("LampData","NeAr","whn16.fits"))
        ]
    },
    "ThAr": {
        "materials": ["Th I", "Th II", "Ar I", "Ar II"],
        "files": [
            os.path.join(act_dir, os.path.join("LampData","ThAr","cdto16.fits")),
            os.path.join(act_dir, os.path.join("LampData","ThAr","cdto21.fits")),
            os.path.join(act_dir, os.path.join("LampData","ThAr","cdto26.fits")),
            os.path.join(act_dir, os.path.join("LampData","ThAr","cdto30.fits")),
            os.path.join(act_dir, os.path.join("LampData","ThAr","cdto30a.fits")),
            os.path.join(act_dir, os.path.join("LampData","ThAr","cdto30b.fits")),
            os.path.join(act_dir, os.path.join("LampData","ThAr","cdto35a.1.fits")),
            os.path.join(act_dir, os.path.join("LampData","ThAr","cdto135.fits")),
            os.path.join(act_dir, os.path.join("LampData","ThAr","cdto185.fits")),
            os.path.join(act_dir, os.path.join("LampData","ThAr","cdto235.fits")),
            os.path.join(act_dir, os.path.join("LampData","ThAr","cdto275.fits")),
            os.path.join(act_dir, os.path.join("LampData","ThAr","cdto325.fits")),
            os.path.join(act_dir, os.path.join("LampData","ThAr","cdto325a.fits")),
            os.path.join(act_dir, os.path.join("LampData","ThAr","thar_sophia.fits")),
            os.path.join(act_dir, os.path.join("LampData","ThAr","toar-DC-R1000.fits")),
            os.path.join(act_dir, os.path.join("LampData","ThAr","toar-DC-R2000.fits")),
            os.path.join(act_dir, os.path.join("LampData","ThAr","toar-DC.fits")),
            os.path.join(act_dir, os.path.join("LampData","ThAr","twrcomp03.fits")),
            os.path.join(act_dir, os.path.join("LampData","ThAr","wfbc03.ec.fits"))
        ]
    }
}

inspect_files_comparison(act_dir=act_dir, files=FILES)
# #file = FILES["CuNeAr"]["files"][0]
# file = FILES["Ar"]["files"][0]
# emp_x, emp_y, emp_head = extract_lamp_info(file, normalize=True)
# emp_real_x = emp_x * emp_head['CD1_1'] + emp_head['CRVAL1']
# print(emp_x)
# print(emp_real_x)
# print(emp_head['CD1_1'])
# print(emp_head['CRVAL1'])

# ---------------------------------------------------------------------
# # Preparar CSV para persistencia de los datos
# CSV_NAME = 'fullDataOutput.csv'
# csv_path = os.path.join(act_dir, CSV_NAME)
# df = pd.DataFrame(
#         columns=['File', 'Cant_Ventanas_Probadas', 'W_STEP', 'W_RANGE', 
#                  'Distance_mejor_ventana', 'IoU_mejor_ventana', 'Error_de_desplazamiento']
#     )
# df.to_csv(csv_path, index=False)

# # Recorre los conjuntos de lamparas de cada material
# for material_set in FILES.keys():
        
#     # Para cada material, obtiene el teorico
#     filter:list=["He I", "Ar I", "Ar II"]
#     teo_x, teo_y = get_Data_NIST(dirpath=os.path.dirname(act_dir), 
#                                  name=CONFIG.TEO_NAME, filter=filter)

#         for file in FILES[material_set]["files"]:


# Aplicar DTW sobre cada archivo

# Guardar resultados:
#   - ¿Hubo match?
#       - No:
#           - Considerar guardar mejor match sin restricciones
#       - Si:
#           - Guardar parametros y metricas