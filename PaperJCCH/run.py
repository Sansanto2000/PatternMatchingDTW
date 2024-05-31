import os
from utils import concateAndExtractTeoricalFiles, process_batch, get_Data_NIST

# Especificar lamparas a analizar
act_dir = os.path.dirname(os.path.abspath(__file__))
files = {   # Usamos las 3 lamparas Full 'Fe' calibradas por Meilan.
    os.path.join(act_dir, os.path.join("Lamparas_Meilan","w1466Blamp1.fits")),
    # os.path.join(act_dir, os.path.join("Lamparas_Meilan","w1466Blamp2.fits")),
    os.path.join(act_dir, os.path.join("Lamparas_Meilan","w1467Blamp1.fits")),
    # os.path.join(act_dir, os.path.join("Lamparas_Meilan","w1467Blamp2.fits")),
    os.path.join(act_dir, os.path.join("Lamparas_Meilan","w1576Blamp1.fits")),
    # os.path.join(act_dir, os.path.join("Lamparas_Meilan","w1576Blamp2.fits")),
    os.path.join(act_dir, os.path.join("Lamparas_Meilan","w1581Blamp1.fits")),
    # os.path.join(act_dir, os.path.join("Lamparas_Meilan","w1581Blamp2.fits")),
    os.path.join(act_dir, os.path.join("Lamparas_Meilan","w1600Blamp1.fits")),
    # os.path.join(act_dir, os.path.join("Lamparas_Meilan","w1600Blamp2.fits"))
}

print('-----------------------------------------')
# Obtencion de datos teoricos (CARTONES)
EPSILON = 0.055
act_dir = os.path.dirname(os.path.abspath(__file__))
teo_files = [
    os.path.join(act_dir, 'Cartones_Info', 'W_Lamp01.fits'),
    os.path.join(act_dir, 'Cartones_Info', 'W_Lamp02.fits'),
    os.path.join(act_dir, 'Cartones_Info', 'W_Lamp03.fits'),
    os.path.join(act_dir, 'Cartones_Info', 'W_Lamp04.fits'),
    os.path.join(act_dir, 'Cartones_Info', 'W_Lamp05.fits'),
    os.path.join(act_dir, 'Cartones_Info', 'W_Lamp06.fits'),
    os.path.join(act_dir, 'Cartones_Info', 'W_Lamp07.fits'),
    os.path.join(act_dir, 'Cartones_Info', 'W_Lamp08.fits'),
    os.path.join(act_dir, 'Cartones_Info', 'W_Lamp09.fits'),
    os.path.join(act_dir, 'Cartones_Info', 'W_Lamp10.fits')
    ]
teo_x, teo_y = concateAndExtractTeoricalFiles(teo_files, EPSILON)
teo_x = teo_x[teo_y >= 0]
teo_y = teo_y[teo_y >= 0]

print(f'Informacion de los cartones:')
process_batch(teo_x, teo_y, files)

print('-----------------------------------------')
# Obtencion de datos teoricos (NIST)
act_dir = os.path.dirname(os.path.abspath(__file__))
teo_x, teo_y = get_Data_NIST(
    dirpath=act_dir, 
    name='Tabla(NIST)_Int_Long_Mat_Ref.csv', 
    filter=["Fe I", "Fe II"], 
    normalize=True)

print(f'Informacion del NIST:')
process_batch(teo_x, teo_y, files)

print('-----------------------------------------')
print(f'Informacion de Iraft:')
#process_batch(teo_x, teo_y, files)

# print("INICIO")
# run_calibrations(
#     teo_x=teo_x, 
#     teo_y=teo_y, 
#     files=config.FILES,
#     window_length=config.WINDOW_LENGTH,
#     window_step=config.WINDOW_STEP,
#     detect_teorical_peaks=False,
#     detect_empirical_peaks=False,
#     zero_padding_bool=True,
#     normalize_windows=True,
#     save_dir=save_dir,
#     graph=config.GRAPH,
#     output_csv_path=output_csv_path,
#     teorical_weed_out=True, 
#     empirical_weed_out=True, 
#     minimal_data_for_weed_out=20
#     )