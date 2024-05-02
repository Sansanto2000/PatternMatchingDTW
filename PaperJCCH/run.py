import os
import pandas as pd
from config import Config
from utils import run_calibrations

# Especificar lamparas a analizar
act_dir = os.path.dirname(os.path.abspath(__file__))
files = {   # Usamos las 3 lamparas Full 'Fe' calibradas por Meilan.
    "materials": ["Fe I", "Fe II"],
    "files": [
        os.path.join(act_dir, os.path.join("LampData","Fe","...")),
        os.path.join(act_dir, os.path.join("LampData","Fe","...")),
        os.path.join(act_dir, os.path.join("LampData","Fe","..."))
    ]
}

teorical_path = os.path.join(act_dir, "TeoricalData", "...")
save_dir = os.path.join(act_dir, 'output')
config = Config(FILES=files, TEORICAL_PATH=teorical_path, SAVE_DIR=save_dir,
                WINDOW_STEP=25, WINDOW_LENGTH=3000, GRAPH=True,
                OUTPUT_CSV_NAME="output.csv")

# Preparar CSV para persistencia de resultados
output_csv_path = os.path.join(config.SAVE_DIR, config.OUTPUT_CSV_NAME)
try: # Si existe lo lee
    df = pd.read_csv(output_csv_path)
except FileNotFoundError: # Si no existe lo crea
    df = pd.DataFrame(columns=[
        'Teorical peaks',
        'Empirical peaks', 
        'Normalized', 
        'Zero Padding',
        'W_STEP', 
        'W_LENGTH', 
        'Count of Windows', 
        'Distance', 
        'IoU', 
        'Scroll Error',
        'Time'
        ])
    df.to_csv(output_csv_path, index=False)

# Como teorico usamos los datos del carton que mando Yael.
csvpath=os.path.join(act_dir, "TeoricalData", "...")
teo_x, teo_y = get_Data(csvpath=config.TEORICAL_PATH, normalize=True)

# Ejecutar calibraciones para todas las combinaciones de interes
total_iteraciones = 2*2*2*2
iteracion_actual = 1
for detect_teorical_peaks in [True, False]:
    for detect_empirical_peaks in [True, False]:
        for normalize_windows in [True, False]:
            for zero_padding_bool in [True, False]:
                print(f"{iteracion_actual}/{total_iteraciones}")
                run_calibrations(
                    teo_x=teo_x, 
                    teo_y=teo_y, 
                    files=config.FILES,
                    window_length=config.WINDOW_LENGTH,
                    window_step=config.WINDOW_STEP,
                    detect_teorical_peaks=detect_teorical_peaks,
                    detect_empirical_peaks=detect_empirical_peaks,
                    zero_padding=zero_padding_bool,
                    normalize_windows=normalize_windows,
                    save_dir=save_dir,
                    graph=config.GRAPH,
                    output_csv_path=output_csv_path
                    )
                iteracion_actual += 1