class Config:
    
    def __init__(self, FILES:list, TEORICAL_PATH:str, SAVE_DIR:str, WINDOW_STEP:int, 
                 WINDOW_LENGTH:int, NORMALIZE_WINDOWS:bool, ZERO_PADDING:bool, 
                 DETECT_TEORICAL_PEAKS:bool, DETECT_EMPIRICAL_PEAKS:bool, 
                 GRAPH:bool=False, OUTPUT_CSV_NAME:str="output.csv"):
        """Funcion de inicializacion de la clase Config

        Args:
            FILES (list): Arreglo con los paths de los archivos a calibrar.
            TEORICAL_PATH (str): Path del archivo de datos teoricos a emplear
            SAVE_DIR (str): Path de la carpeta para almacenar los resultados.
            WINDOW_STEP (int): Cantidad de longitudes de onda entre inicios de ventanas.
            WINDOW_LENGTH (int): Rango de longitudes de onda que una ventana cubre.
            NORMALIZE_WINDOWS (bool): Condicion boleana par saber si se normalizaran las 
            ventanas teoricas posterior a su segmentado.
            ZERO_PADDING (bool): Condicion booleana para saber si se rellenaran de ceros 
            los espacios vacios de longitudes de onda o no.
            DETECT_TEORICAL_PEAKS (bool): Condicion boleana par saber si se recortara el 
            teorico segun los picos disponibles o no.
            DETECT_EMPIRICAL_PEAKS (bool): Condicion boleana par saber si se recortara el 
            empirico segun los picos disponibles o no.
            GRAPH (bool, optional): Condicion boleana par saber si se generaran y 
            guardaran graficos de referencia de las mejores calibraciones o no. Default 
            False.
            OUTPUT_CSV_NAME (str, optional): Nombre del archivo csv que se generara con 
            los resultados. Default "output.csv".
        """
        
        self.FILES = FILES
        self.TEORICAL_PATH = TEORICAL_PATH
        self.SAVE_DIR = SAVE_DIR
        self.WINDOW_STEP = WINDOW_STEP
        self.WINDOW_LENGTH = WINDOW_LENGTH
        self.NORMALIZE_WINDOWS = NORMALIZE_WINDOWS
        self.ZERO_PADDING = ZERO_PADDING
        self.DETECT_TEORICAL_PEAKS = DETECT_TEORICAL_PEAKS
        self.DETECT_EMPIRICAL_PEAKS = DETECT_EMPIRICAL_PEAKS
        self.GRAPH = GRAPH
        self.OUTPUT_CSV_NAME = OUTPUT_CSV_NAME
