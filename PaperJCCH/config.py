class Config:
    
    def __init__(self, FILES:list, TEORICAL_PATH:str, SAVE_DIR:str, WINDOW_STEP:int, 
                 WINDOW_LENGTH:int, GRAPH:bool=False, OUTPUT_CSV_NAME:str="output.csv"):
        """Funcion de inicializacion de la clase Config

        Args:
            FILES (list): Arreglo con los paths de los archivos a calibrar.
            TEORICAL_PATH (str): Path del archivo de datos teoricos a emplear
            SAVE_DIR (str): Path de la carpeta para almacenar los resultados.
            WINDOW_STEP (int): Cantidad de longitudes de onda entre inicios de ventanas.
            WINDOW_LENGTH (int): Rango de longitudes de onda que una ventana cubre.
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
        self.GRAPH = GRAPH
        self.OUTPUT_CSV_NAME = OUTPUT_CSV_NAME
