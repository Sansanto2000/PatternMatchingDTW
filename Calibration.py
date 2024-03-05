class Calibration:
    """Clase que agrupa datos y metricas relevantes de una calibración
    """
    
    def __init__(self, arr_X:list, arr_Y:list, IoU:float, NaC:float, ElM:float = None, ErM:float = None):
        """Metodo de inicialización

        Args:
            arr_X (list): Arreglo de valores que la calibración toma en el eje X
            arr_Y (list): Arreglo de valores que la calibración toma en el eje Y
            IoU (float): Metrica 'IoU'
            NaC (float): Metrica 'Normalized alignment cost'. Un NAC cercano a cero indica una mejor similitud 
            entre las series temporales, tambien sugiere que las dos secuencias están alineadas de manera óptima 
            y son más similares. Un valor lejano a 0 indica que las series temporales son distintas y que no se 
            alinean de forma optima
            ElM (float): Metrica 'Error Lineal Medio'. Un ELM cercano a cero indica una mejor similitud entre las series 
            temporales, tambien sugiere que las dos secuencias están mejor alineadas. Un valor lejano a 0 indica 
            que las series temporales son distintas y que no se alinean de forma optima
            ErM (float): Metrica 'Error Racional Medio'. Un ERM cercano a cero indica una mejor similitud 
            entre las series temporales, tambien sugiere que las dos secuencias están mejor alineadas. Un 
            valor lejano a 0 indica que las series temporales son distintas y que no se alinean de forma optima
        """
        self.arr_X:list = arr_X
        self.arr_Y:list = arr_Y
        self.IoU:float = IoU
        self.NaC:float = NaC
        self.ElM:float = ElM
        self.ErM:float = ErM
    
    def __str__(self):
        """Metodo de especificacion de cadena (str) representativa de instancia de la clase

        Returns:
            str: Cadena representativa de la instancia de la clase
        """
        return f"Calibration:[arr_X(len)={len(self.arr_X)}, \
            arr_Y(len)={len(self.arr_Y)}, \
            IOU={self.IoU}, \
            NAC={self.NaC}]"