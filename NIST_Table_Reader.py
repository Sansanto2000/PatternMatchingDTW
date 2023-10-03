import pandas as pd

class NIST_Table_Reader:
    
    df = None
    
    def __init__(self, csv_filename):
        df = pd.read_csv(csv_filename)
    
    def get_gaussianized_data() -> list:
        return []