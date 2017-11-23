import pandas as pd

def Leaf:
    def setup(self):
        # This is the leaf file for output (see function "output to leaffile"
        self.columns = [
            'date', # 10
            'jday', # 6
            'time', # 7
            'Lvs_Init', # 9
            'Lvs_Apr', # 9
            'Leaf_#', # 9
            'area', # 7
            'mass', # 10
            'Sen_Area', # 10
            'Pntl_Area', # 10
            'Longev', # 9
            'CarbRat', # 9
            'SLA', # 9
            'dropped', # 9
            'state', # 9
            'GDD Sum', # 9
        ]
