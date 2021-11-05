import csv
import numpy as np
import yaml

class loadSettings:
    def __init__(self, file=""):
        ''' Constructor for this class. '''
        self.inFile = file
        
    def loadFileFunc(self, inFile):
        config = yaml.safe_load(open("settings.yml"))