#Calculates all primers

import PrimerFunctions as p

#Chromosome X results
import pandas as pd
import time

TFsdf = pd.read_excel("inputfiles/mockMaterials/TFsdf1392-3503.xlsx", index_col = 0)
outputFile = "inputfiles/mockMaterials/TFsoutputcontinued.xlsx"

faultyPrimers= p.finalPrimersdf(TFsdf, outputFile, returnParameters= True)
