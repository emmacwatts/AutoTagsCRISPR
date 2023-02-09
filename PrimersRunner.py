#Calculates all primers

import NprimerFunctions as n

#Chromosome X results
import pandas as pd
import time

TFsXonly = pd.read_excel("inputfiles/mockMaterials/TFsdf.xlsx", index_col = 0)
outputFile = "inputfiles/mockMaterials/finalTFsdfwithPrimers.xlsx"

faultyPrimers= n.finalPrimersdf(TFsXonly, outputFile, returnParameters= True)

