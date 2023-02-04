import primersFunctions as u

#Chromosome X results
import pandas as pd
import time

TFsXonly = pd.read_excel("inputfiles/mockMaterials/TFsdf_Xonly.xlsx", index_col = 0)
outputFile = "inputfiles/mockMaterials/XonlyPrimers.xlsx"

start = time.perf_counter()
faultyPrimers= u.finalPrimersdf(TFsXonly, outputFile, returnParameters= True)
stop = time.perf_counter()
print(f"ran function for basic case in {stop - start:0.4f} seconds.")

