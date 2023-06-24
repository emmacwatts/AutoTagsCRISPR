#Calculates all primers

import PrimerFunctions as p

import pandas as pd

TFs = 'inputfiles/TFs.xlsx'
refgenome = 'inputfiles/dmel-all-chromosome-r6.48.fasta'
annotation = 'inputfiles/dmel-all-r6.48.gtf'

TFsdf, TFsdict_of_dict = p.make_dataframe_from_TFs_list(TFs, refgenome, annotation)

print (TFsdf)

outputFile = "inputfiles/mockMaterials/TFsfullOutput.xlsx"

p.finalPrimersdf(TFsdf, outputFile, returnParameters= True)
