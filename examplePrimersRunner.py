#Import TFsdf
import pandas as pd
TFsdf = pd.read_excel("inputfiles/mockInputs/TFsdf.xlsx") #note that this is double indexed
TFsdf = TFsdf.drop(columns = ["Unnamed: 0"])
#An example TFsdf - only the first 5 rows
truncTFsdf = TFsdf.truncate(after = 5)
exampleTemplate = truncTFsdf.at[0, "Reference_Seq"]
#An example primer set for exampleTemplate - for HAL-R
exampleHALRprimers = ['GACGAATGGGATTTAGCTTCGG', 'ACGAATGGGATTTAGCTTCGG', 'CGAATGGGATTTAGCTTCGGATG', 'TGGGATTTAGCTTCGGATGTC', 'ATGGGATTTAGCTTCGGATGTC', 'GACGAATGGGATTTAGCTTCG', 'AATGGGATTTAGCTTCGGATGTC', 'GAATGGGATTTAGCTTCGGATGTC', 'GAATGGGATTTAGCTTCGGATG', 'ACGAATGGGATTTAGCTTCGGATG', 'GACGAATGGGATTTAGCTTCGGATG', 'CGAATGGGATTTAGCTTCGGATGTC']

#Import functions
import primersFunctions as u #for utils

#run each function on example data
success, primer_sequence = u.designPrimers(2, "VAL-F", exampleTemplate)
    #this works! returns True, ['GCTGCGATGGGAGATAGCTT', 'AAATCGCTGCGATGGGAGAT', 'GTGGGCGTGGCAAATAGTTG', 'TTTGGCAAATCGCTGCGATG', 'GTTTGGCAAATCGCTGCGAT', 'GTTGTTTGGCAAATCGCTGC', 'ATCGCTGCGATGGGAGATAG', 'TGGGCGTGGCAAATAGTTGT', 'TGTGGGCGTGGCAAATAGTT', 'TTGTGGGCGTGGCAAATAGT']

primer_sequence = u.stringencyAndExtension(exampleTemplate, "VAL-F", nextprimer = 0)
#This one runs as well

mountedPrimer = u.mountedPrimers("HAL-R", exampleHALRprimers, exampleTemplate)
#This works! Given the example, it seleced primer 'TTGACATCCGAAGCTAAATCCCATTCGTC' which starts at 1299 (going reverse) 

sixPrimers = u.addSixPrimers(exampleTemplate)
#Works! Prints: ['GTGGGCGTGGCAAATAGTTG', 'TCCATCCGTCGGTCATTTCG', 'TTGACATCCGAAGCTAAATCCCATTCGTC', 'AAGGGCTTACAACAAATTCCG', 'GCTGATGCTGATGCTGATGC', 'AAACATTGCCGCGTTCTACG']

allPrimers, TFsdf = u.finalPrimersdf(truncTFsdf)
#This works too - example outputs in mock materials

#Example problem primer:
problemprimerindex = 0
problemprimertype = 'VAL-F'

TFsdf, allPrimers,  = u.primersBLASTrecalc(problemprimerindex, problemprimertype, truncTFsdf, allPrimers)
#Works, returns: "Problem primer VAL-F at index 0 has been fixed in all dataframes."
#But not sure if it's replacing dataframe values properly so need to check this - might be that TFsdf isn't commiting to the original dataframe

#Export all final outputs to excel:
TFsdf.to_excel("inputfiles/mockMaterials/finalTFsdfwithPrimers.xlsx")
allPrimers.to_excel("inputfiles/mockMaterials/allPrimers.xlsx")