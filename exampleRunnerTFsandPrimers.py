import pandas as pd

#Import functions
import primersFunctions as u #for utils
import TFsdf as t

TF_list = "inputfiles/TFs.xlsx"
annotation = "inputfiles/dmel-all-r6.48.gtf"
ref_genome = "inputfiles/dmel-all-chromosome-r6.48.fasta"

#Import TFsdf
TFsdf, TFs_dict_of_dict = t.make_dataframe_from_TFs_list(TF_list, ref_genome, annotation)
TFsdf.to_excel("inputfiles/mockMaterials/TFsdf.xlsx")
#this has been checked, and has ATG/stop codons at 1703 site.

#An example TFsdf - only the first 5 rows
truncTFsdf = TFsdf.truncate(after = 5)
exampleTemplate = truncTFsdf.at[0, "Reference_Seq"]
#This seq is: exampleTemplate = "TTTGTGGGCGTGGCAAATAGTTGTTTGGCAAATCGCTGCGATGGGAGATAGCTTAGTTAACGAAGATAAGGCTGCATAGCTGCGGTGGCGCATAGTTGAACTGTGTCTTCAAAATTCTGGCGGCTGCGAGACGCGACAAAGGAAGTGAAACATTTTTATACTCGTTACTCGAAGAGTAAAAGGGTATACTCGATATAAAGTATATACATTCTTGATCCGGATCAATAGCCGTGTCGAGCTGGCGATTTGGTCCATCCGTCGGTCATTTCGTATGAACGTCTCGATCTCAGGAACTATAGAAGATGTTCTTTTACATTTTTGTTAGTTCTATAGATTTTATCGATTTGCCAAACAACTTTTCGGCACGCCCACCATCCGCCCACTACTGCCACGCTCACACTTTTGAAAGATGTGTAGATTTTTCTTCGTTTTATTGTTTATTTTGTAAATTTCTACCGCCATACCAAAAACACATTGGTTACGCCATTCCTTACACCTACTTTTGAATAATTTGAAACAATTTTTCTCATTCTATTCCCCAATATCTATTCCCATTCCCCCTAACTGAGTAACTGGTATCTGATAGTCGGGTAACTCGACCATAGCATTCTCTCTTGTTTAAACATAAAATTTCGGAACAATTCATAACAATTTTCCCGTTTGAATAACATTTGATTTTAGAATTTTTGACAAATTATATCTATAATATATAATATCATTTATTTCGAATATGGCTTCATTACATATCTTGGTTTTTGAATTTTGGCATGTAAGATGTTATTCCACTTTTACTCCCATTTTTTTATACAGTAGTATACAGTTTTGAAGAGAATTGAAATACATTTTAATCAAGATATAAATGGTTAGTAGTCTCTCTAAAATAACGATACCCAGAGCATGTGTATAGTACTCGTAGTATATATCCAGCTGTATGAGAAATTTGTAAGTGGAAAAACACGACATCTTATTAATACATTAGTCTGAAATTTATTCACATGGGAATGCTTAAAATAAAAGATTATACCGGAAAAGAAATGATCAAATACAACAGAGTAATTTTTGAGAACTTTCCTATTGGTTTTTTATTAATAACAGAAAAACGTATTATTTTCTAGGGAATCGGAATTTCAATATTTTGGAAATATAAAGAACTCCATTTGATTATAAAGACAGCAAATTAGGCATATTTGTAATTATCATCGCGTTTACAGTCAGAAAAGTACCAATATTGCAATTGAGTAGAAATATTCCTAATATGATTATGGTGACGTTGACATCCGAAGCTAAATCCCATTCGTCTTGAAAGGGCTTACAACAAATTCCGTATTTTTCCAGTGTGAAAAGGTATGAAACCCAAAAACGACATCCATGCAGGAGGCGCTAAGGGACGGAGGACGAGGACGATGGCTCCAGGGTCAGGTAGAAATAATCCAATTTAGAAACCGCTAGTCCAGCGCTAATCCGAAAACGTTTGCGACTCCCGTTCCCGTAGATGTGGATCAATGAATGAATCCGAAGGCTTTCATGAATCTGTTTGCGCGGCGATTAAAAGGGGAATGGATCGCTAATCCGCGTCCCAGGATACACGGGTAATCCTCGAGTTCAGCCTCTTCCTCTCGTATAAATAGATCGGCGGCCGGTGGACTGGTCACAGTCGTTTTACGAACTCAAACAGTGACACATATTGCAATATTATTACGATGGATACGACACCAATCTTCCAGTCCAGCTTCTCCATCCGATCGCTGCTCTCCGTGGACAAAAAGGAGGAGTCCCCCATTTCGAAGCACAACTCAGGAAGCAGCTTCAGCTCATGCTCGAGCTCCAGTTCCAACTCCAGCTCGGATTCGATGGCTGCAAAGAGCAATGCAAAGCCGGCTTTCACCTACAGCGCCCTCATAGTGATGGCCATCTGGAGCAGCTCCGAGAAGCGTTTGACCCTAAGCGGGATCTGCAAGTGGATCGCGGACAACTTCCCGTACTATCGCACCCGCAAGAGCGTCTGGCAGAACTCGATCCGGCACAACCTGAGTCTCAATCCGTTCTTCGTTCGGGTGCCGCGAGCTCTCGATGATCCTGGACGTGGCCACTACTGGGCACTCGATCCGTATGCCGAGGATCTGTCCATTGGCGAGACGACGGGACGCCTGCGCCGCAGCAATTGGCAGCAGAATACGGGGGCGAGACCCAAGGTCACAGGTCATCCCTATCAGCGAATGCCATATTATGGGCACGGGCACGGTAACGGCCCATATATCAAGGCACACAGCGCCTACTTCCCCATAATGGACCATCAACACCATGCCGCCATGGTGCAGCACTACCAGGCCATGATGCACAGATACCAGATGATGCCACATCCTCACCATCACCAGCATCAGCATCAGCATCAGCATCCTCATTCTCATTTCATTCAGCAATCAAAGCCCCTTCATATCCAGGAGCCATATCATCATACCCGCTACCATCTTCACCAAGAGTGAATTTCAGCAGTCTTAGAGGCGTGAAAACTCACTCAAAGACTACGAACTGTGAAAATTCCTTCTTGGCGTAGAACGCGGCAATGTTTTTAGCTTTAAATCATGTAGGGTAATATCCATAGAAACTATGCCCAGGCTTTAAAATATATCGTAGTAAACGTTACGTGCATTCTCATACAAAACTTTAGGTTGTTTTCAAATAAAGATGATTCGACATTCGAAGAATAGTGCGAAATTTTTATTGTATTTTGGTTTAGTGATTTATATATTTTGCTTATTTAAAAATTGAATTAATATTAACTGCCATTGAAAAATAAAACTTAATATCCGAAGAACTTTGTTCTTTTTTCTTAAATAAAATAGGGCATTTTCTGAAGACAATAATTTTGCATTTATATTGGACTTAAAACCGATTTGTACGAATATTGATTTAATTCCAGCAAAGTCCATATGCTAGTTTTTGTTTTGTGATATTCATGTGCGTTTCGTCACTACGTGAAACCCTAACCAACCGCAAAGCCCATGTATTGTTACAGCAGGTACATGCTACAACTGCTTTTTTTTTTTGGTTTTTAATTGCATCCCTGAAAACACTTGAATATTTCAACACAATTTAAATAAAGGACTTCCCAGCTGGTCTCGAGGATCCCAACTAGCCGCTATTCCGGGAAGTTTTCAACTGTTGAAAACTGTCGGTTGAAGCTGCATCGCCACGTTGGATACCACATCCAAATGCCGGCTGTAAATGACTTAAGTGCCTGGCATTTATCTGCTAAAACGTGTCCTGCTTTCGTCCTTGTCATCGTCACAGTTATTCACTTGGCTTTACTTGTTCGAGTGACAGCAGCCGACTCGTTGTTGCCGCTGTCACACCTGTTAAGTAAGTGAATCGGGCAGTGGGTGTACTTGGATGCGGCTACGGC"

#An example primer set for exampleTemplate - for HAL-R
exampleHALRprimers = ['AATAATATTGCAATATGTGTCACTGTTT', 'GTAATAATATTGCAATATGTGTCACTGTT', 'GTAATAATATTGCAATATGTGTCACTGTTT']

#run each function on example data
success, potential_primers = u.designPrimers(2, "VAL-F", exampleTemplate)
    #this works! returns True, ['GCTGCGATGGGAGATAGCTT', 'AAATCGCTGCGATGGGAGAT', 'GTGGGCGTGGCAAATAGTTG', 'TTTGGCAAATCGCTGCGATG', 'GTTTGGCAAATCGCTGCGAT', 'GTTGTTTGGCAAATCGCTGC', 'ATCGCTGCGATGGGAGATAG', 'TGGGCGTGGCAAATAGTTGT', 'TGTGGGCGTGGCAAATAGTT', 'TTGTGGGCGTGGCAAATAGT']

stringency, potential_primers = u.stringencyIterate(exampleTemplate, "VAL-F")
#This one runs as well w. output 0, ['GGTTACGCCATTCCTTACACC', 'TGGTTACGCCATTCCTTACACC', 'ACGCCATTCCTTACACCTAC', 'TGGTTACGCCATTCCTTACAC', 'CCGCCATACCAAAAACACATTG', 'TCTACCGCCATACCAAAAACAC', 'CTACCGCCATACCAAAAACAC', 'TTCTACCGCCATACCAAAAACAC', 'ACCGCCATACCAAAAACACATTG', 'TACGCCATTCCTTACACCTAC']
stringency, potential_primers = u.stringencyIterate(exampleTemplate, "VAL-F", blastSkip= {"VAL-F":5}, useExtRegion= True)
#using all conditions now

mountedPrimers = u.mountedPrimers("HAL-R", exampleHALRprimers, exampleTemplate)
#This gives the only primer that starts at 1669, and makeMount prints 25 bp from 1669
mountedPrimers = u.mountedPrimers("HAL-R", exampleHALRprimers, exampleTemplate, makeMount= True)
#Not sure about the makeMount moutedprimer rev comp - isn't it already the revcomp at the template?

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

#Example BLAST check
sixPrimers = ["GTGGGCGTGGCAAATAGTTG", "TCCATCCGTCGGTCATTTCG",
"TTGACATCCGAAGCTAAATCCCATTCGTC", "AAGGGCTTACAACAAATTCCG",
"GCTGATGCTGATGCTGATGC", "AAACATTGCCGCGTTCTACG"]

faultySix = ["GTGGGCGTGGCAAATAGTTG", "TCCATCCGTCGGTCATTTCG",
"TTGACATCCGAAGCTAAATCCCATTCGTC", "AAGGGCTTACAACAAATTCCG",
"GCTGATGCTGATGCTGATGC", "AA"]

faultyPrimer = u.BLASTSixPrimers(faultySix, 'X')

#Export all final outputs to excel:
TFsdf.to_excel("inputfiles/mockMaterials/finalTFsdfwithPrimers.xlsx")
allPrimers.to_excel("inputfiles/mockMaterials/allPrimers.xlsx")