#sixPrimers function needs to have chromosome of the template also in the input for this to work

def BLASTSixPrimers(sixPrimers, chromosome):

    #Upload transgenic reference file
        #This will be chromosome 3 of cas_on_2, and all other chromosomes from cas_on_3.
    
    import pandas as pd
    from Bio import SeqIO

    #This corresponds to the order of the list, sixPrimers
    primerTypes = ["VAL-F", "HAL-F", "HAL-R", "HAR-F", "HAR-R", "VAL-R"]

    #This variable will change is a faulty primer is detected
    if chromosome == 3:
        transgenicRefFile = SeqIO.parse("inputfiles/dmel6-nos-Cas9_on_2.fasta", "fasta")    
    else:
        transgenicRefFile = SeqIO.parse("inputfiles/dmel6-nos-Cas9_on_3.fasta", "fasta")    

    for fasta in transgenicRefFile:
        if fasta.id == chromosome:
            transgenicRefSequence = fasta.seq

    for primer in sixPrimers:
        BLASTcount = transgenicRefSequence.count(primer)
        if BLASTcount != 1:
            faultyPrimer = primerTypes[sixPrimers.index(primer)]
        else:
            faultyPrimer = "none"

    return faultyPrimer

sixPrimers = ["GTGGGCGTGGCAAATAGTTG", "TCCATCCGTCGGTCATTTCG",
"TTGACATCCGAAGCTAAATCCCATTCGTC", "AAGGGCTTACAACAAATTCCG",
"GCTGATGCTGATGCTGATGC", "AAACATTGCCGCGTTCTACG"]

faultySix = ["GTGGGCGTGGCAAATAGTTG", "TCCATCCGTCGGTCATTTCG",
"TTGACATCCGAAGCTAAATCCCATTCGTC", "AAGGGCTTACAACAAATTCCG",
"GCTGATGCTGATGCTGATGC", "AA"]

faultyPrimer = BLASTSixPrimers(faultySix, 'X')


            
