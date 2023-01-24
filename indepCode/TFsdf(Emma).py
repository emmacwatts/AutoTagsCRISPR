def make_dataframe_from_TFs_list(TF_list, ref_genome, annotation):
    '''
    Extracts information and sequence region for genes of interest (TFs) for design of primers per gene.

    Input: 
      TFs_list: excel file of query sequences with Gene_ID and Transcript_ID
      ref_genome: fasta file for reference genome
      annotation: .gtf file for ref_genome
    
    Output
      TFsdf: dataframe of TF information and sequences
      TFsdict_of_dict: TFsdf as a dictionary of dictionaries per index
    '''

    import pandas as pd
    from Bio import SeqIO

    #This is the input file containing the TFs we want to query
    #Imported as a pandas dataframe
    queryTFsdf = pd.read_excel(TF_list)

    #This is the .gtf file with annotations for each gene on the reference genome
    refAnnotationsHeaders = ["Chromosome", "Source", "Gene_Region", "Start", "Stop", "Score", "Strand", "Frame", "Attribute"]
    refGenomeAnnotation = pd.read_csv(annotation, sep = "\t", header = None, index_col = False, names = refAnnotationsHeaders)

    #This is the FASTA file of the reference genome sequence
    refSeqPerChromosome = {}
    for seq in SeqIO.parse(open(ref_genome), 'fasta'):
        refSeqPerChromosome[seq.id] = seq.seq
    
    #This is to reformat the "Attribute" category in refGenomeAnnotation, to extract Gene_ID, Gene_Symbol, and Transcript ID
    index = 0

    #Add new categories to the dataframe
    refGenomeAnnotation = refGenomeAnnotation.assign(Gene_ID = "", Gene_Symbol = "", Transcript_ID = "")

    #For each attribute value, extract the gene ID and symbol and add this to the new categories
    for attribute in refGenomeAnnotation['Attribute']:
        fullatt = (refGenomeAnnotation.loc[index]["Attribute"]).replace(";", "")
        fullatt = fullatt.replace('"', "")
        fullattsplit = fullatt.split(" ")
        refGenomeAnnotation.at[index,"Gene_ID"] = fullattsplit[1]
        refGenomeAnnotation.at[index,"Gene_Symbol"] = fullattsplit[3]
        if len(fullattsplit) == 8:
            refGenomeAnnotation.at[index,"Transcript_ID"] = fullattsplit[5]
        index+=1

    #Delete Attributes category
    del refGenomeAnnotation["Attribute"]

    #Select only rows that TFs are in, and keep only the start and stop codon gene regions

    refGenomeAnnotation = refGenomeAnnotation.loc[refGenomeAnnotation["Gene_Region"].isin(["start_codon", "stop_codon"])]

    TFsdf = refGenomeAnnotation[["Gene_ID", "Gene_Symbol", "Transcript_ID", "Chromosome", "Gene_Region", "Start", "Stop", "Strand"]].loc[refGenomeAnnotation["Gene_ID"].isin(queryTFsdf["Flybase_ID"])]

    #Add reference genome sequence per gene region
    #This will correspond to 1.3kb upstream and downstream of ATG/stop codon 
    TFsdf = TFsdf.assign(Reference_Seq = "")

    for index, rowcontents in TFsdf.iterrows():
        if rowcontents["Strand"] == "+":

            #Define 2.6kb gene region
            regionStart = rowcontents["Start"] - 1701
            regionStop = rowcontents["Stop"] + 1700

            #Add reference sequence
            TFsdf.at[index,"Reference_Seq"] = str(refSeqPerChromosome[rowcontents["Chromosome"]][regionStart:regionStop])

        if rowcontents["Strand"] == "-":

            #Define 2.6kb gene region
            regionStart = rowcontents["Start"] - 1701
            regionStop = rowcontents["Stop"] + 1700

            #Add reference sequence
            refPosStrandSeq = str(refSeqPerChromosome[rowcontents["Chromosome"]][regionStart:regionStop]) #This is the + strand seq, so goes from end to beginning
            TFsdf.at[index,"Reference_Seq"] = revComp(refPosStrandSeq)

    #re-index
    TFsdf = TFsdf.reset_index()
    del TFsdf["index"]

    #Create a dictionary of dictionaries for the df
    TFsdict_of_dict = {}
    for index, rowcontents in TFsdf.iterrows():
        TFsdict_of_dict[index] = rowcontents.to_dict()

    return TFsdf, TFsdict_of_dict