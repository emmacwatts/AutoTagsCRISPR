def revComp(inputSeq):
    """
    This function takes an input sequence and returns the reverse complement.

    params: inputSeq in str format
    returns: revComp in str format

    """
    complement = {'A': 'T', 'C': 'G', 'G': 'C', 'T': 'A'}
    
    revComp = ""
    for base in inputSeq[::-1]:
        revComp += complement[(base.upper())]

    return revComp

def make_dataframe_from_TFs_list(TF_list, ref_genome, annotation):
    '''
    Extracts information and sequence region for genes of interest (TFs) for design of primers per gene.

    params: 
      TFs_list: excel file of query sequences with Gene_ID and Transcript_ID
      ref_genome: fasta file for reference genome
      annotation: .gtf file for ref_genome
    
    returns:
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
    
    #This is to reformat the "Attribute" category in refGenomeAnnotation, to extract Gene_ID, Gene_Symbol, and Transcript_ID
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
    
def filter_gRNA(gRNA_file, TF_dict):
    '''
    gRNA_file: gff file
    TF_dict: dictionary
    '''
    import pandas as pd

    # creating a data frame from the gRNA files
    # makes pandas dataframe from csv file
    gRNAFileAnnotation = pd.read_csv(gRNA_file, sep = "\t", index_col = False)

    #This is to reformat the "Attribute" category in refGenomeAnnotation, to extract Gene_ID, Gene_Symbol, and Transcript ID
    #index = 0

    #Add new categories to the dataframe
    gRNAFileAnnotation = gRNAFileAnnotation.assign(target_site_variation= "")

    #This is to reformat the "Attribute" category in refGenomeAnnotation, to extract Gene_ID, Gene_Symbol, and Transcript ID
    index = 0

    #For each attribute value, extract the gene ID and symbol and add this to the new categories
    for attribute in gRNAFileAnnotation['attributes']:

        fullatt = (gRNAFileAnnotation.loc[index]["attributes"]).split(";")
        gRNAFileAnnotation.at[index,"target_site_variation"] = fullatt[8]
        index+=1
    
    # shorten file to essential information
    GenomeCoordinates = gRNAFileAnnotation.loc[:,["target_site_variation", "fmin", "fmax"]]

    #initialize a list that stores the positions of the gRNAs
    sgRNA_list_pos = []
    sgRNA_list_values = []

    
    for gRNA in range(len(GenomeCoordinates)):
        
        if int(GenomeCoordinates['fmin'].iloc[gRNA]) >= TF_dict["Start"] - 20 and int(GenomeCoordinates['fmax'].iloc[gRNA]) <= TF_dict["Stop"] + 20 \
            and GenomeCoordinates['target_site_variation'].iloc[gRNA] == "target_site_variation=none observed":
            
            sgRNA_list_pos.append([int(GenomeCoordinates['fmin'].iloc[gRNA]), int(GenomeCoordinates['fmax'].iloc[gRNA])])
            
            # distance of gRNA_start and gRNA_stop from start or stop codon
            if TF_dict["Gene_Region"] == "start_codon":

                gRNA_start_diff = TF_dict["Start"] - GenomeCoordinates["fmin"].iloc[gRNA]

                if TF_dict["Start"] <= GenomeCoordinates["fmin"].iloc[gRNA]:

                    gRNA_string = TF_dict["Reference_Seq"][1300 + gRNA_start_diff : 1300 + gRNA_start_diff + 22]

                else:
                    
                    gRNA_string = TF_dict["Reference_Seq"][1300 - gRNA_start_diff : 1300 - gRNA_start_diff + 22]

            else:

                gRNA_start_diff = TF_dict["Stop"] - GenomeCoordinates["fmax"].iloc[gRNA]

                if TF_dict["Stop"] <= GenomeCoordinates["fmax"].iloc[gRNA]:

                    gRNA_string = TF_dict["Reference_Seq"][1300 + gRNA_start_diff : 1300 + gRNA_start_diff + 22]

                else:
                    
                    gRNA_string = TF_dict["Reference_Seq"][1300 - gRNA_start_diff : 1300 - gRNA_start_diff + 22]

            sgRNA_list_values.append(gRNA_string)

    TF_dict["sgRNA_list_pos"] = sgRNA_list_pos
    TF_dict["sgRNA_list_values"] = sgRNA_list_values

    return(TF_dict)

def check_start_stop_NGG(df):
    '''
    Checking whether the best case scenario for sgRNA design holds true.
    The best case scenario would be that the PAM (last 3bp of gRNA; NGG) includes an ATG or a stop codon (TAG,TAA,TGA). 

    params: dictionary with the following format (only mock not real gene)

        df={
            "start/stop":"start",#is it N or C termini -> do we need to look at start or stop codon for teh cut 
            'genome_start_codon_pos':400, 
            'genome_stop_codon_pos':700, # or only 1 of those 
            'strand_type':'+',
            "sgRNA_list_pos":[[401,425],[456,467],[478,489],[395,415]],#those wil be as genome positions -assumptions - the coordinates correspond to the 1st and last bp of the strand to which the gsRNA will be complementary to
            "sgRNA_list_values":["AAGCGACTA","AAAAAAAATAAAAA","ATATATTTTTTTTTTAAAAA","AGCGCGAAATAATA"]
        }

    returns: dictionary with the following format (only mock not real gene)

        df={
            "start/stop":"start",#is it N or C termini -> do we need to look at start or stop codon for teh cut 
            'genome_start_codon_pos':400, 
            'genome_stop_codon_pos':700, # or only 1 of those 
            'strand_type':'+',
            "sgRNA_list_pos":[[401,425],[456,467],[478,489],[395,415]],#those wil be as genome positions -assumptions - the coordinates correspond to the 1st and last bp of the strand to which the gsRNA will be complementary to
            "sgRNA_list_values":["AAGCGACTA","AAAAAAAATAAAAA","ATATATTTTTTTTTTAAAAA","AGCGCGAAATAATA"]
            "best_sgRNA": ["AAGCGACTA", [401,425]]
        }

'''
    df["best_sgRNA"]=[]

    for i in range(0,len(df["sgRNA_list_values"])):
        
        cand=df["sgRNA_list_values"][i] #candidate sgRNA
        cand_len=len(cand) #candidate sgRNA lenght
        
        if df["start/stop"]=="start_codon":

            if cand[(cand_len-4):cand_len]=="ATGG":

                df["best_sgRNA"]=df["sgRNA_list_values"][i]

        if df["start/stop"]=="stop_codon":

            if cand[(cand_len-4):cand_len]=="TAGG":

                df["best_sgRNA"]=(df["sgRNA_list_values"][i], df["sgRNA_list_positions"][i])

        return df

def exclude_over15(df):
    '''
    If the best case scenario does not hold true, we are looking for sgRNAs which span the start/stop codon,
    so that max. 15 bp of the 3´ end (incl. PAM) are on one side.
    
    params: dictionary with the following format (only mock not real gene)

        df={
            "start/stop":"start",#is it N or C termini -> do we need to look at start or stop codon for teh cut 
            'genome_start_codon_pos':400, 
            'genome_stop_codon_pos':700, # or only 1 of those 
            'strand_type':'+',
            "sgRNA_list_pos":[[401,425],[456,467],[478,489],[395,415]],#those wil be as genome positions -assumptions - the coordinates correspond to the 1st and last bp of the strand to which the gsRNA will be complementary to
            "sgRNA_list_values":["AAGCGACTA","AAAAAAAATAAAAA","ATATATTTTTTTTTTAAAAA","AGCGCGAAATAATA"]
        }
    returns: same dictionary as in params but only with sgRNAs listed in "sgRNA_list_pos" and "sgRNA_list_values" that fulfil the condition
    
    '''
    start_pos = df["genome_start_codon_pos"]
    stop_pos = df["genom_stop_codon_pos"]

    list_of_ok=[]
    list_of_ok_pos=[]

    for i in range(0, len(df["sgRNA_list_values"])):

        if df["start/stop"]=="start_codon":

            if 0 < df["sgRNA_list_pos"][i][1]-start_pos < 16: # max. 15 bp max. overhang of sgRNA 3´end from start codon
                
                list_of_ok.append(df["sgRNA_list_values"][i])
                list_of_ok_pos.append(df["sgRNA_list_pos"][i])
        else:

            if 0 < df["sgRNA_list_pos"][i][1]-stop_pos < 16: # max. 15bp overhang of sgRNA 3´end from stop codon

                list_of_ok.append(df["sgRNA_list_values"][i])
                list_of_ok_pos.append(df["sgRNA_list_pos"][i])

    df["sgRNA_list_values"] = list_of_ok
    df["sgRNA_list_pos"] = list_of_ok_pos

    return df

def select_closest(df):
    """
    Pick the sgRNA which has the fewest distance between the cutting site (6bp from the 3´end of the sgRNA) and the stop or start codon. 
    Since there could be multiple sgRNA with the same distance (e.g. left and right form the cutting site) pick the one that does not
    cut within the gene.

    params: dictionary with the following format (only mock not real gene)

        df={
            "start/stop":"start",#is it N or C termini -> do we need to look at start or stop codon for teh cut 
            'genome_start_codon_pos':400, 
            'genome_stop_codon_pos':700, # or only 1 of those 
            'strand_type':'+',
            "sgRNA_list_pos":[[401,425],[456,467],[478,489],[395,415]],#those wil be as genome positions -assumptions - the coordinates correspond to the 1st and last bp of the strand to which the gsRNA will be complementary to
            "sgRNA_list_values":["AAGCGACTA","AAAAAAAATAAAAA","ATATATTTTTTTTTTAAAAA","AGCGCGAAATAATA"]
        }

    returns: closest sgRNA as string

    """
    list_of_pos=len(df["sgRNA_list_values"])*[""]
    list_of_dist=len(df["sgRNA_list_values"])*[""]
    list_of_val=len(df["sgRNA_list_values"])*[""]
    cut_outside_of_gene=len(df["sgRNA_list_values"])*[""]

    for i in range(0,len(df["sgRNA_list_values"])):

        cutting_pos = df["sgRNA_list_pos"][i][1]-6
        
        if df["start/stop"]=="start_codon":

            start_stop_pos=df["genome_start_codon_pos"]

            dist=abs(cutting_pos-start_stop_pos)

            if cutting_pos < start_stop_pos:

                cut_outside_of_gene[i] = "yes"
            
            else: cut_outside_of_gene[i] = "no"

        else:
            
            start_stop_pos=df["genome_stop_codon_pos"]

            dist=abs(cutting_pos-start_stop_pos)

            if cutting_pos < start_stop_pos:

                cut_outside_of_gene[i] = "no"
            
            else: cut_outside_of_gene[i] = "yes"

        list_of_dist[i]=dist
        list_of_val[i]=df["sgRNA_list_values"][i]
        list_of_pos[i]=df["sgRNA_list_pos"][i]
    
    smallest=min(list_of_dist)

    indexes = [i for i, x in enumerate(list_of_dist) if x == smallest]
    
    if len(indexes) == 1:

        df["best_sgRNA"]=(list_of_val[indexes[0]], list_of_pos[indexes[0]])
    
    else: 
        
        for i in indexes: 
            
            if cut_outside_of_gene[i] == "yes" and list_of_dist[i] == smallest:

                df["best_sgRNA"]=(list_of_val[i], list_of_pos[i])

    return df

def find_best_gRNA(df):

    """

    This function is dependent on the functions check_start_stop_NGG(), exclude_over15() and select_closest().

    params: df={
            "start/stop":"start",#is it N or C termini -> do we need to look at start or stop codon for teh cut 
            'genome_start_codon_pos':400, 
            'genome_stop_codon_pos':700, # or only 1 of those 
            'strand_type':'+',
            "sgRNA_list_pos":[[401,425],[456,467],[478,489],[395,415]],#those wil be as genome positions -assumptions - the coordinates correspond to the 1st and last bp of the strand to which the gsRNA will be complementary to
            "sgRNA_list_values":["AAGCGACTA","AAAAAAAATAAAAA","ATATATTTTTTTTTTAAAAA","AGCGCGAAATAATA"]
        }

    returns:

    """

    # check whether best case scenarios hold true
    best_sgRNA=check_start_stop_NGG(df)
  
    if best_sgRNA:

        return best_sgRNA

    else:
        
        allowed_sgRNA = exclude_over15(df) 

        if len(allowed_sgRNA["sgRNA_list_values"])>1: #checking whether multiple candidates meet the condition

            best_sgRNA=select_closest(allowed_sgRNA)

            return best_sgRNA
        
        elif len(allowed_sgRNA["sgRNA_list_values"])==1:

            return allowed_sgRNA
        
        else: return("no sgRNA found for this gene")

  #if best_sgRNA!=[]:
    #get a new batch of sgRNA's and reapet process


def primerSpacing(inputslist):