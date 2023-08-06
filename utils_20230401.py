def revComp(inputSeq):
  """
  This function takes an input sequence and returns the reverse complement.

  Input: inputSeq in str format
  Output: revComp in str format

  """
  complement = {'A': 'T', 'C': 'G', 'G': 'C', 'T': 'A'}
  
  revComp = ""
  for base in inputSeq[::-1]:
    revComp += complement[(base.upper())]

  return revComp

def make_dataframe_from_TFs_list(TF_list, refSeqPerChromosome = "inputfiles/dmel-all-chromosome-r6.48.fasta", annotation = "inputfiles/dmel-all-r6.48.gtf"):
    '''
    Creating a dataframe from sequence information and genes of interest. Depends on function revComp(). 

    params:
            TF_list: xlsx file listing the genes of interest
            ref_genome: fasta file containing the gene sequence information for the reference genome
            annotation: gtf file containing the gene annotation information for the reference genome


    returns: dataframe with dictionaries for every isotype for each gene of interest
        {..., 
        108843:{    'Gene_ID': 'FBgn0004652',
                    'Transcript_ID': 'FBtr0083651', 
                    'Chromosome': '3R', 
                    'Gene_Region': 'stop_codon', 
                    'Start': 18426145, 
                    'Stop': 18426147, 
                    'Strand': '-', 
                    'Reference_Seq': 'TTGATCGTAGGACAC', 
                }
        ...,
        }

    '''
    import pandas as pd
    
    #This is the input file containing the TFs we want to query
    #Imported as a pandas dataframe
    queryTFsdf = pd.read_excel(TF_list)
    #Ths contains 765 TFs total

    #This is the .gtf file with annotations for each gene on the reference genome
    #Note that I'll use these for the info categories of the final pandas df (rather than the transgenic genome annotations)
    refAnnotationsHeaders = ["Chromosome", "Source", "Gene_Region", "Start", "Stop", "Score", "Strand", "Frame", "Attribute"]
    refGenomeAnnotation = pd.read_csv(annotation, sep = "\t", header = None, index_col = False, names = refAnnotationsHeaders)
    
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
        if len(fullattsplit) > 4:
            refGenomeAnnotation.at[index,"Transcript_ID"] = fullattsplit[5]
            refGenomeAnnotation.at[index,"Transcript_Symbol"] = fullattsplit[7]
        index+=1

    #Delete Attributes category
    del refGenomeAnnotation["Attribute"]

    #Select only rows that TFs are in, and keep only the start and stop codon gene regions

    refGenomeAnnotation = refGenomeAnnotation.loc[refGenomeAnnotation["Gene_Region"].isin(["start_codon", "stop_codon"])]

    TFsdf = refGenomeAnnotation[["Gene_ID", "Transcript_ID", "Chromosome", "Gene_Region", "Start", "Stop", "Strand"]].loc[refGenomeAnnotation["Gene_ID"].isin(queryTFsdf["Flybase_ID"])]

    #Add reference genome sequence per gene region
    #This will correspond to 1.6kb upstream and downstream of ATG/stop codon 
    TFsdf = TFsdf.assign(Reference_Seq = "")

    for index, rowcontents in TFsdf.iterrows():
            #Define 3.2kb gene region
            regionStart = rowcontents["Start"] - 1601
            regionStop = rowcontents["Stop"] + 1600

            #Add reference sequence
            refPosStrandSeq = str(refSeqPerChromosome[rowcontents["Chromosome"]][regionStart:regionStop]) #This is the + strand seq, so goes from end to beginning
            TFsdf.at[index,"Reference_Seq"] = revComp(refPosStrandSeq)

    #Create fragments for the HDR-arm
    TFsdf = make_homology_arm_fragments(TFsdf)

    #re-index
    TFsdf = TFsdf.reset_index()
    del TFsdf["index"]
           
    #Create a dictionary of dictionaries for the df
    TFsdict_of_dict = {}
    for index, rowcontents in TFsdf.iterrows():
        TFsdict_of_dict[index] = rowcontents.to_dict()

    return TFsdf, TFsdict_of_dict

def filter_gRNA(gRNA_file, TF_dict, refSeqPerChromosome = "inputfiles/dmel-all-chromosome-r6.48.fasta"):
    '''
    params:
        gRNA_file: gff file
        TF_dict: dataframe with the following format:

            {    'Gene_ID': 'FBgn0004652',
                'Transcript_ID': 'FBtr0083651', 
                'Chromosome': '3R', 
                'Gene_Region': 'stop_codon', 
                'Start': 18426145, 
                'Stop': 18426147, 
                'Strand': '-', 
                'Reference_Seq': 'TTGATCGTAGGACAC', 
                'Transgenic_Seq': 'GACCCTAGGACCGG'
            }
        returns: dictionary with the following format

            df={
            "start/stop":"start",#is it N or C termini -> do we need to look at start or stop codon for teh cut 
            'genome_start_codon_pos':400, 
            'genome_stop_codon_pos':700, # or only 1 of those 
            'strand_type':'+',
            "sgRNA_list_positions":[[401,425],[456,467],[478,489],[395,415]],#those wil be as genome positions -assumptions - the coordinates correspond to the 1st and last bp of the strand to which the gsRNA will be complementary to
            "sgRNA_list_values":["AAGCGACTA","AAAAAAAATAAAAA","ATATATTTTTTTTTTAAAAA","AGCGCGAAATAATA"]
            "sgRNA_strand" = ['-']
        }
    '''
    import pandas as pd

    # create a dataframe from the gRNA files
    gRNAFileAnnotation = pd.read_csv(gRNA_file, sep = "\t", index_col = False)

    # add a new category to the dataframe that provides information about whether the sequence deviates from the transgenic strain 
    gRNAFileAnnotation = gRNAFileAnnotation.assign(target_site_variation= "")

    # reformat the "Attribute" category in refGenomeAnnotation, to extract Gene_ID, Gene_Symbol, and Transcript ID
    index = 0 #TODO@Marina improve this code

    # for each attribute value, extract the gene ID and symbol and add this to the new categories
    for attribute in gRNAFileAnnotation['attributes']:

        fullatt = (attribute).split(";")
        gRNAFileAnnotation.at[index,"target_site_variation"] = fullatt[8]
        index+=1
    
    # shorten file to essential information
    GenomeCoordinates = gRNAFileAnnotation.loc[:,["target_site_variation", "fmin", "fmax", "#chr", "strand"]]

    # initialize a list that stores the positions of the gRNAs
    sgRNA_list_positions = []
    sgRNA_list_values = []

    # check whether the sgRNAs match the transgeneic genome, whether the sgRNAs match to the same chromosome as the transcription factor
    # select sgRNA that are located maximally 20 bp upstream of the start/stop codon of the transcription factors
    # select sgRNA that are located maximally 20 bp downstream of the start/stop codon of the transcription factors

    for gRNA in range(len(GenomeCoordinates)):
        
        if GenomeCoordinates['target_site_variation'].iloc[gRNA] == "target_site_variation=none observed" and \
        GenomeCoordinates['#chr'].iloc[gRNA] == TF_dict["Chromosome"] and \
        int(GenomeCoordinates['fmin'].iloc[gRNA])-1 >= TF_dict["Start"] - 20 and \
        int(GenomeCoordinates['fmax'].iloc[gRNA]) <= TF_dict["Stop"] + 20:
                
            sgRNA_list_positions.append([int(GenomeCoordinates['fmin'].iloc[gRNA])-1, int(GenomeCoordinates['fmax'].iloc[gRNA])])
        
            # to retrieve the nucleotide sequence of the sgRNA, load the fruit fly reference genome and extract the sequence at the genomic coordinates provided in the sgRNA file

            gRNA_string = str(refSeqPerChromosome[GenomeCoordinates['#chr'].iloc[gRNA]][int(GenomeCoordinates['fmin'].iloc[gRNA])-1:int(GenomeCoordinates['fmax'].iloc[gRNA])])

            sgRNA_list_values.append(gRNA_string)

    TF_dict["sgRNA_list_positions"] = sgRNA_list_positions
    TF_dict["sgRNA_list_values"] = sgRNA_list_values
    TF_dict["sgRNA_strand"] = GenomeCoordinates['strand'].iloc[gRNA]


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
            "sgRNA_list_positions":[[401,425],[456,467],[478,489],[395,415]],#those wil be as genome positions -assumptions - the coordinates correspond to the 1st and last bp of the strand to which the gsRNA will be complementary to
            "sgRNA_list_values":["AAGCGACTA","AAAAAAAATAAAAA","ATATATTTTTTTTTTAAAAA","AGCGCGAAATAATA"]
        }

    returns: dictionary with the following format (only mock not real gene)

        df={
            "start/stop":"start",#is it N or C termini -> do we need to look at start or stop codon for teh cut 
            'genome_start_codon_pos':400, 
            'genome_start_codon_pos':700, # or only 1 of those 
            'strand_type':'+',
            "sgRNA_list_positions":[[401,425],[456,467],[478,489],[395,415]],#those wil be as genome positions -assumptions - the coordinates correspond to the 1st and last bp of the strand to which the gsRNA will be complementary to
            "sgRNA_list_values":["AAGCGACTA","AAAAAAAATAAAAA","ATATATTTTTTTTTTAAAAA","AGCGCGAAATAATA"]
        }

    '''
    # loop through list of sgRNAs

    for i in range(0,len(df["sgRNA_list_values"])):
        
        cand=df["sgRNA_list_values"][i] #candidate sgRNA
        cand_len=len(cand) #candidate sgRNA lenght
        
        # if the region is a start codon, then (if possible) select the sgRNA where the ATG is part of the PAM

        if df["start/stop"]=="start_codon":

            if cand[(cand_len-4):cand_len]=="ATGG":

                # return a dictionary with the sequence and the position of the selected sgRNA and the note that the PAM must not be mutated in the HDR plasmid

                df["sgRNA_list_values"]=df["sgRNA_list_values"][i]
                df["sgRNA_list_positions"]=df["sgRNA_list_positions"][i]
                df["must_PAM_be_mutated_in_HDR_plasmid?"] = "no"
                
                return df
        
        # if the region is a stop codon, then (if possible) select the sgRNA where the TAG is part of the PAM

        else:

            if cand[(cand_len-4):cand_len]=="TAGG": 

                # return a dictionary with the sequence and the position of the selected sgRNA and the note that the PAM must not be mutated in the HDR plasmid

                df["sgRNA_list_values"]=df["sgRNA_list_values"][i]
                df["sgRNA_list_positions"]=df["sgRNA_list_positions"][i]
                df["must_PAM_be_mutated_in_HDR_plasmid?"] = "no"
                
                return df
    
    # if looping through all of the given sgRNAs is finished without having encountered a best case sgRNA, return None

    return

def check_over_15(df):
    '''
    If the best case scenario does not hold true, we are looking for sgRNAs which span the start/stop codon,
    so that max. 15 bp of the 3´ end (incl. PAM) are on one side.
    
    params: dictionary with the following format (only mock not real gene)

        df={
            "start/stop":"start",#is it N or C termini -> do we need to look at start or stop codon for teh cut 
            'genome_start_codon_pos':400, 
            'genome_stop_codon_pos':700, # or only 1 of those 
            'strand_type':'+',
            "sgRNA_list_positions":[[401,425],[456,467],[478,489],[395,415]],#those wil be as genome positions -assumptions - the coordinates correspond to the 1st and last bp of the strand to which the gsRNA will be complementary to
            "sgRNA_list_values":["AAGCGACTA","AAAAAAAATAAAAA","ATATATTTTTTTTTTAAAAA","AGCGCGAAATAATA"]
        }
    returns: same dictionary as in params but only with sgRNAs listed in "sgRNA_list_positions" and "sgRNA_list_values" that fulfil the condition
    
    '''
    start_pos = df["genome_start_codon_pos"]
    stop_pos = df["genome_stop_codon_pos"]

    list_of_max_15=[]
    list_of_max_15_pos=[]
    list_of_more_than_15=[]
    list_of_more_than_15_pos=[]

    # loop through the list of sgRNAs provided

    for i in range(0, len(df["sgRNA_list_values"])):

        # check whether the overhang of the sgRNA 3´end from the start codon is max. 15 bp
        # if yes, store the sequence and position of the sgRNA in the respective list
        # if no, store the sequence and position of the sgRNA in the respective list

        if df["start/stop"]=="start_codon":

            if 0 < df["sgRNA_list_positions"][i][1] - 2 - start_pos < 16:
                
                list_of_max_15.append(df["sgRNA_list_values"][i])
                list_of_max_15_pos.append(df["sgRNA_list_positions"][i])

            else:
                
                list_of_more_than_15.append(df["sgRNA_list_values"][i])
                list_of_more_than_15_pos.append(df["sgRNA_list_positions"][i])
        
        # check whether the overhang of the sgRNA 3´end from the stop codon is max. 15 bp
        # if yes, store the sequence and position of the sgRNA in the respective list
        # if no, store the sequence and position of the sgRNA in the respective list

        else:

            if 0 < df["sgRNA_list_positions"][i][1] - 2 - stop_pos < 16:

                list_of_max_15.append(df["sgRNA_list_values"][i])
                list_of_max_15_pos.append(df["sgRNA_list_positions"][i])

            else:
                
                list_of_more_than_15.append(df["sgRNA_list_values"][i])
                list_of_more_than_15_pos.append(df["sgRNA_list_positions"][i])
    
    # after looping through all of the provided sgRNAs, (if possible,) return the dictionary with the sgRNAs with max. 15 bp overhang on the 3´prime end of the start/ stop codon
    # if possible, add a note that the PAM does not need to be mutated in the HDR plasmid
    
    if list_of_max_15:

        df["sgRNA_list_values"] = list_of_max_15
        df["sgRNA_list_positions"] = list_of_max_15_pos
        df["must_PAM_be_mutated_in_HDR_plasmid?"] = "no"
    
    # if no sgRNA with a max. overhang of 15 bp from the 3`end of the start/ stop codon was found, return a dictionary with the sgRNAs with a larger overhang
    else: 

        df["sgRNA_list_values"] = list_of_more_than_15
        df["sgRNA_list_positions"] = list_of_more_than_15_pos
        df["must_PAM_be_mutated_in_HDR_plasmid?"] = "yes"

    return df

def select_closest(df):
    """
    Pick the sgRNA which has the shortest distance between the cutting site (6bp from the 3´end of the sgRNA) and the stop or start codon.

    params: dictionary with the following format (only mock not real gene)

        df={
            "start/stop":"start",#is it N or C termini -> do we need to look at start or stop codon for the cut 
            'genome_start_codon_pos':400, 
            'genome_stop_codon_pos':700, # or only 1 of those 
            'strand_type':'+',
            "sgRNA_list_positions":[[401,425],[456,467],[478,489],[395,415]],#those wil be as genome positions -assumptions - the coordinates correspond to the 1st and last bp of the strand to which the gsRNA will be complementary to
            "sgRNA_list_values":["AAGCGACTA","AAAAAAAATAAAAA","ATATATTTTTTTTTTAAAAA","AGCGCGAAATAATA"]
        }

    returns: same dictionary as above storing only the sequence and position of the sgRNA with the shortest distance between the cutting site and the start/ stop codon.

    """
    from os import sys

    # make some empty lists of strings of the right length without any values in them

    list_of_pos=len(df["sgRNA_list_values"])*[""]
    list_of_dist=len(df["sgRNA_list_values"])*[""]
    list_of_val=len(df["sgRNA_list_values"])*[""]

    # loop through the list of sgRNAs, calculate and store the distance from the start/stop codon

    for i in range(0,len(df["sgRNA_list_values"])):

        cutting_pos = df["sgRNA_list_positions"][i][1]-6
        
        if df["start/stop"]=="start_codon":

            start_codon_pos=df["genome_start_codon_pos"]

            dist=abs(cutting_pos-start_codon_pos)

        else:
            
            stop_codon_pos=df["genome_stop_codon_pos"]

            dist=abs(cutting_pos-stop_codon_pos)

        list_of_dist[i]=dist
        list_of_val[i]=df["sgRNA_list_values"][i]
        list_of_pos[i]=df["sgRNA_list_positions"][i]
    
    # determine the shortest distance of cutting site from start/ stop codon

    smallest=min(list_of_dist)

    # make list of sgRNAs that have the shortest distance from start/ stop

    indexes_smallest_distance = [i for i, x in enumerate(list_of_dist) if x == smallest]

    # if there is only one sgRNA that has the shortest distance from start/ stop, store the sequence and position of this sgRNA
    # return a dictionary for this sgRNA with the shortest distance
    
    if 0 < len(indexes_smallest_distance) <= 2:

        df["sgRNA_list_values"]=list_of_val[indexes_smallest_distance[0]]
        df["sgRNA_list_positions"]=list_of_pos[indexes_smallest_distance[0]]

    elif len(indexes_smallest_distance) == 0:

        print("Error! There was no sgRNA found that is the closest to the cutting site.")

        sys.exit()

        df = {}

    elif len(indexes_smallest_distance) > 2:

        print("Error! More than two sgRNA were found that are the closest to the cutting site.")

        sys.exit()

        df = {}

    # if there is two sgRNAs with the shortest distance between the start/ stop codon and the cutting site, print an error
    # before running the function select_closest(), the correct side of the cutting site relative to the start/ stop codon should be chosen

    return df

def check_cutting_site_inside_CDS(df):

    """
    Pick the sgRNA which whose cutting site (6bp from the 3´end of the sgRNA) lies within the gene.
    This is desirable because it allows us to design a synonymus mutation in the PAM of the HDR plasmid knowing, without changing the peptide sequence,
    while we are not sure what consequences a mutation in the UTR would have. 
    
    params: dictionary with the following format (only mock not real gene)

        df={
            "start/stop":"start",#is it N or C termini -> do we need to look at start or stop codon for teh cut 
            'genome_start_codon_pos':400, 
            'genome_stop_codon_pos':700, # or only 1 of those 
            'strand_type':'+',
            "sgRNA_list_positions":[[401,425],[456,467],[478,489],[395,415]],#those wil be as genome positions -assumptions - the coordinates correspond to the 1st and last bp of the strand to which the gsRNA will be complementary to
            "sgRNA_list_values":["AAGCGACTA","AAAAAAAATAAAAA","ATATATTTTTTTTTTAAAAA","AGCGCGAAATAATA"]
        }

    returns: dictionary with the same format as the params with the sequence and the positions of the most suitable sgRNA in "sgRNA_list_values" and "sgRNA_list_positions" respectively.

    """
    # make some empty list of strings of the right length without any values in them

    cut_outside_of_CDS=len(df["sgRNA_list_values"])*[""]

    # loop through the list of sgRNAs, calculate whether the cutting site of the sgRNAs would be inside or outside the CDS
    # store whether the sgRNA cuts outside of the CDS or not

    for i in range(0,len(df["sgRNA_list_values"])):

        cutting_pos = df["sgRNA_list_positions"][i][1]-6
        
        if df["start/stop"]=="start_codon":

            start_codon_pos=df["genome_start_codon_pos"]

            if cutting_pos < start_codon_pos:

                cut_outside_of_CDS[i] = "yes"
            
            else: cut_outside_of_CDS[i] = "no"

        else:
            
            stop_codon_pos=df["genome_stop_codon_pos"]

            if cutting_pos < stop_codon_pos:

                cut_outside_of_CDS[i] = "no"
            
            else: cut_outside_of_CDS[i] = "yes"
    
    # check whether one of the sgRNAs cuts inside the CDS

    if any(x == "no" for x in cut_outside_of_CDS):
    
        df["sgRNA_list_values"]=[df["sgRNA_list_values"][i] for i, x in enumerate(cut_outside_of_CDS) if x == "no"]
        df["sgRNA_list_positions"]=[df["sgRNA_list_positions"][i] for i, x in enumerate(cut_outside_of_CDS) if x == "no"]
        df["cut_outside_of_CDS"]="no"

    else:

        df["cut_outside_of_CDS"]="yes"
    
    return df

def make_homology_arm_fragments(tfsDF):
    """
    Takes in the dataframe of information about start/stop codon regions, 
    and appends with columns for the 225 bp upstream and downstream of
    this site (the homlogy arm fragments).

    input: tfsDF - as defined in previous function and on GitHub README, with Reference_Seq 1700bp either side of the start/stop.
    output: tfsDF appended with upstreamHA and downstreamHA

    """

    tfsDF["upstreamHA"] = tfsDF.Reference_Seq.str[1375:1600]
    tfsDF["downstreamHA"] = tfsDF.Reference_Seq.str[1604:1829]

    return tfsDF

def find_synonymous_codons(query_codon, codon_table_excel):

    '''
    Uses the amino acids table to select codons that encode for the same amino acid as the query codon.

    Params:
        codon: string, codon to select synonymous codons for
        codon_table_excel: string, path to an excel file that lists per codon which amino acid that codon encodes.
    
    Returns:
        synonymous_codons: list of strings, each string is a codon that encodes for the same amino acid as the query codon.
    
    '''

    import pandas as pd

    synonymous_codons = []

    codon_table = pd.read_excel(codon_table_excel)

    amino_acid = codon_table.query(f"codon =='{query_codon}'")

    amino_acid_query = amino_acid.iloc[0]['amino_acid']

    codons = codon_table.query(f"amino_acid=='{amino_acid_query}' & codon !='{query_codon}'")

    for codon in range(0,len(codons.index)):

        synonymous_codon = codons.iloc[codon]['codon']
        synonymous_codons.append(synonymous_codon)

    return synonymous_codons

def mutate_PAM_in_codon(query_codon, synonymous_codons, sgRNA_strand):

    '''
    Mutates one of the G/C`s from the PAM in the query codon to another nucleotide based on the codon list provided. 
    Returns an empty list if mutation was not possible.

    params: 
        query_codon: string, codon containing the first G of the PAM.
        synonymous_codons: list of strings, codons that incode for the same amino acid as the query codon.

    returns:
        selected_codon: string, codon that is the same as the query codon apart from one G, which also encodes for the same amino acid as the query codon.
    '''
    if synoynmous_codons:

    if sgRNA_strand == '+' # in the case that the sgRNA is lying on the + strand, the part of the fragment corresponding to the PAM will be NGG

        list_query_codon = list(query_codon)
    
        if list_query_codon[2] == 'G':

            for synonymous_codon in synonymous_codons:

                list_synonymous_codon = list(synonymous_codon)

                if list_synonymous_codon[2] != 'G' and listsynonymous_codon != 'A': # this is because the PAM should not be mutated into NGA

                    selected_codon = synonymous_codon

                    break

                else: selected_codon = ''
        
        elif list_query_codon[0] == 'G':

            for synonymous_codon in synonymous_codons:

                list_synonymous_codon = list(synonymous_codon)

                if list_synonymous_codon[0] != 'G':

                    selected_codon = synonymous_codon

                    break

                else: selected_codon = ''

        else: 

            print('There is no synonymous codon that can be used to mutate the PAM. Searching for other sgRNAs.')

            selected_codon = ''

    else: 

        print('There is no synonymous codon that can be used to mutate the PAM. Searching for other sgRNAs.')

        selected_codon = ''
        
    return selected_codon

def make_synonymous_mutation(sequence, position_of_mutation, codon_table_excel = 'inputfiles/codon_table.xlsx', sgRNA_strand):

    '''
    Mutates a G/C from the PAM in a sequence in a way, that the codon containing the G still encodes for the same amino acid.
    The genetic code does not allow for mutating away middle G´s, it does only allow for mutating away G´s located in the wobble position.

    params:
        sequence: string, containing the sequence that is supposed to be mutated.
        positions_of_mutation: int, position of nucleotide that is supposed to be mutated in sequence starting with 0.
    
    returns: string, mutated sequence.
    '''

    # a codon is 3 nucleotides long
    x = 3

    list_of_codons=[sequence[y-x:y] for y in range(x, len(sequence)+x,x)]

    codon_position = translate_nucleotide_position_into_codon_position(sequence, position_of_mutation)

    codon_to_mutate = list_of_codons[codon_position]

    synonymous_codons = find_synonymous_codons(codon_to_mutate, codon_table_excel)

    selected_codon = mutate_PAM_in_codon(codon_to_mutate, synonymous_codons, sgRNA_strand)

    if selected_codon: 

        list_of_codons[codon_position] = selected_codon

        mutated_sequence = "".join(list_of_codons)

    else: mutated_sequence = ""

    return mutated_sequence

def translate_nucleotide_position_into_codon_position(sequence, nucleotide_position):

    '''
    params: 
        sequence: string, nucleotide sequence
        nucleotide_position: integer, position of nucleotide that you want to translate into the codon position

    returns: 
        count: integer, position of codon containing nucleotide from nucleotide position
    '''

    count = 0

    # since position of first letter in string is zero, function has to be adapted accordingly

    for n in range(0,len(sequence)):

        if nucleotide_position - 3 >= 0:

            count = count +1

            nucleotide_position = nucleotide_position - 3

        else:

            break

    return(count)

def check_PAM_inside_CDS(df):

    """
    Check whether PAM (3bp from the 3´end of the sgRNA) lies within the CDS.
    This needs to be known to figure out whether it would be best to mutate the PAM or the recognition site.
    We would rather mutate something in the CDS because we know that synonymous mutations should not destroy
    protein function wheareas we are not sure if mutations might destroy relevant motifs in the UTRs. 
    
    params: dictionary with the following format (only mock not real gene)

        df={
            'gene_ID': 'FBgn0004652', 
            'transcript_ID': 'FBtr0083648', 
            'chromosome': '3R', 
            'start/stop': 'start_codon', 
            'strand_type': '-', 
            'genome_start_codon_pos': 18506517, 
            'genome_stop_codon_pos': 'N/A', 
            'sgRNA_list_positions': [18506502, 18506525], 
            'sgRNA_list_values': 'GAGTTCTTTTGCAGCATTTATGG', 
            'must_PAM_be_mutated_in_HDR_plasmid?': 'no'
            }

    returns: bolean, 'yes' or 'no'

    """
    if df['sgRNA_strand'] == '+':

        if df['strand_type'] == '+' and df["start/stop"]=="start_codon":

            start_codon_pos = df["genome_start_codon_pos"]

            pos_of_interest = df["sgRNA_list_positions"][1]-2

            if start_codon_pos <= pos_of_interest:

                pos_inside_CDS = "yes"

            elif start_codon_pos > pos_of_interest:

                pos_inside_CDS = "no"

        elif df['strand_type'] == '+' and df["start/stop"]=="stop_codon":

            stop_codon_pos = df["genome_stop_codon_pos"]

            pos_of_interest = df["sgRNA_list_positions"][1]

            if stop_codon_pos >= pos_of_interest:

                pos_inside_CDS = "yes"

            elif stop_codon_pos < pos_of_interest:

                pos_inside_CDS = "no"
        
        elif df['strand_type'] == '-' and df['start/stop'] == 'start_codon':

            start_codon_pos = df['genome_start_codon_pos']

            pos_of_interest = df['sgRNA_list_position'][1]

            if start_codon >= pos_of_interest:

                pos_inside_CDS = "yes"

            elif start_codon_pos < pos_of_interest:

                pos_inside_CDS = "no"

        elif df['strand_type'] == '-' and df['start/stop'] == 'stop_codon':

            stop_codon_pos = df['genome_stop_codon_pos']

            pos_of_interest = df['sgRNA_list_position'][1]-2

            if stop_codon <= pos_of_interest:

                pos_inside_CDS = "yes"

            elif stop_codon > pos_of_interest:

                pos_inside_CDS = "no"
        
    elif df['sgRNA_strand'] == '-':

        if df['strand_type'] == '+' and df["start/stop"]=="start_codon":

            start_codon_pos = df["genome_start_codon_pos"]

            pos_of_interest = df["sgRNA_list_positions"][0]

            if start_codon_pos <= pos_of_interest:

                pos_inside_CDS = "yes"

            elif start_codon_pos > pos_of_interest:

                pos_inside_CDS = "no"

        elif df['strand_type'] == '+' and df["start/stop"]=="stop_codon":

            stop_codon_pos = df["genome_stop_codon_pos"]

            pos_of_interest = df["sgRNA_list_positions"][0]+2

            if stop_codon_pos >= pos_of_interest:

                pos_inside_CDS = "yes"

            elif stop_codon_pos < pos_of_interest:

                pos_inside_CDS = "no"
        
        elif df['strand_type'] == '-' and df['start/stop'] == 'start_codon':

            start_codon_pos = df['genome_start_codon_pos']

            pos_of_interest = df['sgRNA_list_position'][0]+2

            if start_codon >= pos_of_interest:

                pos_inside_CDS = "yes"

            elif start_codon_pos < pos_of_interest:

                pos_inside_CDS = "no"

        elif df['strand_type'] == '-' and df['start/stop'] == 'stop_codon':

            stop_codon_pos = df['genome_stop_codon_pos']

            pos_of_interest = df['sgRNA_list_position'][0]

            if stop_codon <= pos_of_interest:

                pos_inside_CDS = "yes"

            elif stop_codon > pos_of_interest:

                pos_inside_CDS = "no"


    return pos_inside_CDS

def check_sgRNA_recognition_site_inside_CDS(df):

    """
    Check whether sgRNA recognition site (9bp from the 3´end of the sgRNA) lies within the CDS.
    This needs to be known to figure out whether it would be best to mutate the PAM or the recognition site.
    We would rather mutate something in the CDS because we know that synonymous mutations should not destroy
    protein function wheareas we are not sure if mutations might destroy relevant motifs in the UTRs. 
    
    params: dictionary with the following format (only mock not real gene)

        df={
            "start/stop":"start",#is it N or C termini -> do we need to look at start or stop codon for teh cut 
            'genome_start_codon_pos':400, 
            'genome_stop_codon_pos':700, # or only 1 of those 
            'strand_type':'+',
            "sgRNA_list_positions":[[401,425],[456,467],[478,489],[395,415]],#those wil be as genome positions -assumptions - the coordinates correspond to the 1st and last bp of the strand to which the gsRNA will be complementary to
            "sgRNA_list_values":["AAGCGACTA","AAAAAAAATAAAAA","ATATATTTTTTTTTTAAAAA","AGCGCGAAATAATA"]
        }

    returns: bolean, 'yes' or 'no'

    """
    if df['sgRNA_strand'] == '+':

        if df['strand_type'] == '+' and df["start/stop"]=="start_codon":

            start_codon_pos = df["genome_start_codon_pos"]

            pos_of_interest = df["sgRNA_list_positions"][1]-2-6

            if start_codon_pos <= pos_of_interest:

                pos_inside_CDS = "yes"

            elif start_codon_pos > pos_of_interest:

                pos_inside_CDS = "no"

        elif df['strand_type'] == '+' and df["start/stop"]=="stop_codon":

            stop_codon_pos = df["genome_stop_codon_pos"]

            pos_of_interest = df["sgRNA_list_positions"][1]-2

            if stop_codon_pos >= pos_of_interest:

                pos_inside_CDS = "yes"

            elif stop_codon_pos < pos_of_interest:

                pos_inside_CDS = "no"
        
        elif df['strand_type'] == '-' and df['start/stop'] == 'start_codon':

            start_codon_pos = df['genome_start_codon_pos']

            pos_of_interest = df['sgRNA_list_position'][1]-2

            if start_codon >= pos_of_interest:

                pos_inside_CDS = "yes"

            elif start_codon_pos < pos_of_interest:

                pos_inside_CDS = "no"

        elif df['strand_type'] == '-' and df['start/stop'] == 'stop_codon':

            stop_codon_pos = df['genome_stop_codon_pos']

            pos_of_interest = df['sgRNA_list_position'][1]-2-6

            if stop_codon <= pos_of_interest:

                pos_inside_CDS = "yes"

            elif stop_codon > pos_of_interest:

                pos_inside_CDS = "no"
        
    elif df['sgRNA_strand'] == '-':

        if df['strand_type'] == '+' and df["start/stop"]=="start_codon":

            start_codon_pos = df["genome_start_codon_pos"]

            pos_of_interest = df["sgRNA_list_positions"][0]+2

            if start_codon_pos <= pos_of_interest:

                pos_inside_CDS = "yes"

            elif start_codon_pos > pos_of_interest:

                pos_inside_CDS = "no"

        elif df['strand_type'] == '+' and df["start/stop"]=="stop_codon":

            stop_codon_pos = df["genome_stop_codon_pos"]

            pos_of_interest = df["sgRNA_list_positions"][0]+2+6

            if stop_codon_pos >= pos_of_interest:

                pos_inside_CDS = "yes"

            elif stop_codon_pos < pos_of_interest:

                pos_inside_CDS = "no"
        
        elif df['strand_type'] == '-' and df['start/stop'] == 'start_codon':

            start_codon_pos = df['genome_start_codon_pos']

            pos_of_interest = df['sgRNA_list_position'][0]+2+6

            if start_codon >= pos_of_interest:

                pos_inside_CDS = "yes"

            elif start_codon_pos < pos_of_interest:

                pos_inside_CDS = "no"

        elif df['strand_type'] == '-' and df['start/stop'] == 'stop_codon':

            stop_codon_pos = df['genome_stop_codon_pos']

            pos_of_interest = df['sgRNA_list_position'][0]+2

            if stop_codon <= pos_of_interest:

                pos_inside_CDS = "yes"

            elif stop_codon > pos_of_interest:

                pos_inside_CDS = "no"

    return pos_inside_CDS

def mutate_PAM_in_HDR_plasmid_old(HAL_R, HAR_F, df):

    # TODO: delete this as soon as Emma has copied the old parts that are still relevant into her new function

    from os import sys

    codon_table_excel = "inputfiles/codon_table.xlsx"

    if df["start/stop"] == "start_codon":

        pos_of_interest = df["genome_start_codon_pos"]

    if df["start/stop"] == "stop_codon":

        pos_of_interest = df["genome_stop_codon_pos"]

    if df['sgRNA_strand'] == df['strand_type']:

        PAM_pos = df["sgRNA_list_positions"][1]

        # check whether PAM in HAL-R

        if PAM_pos < pos_of_interest:

            # calculate distance between end of HAL-R and end of PAM

            PAM_pos_in_HAL_R = PAM_pos - pos_of_interest
        
            # sanity check whether last nucleotide of PAM is a G

            if HAL_R[PAM_pos_in_HAL_R] != "G" :

                # if there was no PAM found there is a bug and the program should stop running 
                
                print(HAL_R, df)

                print("Error! The PAM position in HAL-R is not correct")

                sys.exit()

            # try mutating the PAM in the HDR fragment
            mutated_HAL_R = make_synonymous_mutation(HAL_R, PAM_pos, codon_table_excel)

        elif PAM_pos-2 > pos_of_interest + 2:



            # if you were able to mutate the PAM in the primer/ fragment sequence, add mutated primer 
            # and information about which primer was mutated to the dataframe
            if mutated_HAL_R:

                df["HAL-R"] = mutated_HAL_R
                df["HAR-F"] = HAR_F
                df["mutated?"] = "HAL_R"

            # if you were not successful in mutating the PAM in the primer/ fragment sequence,
            # try mutating the sgRNA recognition sequence in the primer/ fragement sequence
            else: 

                # calculate the distance between the start/ stop condon and the end of the PAM 
                # to determine where the sgRNA recognition site is located in the primer/ fragement sequence
                distance = pos_of_interest - df["sgRNA_list_positions"][1] - 2

                mutated_HAL_R = mutate_sgRNA_recognition_site_in_HDR_plasmid('HAL_R', HAL_R, distance)
                
                if mutated_HAL_R:
                    
                    #if mutation was successful,  update the sgRNA dictionary with mutated primers 
                    # and add information about mutation status
                    df["HAR-F"] = HAR_F
                    df["HAL-R"] = mutated_HAL_R
                    df["mutated?"] = "HAL_R"

                else:

                    # if mutation was not successful, add information about mutation status 
                    # and keep old unmutated primer
                    df["HAR-F"] = HAR_F
                    df["HAL-R"] = HAL_R
                    df["mutated?"] = "no"

    # check whether PAM in HAR-F
    
    elif df["sgRNA_list_positions"][1] - 2 - pos_of_interest >= 16:

        # sanity check since funtion x.find() returns -1 if string is not found

        PAM_pos = HAR_F.find("GG")

        if PAM_pos == -1:

            print("Error! No GG found in HAR-F even though PAM should be located in HAR-F")

            sys.exit()

        else: 
            
            # try mutating the PAM in the primer/ fragment sequence
            mutated_HAR_F = make_synonymous_mutation(HAR_F, PAM_pos, codon_table_excel)

            # if you were able to mutate the PAM in the primer sequence, add mutated primer 
            # and information about which primer was mutated to the dataframe

            if mutated_HAR_F:

                df["HAR-F"] = mutated_HAR_F
                df["HAL-R"] = HAL_R
                df["mutated?"] = "HAR_F"

            # if there was no synonymous mutation for the PAM, mutate the sgRNA recognition site

            else:

                # calculate distance from the beginning of the PAM to the beginning of the start/ stop codon
                distance = df["sgRNA_list_positions"][1] - 2 - pos_of_interest - 3

                # try to mutate sgRNA recognition site in primer
                mutated_HAR_F = mutate_sgRNA_recognition_site_in_HDR_plasmid('HAR_F', HAR_F, distance)

                if mutated_HAR_F:
                    
                    #if mutation was successful,  update the sgRNA dictionary with mutated primers 
                    # and add information about mutation status
                    df["HAR-F"] = mutated_HAR_F
                    df["HAL-R"] = HAL_R
                    df["mutated?"] = "HAR_F"
                
                else:

                    # if mutation was not successful, add information about mutation status 
                    # and keep old unmutated primer
                    df["HAR-F"] = HAR_F
                    df["HAL-R"] = HAL_R
                    df["mutated?"] = "no"
                
    return df

def mutate_PAM_in_HDR_plasmid_new(HAL_R, HAR_F, df):

    from os import sys

    codon_table_excel = "inputfiles/codon_table.xlsx"

    if df["start/stop"] == "start_codon":

        pos_of_interest = df["genome_start_codon_pos"]

    if df["start/stop"] == "stop_codon":

        pos_of_interest = df["genome_stop_codon_pos"]

    if df['strand_type'] == "+" and df['sgRNA_strand'] == '+':

        PAM_pos = df["sgRNA_list_positions"][1]

        # check whether PAM in HAL-R

        if PAM_pos < pos_of_interest:

            # calculate distance between end of HAL-R and end of PAM

            PAM_pos_in_HAL_R = PAM_pos - pos_of_interest
        
            # sanity check whether last nucleotide of PAM is a G

            if HAL_R[PAM_pos_in_HAL_R] != "G" :

                # if there was no G found in the last nucleotide of the PAM, there is a bug somwhere and the program should stop running 
                
                print(HAL_R, df)

                print("Error! The PAM position in HAL-R is not correct")

                sys.exit()

            # try mutating the PAM in the HDR fragment
            mutated_HAL_R = make_synonymous_mutation(HAL_R, PAM_pos, codon_table_excel, df['sgRNA_strand'])

        # check whether PAM is in HAR-F

        elif PAM_pos - 2 > pos_of_interest + 2:

            # calculate distance between end of HAR-F and end of PAM

            PAM_pos_in_HAR_F = PAM_pos - (pos_of_interest + 2)
        
            # sanity check whether last nucleotide of PAM is a G

            if HAR_F[PAM_pos_in_HAR_F - 1] != "G" : # the minus 1 was added here because the first letter of a string is indexed as zero in python

                # if there was no G found in the last nucleotide of the PAM, there is a bug somwhere and the program should stop running 
                
                print(HAR_F, df)

                print("Error! The PAM position in HAR-F is not correct")

                sys.exit()

            # try mutating the PAM in the HDR fragment
            mutated_HAR_F = make_synonymous_mutation(HAR_F, PAM_pos, codon_table_excel, df['sgRNA_strand'])
    
    elif df['strand_type'] == '+' and df['sgRNA_strand'] == '-':
        
        PAM_pos = df["sgRNA_list_positions"][0]

        # check whether PAM in HAL-R

        if PAM_pos + 2 < pos_of_interest:

            # calculate distance between end of HAL-R and end of PAM

            PAM_pos_in_HAL_R = PAM_pos - pos_of_interest
        
            # sanity check whether last nucleotide of PAM is a G

            if HAL_R[PAM_pos_in_HAL_R] != "G" :

                # if there was no G found in the last nucleotide of the PAM, there is a bug and the program should stop running 
                
                print(HAL_R, df)

                print("Error! The PAM position in HAL-R is not correct")

                sys.exit()

            # try mutating the PAM in the HDR fragment
            mutated_HAL_R = make_synonymous_mutation(HAL_R, PAM_pos, codon_table_excel, df['sgRNA_strand'])

        # check whether PAM is in HAR-F

        elif PAM_pos > pos_of_interest + 2:

            # calculate distance between end of HAR-F and end of PAM

            PAM_pos_in_HAR_F = PAM_pos - (pos_of_interest + 2)
        
            # sanity check whether last nucleotide of PAM is a G

            if HAR_F[PAM_pos_in_HAR_F - 1] != "G" : # the minus 1 was added here because the first letter of a string is indexed as zero in python

                # if there was no G found in the last nucleotide of the PAM, there is a bug somwhere and the program should stop running 
                
                print(HAR_F, df)

                print("Error! The PAM position in HAR-F is not correct")

                sys.exit()

            # try mutating the PAM in the HDR fragment
            mutated_HAR_F = make_synonymous_mutation(HAR_F, PAM_pos, codon_table_excel, df['sgRNA_strand'])

    elif df['strand_type'] == '-' and df['sgRNA_strand'] == '+':

        PAM_pos = df["sgRNA_list_positions"][1]

        # check whether PAM in HAL-R

        if PAM_pos - 2 > pos_of_interest + 2:

            # calculate distance between end of HAL-R and end of PAM

            PAM_pos_in_HAL_R = PAM_pos - (pos_of_interest + 2)
        
            # sanity check whether last nucleotide of PAM is a G

            if HAL_R[PAM_pos_in_HAL_R - 1] != "G" : # the minus 1 was added here because the first letter of a string is indexed as zero in python

                # if there was no PAM found there is a bug and the program should stop running 
                
                print(HAL_R, df)

                print("Error! The PAM position in HAL-R is not correct")

                sys.exit()

            # try mutating the PAM in the HDR fragment
            mutated_HAL_R = make_synonymous_mutation(HAL_R, PAM_pos, codon_table_excel, df['sgRNA_strand'])

        # check whether PAM is in HAR-F

        elif PAM_pos < pos_of_interest:

            # calculate distance between end of HAL-R and end of PAM

            PAM_pos_in_HAR_F = PAM_pos - pos_of_interest
        
            # sanity check whether last nucleotide of PAM is a G

            if HAR_F[PAM_pos_in_HAR_F] != "G" :

                # if there was no G found in the last nucleotide of the PAM, there is a bug and the program should stop running 
                
                print(HAR_F, df)

                print("Error! The PAM position in HAL-R is not correct")

                sys.exit()

            # try mutating the PAM in the HDR fragment
            mutated_HAR_F = make_synonymous_mutation(HAR_F, PAM_pos, codon_table_excel, df['sgRNA_strand'])    

    elif df['strand_type'] == '-' and df['sgRNA_strand'] == '-':

        PAM_pos = df["sgRNA_list_positions"][0]

        # check whether PAM in HAL-R

        if PAM_pos > pos_of_interest + 2:

            # calculate distance between end of HAL-R and end of PAM

            PAM_pos_in_HAL_R = PAM_pos - pos_of_interest + 2
        
            # sanity check whether last nucleotide of PAM is a G

            if HAL_R[PAM_pos_in_HAL_R - 1] != "G" : # the minus 1 was added here because the first letter of a string is indexed as zero in python

                # if there was no PAM found there is a bug and the program should stop running 
                
                print(HAL_R, df)

                print("Error! The PAM position in HAL-R is not correct")

                sys.exit()

            # try mutating the PAM in the HDR fragment
            mutated_HAL_R = make_synonymous_mutation(HAL_R, PAM_pos, codon_table_excel, df['sgRNA_strand'])

        # check whether PAM is in HAR-F

        elif PAM_pos + 2 < pos_of_interest:

            # calculate distance between end of HAL-R and end of PAM

            PAM_pos_in_HAR_F = PAM_pos - pos_of_interest
        
            # sanity check whether last nucleotide of PAM is a G

            if HAR_F[PAM_pos_in_HAR_F] != "G" :

                # if there was no G found in the last nucleotide of the PAM, there is a bug and the program should stop running 
                
                print(HAR_F, df)

                print("Error! The PAM position in HAL-R is not correct")

                sys.exit()

            # try mutating the PAM in the HDR fragment
            mutated_HAR_F = make_synonymous_mutation(HAR_F, PAM_pos, codon_table_excel, df['sgRNA_strand'])    

    # if you were able to mutate the PAM in the fragment sequence, add mutated primer 
    # and information about which primer was mutated to the dataframe
    if mutated_HAL_R:

        df["HAL-R"] = mutated_HAL_R
        df["HAR-F"] = HAR_F
        df["mutated?"] = "HAL_R"

    elif mutated_HAR_F:

        df["HAR-F"] = mutated_HAR_F
        df["HAL-R"] = HAL_R
        df["mutated?"] = "HAR_F"

    # if you were not successful in mutating the PAM in the primer/ fragment sequence,
    # return the primer/ fragments in their unmutated forms and flag them as unmutated
    
    else:
        
        df["HAR-F"] = HAR_F
        df["HAL-R"] = HAL_R
        df["mutated?"] = "no"
                
    return df

def mutate_sgRNA_recognition_site_in_HDR_plasmid(sequenceType, sequenceToMutate, positionCount):
    """
    Mutates HDR arm fragments in the sgRNA recognition site in the 2-3 codons closest to PAM in the case where PAM itself cannot be mutated. 1-2 mutations are made in separate amino acids.

    params:
        sequenceType: 'HAL-R' or 'HAR-F' (corresponding to left or right homology arm)
        sequenceToMutate: the primer/fragment sequence. If reverse primer, this should be the + strand (revComp)
        positionCount: Distance of start/stop from PAM (excluding NGG).

    returns:
        newSequence: The mutated HDR arm. If no mutations were not possible, an empty string will be returned.

    """
    #Define all codons in the homology arm fragment
    allCodons = []
    for codonStart in range(0, len(sequenceToMutate), 3):
        allCodons.append(sequenceToMutate[codonStart:codonStart+3])

    if sequenceType == 'HAL-R':
        allCodons.reverse() #For HAL-R, order codons moving away from te STOP site
    
    #Take the 1-3 codons closest to PAM in order
    if positionCount <= 3: #Can only attempt to mutate one codon if positionCount is less than three
        orderedMutableCodons = {0: allCodons[0]}
    elif positionCount <= 6: #Can only attempt to mutate two codon if positionCount is less than six
        orderedMutableCodons = {0: allCodons[0], 1: allCodons[1]}
    else:
        PAMcodon = int((positionCount-1)/3) #This is the codon number that PAM is in
        if positionCount%3 !=0: #If PAM is not in frame, take the codon that it is in and the two after that
            orderedMutableCodons = {PAMcodon:allCodons[PAMcodon], PAMcodon-1:allCodons[PAMcodon-1], PAMcodon-2:allCodons[PAMcodon-2]} 
        else: #If in frame, just take the 2 codons after PAM
            orderedMutableCodons = {PAMcodon-1:allCodons[PAMcodon-1], PAMcodon-2:allCodons[PAMcodon-2]} 
        
    #This is the count of mutations completed in the sequence.
    mutatedCount = 0
    newSequence = sequenceToMutate #Make copy to avoid mutating original sequence

    #loop through the codons and mutate if possible to a maximum of two mutations.
    for key, codon in orderedMutableCodons.items():
        #find synonymous codons - will provide list of synonymous codons
        codonOptions = find_synonymous_codons(codon)

        #if the list is not empty, select the first synonymous codon and replace this within the mutable region. Add 1 to mutatedCount.
        if codonOptions != []:
            newcodon = codonOptions[0]
            newSequence = newSequence[:(len(allCodons)-key-1)*3] + newcodon + newSequence[(len(allCodons)-key)*3:] #replace the appropriate codon in sequenceToMutate
            mutatedCount +=1
        #Once you reach 2, stop mutating
        if mutatedCount == 2:
            print(mutatedCount)
            break

    if mutatedCount == 0: #If nothing was mutated, return an empty string
        return ""
    else:
        return newSequence #Otherwise, return the mutated sequence

def codonFragmenter(sequence, direction = 'HAL', geneStrand = '+'):
    """
    Fragments codons in order moving away from the start or stop site. If gene on - strand, will take the reverse complement codon.

    params:
        sequence: sequence to fragment into codons as string.
        type: one of "homologyArm" or "primer". If primer, will truncate at far end into a fragment divisible by 3 before fragmenting into codons.
        direction: 'HAL' or 'HAR'. If on left homology arm, codon order is from right to left. If on right, codon order is from left to right.
        geneStrand: if +, codons are the original + strand sequence. If -, each individual codon has been reverse complemented.
    """

    orderedCodons = []

    #Define codons in threes along to sequence.
    for codonBase1 in range(0, len(sequence), 3):
        orderedCodons.append(sequence[codonBase1:codonBase1+3]) #Append this to orderedCodons
    if direction == 'HAL': #If on the left, reverse the ordered codons list to be moving away from the start/stop
        orderedCodons.reverse()
    if geneStrand == '-': #If on the minus strand, take revComp per codon
        for ind, codon in enumerate(orderedCodons):
            orderedCodons[ind] = revComp(codon)
    
    return orderedCodons

def positionScore(df):
    """
    Boolean evaluations of whether sgRNAs for a site meets required best-pick conditions and positional information about these sgRNAs.

    params:
        df: a dataframe for one start/stop site of the format: 
            df= {"start/stop":"start", #is it N or C termini -> do we need to look at start or stop codon for teh cut 
            'genome_start_codon_pos':403, 
            'genome_stop_codon_pos':406, # only one of these - the other is n/a
            'strand_type':'+',
            "sgRNA_list_positions":[[401,425],[403, 427]],#those wil be as genome positions -assumptions - the coordinates correspond to the 1st and last bp of the strand to which the gsRNA will be complementary to
            "sgRNA_list_values":["AAGCGACTA", "CCTGTAA"],
            "sgRNA_strand" : ['-', '+']}
    """

    import pandas as pd

    #Set up a new dataframe where boolean values will be stored per sgRNA in the input df
    sgRNACatalogue = pd.DataFrame(columns = ["sgRNA_sequence", "sgRNA_stop", "sgRNA_strand"], index = range(1, len(df["sgRNA_list_values"])))

    #Iterating through sgRNAs in df, add to the sgRNACatalogue
    for ind, sgRNA in enumerate(df["sgRNA_list_values"]):
        newCatalogueRow = [sgRNA, df["sgRNA_list_positions"][ind][1], df["sgRNA_strand"][ind]]
        sgRNACatalogue.loc[ind] = newCatalogueRow

    #Re-order the index as we've been appending to the top
    sgRNACatalogue = sgRNACatalogue.sort_index()

    #Calculate fmax of sgRNA - stop of genome start/stop position - this will be used for position scoring throughout
    #only have fmin in case of start codon - recalc fmax from this for homogeneity:
    if df["start/stop"] == "start":
        df["genome_stop_codon_pos"] == df["genome_stop_codon_pos"] + 2
    
    sgRNACatalogue["positionScore"] = sgRNACatalogue["sgRNA_stop"] - df["genome_stop_codon_pos"] #This is fmax - stop

    #Dataframe containing parameter ranges to interpret the positon score, based on gene strand, sgRNA strand, and start/stop
    positionScoreParameters = pd.read_excel("inputfiles/fmaxStopScore.xlsx")

    #Add the boolean columns to our output dataframe
    booleanColumns = ["P6inCDS", "Mutate2", "Mutate1", "OutsideCDS", "CutSiteInCDS", "PAMinCDS", "15bpOverhangL", "15bpOverhangR", "PAMoutsideCDS"]
    PAMColumns = ["PAMRelativeEnd", "PAMinStartStop"]
    sgRNACatalogue = sgRNACatalogue.reindex(columns = sgRNACatalogue.columns.tolist() + booleanColumns + PAMColumns)

    for ind, sgRNA in sgRNACatalogue.iterrows():

        #Extract the appropriate parameter row per sgRNA
        conditions = positionScoreParameters.loc[(positionScoreParameters["start/stop"] == df["start/stop"]) & (positionScoreParameters["strand_type"] == df["strand_type"]) & (positionScoreParameters["sgRNA_strand"] == sgRNACatalogue.at[ind, 'sgRNA_strand'])]
        conditions = conditions.reset_index(drop = True)

        #Per column, input true/false as to whether the position score meets that condition
        for col in booleanColumns:
            colValue = conditions.at[0,col] #extract parameter range values from dataframe
            #Process the value into a range (in format list [min, max])
            #If the values should be 'more than' or 'less than', 250 is used as a max or -250 as min because fragments are 225bp, so distances cannot be more than this
            if ">" in colValue: 
                minMax = [int(colValue[1:]), 250]
            elif "<" in colValue:
                minMax = [-250,int(colValue[1:])]
            elif ":" in colValue:
                min, max = colValue.split(":")
                minMax = [int(min), int(max)]
            else:
                print("Incorrect format of range value. Verify inputs.")
            
            #Into the output dataframe, print true/false as to whether the positionScore has met the condition for that column
            sgRNACatalogue.at[ind, col] = bool(sgRNACatalogue.at[ind, "positionScore"] in range(minMax[0], minMax[1]))

        #Calculate PAM position as a relative distance to the left end of start/stop if PAM is on the left, or right if on right (this can be calculated from the position score with the same equation for all cases)
        sgRNACatalogue.at[ind, "PAMRelativeEnd"] = int(abs(sgRNACatalogue.at[ind, "positionScore"]))

        #Calculate if PAM is in start/stop - where positionScore is 2 or less
        sgRNACatalogue.at[ind, "PAMinStartStop"] = bool(sgRNACatalogue.at[ind,"PAMRelativeEnd"] < 3)

    return sgRNACatalogue

def mutate_HDR_plasmid(HAL_R, HAR_F, df, sequenceType = "homologyArm"):
    """
    In the case where a fragment or primer needs to be mutated, will mutate in CDS (preferably PAM, if not in the sgRNA). If not possible, will mutate PAM outside of CDS to NGT.

    params:
        sequenceToMutate: the homology arm or primer fragment that will be mutated.
        winner_sgRNACatalogue: sgRNACatalogue in format as above, with only the row for the winner sgRNA selected.
        sequenceType: one of 'homologyArm' or 'primer'.
    """
    sgRNAScoreCatalogue = positionScore(df)

    #Check mutable conditions and posiions

    #1. Mutate PAM in CDS
    if sgRNAScoreCatalogue.at[1,"PAMinCDS"] is True:
        df = mutate_PAM_in_HDR_plasmid(HAL_R, HAR_F, df)

    #2. Mutate sgRNA in CDS
    if df["mutated?"] == "no": #Only proceed with this step if previous mutation was not possible
        df = mutate_sgRNA_recognition_site_in_HDR_plasmid(HAL_R, HAR_F, df, sgRNAScoreCatalogue) #will mutate 1 or 2 codons
    
    #3. Mutate PAM outside of CDS
    if df["mutated?"] == "no":
        df = mutatePamOutsideCDS(HAL_R, HAR_F, df, sgRNAScoreCatalogue)
    
    return df
    
def mutatePamOutsideCDS(HAL_R, HAR_F, df, sgRNAScoreCatalogue):
    """
    mutates PAM outside CDS from NGG to NGT
    """

    #Determine direction to go outside CDS, plus PAM position
    df[]




    #after mutating, revComp back

    return df

##Final runner

def retrieve_HDR_arm(df):

    """
    Accesses the HDR arm fragment from the dataframe with the start/ stop information about the transcription factors
    
    params: 
            df: dataframe with start/stop information about the transcription factors
        
    returns: 
            fragment: fragment to synthesize for obtaining the HDR arm
    """

    HAL = df["upstreamHA"]
    HAR = df["downstreamHA"] 

    return HAL, HAR

def retrieve_primer(primer_file, df):
    return

def find_best_gRNA(df):

    """

    This function is dependent on the functions check_start_stop_NGG(), exclude_over15() and select_closest().

    params: df={
            "start/stop":"start",#is it N or C termini -> do we need to look at start or stop codon for teh cut 
            'genome_start_codon_pos':400, 
            'genome_stop_codon_pos':700, # or only 1 of those 
            'strand_type':'+',
            "sgRNA_list_positions":[[401,425],[456,467],[478,489],[395,415]],#those wil be as genome positions -assumptions - the coordinates correspond to the 1st and last bp of the strand to which the gsRNA will be complementary to
            "sgRNA_list_values":["AAGCGACTA","AAAAAAAATAAAAA","ATATATTTTTTTTTTAAAAA","AGCGCGAAATAATA"]
        }

    returns:

    """
    from os import sys

    # check whether best case scenario holds true
    # if best case scenario holds true, return dictionary for best case scenario sgRNA

    best_sgRNA=check_start_stop_NGG(df)
  
    if best_sgRNA:

        ("best case scenario for sgRNA has occurred")

        winner_sgRNA = best_sgRNA #@TODO: This was a string before. Weird because that forbids the runner to assign an extra idem to the dictionary. Find out why does this have to be a string?

    # if best case scenario does not hold true check whether the 3´prime overhang of the sgRNA from the start/ stop is max. 15 bp and min. 1 bp
    # tag sgRNAs with a note about whether the PAM in the HDR needs to be mutated or not

    else:
        
        tagged_sgRNA = check_over_15(df)

        # if 3´prime overhang from start/ stop codon is max. 15 bp and min. 1 bp, check whether multiple sgRNA fulfil this condition

        if tagged_sgRNA["must_PAM_be_mutated_in_HDR_plasmid?"] == "no":

            print("PAM must not be mutated in HDR plasmid. Max 15 bp overhang on 3´prime end")

            max_15_sgRNA = tagged_sgRNA

            print(max_15_sgRNA)

            # if only one sgRNA fulfils this condition, return this sgRNA

            if len(max_15_sgRNA["sgRNA_list_values"]) == 1:

                winner_sgRNA == tagged_sgRNA

                winner_sgRNA["sgRNA_list_values"]=winner_sgRNA["sgRNA_list_values"][0]
                winner_sgRNA["sgRNA_list_positions"]=winner_sgRNA["sgRNA_list_positions"][0]
            
            # if multiple sgRNAs fulfil this condition, check if sgRNAs cut inside the CDS and whose cutting sites are closest to the start/ stop
            # return the most ideal sgRNA

            elif len(max_15_sgRNA["sgRNA_list_values"]) > 1:

                checked_cutting_inside_CDS_sgRNA = check_cutting_site_inside_CDS(max_15_sgRNA)

                winner_sgRNA = select_closest(checked_cutting_inside_CDS_sgRNA)

            # if neither of both cases is true, something is wrong with the pipeline, return error

            elif len(max_15_sgRNA["sgRNA_list_values"]) < 1: 

                print("Error! Empty sgRNA_list_values was returned even though the dictionary was tagged as max. 15 bp overhang!")

                sys.exit()

        # if tagged as 3´prime overhang more than 15 bp or less than 1 bp, check if there are multiple sgRNAs to select from

        elif tagged_sgRNA["must_PAM_be_mutated_in_HDR_plasmid?"] == "yes":

            print("PAM must be mutated in HDR plasmid")

            more_than_15_sgRNA = tagged_sgRNA

            # if there is only one sgRNA, return this one

            if len(more_than_15_sgRNA) == 1:

                winner_sgRNA = more_than_15_sgRNA

                winner_sgRNA["sgRNA_list_values"]=winner_sgRNA["sgRNA_list_values"][0]
                winner_sgRNA["sgRNA_list_positions"]=winner_sgRNA["sgRNA_list_positions"][0]

                HAL, HAR = winner_sgRNA["downstreamHA"], winner_sgRNA["upstreamHA"]

                mutated_winner = mutate_HDR_plasmid(HAL, HAR, winner_sgRNA)

                if mutated_winner:

                # if the primers output delivers an HDR arm, then mutate the HDR arm

                   winner_sgRNA = mutated_winner

                else:

                # if the primers output does not deliver an HDR arm, then mutate the primer for the HDR arm if possible

                    #HAL_R, HAR_F = retrieve_primer(primer_file, TF)

                    #mutated_winner = mutate_PAM_in_HDR_plasmid(HAL_R, HAR_F, df)

                    #if not mutated_winner:

                        #print("Neither the HDR-arm nor the primers could be mutated!")

                        #winner_sgRNA = {}

                        #break

                    #else: 

                        #winner_sgRNA = mutated_winner
            
            # if there are multiple sgRNAs, check whether they all cut inside the CDS
            
            elif len(more_than_15_sgRNA) > 1: 
                
                checked_cutting_inside_CDS_sgRNA = check_cutting_site_inside_CDS(more_than_15_sgRNA)

                # if there is only one sgRNA that cuts inside the CDS return this one

                if len(checked_cutting_inside_CDS_sgRNA)== 1:

                    winner_sgRNA = checked_cutting_inside_CDS_sgRNA

                    winner_sgRNA["sgRNA_list_values"]=winner_sgRNA["sgRNA_list_values"][0]
                    winner_sgRNA["sgRNA_list_positions"]=winner_sgRNA["sgRNA_list_positions"][0]
            

                    #@TODO: write the retrieve_HDR_arm() function

                    # check whether primers output delivers HDR arm

                    #if retrieve_HDR_arm():

                        # if the primers output delivers an HDR arm, then mutate the HDR arm

                        #mutate_PAM_in_HDR_arm()

                    #else:

                        # if the primers output does not deliver an HDR arm, then mutate the primer in the HDR arm if possible

                        #mutate_PAM_in_HDR_plasmid()

                # if there are multiple sgRNAs that cut insite the CDS, return the one whose cutting site is closest to start/ stop

                elif len(checked_cutting_inside_CDS_sgRNA) > 1: 

                    winner_sgRNA = select_closest(checked_cutting_inside_CDS_sgRNA)

                    #@TODO: write the retrieve_HDR_arm() function

                    # check whether primers output delivers HDR arm

                    #if retrieve_HDR_arm():

                        # if the primers output delivers an HDR arm, then mutate the HDR arm

                        #mutate_PAM_in_HDR_arm()

                    #else:

                        # if the primers output does not deliver an HDR arm, then mutate the primer in the HDR arm if possible

                        #mutate_PAM_in_HDR_plasmid()

            elif len(more_than_15_sgRNA) < 1: 

                print("Error! It should have checked earlier whether there is a sgRNA! Something is wrong with the pipeline!")

                sys.exit()

    return winner_sgRNA
