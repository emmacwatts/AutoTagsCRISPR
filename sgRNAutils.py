#Complete functions
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

def refSeq():
    """
    Creates a Bio SeqIO element from the Drosophila melanogaster reference genome.

    returns: refSeqPerChromosome - a dictionary of seq ID: sequence (in SeqIO fasta format) for the full reference genome.
    """
    from Bio import SeqIO

    refSeqPerChromosome = {}
    for seq in SeqIO.parse(open("inputfiles/dmel-all-chromosome-r6.48.fasta"), 'fasta'): #This is the FASTA file of the reference genome sequence
        refSeqPerChromosome[seq.id] = seq.seq 

    return refSeqPerChromosome

def make_dataframe_from_TFs_list(TF_list, refSeqPerChromosome, annotation = "inputfiles/dmel-all-r6.48.gtf"):
    '''
    Creating a dataframe from sequence information and genes of interest. Depends on functions revComp(). 

    params:
        TF_list: xlsx file listing the genes of interest
        refSeqPerChromosome: reference chromosome stored in SeqIO sequence format, indexed by refSeqPerChromosome[seq.id] = seq.seq
        annotation: gtf file containing the gene annotation information for the reference genome

    returns: dataframe of information about each start/stop site, in format:
        'Gene_ID': 'FBgn0004652',
        'Transcript_ID': 'FBtr0083651', 
        'Chromosome': '3R', 
        'Gene_Region': 'stop_codon', 
        'Start': 18426145, 
        'Stop': 18426147, 
        'Strand': '-', 
        'Reference_Seq': 'TTGATCGTAGGACAC...', 
    '''

    import pandas as pd
    
    #This is the input file containing the TFs we want to query
    queryTFsdf = pd.read_excel(TF_list)

    #This is the .gtf file with annotations for each gene on the reference genome
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

    #re-index
    TFsdf = TFsdf.reset_index()
    del TFsdf["index"]

    return TFsdf

def make_homology_arm_fragments(TFsdf, refSeqPerChromosome):
    """
    Takes in the dataframe of information about start/stop codon regions, 
    and appends with columns for the 225 bp upstream and downstream of
    this site (the homlogy arm fragments).

    input: TFsDF - as defined in previous function and on GitHub README, with Reference_Seq 1600bp either side of the start/stop.
           refSeqPerChromosome: reference chromosome stored in SeqIO sequence format, indexed by refSeqPerChromosome[seq.id] = seq.seq
    output: TFsDF appended with HAL and HAR

    """
    for index, rowcontents in TFsdf.iterrows():

            #Define 225 bp region left and right from Start and Stop of gene region
            leftStart = rowcontents["Start"] - 226 #SeqIO takes the start coordinate as exclusive, which is why we include one extra bp
            leftStop = rowcontents["Start"]  - 1 #SeqIO takes the end coordinate as inclusive, so we subtract one bp because we want to exlculde the first position of the gene region
            rightStart = rowcontents["Stop"] #SeqIO takes the start coordinate as exclusive, which is fine because we want to exclude the last position of the gene region
            rightStop = rowcontents["Stop"] + 225 #SeqIO takes the end coordinate as inclusive, which is fine
            
            #Add left and right homology arm (HA)
            HAL = revComp(str(refSeqPerChromosome[rowcontents["Chromosome"]][leftStart:leftStop]))
            HAR = revComp(str(refSeqPerChromosome[rowcontents["Chromosome"]][rightStart:rightStop]))
            TFsdf.at[index,"HAL"] = HAL
            TFsdf.at[index, "HAR"] = HAR

    return TFsdf

def find_synonymous_codons(query_codon, base_to_change, codon_table_excel = "inputfiles/codon_table.xlsx"):

    '''
    Uses the amino acids table to select codons that encode for the same amino acid as the query codon. Will ensure the specified base has been mutated.

    Params:
        query_codon: string, codon to select synonymous codons for
        base_to_change: the base within the codon that needs to change. Numeric value from 1-3.
        codon_table_excel: string, path to an excel file that lists per codon which amino acid that codon encodes.
    
    Returns:
        synonymous_codons: list of strings, each string is a codon that encodes for the same amino acid as the query codon.
    
    '''

    import pandas as pd

    codon_table = pd.read_excel(codon_table_excel)

    #Extract amino acid given codon
    amino_acid_query = codon_table[codon_table['codon'] == query_codon].iloc[0]["amino_acid"]

    #Subset df for other rows corresponding to this amino acid
    same_aa_df=codon_table[codon_table["amino_acid"] == amino_acid_query]

    #Convert the codons to a list
    codon_list = same_aa_df["codon"].values.tolist()

    #Keep only codons where the indicated base has changed (this will also remove original codon)
    changedCodons = [codon for codon in codon_list if codon[base_to_change-1] != query_codon[base_to_change-1]]
            
    return changedCodons

def filter_gRNA(gRNA_file, window, tfSingleRow, refSeqPerChromosome):
    """
    Selects gRNAs within window of the start/stop site and returns these in a dataframe with information about their coordinates and the start/stop site.
    params:
        gRNA_file: gff file
        tfSingleRow: a pandas series for one row of the original tfsDF with the following format
            Gene_ID                                                    FBgn0000022
            Transcript_ID                                              FBtr0070072
            Chromosome                                                           X
            Gene_Region                                                start_codon
            Start                                                           370094
            Stop                                                            370096
            Strand                                                               +
            Reference_Seq        CTAATGAATAGATTGGTGTGTGATGTAGTGATCTAATATGGTGAAG...
            upstreamHA           TTTGCTCAGTTTTTTATTGGCGCCGGGACCAATTCCCCGGCGACCA...
            downstreamHA         TTAAGAGATAGTATAACGTTATTGTGTGACGATGCTCCTTGCTTCA...
        
    returns: filtered_gRNAs - a dataframe of sgRNAs for one start/stop site of the format:
            fmin	fmax	#chr	strand	sgRNA_sequence	        Gene_ID	    Transcript_ID	Chromosome	Gene_Region	Start	...	downstreamHA
        0	370082	370104	X	    +	    CTATCTCTTAAAATGGCTTTGGG	FBgn0000022	FBtr0070072	    X	        start_codon	370094	...	TTAAGAGATAGTATAACGTTATTGTGTGACGATGCTCCTTGCTTCA...
        1	370693	370715	X	    -	    CCTGTAAAAAAACAGATCAAATC	FBgn0000022	FBtr0070072	    X	        start_codon	370694	...	TTAAGAGATAGTATAACGTTATTGTGTGACGATGCTCCTTGCTTCA...

    """
    import pandas as pd

    # create a dataframe from the gRNA file
    gRNAFileAnnotation = pd.read_csv(gRNA_file, sep = "\t", index_col = False)

    # add a new category to the dataframe that provides information about whether the sequence deviates from the transgenic strain 
    gRNAFileAnnotation = gRNAFileAnnotation.assign(target_site_variation= "")

    # reformat the "Attribute" category in refGenomeAnnotation, to extract Gene_ID, Gene_Symbol, and Transcript ID
    # for each attribute value, extract the gene ID and symbol and add this to the new categories
    for index, attribute in enumerate(gRNAFileAnnotation['attributes']):
        fullatt = (attribute).split(";")
        gRNAFileAnnotation.at[index,"target_site_variation"] = fullatt[8]
    
    # shorten file to essential information
    GenomeCoordinates = gRNAFileAnnotation[["target_site_variation", "fmin", "fmax", "#chr", "strand"]]

    #Filter for the appropriate conditions:
        # check whether the sgRNAs match the transgeneic genome, whether the sgRNAs match to the same chromosome as the transcription factor
        # select sgRNA that are located maximally 20 bp upstream of the start/stop codon of the transcription factors
        # select sgRNA that are located maximally 20 bp downstream of the start/stop codon of the transcription factors
    filtered_gRNAs = GenomeCoordinates[(GenomeCoordinates['target_site_variation'] == "target_site_variation=none observed") & 
                                                  (GenomeCoordinates['#chr'] == tfSingleRow["Chromosome"]) & 
                                                  (GenomeCoordinates['fmin'] >= tfSingleRow["Start"] - window) & 
                                                  (GenomeCoordinates['fmax'] <= tfSingleRow["Stop"] + window)]

    #add sgRNA sequence using coordinates
    filtered_gRNAs = filtered_gRNAs.assign(sgRNA_sequence= "")
    for index, row in filtered_gRNAs.iterrows():
        filtered_gRNAs.loc[index,'sgRNA_sequence'] = str(refSeqPerChromosome[row['#chr']][int(row["fmin"])-1:int(row["fmax"])])
    
    #Add columns from tfSingleRow containing start/stop site info
    for (columnName, columnData) in tfSingleRow.items():
        filtered_gRNAs = filtered_gRNAs.assign(**{columnName: columnData})

    #drop target site variation
    filtered_gRNAs = filtered_gRNAs.drop(columns=["target_site_variation"])

    #reset index
    filtered_gRNAs = filtered_gRNAs.reset_index(drop = True)

    return filtered_gRNAs

def gRNA_stringencyIterator(tfSingleRow, window, refSeqPerChromosome):
    """
    Iterates trough gRNA stringency files until at least one gRNA is returned from filter_gRNAs function.
    Function depends on filter_gRNA() functioon.

    params:
        tfSingleRow: a pandas Series for one row of the tfsDF dataframe produced by make_dataframe_from_TFs_list
        refSeqPerChromosome: a Bio SeqIO dictionary for the D. melanogaster genome
        window: the window size to search for gRNAs around the start/stop site

    returns:
        filtered_gRNAs: a df of gRNAs that meet the conditions at this site
    """
    import pandas as pd

    gRNAfiles = ["inputfiles/NoOffTarget_high_stringency.gff", "inputfiles/NoOffTarget_med_stringency.gff", 
    "inputfiles/NoOffTarget_low_stringency.gff", "inputfiles/1to3NonCdsOffTarget_low_stringency.gff", 
    "inputfiles/ManyOffTarget_low_stringency.gff"]

    #Set up sgRNA variable
    sgRNA = pd.DataFrame()
    file_ind = 0 

    #Loop through files and retain the filtered sgRNA df if at least one row is present
    while sgRNA.empty and file_ind <=4:
        sgRNA = filter_gRNA(gRNAfiles[file_ind], window, tfSingleRow, refSeqPerChromosome)
        print(f"stringency {file_ind}")
        file_ind +=1
    
    return sgRNA

def sgRNApositionCheck(sgRNAdf, minDistance, maxDistance):
    """
    Given filtered sgRNAs for a start/stop site in the 'filtered_gRNAs' format, will create dataframe containing positional information and condition checks
    for each sgRNA.
    minDistance and maxDistance must be divisible by 3.

    params: minDistance, maxDistance: minimum and maximum distances from the start/stop site that the sgRNA can be in
            df: a pandas df of sgRNAs for one start/stop site:
            fmin	fmax	#chr	strand	sgRNA_sequence	        Gene_ID	    Transcript_ID	Chromosome	Gene_Region	Start	... downstreamHA
        0	370082	370104	X	    +	    CTATCTCTTAAAATGGCTTTGGG	FBgn0000022	FBtr0070072	    X	        start_codon	370094	...	TTAAGAGATAGTATAACGTTATTGTGTGACGATGCTCCTTGCTTCA...
        1	370693	370715	X	    -	    CCTGTAAAAAAACAGATCAAATC	FBgn0000022	FBtr0070072	    X	        start_codon	370694	...	TTAAGAGATAGTATAACGTTATTGTGTGACGATGCTCCTTGCTTCA...
        
    output: the same dataframe with added columns:
            fmin	fmax	#chr	strand	sgRNA_sequence	        Gene_ID	    Transcript_ID	Chromosome	Gene_Region	Start	...	downstreamHA	                                    positionScore	PAM_in_start/stop	<15_bp3’_overhang	PAM_in_CDS	PAM_outside_CDS	CDS_boundary	lastG	cutSite	
        0	370082	370104	X	    +	    CTATCTCTTAAAATGGCTTTGGG	FBgn0000022	FBtr0070072	    X	        start_codon	370094	...	TTAAGAGATAGTATAACGTTATTGTGTGACGATGCTCCTTGCTTCA...	8	            False	            True	            True	    False	        >22	            31.0	26.0	
        1	370693	370715	X	    -	    CCTGTAAAAAAACAGATCAAATC	FBgn0000022	FBtr0070072	    X	        start_codon	370694	...	TTAAGAGATAGTATAACGTTATTGTGTGACGATGCTCCTTGCTTCA...	16	            False	            False	            True	    False	        <20	            16.0	21.0	
                    
    """
    import pandas as pd
    import sys
    import numpy as np

    # Sanity check for divisibility by 3
    if minDistance % 3 != 0 or maxDistance % 3 != 0:
        raise ValueError("minDistance and maxDistance must be divisible by 3.")

    #Adding position scores (fmax - stop)
    sgRNAdf["positionScore"] = sgRNAdf["fmax"] - sgRNAdf["Stop"]
    
    #Dataframe containing parameter ranges to interpret the positon score, based on gene strand, sgRNA strand, and start/stop
    positionScoreParameters = pd.read_excel("inputfiles/fmaxStopScoreML.xlsx")

    #Per parameter, append the sgRNAdf with a TRUE/FALSE value per condition.
    #Also add positional information about the CDS boundary, position of the last G/C, the cut site and the sgRNA recognition site
    #from the conditions file.

    #Start by adding columns for these to sgRNAdf
    booleanColumns = ["PAM_in_start/stop", "max_15_bp_3’_overhang", "PAM_in_CDS", "PAM_outside_CDS", "SRS_in_CDS"]
    positionColumns = ["CDS_boundary", "non_CDS_boundary", "SRS_boundary", "mutable_PAM"]

    # Use assign for adding new columns
    sgRNAdf = sgRNAdf.assign(**{col: None for col in ["CDS_side"] + ["non_CDS_side"]+ booleanColumns + positionColumns + ["cut_site"]})
    
    #Iterating through sgRNAdf, extract the appropriate row of the conditions file and add new information to sgRNAdf.
    for ind, sgRNA in sgRNAdf.iterrows():
        conditions = positionScoreParameters.loc[
            (positionScoreParameters["start/stop"] == sgRNA['Gene_Region']) &
            (positionScoreParameters["strand_type"] == sgRNA['Strand']) &
            (positionScoreParameters["sgRNA_strand"] == sgRNA['strand'])
        ].reset_index(drop=True)

        for col in booleanColumns + positionColumns:
            colValue = conditions.at[0, col]
            
            if ":" in colValue:
                min, max = colValue.split(":")
                minMax = [min, max]
                newMinMax = []
                
                for value in minMax:
                    if value == "minDistance":
                        newMinMax.append(minDistance)
                    elif value == "maxDistance":
                        newMinMax.append(maxDistance+1)
                    else:
                        newMinMax.append(int(value))

            if col in booleanColumns:
                sgRNAdf.at[ind, col] = bool(sgRNA["positionScore"] in range(newMinMax[0], newMinMax[1]))

            if col in positionColumns:
                sgRNAdf.at[ind, col] = list(range(newMinMax[0], newMinMax[1]))

        #Calculate cutsite coordinate using positionScore. Add this to the sgRNA catalogue.
        sgRNAdf.loc[ind, "cut_site"] = int(conditions.at[0, "cut_site"]) - sgRNA["positionScore"]

        # Add CDS side to sgRNA catalogue
        sgRNAdf.loc[ind, "CDS_side"] = conditions.at[0, "CDS_side"]

        # Add non CDS side to sgRNA catalogue
        sgRNAdf.loc[ind, "non_CDS_side"] = conditions.at[0, "non_CDS_side"]

    return sgRNAdf
    
def checkCDSCutandOrder(sgRNAdf):
    """
    Given an sgRNAdf, will calculate firstly rows where the sgRNA cuts inside CDS. If multiple, selects that which cuts closest to start/stop.
    If none, selects closest cut sgRNA that is outside CDS.
    
    params: sgRNAdf in format:
            fmin	fmax	#chr	strand	sgRNA_sequence	        Gene_ID	    Transcript_ID	Chromosome	Gene_Region	Start	...	downstreamHA	                                    positionScore	PAM_in_start/stop	<15_bp3’_overhang	PAM_in_CDS	PAM_outside_CDS	CDS_boundary	lastG	cutSite mutated	
        0	370082	370104	X	    +	    CTATCTCTTAAAATGGCTTTGGG	FBgn0000022	FBtr0070072	    X	        start_codon	370094	...	TTAAGAGATAGTATAACGTTATTGTGTGACGATGCTCCTTGCTTCA...	8	            False	            True	            True	    False	        >22	            31.0	26.0	False
        1	370693	370715	X	    -	    CCTGTAAAAAAACAGATCAAATC	FBgn0000022	FBtr0070072	    X	        start_codon	370694	...	TTAAGAGATAGTATAACGTTATTGTGTGACGATGCTCCTTGCTTCA...	16	            False	            False	            True	    False	        <20	            16.0	21.0	False
     
    output: pandas series of the winning sgRNA in format:
    fmin                                                            409794
    fmax                                                            409816
    #chr                                                                 X
    strand                                                               -
    sgRNA_sequence                                 CCATGACGAGCATTTGCAGCAGC
    Gene_ID                                                    FBgn0002561
    Transcript_ID                                              FBtr0070074
    Chromosome                                                           X
    Gene_Region                                                start_codon
    Start                                                           409796
    Stop                                                            409798
    Strand                                                               +
    Reference_Seq        CGACGTCGTCCAGCTTGATGAAAAACTTTGTGGCGGTGTGCGCCAG...
    upstreamHA           CATGTTATTGTAGTTGAACTTCTTCTTCTGGGAGGACAACATGGGT...
    downstreamHA         GTAATCCTTGCGAGAGTTTTCTAAGATTTAGTTTACAGATGTTGAC...
    positionScore                                                       18
    PAM_in_start/stop                                                 True
    <15_bp3’_overhang                                                 True
    PAM_in_CDS                                                       False
    PAM_outside_CDS                                                  False
    CDS_boundary                                                       >22
    lastG                                                             18.0
    cutSite                                                           23.0

    """

    #extract CDS boundary value
    sgRNAdf.reset_index(inplace= True, drop= True)
    #Calculate sgRNA/cut site difference
    sgRNAdf["cutsite-CDSbound"] = abs(sgRNAdf["cut_site"])
    #C. Check cut site in CDS
    #filter for only sgRNAs that cut in CDS
    conditionC = sgRNAdf[sgRNAdf["cut_site"].isin(sgRNAdf["CDS_boundary"])]

    #if there are sgRNAs that cut in CDS, sort by distance of cutsite and CDS start. Select the one that cuts closest.
    if len(conditionC) > 0: #cuts in CDS, closest cut (C1, C2)
        conditionCclosestCut = conditionC.sort_values('cutsite-CDSbound')
        conditionCclosestCut = conditionCclosestCut.reset_index(drop = True) #reset index
        winnersgRNA = conditionCclosestCut.at[0, "sgRNA_sequence"]
    else: #no sgRNAs cut in CDS, select closest that still met condition B (C3)
        nonCDSclosestCut = sgRNAdf.sort_values('cutsite-CDSbound')
        winnersgRNA = nonCDSclosestCut.at[0, "sgRNA_sequence"]

    return winnersgRNA

def find_best_gRNA(sgRNAdf):
    """
    Select the best guide RNA from a dataframe of guideRNAs that are ±20bp from the start/stop site.

    params:
        df: the dataframe of potential guideRNAs in format:
            fmin	fmax	#chr	strand	sgRNA_sequence	        Gene_ID	    Transcript_ID	Chromosome	Gene_Region	Start	...	downstreamHA	                                    positionScore	PAM_in_start/stop	<15_bp3’_overhang	PAM_in_CDS	PAM_outside_CDS	CDS_boundary	lastG	cutSite	mutated
        0	370082	370104	X	    +	    CTATCTCTTAAAATGGCTTTGGG	FBgn0000022	FBtr0070072	    X	        start_codon	370094	...	TTAAGAGATAGTATAACGTTATTGTGTGACGATGCTCCTTGCTTCA...	8	            False	            True	            True	    False	        >22	            31.0	26.0	False
        1	370693	370715	X	    -	    CCTGTAAAAAAACAGATCAAATC	FBgn0000022	FBtr0070072	    X	        start_codon	370694	...	TTAAGAGATAGTATAACGTTATTGTGTGACGATGCTCCTTGCTTCA...	16	            False	            False	            True	    False	        <20	            16.0	21.0	False
   
    output:
        winnersgRNA: the sequence of the best guideRNA (str)
        mutationNeeded: True/False of whether this sgRNA needs to be mutated.

    """
    #Set up the winning guide RNA
    winnerFound = False
    mutationNeeded = False
    winnersgRNA = ""

    #A. Ideal condition - PAM in start/stop
    conditionA = sgRNAdf[sgRNAdf["PAM_in_start/stop"] == True] #This is the subset df for which PAM is in the start/stop
    if len(conditionA) > 0: #if there is one or more sgRNAs for this condition, select the first as the winner
        conditionA = conditionA.reset_index(drop = True) #reset index
        winnersgRNA = conditionA.at[0, "sgRNA_sequence"]
        winnerFound = True

    #B. sgRNA overhang is max 15bp
    if winnerFound == False:
        conditionB = sgRNAdf[sgRNAdf["<15_bp3’_overhang"] == True]
        if len(conditionB) == 1:
            conditionB = conditionB.reset_index(drop = True)
            winnersgRNA = conditionB.at[0, "sgRNA_sequence"]
            winnerFound = True
        elif len(conditionB) > 1: #if more than one, check condition C (select for cut site in CDS preferencially, and the sgRNA where cut site is closest to CDS)
            winnersgRNA = checkCDSCutandOrder(conditionB) 
            winnerFound = True   

    #D. sgRNA overhang is more than 15bp, need to mutate
    if winnerFound == False:
        if len(sgRNAdf) == 1: #if just 1 sgRNA, select this one
            winnersgRNA = sgRNAdf.at[0, "sgRNA_sequence"]
        else: #more than 1 sgRNA
            winnersgRNA = checkCDSCutandOrder(sgRNAdf) #Again, iterate through condition C to select an sgRNA
            winnerFound = True         

        mutationNeeded = True

    return winnersgRNA, mutationNeeded

def codonFragmenter(winnerdf, side, side_boundary):
    """
    For the CDS of each start/stop, will create a list of codons.
    If gene is on - strand, these will be revComp codons.

    params:
        side: "HAL" or "HAR" depending on where base to mutate
        winnerdf: pandas Series for the row of the sgRNA dataframe containing information about the winner sgRNA:
            fmin                                                            396165
            fmax                                                            396187
            #chr                                                                 X
            strand                                                               -
            sgRNA_sequence                                 CCGAGTGTGTTAATGAAAAATAA
            Gene_ID                                                    FBgn0004170
            Transcript_ID                                              FBtr0070073
            Chromosome                                                           X
            Gene_Region                                                start_codon
            Start                                                           396177
            Stop                                                            396179
            Strand                                                               +
            Reference_Seq        ACCTGCGATAATTTGACATTCTTAGAAACTACCTGAAGAAATAGGA...
            upstreamHA           TCTGGTCAGTGCCATACCCCTTGGTGTATACTTGCGAGTCTTAATT...
            downstreamHA         AACACACTCGGAGCTTTCTTTAACTTTCCGGATAACGATCAACAGA...
            positionScore                                                        8
            PAM_in_start/stop                                                False
            <15_bp3’_overhang                                                 True
            PAM_in_CDS                                                       False
            PAM_outside_CDS                                                   True
            CDS_boundary                                                       >22
            lastG                                                              8.0
            cutSite                                                           13.0
            mutated                                                          False

    returns:
        codonList: codons from CDS respective to the region of interest in mutable format (revComp if gene is -)
    """
    import numpy as np
    
    #Extract HA sequence where CDS
    HA = winnerdf[side]

    if side == "HAR":
        winnerdf[f"{side_boundary}_HA"] = np.array(winnerdf[side_boundary]) - 1 #HA is 1bp to the right of stop
        mutable_region = HA[:winnerdf["side_boundary_HA"][-1]+1]
    
    if side == "HAL":
        winnerdf[f"{side_boundary}_HA"] = np.array(winnerdf[side_boundary]) + 2 #HA is 3pb to the left of stop but because we index backwards we calculate - 1 
        mutable_region = HA[winnerdf["side_boundary_HA"][0]+1:]
    
    #Start codon list
    codonList = []

    #Chop up the CDS into codons
    for codonBase1 in range(0, len(mutable_region), 3):
        codon = mutable_region[codonBase1:codonBase1+3]    
        if winnerdf["Strand"] == '-': #Create reverse complement of every codon if gene is on - strand
            codonList.append(revComp(codon))
        else: codonList.append(codon)

    return codonList, winnerdf

def codonReverseFragmenter(codonsList, winnerdf, side, side_boundary):
    """
    Will take a fragmented list of codons of CDS around the region of interest and replace the mutated HAL/HAR arms in tfsDF.
    If gene is on - strand, revComp codons will be reversed back to the + strand sequences.

    params:
        codonList: codons from the ±21bp region around the start/stop site in mutable format (revComp if gene is -)
        winnerdf: pandas series of the winning df in format:

        fmin                                                            396165
        fmax                                                            396187
        #chr                                                                 X
        strand                                                               -
        sgRNA_sequence                                 CCGAGTGTGTTAATGAAAAATAA
        Gene_ID                                                    FBgn0004170
        Transcript_ID                                              FBtr0070073
        Chromosome                                                           X
        Gene_Region                                                start_codon
        Start                                                           396177
        Stop                                                            396179
        Strand                                                               +
        Reference_Seq        ACCTGCGATAATTTGACATTCTTAGAAACTACCTGAAGAAATAGGA...
        upstreamHA           TCTGGTCAGTGCCATACCCCTTGGTGTATACTTGCGAGTCTTAATT...
        downstreamHA         AACACACTCGGAGCTTTCTTTAACTTTCCGGATAACGATCAACAGA...
        positionScore                                                        8
        PAM_in_start/stop                                                False
        <15_bp3’_overhang                                                 True
        PAM_in_CDS                                                       False
        PAM_outside_CDS                                                   True
        CDS_boundary                                                       >22
        lastG                                                              8.0
        cutSite                                                           13.0
        mutated                                                          False
    
    output: the winningdf as previously, now with mutated upstreamHA or downstreamHA.
        fmin                                                            396165
        fmax                                                            396187
        #chr                                                                 X
        strand                                                               -
        sgRNA_sequence                                 CCGAGTGTGTTAATGAAAAATAA
        Gene_ID                                                    FBgn0004170
        Transcript_ID                                              FBtr0070073
        Chromosome                                                           X
        Gene_Region                                                start_codon
        Start                                                           396177
        Stop                                                            396179
        Strand                                                               +
        Reference_Seq        ACCTGCGATAATTTGACATTCTTAGAAACTACCTGAAGAAATAGGA...
        upstreamHA           TCTGGTCAGTGCCATACCCCTTGGTGTATACTTGCGAGTCTTAATT...
        downstreamHA         AACACACTCGGAGCTTTCTTTAACTTTCCGGATAACGATCAACAGA...
        positionScore                                                        8
        PAM_in_start/stop                                                False
        <15_bp3’_overhang                                                 True
        PAM_in_CDS                                                       False
        PAM_outside_CDS                                                   True
        CDS_boundary                                                       >22
        lastG                                                              8.0
        cutSite                                                           13.0
        mutated                                                          False

    """
    #return to + strand if the gene is on -
    if winnerdf["Strand"] == '-':
        for ind, codon in enumerate(codonsList):
            codonsList[ind] = revComp(codon)
    
    #recombine the codons of the CDS
    mutable_region = ''.join(codonsList)

    #replace mutated HA
    if winnerdf[side] == "HAR":
        winnerdf.at["HAR"] = mutable_region + winnerdf["HAR"][winnerdf[f"{side_boundary}_HA"][-1]+1:]
    if winnerdf[side] == "HAL":
        winnerdf.at["HAL"] =  winnerdf["HAL"][:winnerdf[f"{side_boundary}_HA"][0]+1] + mutable_region

    return winnerdf

def mutator_Marina(basePosition, fragmentedCDS, winnerdf, codonCoordinates, PAM = False):

    """
    (1) converts the position of the base that we would like to try and mutate into a codon position
    (2) finds synonymous mutations for the codon in the the codon position of the fragmented CDS
    (3) picks a synonymous mutation that mutates the base that we would like to try and mutate if possible
    (4) returns the mutated fragemented CDS if mutation was possible or throws exemption

    params:
        basePosition: integer, base position within sgRNA relative to fmax
        fragmentedCDS: list, codons in direction of translation between gene region of interest and maxDistance or minDistance 
    """
    pos_relative_to_CDS = basePosition + winnerdf["positionScore"] - winnerdf["CDS_boundary"][0] # explained in Positions section of MutationLogic.pptx
    codonNumber = codonCoordinates.at[pos_relative_to_CDS, 'codon'] #this is the codon number (as an index in mutableCodons)
    codon = fragmentedCDS[codonNumber] #This is the codon we want to mutate
    base = codonCoordinates.at[pos_relative_to_CDS, 'base'] #This is the base within that codon (1-3)
    potentialCodons = find_synonymous_codons(query_codon =codon, base_to_change= base)
        
    if potentialCodons:
        if PAM == True:
            #check for NGA and remove if present
            potentialCodons = [codons for codons in potentialCodons if codons[-2:] != "GA"]
        fragmentedCDS[codonNumber] = potentialCodons[0]
        
    else: 
        raise ValueError("No synonymous codons found.")
    
    return fragmentedCDS

def find_best_mutation(winnerdf, maxDistance):
    """
        In the case where a fragment needs to be mutated, will mutate in CDS (preferably PAM, if not then in the sgRNA). 
        If not possible, will mutate PAM outside of CDS to NTG/CTN.
        
        params:
            winnerdf: the pandas series of the winning sgRNA in format:
                fmin                                                            396165
                fmax                                                            396187
                #chr                                                                 X
                strand                                                               -
                sgRNA_sequence                                 CCGAGTGTGTTAATGAAAAATAA
                Gene_ID                                                    FBgn0004170
                Transcript_ID                                              FBtr0070073
                Chromosome                                                           X
                Gene_Region                                                start_codon
                Start                                                           396177
                Stop                                                            396179
                Strand                                                               +
                Reference_Seq        ACCTGCGATAATTTGACATTCTTAGAAACTACCTGAAGAAATAGGA...
                upstreamHA           TCTGGTCAGTGCCATACCCCTTGGTGTATACTTGCGAGTCTTAATT...
                downstreamHA         AACACACTCGGAGCTTTCTTTAACTTTCCGGATAACGATCAACAGA...
                positionScore                                                        8
                PAM_in_start/stop                                                False
                <15_bp3’_overhang                                                 True
                PAM_in_CDS                                                       False
                PAM_outside_CDS                                                   True
                CDS_boundary                                                       >22
                lastG                                                              8.0
                cutSite                                                           13.0
                mutated                                                          False

        """
    import pandas as pd

    if winnerdf["Strand"] == '+':
        codonCoordinates = pd.read_excel("inputfiles/fmaxStopScore.xlsx", sheet_name="CodonCoordinatePlus", index_col= 0)
    elif winnerdf["Strand"] == '-':
        codonCoordinates = pd.read_excel("inputfiles/fmaxStopScore.xlsx", sheet_name="CodonCoordinateMinus", index_col= 0)

    # Condition A: PAM in CDS
    if winnerdf["PAM_in_CDS"] == True and winnerdf["mutated?"] == False:
        fragmentedCDS, winnerdf = codonFragmenter(winnerdf, winnerdf["CDS_side"], winnerdf["CDS_boundary"])

        for basePosition in winnerdf["mutable_PAM"]:
            
            try: 
                fragmentedCDS = mutator(basePosition, fragmentedCDS, winnerdf, PAM=True)
                winnerdf["mutated?"] = True
                break

            except ValueError: continue

        winnerdf = codonReverseFragmenter(fragmentedCDS, winnerdf, winnerdf["CDS_side"], winnerdf["CDS_boundary"])

    # Condition B: SRS in CDS
    if winnerdf["sgRNA_recognition_site_in_CDS"] == True and winnerdf["mutated?"] == False:
        number_of_mutations = 0
        position = 0
        fragmentedCDS, winnerdf = codonFragmenter(winnerdf, winnerdf["CDS_side"], winnerdf["CDS_boundary"])

        while number_of_mutations < 2 and position in range(0, len(winnerdf["sgRNA_recognition_site_in_CDS"])):
                
            try: 
                fragmentedCDS = mutator(basePosition, fragmentedCDS, winnerdf)
                winnerdf["mutated?"] = True
                position += 1
                number_of_mutations += 1

            except ValueError: continue
        
        winnerdf = codonReverseFragmenter(fragmentedCDS, winnerdf, winnerdf["CDS_side"], winnerdf["CDS_boundary"])
    
    # Condition C: PAM outside of CDS
    if winnerdf["mutated?"] == False:
        fragmented_non_CDS_side, winnerdf = codonFragmenter(winnerdf, winnerdf["non_CDS_side"], winnerdf["non_CDS_boundary"])
        codonCoordinates = pd.read_excel("inputfiles/fmaxStopScoreML.xlsx", sheet_name=winnerdf["coordinates"], index_col= 0)
        pos_relative_to_gene_region = winnerdf["mutable_PAM"][0] + winnerdf["positionScore"]
        base = codonCoordinates.at[pos_relative_to_gene_region, 'relative_coordinate_CDS'] #This is one base of the PAM within HA 
        newHA = winnerdf[f"{PAM_side}"][:base-1] + 'T' + winnerdf[f"{PAM_side}"][base:] #add a 'T' where the old base was in HA
        winnerdf[f"{PAM_side}"] = newHA
        winnerdf["mutated?"] = True

    return winnerdf

def mutator(winnerdf):
    """
    In the case where a fragment needs to be mutated, will mutate in CDS (preferably PAM, if not then in the sgRNA). 
    If not possible, will mutate PAM outside of CDS to NTG/CTN.
    
    params:
        winnerdf: the pandas series of the winning sgRNA in format:
            fmin                                                            396165
            fmax                                                            396187
            #chr                                                                 X
            strand                                                               -
            sgRNA_sequence                                 CCGAGTGTGTTAATGAAAAATAA
            Gene_ID                                                    FBgn0004170
            Transcript_ID                                              FBtr0070073
            Chromosome                                                           X
            Gene_Region                                                start_codon
            Start                                                           396177
            Stop                                                            396179
            Strand                                                               +
            Reference_Seq        ACCTGCGATAATTTGACATTCTTAGAAACTACCTGAAGAAATAGGA...
            upstreamHA           TCTGGTCAGTGCCATACCCCTTGGTGTATACTTGCGAGTCTTAATT...
            downstreamHA         AACACACTCGGAGCTTTCTTTAACTTTCCGGATAACGATCAACAGA...
            positionScore                                                        8
            PAM_in_start/stop                                                False
            <15_bp3’_overhang                                                 True
            PAM_in_CDS                                                       False
            PAM_outside_CDS                                                   True
            CDS_boundary                                                       >22
            lastG                                                              8.0
            cutSite                                                           13.0
            mutated                                                          False

    """
    import pandas as pd

    #Mutated starts as false
    mutated = False

    #Codon fragmenter - codons are now in mutable format in a list, from last 21bp of HAL, the start/stop site, then the first 21bp of HAR. RevComp if gene is -
    mutableCodons = codonFragmenter(winnerdf)

    ## Coordinate information per codon of mutableCodons. This accounts for gene strand.
    if winnerdf["Strand"] == '+':
        codonCoordinates = pd.read_excel("inputfiles/fmaxStopScore.xlsx", sheet_name="CodonCoordinatePlus", index_col= 0)
    elif winnerdf["Strand"] == '-':
        codonCoordinates = pd.read_excel("inputfiles/fmaxStopScore.xlsx", sheet_name="CodonCoordinateMinus", index_col= 0)

    #Define CDS boundary range
    if ">" in winnerdf["CDS_boundary"]:
        CDSboundary = [int(winnerdf["CDS_boundary"][1:])+1, 43]
    elif "<" in winnerdf["CDS_boundary"]:
        CDSboundary = [0,int(winnerdf["CDS_boundary"][1:])]

    #Extract PAM last G coordinate
    lastG = winnerdf['lastG'] #this is a C if on minus strand

    #Coordinates for sgRNA bases sgRNA1, sgRNA2, sgRNA3, sgRNA4, sgRNA5, sgRNA6
    #Any of these coordinates that are outside of the CDS will be removed, leaving only sgRNA coordinates in CDS.
    #We also define the second PAM base based on sgRNA +/- here - secondG (again, might be c)
    sgRNAcoordinates = []
    if winnerdf['strand'] == '+': #sgRNA on + strand
        secondG = lastG-1
        for shift in range(3, 9):
            sgRNAcoordinates.append(lastG - shift)
    else:
        secondG = lastG+1
        for shift in range(3, 9): #sgRNA on - strand
            sgRNAcoordinates.append(lastG + shift)

    sgRNAcoordinates = [x for x in sgRNAcoordinates if x in range(CDSboundary[0], CDSboundary[1])]

    #A. Mutate PAM in CDS (1 mutation)
    if winnerdf["PAM_in_CDS"] == True: #check condition
        #mutate last G if possible using relative coordinates
        codonNumber = codonCoordinates.at[lastG, 'codon'] #this is the codon number (as an index in mutableCodons)
        codon = mutableCodons[codonNumber] #This is the codon we want to mutate
        base = codonCoordinates.at[lastG, 'base'] #This is the base within that codon (1-3)
        potentialCodons = find_synonymous_codons(query_codon =codon, base_to_change= base)
        if len(potentialCodons) >0:
        #check for NGA and remove this
            for cod in potentialCodons:
                if cod[:-1] == "A":
                    potentialCodons.remove(cod)
            #If we still have replacement codons left, swap the first one in:
            if len(potentialCodons) >0:
                mutableCodons[codonNumber] = potentialCodons[0]
                mutated = True
        else: #if no mutation is found, try the second G.
            codonNumber = codonCoordinates.at[secondG, 'codon']
            codon = mutableCodons[codonNumber]
            base = codonCoordinates.at[secondG, 'base']
            potentialCodons = find_synonymous_codons(query_codon = codon, base_to_change = base)
             #If we have replacement codons, swap the first one in:
            if len(potentialCodons) >0:
                mutableCodons[codonNumber] = potentialCodons[0]
                mutated = True
    
    #B. Mutate sgRNA in CDS (up to 2 mutations)
    #Check if mutation was successful in A, then check if there are mutable bases of the sgRNA in CDS 
    if mutated == False and len(sgRNAcoordinates) != 0:
        print("Condition B")
        sgMutatedBases = 0 #count of bases mutated
        excludeCodons = [] #this is the codon number of codons we've already mutated. This ensures that we won't re-mutate already mutated codons.
        index = 0 #this is the index of the sgRNA coordinate we are trying to mutate
        while sgMutatedBases < 2 and index <= len(sgRNAcoordinates): #while 2 mutations have not yet been made and the we've not uet used all coordinates in the sgRNA coordinate list
            print(f"mutating base {index}, mutations so far are {sgMutatedBases}")
            coordinate = sgRNAcoordinates[index] #This is the coordinate of the sgRNA base
            codonNumber = codonCoordinates.at[coordinate, 'codon']
            codon = mutableCodons[codonNumber] #this is the codon number that this coordinate is in
            base = codonCoordinates.at[coordinate, 'base'] #this is the base within that codon that the coordinate corresponds to (1-3)
            
            #check if we've already mutated this codon
            if codonNumber not in excludeCodons: #if not, attempt to mutate this base
                potentialCodons = find_synonymous_codons(query_codon = codon, base_to_change = base)
                if len(potentialCodons) > 0: #if mutation is possible, accept this mutation
                    mutableCodons[codonNumber] = potentialCodons[0] #replace this codon in the fragmenter list
                    excludeCodons.append(codonNumber) #add this to already mutated codons
                    sgMutatedBases +=1 #add 1 to count of mutated bases
            
            index += 1 #move to next sgRNA base

        if sgMutatedBases > 0: #if we've mutated at least 1 sgRNA base, set mutated to True
            mutated = True

    #C. Mutate PAM outside CDS - to NTG or CTN
    if mutated == False:
        codonNumber = codonCoordinates.at[secondG, 'codon']
        codon = mutableCodons[codonNumber]
        base = codonCoordinates.at[secondG, 'base']
        newCodon = codon[:base-1] + 'T' + codon[base:] #add a 'T' where the old base was in this codon
        mutableCodons[codonNumber] = newCodon
        mutated = True

    #Replace the potentially mutated homology arms back in the upstreamHA and downstreamHA of our winnerdf
    winnerdfmutated = codonReverseFragmenter(mutableCodons, winnerdf)

    #Indicate if mutation has occurred within winner df
    winnerdfmutated['mutated'] = mutated

    return winnerdfmutated

def sgRNArunner(inputfile, window = 21):
    """
    Check functions so far work with the proposed changes.
    params:
        inputfile:
        minDistance: window size left of region of interest to search for sgRNAs in, has to be divisible by 3, minimum -42
        maxDistance window size right of region of interest to search for sgRNAs in, has to be divisible by 3, maximum 42
    """
    import pandas as pd
    #set up reference sequence Bio SeqIO element
    refSeqPerChromosome = refSeq()

    #Make the dataframe for all transcription factor start/stop sites
    TFsdf = make_dataframe_from_TFs_list(inputfile, refSeqPerChromosome)

    #Create fragments for the HDR-arm and append as new column to dataframe
    TFsdf = make_homology_arm_fragments(TFsdf, refSeqPerChromosome)
    
    #set up the output DF - will contain a winning sgRNA per site in TFsdf
    TFsdfWinnersandMutated = pd.DataFrame(columns=["fmin", "fmax", "strand", "sgRNA_sequence", "Gene_ID",
                                                   "Transcript_ID", "Chromosome", "Gene_Region", "Start", "Stop",
                                                   "Strand", "Reference_Seq", "upstreamHA", "downstreamHA", "positionScore",
                                                    "PAM_in_start/stop", "<15_bp3’_overhang", "PAM_in_CDS", "PAM_outside_CDS",
                                                    "CDS_boundary", "lastG", "cutSite", "mutated"])
    
    #Per row of transcription factor start/stop site dataframe, select a guideRNA
    for ind, row in TFsdf.iterrows():
        print(f"Selecting guide RNA for TFsdf row {ind}")

        filtered_sgRNA = gRNA_stringencyIterator(row, window, refSeqPerChromosome)

        #Derive min and max distance from stop position
        minDistance = (window + 3) * -1 # leave a 3bp gap to account for start/stop site
        maxDistance = window
        #Score the sgRNAs for this site
        sgRNAdf = sgRNApositionCheck(filtered_sgRNA, minDistance, maxDistance)

        #Add a column to establish whether mutation has occurred. This will be set to 'True' in the mutator function and starts as False by default.
        sgRNAdf['mutated'] = False

        #If no sgRNAs were found at any stringency
        if not sgRNAdf:
            print(f"No sgRNAs found at all stringencies for {ind}.")

            #Just input the information about this site into the final DF. The columns about the guideRNA will be filled with 'NaN', indicating no guideRNA could be found.
            for col in list(row.index):
                TFsdfWinnersandMutated.at[ind, col] = row[col]

        else: #If we have at least one sgRNA
            #Select winner
            winnersgRNA, mutationNeeded = find_best_gRNA(sgRNAdf) #winnersgRNA is a string of the winning sequence, mutationNeeded is a bool variable
            print(winnersgRNA)
            winnerdf = sgRNAdf[sgRNAdf["sgRNA_sequence"] == winnersgRNA] #This is the dataframe row for the winning sgRNA from the original sgRNA dataframe
            winnerdf.reset_index(inplace= True, drop= True)
            winnerdf = winnerdf.loc[0] #This is the pandas series for the same information
            if mutationNeeded == True: #run mutator if indicated
                winnerdf = find_best_mutation(winnerdf, maxDistance) #this will mutate HAL or HAR as needed and return the original DF with mutated sequences, and the mutated column set to True
        
            #Add the winning sgRNA into the output df for this TF start/stop site
            TFsdfWinnersandMutated.loc[ind] = winnerdf
    
        #Return dataframe as an excel file (This should be unindented one, but while testing I'd like to see the output file updated after each row of the TFsdf)
        TFsdfWinnersandMutated.to_excel("outputFiles/winningsgRNAs.xlsx")

    return TFsdfWinnersandMutated