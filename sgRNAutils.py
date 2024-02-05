"""
Author: Marina Luchner, Emma Watts, Giulia Biasi
Email: marina.luchner@eng.ox.ac.uk, emma.watts@stx.ox.ac.uk, giulia.biasi@linacre.ox.ac.uk
"""

#Complete functions
def revComp(inputSeq):
  """
  Creates the reverse complement of a DNA sequence.

  Input: 
        inputSeq: DNA sequence, str
  Output: 
        revComp: reverse complement of inputSeq, str

  """
  complement = {'A': 'T', 'C': 'G', 'G': 'C', 'T': 'A'}
  
  revComp = ""
  for base in inputSeq[::-1]:
    revComp += complement[(base.upper())]

  return revComp

def refSeq(fastaFile):
    """
    Creates a dictonary listing the sequence of every chromosome.
    
    Input: 
        fastaFile: the file path to the reference genome in fasta format, str

    Output:
        refSeqPerChromosome: list of chromosomes with DNA sequence, dictionary of format:
        {   '2L': Seq('CGACAATGCACGACAGAGGAAGCAGAACAGATATTTAGATTGCCTCTCATTTTC...GAG'), 
            '2R': Seq('CTCAAGATACCTTCTACAGATTATTTAAAGCTAGTGCACAACAACAATAAATTG...TTC'), 
            ...
        }
    """
    from Bio import SeqIO

    refSeqPerChromosome = {}
    for seq in SeqIO.parse(open(fastaFile), 'fasta'): #This is the FASTA file of the reference genome sequence
        refSeqPerChromosome[seq.id] = seq.seq 

    return refSeqPerChromosome

def make_dataframe_from_TFs_list(TF_list, refSeqPerChromosome, annotation):
    '''
    Creating a dataframe from sequence information and genes of interest. 

    params:
        TF_list: excel file listing the genes of interest
        refSeqPerChromosome: reference chromosome stored in SeqIO sequence format, indexed by refSeqPerChromosome[seq.id] = seq.seq
        annotation: gtf file containing the gene annotation information for the reference genome

    returns: pandas dataframe of information about each annotated termini (start/ stop sites) of genes of interest, in format:
        Gene_ID         Transcript_ID   Chromosome  Gene_Region     Start      Stop     Strand                                      Reference_Seq 
        FBgn0000022     FBtr0070072          X      start_codon     370094     370096   +           CTAATGAATAGATTGGTGTGTGATGTAGTGATCTAATATGGTGAAG... 
        FBgn0000022     FBtr0070072          X      stop_codon      370697     370699   +           GACATTGTGTCGTTCGTATGTCGCCCATTGAGACCGCCAATGGAGG... 
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
    and appends homology arm sequences. Homology arms are defined as 225 bp left and right from the start/stop site on the (+) strand.

    input: TFsdf: pandas dataframe of information about start/stop codon.
           refSeqPerChromosome: reference chromosome stored in SeqIO sequence format, indexed by refSeqPerChromosome[seq.id] = seq.seq
    output: TFsdf pandas dataframe appended with HAL and HAR of the format:
        Gene_ID Transcript_ID Chromosome  Gene_Region   Start    Stop Strand  \
        0  FBgn0000022   FBtr0070072          X  start_codon  370094  370096      +   
        1  FBgn0000022   FBtr0070072          X   stop_codon  370697  370699      +   

                                            Reference_Seq  \
        0  CTAATGAATAGATTGGTGTGTGATGTAGTGATCTAATATGGTGAAG...   
        1  GACATTGTGTCGTTCGTATGTCGCCCATTGAGACCGCCAATGGAGG...   

                                                        HAL  \
        0  TACTACCTCTCTATTAAAATCAGAGAAAACACTCATCTCAAGAGAC...   
        1  CACCAAGAGTTGCAGTTGCAATCTCCAACTGGCAGCACAAGTTCCT...   

                                                        HAR  
        0  GCTTTGGGCAGCGAAAATCACTCTGTTTTCAACGACGACGAGGAGT...  
        1  AAAAACAGATCAAATCTTCAGCTATTGCTAGTCGCACCCAACCATA...  

    """
    for index, row in TFsdf.iterrows():

            #Define 225 bp region left and right from Start and Stop of gene region
            leftStart = row["Start"] - 224 -1
            leftStop = row["Start"]  - 1 
            rightStart = row["Stop"] + 1 
            rightStop = row["Stop"] + 224 + 1
            
            #Add left and right homology arm (HA)
            HAL = str(refSeqPerChromosome[row["Chromosome"]][leftStart-1:leftStop])
            HAR = str(refSeqPerChromosome[row["Chromosome"]][rightStart-1:rightStop])
            TFsdf.at[index,"HAL"] = HAL
            TFsdf.at[index, "HAR"] = HAR

    return TFsdf

def filter_gRNA(gRNA_file, window, tfSingleRow, refSeqPerChromosome):
    """
    Selects sgRNAs within given window around the start/stop site and returns these in a pandas dataframe with information about their coordinates and the start/stop site.
    
    params:
        gRNA_file: gff file
        window: number of bp left and right from the start/stop site to search for gRNAs in
        tfSingleRow: a pandas Series for one row of the TFsdf dataframe produced by make_dataframe_from_TFs_list() function
        refSeqPerChromosome: reference chromosome stored in SeqIO sequence format, indexed by refSeqPerChromosome[seq.id] = seq.seq
        
    returns: filtered_gRNAs: pandas dataframe of sgRNAs for one start/stop site of the format:
       	fmin	fmax	#chr	strand	sgRNA_sequence	        Gene_ID	    Transcript_ID	Chromosome	Gene_Region	Start	Stop	\
        370082	370104	X	    +	    CTATCTCTTAAAATGGCTTTGGG	FBgn0000022	FBtr0070072	    X	        start_codon	370094	370096	\
        
        Strand	Reference_Seq	                                    HAL	                                                HAR
        +	    CTAATGAATAGATTGGTGTGTGATGTAGTGATCTAATATGGTGAAG...	TACTACCTCTCTATTAAAATCAGAGAAAACACTCATCTCAAGAGAC...	GCTTTGGGCAGCGAAAATCACTCTGTTTTCAACGACGACGAGGAGT...
    """
    import pandas as pd

    # create a dataframe from the gRNA file
    gRNAFileAnnotation = pd.read_csv(gRNA_file, sep = "\t", index_col = False)

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

def gRNA_stringencyIterator(tfSingleRow, window, refSeqPerChromosome, sgRNAFolder):
    """
    Iterates trough gRNA stringency files until at least one gRNA is returned from filter_gRNAs() function.
    Function depends on filter_gRNA() function.

    params:
        tfSingleRow: a pandas Series for one row of the TFsdf panas dataframe produced by make_dataframe_from_TFs_list
        refSeqPerChromosome: a Bio SeqIO dictionary for the D. melanogaster genome
        window: number of bp left and right of the start/stop codon to search for gRNAs in
        sgRNAFolder: folder containing gRNA files for different stringencies

    returns:
        sgRNA: pandas dataframe of gRNAs of a particular stringecy for a given start/stop codon in the format:
        fmin	fmax	#chr	strand	sgRNA_sequence	        Gene_ID	    Transcript_ID	Chromosome	Gene_Region	Start	Stop	\
        370082	370104	X	    +	    CTATCTCTTAAAATGGCTTTGGG	FBgn0000022	FBtr0070072	    X	        start_codon	370094	370096	\
        
        Strand	Reference_Seq	                                    HAL	                                                HAR
        +	    CTAATGAATAGATTGGTGTGTGATGTAGTGATCTAATATGGTGAAG...	TACTACCTCTCTATTAAAATCAGAGAAAACACTCATCTCAAGAGAC...	GCTTTGGGCAGCGAAAATCACTCTGTTTTCAACGACGACGAGGAGT...
    """
    import pandas as pd
    import os
    
    gRNAfiles = [os.path.join(sgRNAFolder, f) for f in os.listdir(sgRNAFolder)]

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
    Given filtered sgRNAs for a start/stop site, creates pandas dataframe containing positional information and condition checks
    for each sgRNA.

    params: 
        minDistance, maxDistance: minimum and maximum distances from the start/stop site that the sgRNA can be in, must be divisible by 3
        sgRNAdf: a pandas dataframe of sgRNAs for one start/stop site
        
    output: 
        sgRNAdf: pandas dataframe of sgRNAs for one start/stop site, appended with positional information and condition checks in format:
        	fmin	fmax	#chr	strand	sgRNA_sequence	        Gene_ID	    Transcript_ID	Chromosome	Gene_Region	Start	Stop	Strand \   
            9792987	9793009	2L	    +	    CACAGTTTACGCAGGTCCATGGG	FBgn0032150	FBtr0079917	    2L	        start_codon	9793004	9793006	-	   \
            9792986	9793008	2L	    +	    GCACAGTTTACGCAGGTCCATGG	FBgn0032150	FBtr0079917	    2L	        start_codon	9793004	9793006	-	   \
           
            Reference_Seq	                                    HAL	                                                HAR	\                                                                                                
            AATTTAATTTTTTTTATAATTAATTTAGTGCTAATCTTTGAGCAGC...	AATCATTACAGATCATGGGCAGCTCCTCAGTAAGATTAAGTGCTAT...	GGGTTTTATTATTTAATTAATGTAAATAAACTGTAATGTTAATGTT...\	
            AATTTAATTTTTTTTATAATTAATTTAGTGCTAATCTTTGAGCAGC...	AATCATTACAGATCATGGGCAGCTCCTCAGTAAGATTAAGTGCTAT...	GGGTTTTATTATTTAATTAATGTAAATAAACTGTAATGTTAATGTT...\
            
            positionScore	CDS_side	non_CDS_side	PAM_in_start/stop	max_15_bp_3’_overhang	PAM_in_CDS	PAM_outside_CDS	SRS_in_CDS\
            3	            HAL	        HAR	            False	            True	                False	    True	        True	\  
            2	            HAL	        HAR	            True	            False	                False	    False	        True	\   
    
            CDS_boundary	                                    non_CDS_boundary	                                SRS_boundary	            mutable_PAM	    cut_site         
            [-24, -23, -22, -21, -20, -19, -18, -17, -16, ...	[1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14...	[-8, -7, -6, -5, -4, -3]	[-1, 0]	        -2
            [-24, -23, -22, -21, -20, -19, -18, -17, -16, ...	[1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14...	[-8, -7, -6, -5, -4, -3]	[-1, 0]	        -3
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
        sgRNAdf.loc[ind, "cut_site"] = int(conditions.at[0, "cut_site"]) + sgRNA["positionScore"]

        # Add CDS side to sgRNA catalogue
        sgRNAdf.loc[ind, "CDS_side"] = conditions.at[0, "CDS_side"]

        # Add non CDS side to sgRNA catalogue
        sgRNAdf.loc[ind, "non_CDS_side"] = conditions.at[0, "non_CDS_side"]

    return sgRNAdf
    
def checkCDSCutandOrder(sgRNAdf):
    """
    Given an sgRNAdf, will calculate firstly rows where the sgRNA cuts inside CDS. If multiple, selects that which cuts closest to start/stop.
    If none, selects closest cut sgRNA that is outside CDS.
    
    params:
        sgRNAdf: pandas dataframe of sgRNAs for one start/stop site as produced by sgRNApositionCheck() function
       
    output: sgRNA sequence, str
    """
    import ast
    #extract CDS boundary value
    sgRNAdf.reset_index(inplace= True, drop= True)
    #Calculate sgRNA/cut site difference
    sgRNAdf["cutsite-CDSbound"] = abs(sgRNAdf["cut_site"])

    #C. Check cut site in CDS
    #filter for only sgRNAs that cut in CDS
    conditionC = sgRNAdf[sgRNAdf.apply(lambda row: row['cut_site'] in row['CDS_boundary'], axis=1)]

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
    Selects the best sgRNA from a dataframe of sgRNAs based on condition A-C (see PipelineLogic.pptx).

    params:
        sgRNAdf: pandas dataframe of sgRNAs for one start/stop site as produced by sgRNApositionCheck() function
        
    output:
        winnersgRNA: the sequence of the best guideRNA (str)
        mutationNeeded: True/False of whether this sgRNA needs to be mutated (bool)
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
        conditionB = sgRNAdf[sgRNAdf["max_15_bp_3’_overhang"] == True]
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
    For the CDS of each start/stop, will create a list of codons. If gene is on - strand, will reverse complement each codon individually.

    params:
        side: "HAL" or "HAR" depending on where base to mutate, str
        side_boundary: list of integers, positionScores for bp that are located in the CDS
        winnerdf: pandas Series for the row of the sgRNA dataframe containing information about the winner sgRNA. Format of sgRNAdf
        as produced by sgRNApositionCheck() function.

    returns:
        codonList: list, codons from CDS respective to the region of interest in mutable format. If gene is on - strand, codons are individually reverse complemented.
        winnerdf_copy: same as input pandas series, now with HAR_HA or HAL_HA column added. If CDS is in HAR, the HAR_HA column stores a list of bp positions relative 
        to the start of the CDS. If CDS is in HAL, the HAL_HA column stores a list of bp positions relative to the end of the CDS.
        
        fmin                                                               9792986
        fmax                                                               9793008
        #chr                                                                    2L
        strand                                                                   +
        sgRNA_sequence                                     GCACAGTTTACGCAGGTCCATGG
        Gene_ID                                                        FBgn0032150
        Transcript_ID                                                  FBtr0079917
        Chromosome                                                              2L
        Gene_Region                                                    start_codon
        Start                                                              9793004
        Stop                                                               9793006
        Strand                                                                   -
        Reference_Seq            AATTTAATTTTTTTTATAATTAATTTAGTGCTAATCTTTGAGCAGC...
        HAL                      AATCATTACAGATCATGGGCAGCTCCTCAGTAAGATTAAGTGCTAT...
        HAR                      GGGTTTTATTATTTAATTAATGTAAATAAACTGTAATGTTAATGTT...
        positionScore                                                            2
        CDS_side                                                               HAL
        non_CDS_side                                                           HAR
        PAM_in_start/stop                                                     True
        max_15_bp_3’_overhang                                                False
        PAM_in_CDS                                                           False
        PAM_outside_CDS                                                      False
        SRS_in_CDS                                                            True
        CDS_boundary             [-24, -23, -22, -21, -20, -19, -18, -17, -16, ...
        non_CDS_boundary         [1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14...
        SRS_boundary                                      [-8, -7, -6, -5, -4, -3]
        mutable_PAM                                                        [-1, 0]
        cut_site                                                                -3
        HAL_HA                   [-22, -21, -20, -19, -18, -17, -16, -15, -14, ...
    """
    import numpy as np
    
    #Extract HA sequence where CDS
    HA = winnerdf[side]

    # create a copy of winnerdf to not trigger the "SettingWithCopyWarning"
    winnerdf_copy = winnerdf.copy()

    if side == "HAR":
        winnerdf_copy[f"{side}_HA"] = np.array(side_boundary) - 1 #HA is 1bp to the right of stop
        mutable_region = HA[:winnerdf_copy[f"{side}_HA"][-1]+1]
    
    if side == "HAL":
        winnerdf_copy[f"{side}_HA"] = np.array(side_boundary) + 2 #HA is 3pb to the left of stop but because we index backwards we calculate - 1 
        mutable_region = HA[winnerdf_copy[f"{side}_HA"][0]+1:]
    
    #Start codon list
    codonList = []

    #Chop up the CDS into codons
    for codonBase1 in range(0, len(mutable_region), 3):
        codon = mutable_region[codonBase1:codonBase1+3]    
        if winnerdf["Strand"] == '-': #Create reverse complement of every codon if gene is on - strand
            codonList.append(revComp(codon))
        else: codonList.append(codon)

    return codonList, winnerdf_copy

def codonReverseFragmenter(codonsList, winnerdf, side, side_boundary):
    """
    Will recombine a fragmented list of codons of CDS around the region of interest to a string. Will then replace the mutated HAL/HAR arms in TFsdf.
    If gene is on - strand, reverse complemented codons will be translated back to the + strand sequences.

    params:
        codonList: codons window region around the start/stop site in mutable format (revComp if gene is -), list
        winnerdf: pandas series of the winning sgRNA taken from one row of the sgRNAdf produced by condonFragmenter() function.
        side: "HAL" or "HAR" depending on where base to mutate, str
        side_boundary: list of integers, positionScores for bp that are located in the CDS
    
    output: 
        winningdf_copy: same as input pandas series, now with mutated HAL and HAR in the format:
        
        fmin	fmax	#chr	strand	sgRNA_sequence	        Gene_ID	    Transcript_ID	Chromosome	Gene_Region	Start	Stop	Strand \   
        9792987	9793009	2L	    +	    CACAGTTTACGCAGGTCCATGGG	FBgn0032150	FBtr0079917	    2L	        start_codon	9793004	9793006	-	   \
        
        Reference_Seq	                                    HAL	                                                HAR	\                                                                                                
        AATTTAATTTTTTTTATAATTAATTTAGTGCTAATCTTTGAGCAGC...	AATCATTACAGATCATGGGCAGCTCCTCAGTAAGATTAAGTGCTAT...	GGGTTTTATTATTTAATTAATGTAAATAAACTGTAATGTTAATGTT...\	
        
        positionScore	CDS_side	non_CDS_side	PAM_in_start/stop	max_15_bp_3’_overhang	PAM_in_CDS	PAM_outside_CDS	SRS_in_CDS \
        3	            HAL	        HAR	            False	            True	                False	    True	        True	\  
    
        CDS_boundary	                                    non_CDS_boundary	                                SRS_boundary	            mutable_PAM	    cut_site         
        [-24, -23, -22, -21, -20, -19, -18, -17, -16, ...	[1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14...	[-8, -7, -6, -5, -4, -3]	[-1, 0]	        -2
    """
    #return to + strand if the gene is on -
    if winnerdf["Strand"] == '-':
        for ind, codon in enumerate(codonsList):
            codonsList[ind] = revComp(codon)
    
    #recombine the codons of the CDS
    mutable_region = ''.join(codonsList)

    # make a copy of winnerdf to avoid the SettingWithCopyWarning
    winnerdf_copy = winnerdf.copy()

    #replace mutated HA
    if side == "HAR":
        winnerdf_copy.at["HAR"] = mutable_region + winnerdf["HAR"][winnerdf[f"{side}_HA"][-1]+1:]
    if side == "HAL":
        winnerdf_copy.at["HAL"] =  winnerdf["HAL"][:winnerdf[f"{side}_HA"][0]+1] + mutable_region
    return winnerdf_copy

def find_synonymous_codons(query_codon, base_to_change, codon_table_excel = "inputfiles/codon_table.xlsx"):

    '''
    Uses the codon table to select codons that encode for the same amino acid as the query codon. Will ensure the specified base has been mutated.

    Params:
        query_codon: string, codon to select synonymous codons for
        base_to_change: the base within the codon that needs to change. Numeric value from 1-3.
        codon_table_excel: string, path to an excel file that lists per codon which amino acid that codon encodes.
    
    Returns:
        synonymous_codons: list of strings, each string is a codon that encodes for the same amino acid as the query codon but has a change 
        in the base position specified.
    
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

def mutator(basePosition, fragmentedCDS, winnerdf, codonCoordinates, PAM = False):

    """
    (1) converts the position of the base that we would like to try and mutate into a codon position
    (2) finds synonymous mutations for the codon in the the codon position of the fragmented CDS calling synonymous_codons() function
    (3) picks a synonymous mutation that mutates the base that we would like to try and mutate if possible
    (4) returns the mutated fragemented CDS if mutation was possible or throws exemption

    params:
        basePosition: integer, base position within sgRNA relative to fmax
        fragmentedCDS: list, codons in direction of translation between gene region of interest and maxDistance or minDistance 
        winnwerdf: pandas series of the winning sgRNA as extracted from one row of the sgRNAsf. sgRNAdf format as produced by codonFragmenter function
        codonCoordinates: str, pandas dataframe relating the positionScore to codon position and base position within codon
        PAM: boolean, specifying we are trying to mutate the PAM or not
    
    returns:
        fragmentedCDS: list, codons in direction of translation between gene region of interest and maxDistance or minDistance, 
        with the codon at basePosition mutated if possible
    """
    from sgRNAutils import find_synonymous_codons
    pos_relative_to_CDS = basePosition + winnerdf["positionScore"] - winnerdf["CDS_boundary"][0] # explained in Positions section of MutationLogic.pptx
    codonNumber = codonCoordinates.at[pos_relative_to_CDS, 'codon'] #this is the codon number (as an index in mutableCodons)
    codon = fragmentedCDS[codonNumber] #This is the codon we want to mutate
    base = codonCoordinates.at[pos_relative_to_CDS, 'base'] #This is the base within that codon (1-3)
    potentialCodons = find_synonymous_codons(query_codon =codon, base_to_change= base)

    if PAM == True:
        #check for NGA and remove if present
        potentialCodons = [codons for codons in potentialCodons if codons[-2:] != "GA"]
        
    if potentialCodons:
        fragmentedCDS[codonNumber] = potentialCodons[0]
        
    else: 
        raise ValueError("No synonymous codons found.")
    
    return fragmentedCDS

def find_best_mutation(winnerdf):
    """
        In the case where a fragment needs to be mutated, will mutate in CDS (preferably PAM, if not then in the sgRNA). 
        If not possible, will mutate PAM outside of CDS.
        
        params:
            winnerdf: pandas series of the winning sgRNA as extracted from one row of the sgRNAsf. sgRNAdf format as produced by codonFragmenter() function.
        
        returns:
            winnerdf: same as input but with the "mutated?" column set to True and the mutated HA if a mutation was possible. Dataframe of the following format:
            
            fmin                                                              11278450
            fmax                                                              11278472
            #chr                                                                    2L
            strand                                                                   +
            sgRNA_sequence                                     CATGCAGCTGCATCCCAATGCGG
            Gene_ID                                                        FBgn0052831
            Transcript_ID                                                  FBtr0091683
            Chromosome                                                              2L
            Gene_Region                                                    start_codon
            Start                                                             11278451
            Stop                                                              11278453
            Strand                                                                   +
            Reference_Seq            TTAGCGGACCATTAGAAACACAATTGAGTTTGCATTGGGTGCTTTA...
            HAL                      GCGGGAGAGAGTCGAGAGTCATTCGTTTGTCTCACTCCTGCGAGTG...
            HAR                      CAGCTGCATCCCAATGCTGAGTCGCCATCGGGGTAAGTTGATTCTA...
            positionScore                                                           19
            CDS_side                                                               HAR
            non_CDS_side                                                           HAL
            PAM_in_start/stop                                                    False
            max_15_bp_3’_overhang                                                False
            PAM_in_CDS                                                            True
            PAM_outside_CDS                                                      False
            SRS_in_CDS                                                            True
            CDS_boundary             [1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14...
            non_CDS_boundary         [-24, -23, -22, -21, -20, -19, -18, -17, -16, ...
            SRS_boundary                                      [-8, -7, -6, -5, -4, -3]
            mutable_PAM                                                        [-1, 0]
            cut_site                                                                14
            mutated?                                                              True
            HAR_HA                   [0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13,...
        """
    import pandas as pd
    from sgRNAutils import codonFragmenter, codonReverseFragmenter, mutator

    if winnerdf["Strand"] == '+':
        codonCoordinates = pd.read_excel("inputfiles/fmaxStopScoreML.xlsx", sheet_name="CodonCoordinatePlus", index_col= 0)
    
    elif winnerdf["Strand"] == '-':
        codonCoordinates = pd.read_excel("inputfiles/fmaxStopScoreML.xlsx", sheet_name="CodonCoordinateMinus", index_col= 0)

    # Condition A: PAM in CDS
    if winnerdf["PAM_in_CDS"] == True and winnerdf["mutated?"] == False:
        fragmentedCDS, winnerdf = codonFragmenter(winnerdf, winnerdf["CDS_side"], winnerdf["CDS_boundary"])

        for basePosition in winnerdf["mutable_PAM"]:
            
            try: 
                fragmentedCDS = mutator(basePosition, fragmentedCDS, winnerdf, codonCoordinates, PAM=True)
                # make a copy of winner df to avoid "SettingWithCopyWarning"
                winnerdf_copy = winnerdf.copy()
                winnerdf_copy["mutated?"] = True
                winnerdf = winnerdf_copy
                break

            except ValueError: continue

        winnerdf = codonReverseFragmenter(fragmentedCDS, winnerdf, winnerdf["CDS_side"], winnerdf["CDS_boundary"])

    # Condition B: sgRNA in CDS
    if winnerdf["SRS_in_CDS"] == True and winnerdf["mutated?"] == False:
        fragmentedCDS, winnerdf = codonFragmenter(winnerdf, winnerdf["CDS_side"], winnerdf["CDS_boundary"])
        number_of_mutations = 0

        for basePosition in winnerdf["SRS_boundary"]:
            
            try: 
                fragmentedCDS = mutator(basePosition, fragmentedCDS, winnerdf, codonCoordinates)
                number_of_mutations += 1
                # make a copy of winner df to avoid "SettingWithCopyWarning"
                winnerdf_copy = winnerdf.copy()
                winnerdf_copy["mutated?"] = True
                winnerdf = winnerdf_copy
                if number_of_mutations == 2:
                    break
                else: ValueError("for mutating the SRS, two mutations are needed")

            except ValueError: continue
        
        winnerdf = codonReverseFragmenter(fragmentedCDS, winnerdf, winnerdf["CDS_side"], winnerdf["CDS_boundary"])
    
    # Condition C: PAM outside of CDS
    if winnerdf["mutated?"] == False:
        fragmentedRegion, winnerdf = codonFragmenter(winnerdf, winnerdf["non_CDS_side"], winnerdf["non_CDS_boundary"])
        pos_relative_to_region = winnerdf["mutable_PAM"][0] + winnerdf["positionScore"] - winnerdf["non_CDS_boundary"][0] # formula explained in MutationLogicML.pptx
        codonNumber = codonCoordinates.at[pos_relative_to_region, 'codon'] #this is the codon number (as an index in mutableCodons)
        codon = fragmentedRegion[codonNumber] #This is the codon we want to mutate
        base = codonCoordinates.at[pos_relative_to_region, 'base'] #This is the base within that codon (1-3)
        base_to_mutate_to = [i for i in ['A', 'C', 'G', 'T'] if i != codon[base-1]][0]
        newCodon = codon[:base-1] + base_to_mutate_to + codon[base:] #add a 'T' where the old base was in this codon
        fragmentedRegion[codonNumber] = newCodon
        winnerdf = codonReverseFragmenter(fragmentedRegion, winnerdf, winnerdf["non_CDS_side"], winnerdf["non_CDS_boundary"])
        winnerdf["mutated?"] = True

    return winnerdf

def sgRNArunner(inputfile, fastaFile, annotationFile, sgRNAFolder, window):
    """
    Will design sgRNAs and Homology Arms to tag every annotated termini (every start/stop codon) of given genes. Saves output as an excel file in ./outputFiles.
    
    params:
        inputfile: path to excel file listing gene IDs.
        fastaFile: path to the fasta file containing the reference genome.
        annotationFile: path to the gff file containing the annotation of the reference genome.
        sgRNAFolder: path to the folder containing the files listing sgRNAs at different stringencies.
        window: window size left and right of region of interest to search for sgRNAs in, has to be divisible by 3, maximum 42
        
    output:
        TFsdfWinnersandMutated: pandas dataframe of the designed sgRNAs for each annotated termini of given genes of the format as produced by the find_best_mutation() function.
    """
    import pandas as pd
    from datetime import datetime

    #set up reference sequence Bio SeqIO element
    refSeqPerChromosome = refSeq(fastaFile)

    #Make the dataframe for all transcription factor start/stop sites
    TFsdf = make_dataframe_from_TFs_list(inputfile, refSeqPerChromosome, annotationFile)

    #Create fragments for the HDR-arm and append as new column to dataframe
    TFsdf = make_homology_arm_fragments(TFsdf, refSeqPerChromosome)
    
    #set up the output DF - will contain a winning sgRNA per site in TFsdf
    TFsdfWinnersandMutated = pd.DataFrame(columns=["fmin", "fmax", "#chr", "strand", "sgRNA_sequence", "Gene_ID",
                                                   "Transcript_ID", "Chromosome", "Gene_Region", "Start", "Stop",
                                                   "Strand", "Reference_Seq", "HAL", "HAR", "positionScore", "CDS_side", 
                                                   "non_CDS_side", "PAM_in_start/stop", "max_15_bp_3’_overhang", "PAM_in_CDS", 
                                                   "PAM_outside_CDS", "SRS_in_CDS", "CDS_boundary", "non_CDS_boundary", 
                                                   "SRS_boundary", "mutable_PAM", "cut_site", "mutated?"])
    
    #Per row of transcription factor start/stop site dataframe, select a guideRNA
    for ind, row in TFsdf.iterrows():
        print(f"Selecting guide RNA for TFsdf row {ind}")

        filtered_sgRNA = gRNA_stringencyIterator(row, window, refSeqPerChromosome, sgRNAFolder)

        #Derive min and max distance from stop position
        minDistance = (window + 3) * -1 # leave a 3bp gap to account for start/stop site
        maxDistance = window
        #Score the sgRNAs for this site
        sgRNAdf = sgRNApositionCheck(filtered_sgRNA, minDistance, maxDistance)

        #Add a column to establish whether mutation has occurred. This will be set to 'True' in the mutator function and starts as False by default.
        sgRNAdf["mutated?"] = False

        #If no sgRNAs were found at any stringency
        if sgRNAdf.empty:
            print(f"No sgRNAs found at all stringencies for {ind}.")

            #Just input the information about this site into the final DF. The columns about the guideRNA will be filled with 'NaN', indicating no guideRNA could be found.
            for col in list(row.index):
                TFsdfWinnersandMutated.at[ind, col] = row[col]

        else: #If we have at least one sgRNA
            #Select winner
            winnersgRNA, mutationNeeded = find_best_gRNA(sgRNAdf) #winnersgRNA is a string of the winning sequence, mutationNeeded is a bool variable
            print("best sgRNA:", winnersgRNA)
            winnerdf = sgRNAdf[sgRNAdf["sgRNA_sequence"] == winnersgRNA] #This is the dataframe row for the winning sgRNA from the original sgRNA dataframe
            winnerdf.reset_index(inplace= True, drop= True)
            winnerdf = winnerdf.loc[0] #This is the pandas series for the same information
            if mutationNeeded == True: #run mutator if indicated
                winnerdf = find_best_mutation(winnerdf) #this will mutate HAL or HAR as needed and return the original DF with mutated sequences, and the mutated column set to True
        
            #Add the winning sgRNA into the output df for this TF start/stop site
            TFsdfWinnersandMutated.loc[ind] = winnerdf
    
    # Get the current date
    current_date = datetime.now()

    # Format the date as "year month day" to name the outputfile
    formatted_date = current_date.strftime("%Y %m %d")

    #Return dataframe as an excel file (This should be unindented one, but while testing I'd like to see the output file updated after each row of the TFsdf)
    TFsdfWinnersandMutated.to_excel(f"outputFiles/winning_sgRNAs_window_{window}_bp_{formatted_date}.xlsx")

    return TFsdfWinnersandMutated
