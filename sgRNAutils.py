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

def make_homology_arm_fragments(tfsDF):
    """
    Takes in the dataframe of information about start/stop codon regions, 
    and appends with columns for the 225 bp upstream and downstream of
    this site (the homlogy arm fragments).

    input: tfsDF - as defined in previous function and on GitHub README, with Reference_Seq 1600bp either side of the start/stop.
    output: tfsDF appended with upstreamHA and downstreamHA

    """

    tfsDF["upstreamHA"] = tfsDF.Reference_Seq.str[1375:1600]
    tfsDF["downstreamHA"] = tfsDF.Reference_Seq.str[1604:1829]

    return tfsDF

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
    Creating a dataframe from sequence information and genes of interest. Depends on functions revComp() and make_homology_arm_fragments(). 

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
        'Reference_Seq': 'TTGATCGTAGGACAC', 
        'upstreamHA': 'ATGCCTG', 
        'downstreamHA': 'CTGGATC'

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

    #Create fragments for the HDR-arm
    TFsdf = make_homology_arm_fragments(TFsdf)

    #re-index
    TFsdf = TFsdf.reset_index()
    del TFsdf["index"]

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
    changedCodons = []
    for codon in codon_list:
        if codon[base_to_change-1] != query_codon[base_to_change-1]:
            changedCodons.append(codon)
            
    return changedCodons

def filter_gRNA(gRNA_file, tfSingleRow, refSeqPerChromosome):
    """
    Selects gRNAs within 20pb of the start/stop site and returns these in a dataframe with information about their coordinates and the start/stop site.
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
    index = 0 #TODO@Marina improve this code

    # for each attribute value, extract the gene ID and symbol and add this to the new categories
    for attribute in gRNAFileAnnotation['attributes']:

        fullatt = (attribute).split(";")
        gRNAFileAnnotation.at[index,"target_site_variation"] = fullatt[8]
        index+=1
    
    # shorten file to essential information
    GenomeCoordinates = gRNAFileAnnotation[["target_site_variation", "fmin", "fmax", "#chr", "strand"]]

    #Filter for the appropriate conditions:
        # check whether the sgRNAs match the transgeneic genome, whether the sgRNAs match to the same chromosome as the transcription factor
        # select sgRNA that are located maximally 20 bp upstream of the start/stop codon of the transcription factors
        # select sgRNA that are located maximally 20 bp downstream of the start/stop codon of the transcription factors
    filtered_gRNAs = GenomeCoordinates[(GenomeCoordinates['target_site_variation'] == "target_site_variation=none observed") & 
                                                  (GenomeCoordinates['#chr'] == tfSingleRow["Chromosome"]) & 
                                                  (GenomeCoordinates['fmin']-1 >= tfSingleRow["Start"] - 20) & 
                                                  (GenomeCoordinates['fmax'] <= tfSingleRow["Stop"] + 20)]

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

def gRNA_stringencyIterator(tfSingleRow, refSeqPerChromosome):
    """
    Iterates trough gRNA stringency files until at least one gRNA is returned from filter_gRNAs function.

    params:
        tfSingleRow: a pandas Series for one row of the tfsDF dataframe produced by make_dataframe_from_TFs_list
        refSeqPerChromosome: a Bio SeqIO dictionary for the D. melanogaster genome

    returns:
        filtered_gRNAs: a df of gRNAs that meet the conditions at this site
    """

    gRNAfiles = ["inputfiles/Hu.2019.8.28.sgRNA_designs/NoOffTarget_high_stringency.gff", "inputfiles/Hu.2019.8.28.sgRNA_designs/NoOffTarget_med_stringency.gff", "inputfiles/Hu.2019.8.28.sgRNA_designs/NoOffTarget_low_stringency.gff", "inputfiles/Hu.2019.8.28.sgRNA_designs/1to3NonCdsOffTarget_low_stringency.gff", "inputfiles/Hu.2019.8.28.sgRNA_designs/ManyOffTarget_low_stringency.gff"]

    #Set up sgRNA variable
    sgRNA = ""
    file_ind = 0

    #Loop through files and retain the filtered sgRNA df if at least one row is present
    while len(sgRNA) == 0:
        sgRNA = filter_gRNA(gRNAfiles[file_ind], tfSingleRow, refSeqPerChromosome)
        print(f"stringency {file_ind}")
        file_ind +=1
        if file_ind == 4: #If we've tried all files and sgRNA is still 0, break out of loop
            break
    
    return sgRNA

def sgRNApositionCheck(sgRNAdf):
    """
    Given filtered sgRNAs for a start/stop site in the 'filtered_gRNAs' format, will create dataframe containing positional information and condition checks
    for each sgRNA.

    params: df: a pandas df of sgRNAs for one start/stop site:
            fmin	fmax	#chr	strand	sgRNA_sequence	        Gene_ID	    Transcript_ID	Chromosome	Gene_Region	Start	... downstreamHA
        0	370082	370104	X	    +	    CTATCTCTTAAAATGGCTTTGGG	FBgn0000022	FBtr0070072	    X	        start_codon	370094	...	TTAAGAGATAGTATAACGTTATTGTGTGACGATGCTCCTTGCTTCA...
        1	370693	370715	X	    -	    CCTGTAAAAAAACAGATCAAATC	FBgn0000022	FBtr0070072	    X	        start_codon	370694	...	TTAAGAGATAGTATAACGTTATTGTGTGACGATGCTCCTTGCTTCA...
        
    output: the same dataframe with added columns:
            fmin	fmax	#chr	strand	sgRNA_sequence	        Gene_ID	    Transcript_ID	Chromosome	Gene_Region	Start	...	downstreamHA	                                    positionScore	PAM_in_start/stop	<15_bp3’_overhang	PAM_in_CDS	PAM_outside_CDS	CDS_boundary	lastG	cutSite	
        0	370082	370104	X	    +	    CTATCTCTTAAAATGGCTTTGGG	FBgn0000022	FBtr0070072	    X	        start_codon	370094	...	TTAAGAGATAGTATAACGTTATTGTGTGACGATGCTCCTTGCTTCA...	8	            False	            True	            True	    False	        >22	            31.0	26.0	
        1	370693	370715	X	    -	    CCTGTAAAAAAACAGATCAAATC	FBgn0000022	FBtr0070072	    X	        start_codon	370694	...	TTAAGAGATAGTATAACGTTATTGTGTGACGATGCTCCTTGCTTCA...	16	            False	            False	            True	    False	        <20	            16.0	21.0	
                    
    """
    import pandas as pd

    #Adding position scores (fmax - stop)
    sgRNAdf["positionScore"] = sgRNAdf["fmax"] - sgRNAdf["Stop"]
    
    #Dataframe containing parameter ranges to interpret the positon score, based on gene strand, sgRNA strand, and start/stop
    positionScoreParameters = pd.read_excel("inputfiles/fmaxStopScore.xlsx")

    #Per parameter, append the sgRNAdf with a TRUE/FALSE value per condition.
    #Also add positional information about the CDS boundary, position of the last G/C, and the cut site from the conditions file.

    #Start by adding columns for these to sgRNAdf
    booleanColumns = ["PAM_in_start/stop", "<15_bp3’_overhang", "PAM_in_CDS", "PAM_outside_CDS"]
    positionColumns = ["CDS_boundary", "lastG", "cutSite"]
    sgRNAdf = sgRNAdf.reindex(columns = sgRNAdf.columns.tolist() + booleanColumns + positionColumns)
    
    #Iterating through sgRNAdf, extract the appropriate row of the conditions file and add new information to sgRNAdf.
    for ind, sgRNA in sgRNAdf.iterrows():

        #Extract the appropriate parameter row per sgRNA
        conditions = positionScoreParameters.loc[(positionScoreParameters["start/stop"] == sgRNA['Gene_Region']) & (positionScoreParameters["strand_type"] == sgRNA['Strand']) & (positionScoreParameters["sgRNA_strand"] == sgRNA['strand'])]
        conditions = conditions.reset_index(drop = True)
        
        #Per column, input true/false as to whether the position score meets that condition
        for col in booleanColumns:
            colValue = conditions.at[0,col] #extract parameter range values from dataframe
            #Process the value into a range (in format [min, max])
            #If the values should be 'more than' or 'less than', 25 is used as a max or -25 as min because distances cannot be more than 20
            if ">" in colValue: #could simplify this further by just defining all as ranges in excel
                minMax = [int(colValue[1:])+1, 25]
            elif "<" in colValue:
                minMax = [-25,int(colValue[1:])]
            elif ":" in colValue:
                min, max = colValue.split(":")
                minMax = [int(min), int(max)]
            else:
                print("Incorrect format of range value. Verify inputs.")

            #Into the output dataframe, print true/false as to whether the positionScore has met the condition for that column
            sgRNAdf.at[ind, col] = bool(sgRNA["positionScore"] in range(minMax[0], minMax[1]))
       
        #Calculate PAM last G or C position using positionScore. Add this to the sgRNA catalogue.
        sgRNAdf.at[ind, "lastG"] = sgRNAdf.at[ind, "positionScore"] + int(conditions.at[0, "lastG"])

        #Calculate cutsite coordinate using positionScore. Add this to the sgRNA catalogue.
        sgRNAdf.at[ind, "cutSite"] = sgRNAdf.at[ind, "positionScore"] + int(conditions.at[0, "cutSite"])

        #Keep CDS boundary position in output df
        sgRNAdf.at[ind, "CDS_boundary"] = conditions.at[0, "CDS_boundary"]
    
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
    CDSBoundary = sgRNAdf.at[0, 'CDS_boundary']

    #Calculate sgRNA/cut site difference
    sgRNAdf["cutsite-CDSbound"] = abs(sgRNAdf["cutSite"] - 20)
    
    #C. Check cut site in CDS
    
    #Define CDS boundary range
    if ">" in CDSBoundary: #could simplify this further by just defining all as ranges in excel
        CDSBoundaryRange = [int(CDSBoundary[1:])+1, 43]
    elif "<" in CDSBoundary:
        CDSBoundaryRange = [0,int(CDSBoundary[1:])]
    else:
        print("CDS boundary value not found/invalid")

    #filter for only sgRNAs that cut in CDS
    conditionC = sgRNAdf[sgRNAdf["cutSite"].between(CDSBoundaryRange[0], CDSBoundaryRange[1])]

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

    #B. sgRNA overhang is less than 15bp
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

def codonFragmenter(winnerdf):
    """
    For the start/stop site, will create a list of codons in the appropriate range where sgRNAs might be found (start/stop ±20 bp).
    If gene is on - strand, these will be revComp codons.

    params:
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
        codonList: codons from the ±21bp region around the start/stop site in mutable format (revComp if gene is -)
    """
    
    #Extract HA sequences - just the 21bp near start/stop
    HAL = winnerdf["upstreamHA"][-21:]
    HAR = winnerdf["downstreamHA"][0:21]

    #Start codon list
    codonList = []

    #Add HAL codons
    for codonBase1 in range(0, len(HAL), 3):
        codonList.append(HAL[codonBase1:codonBase1+3])

    #add start/stop - this is added as 'ATG' (even if stop) but will never be mutated. As HAL and HAR only will be replaced into tfsDF,
    #this will not affect original sequences.
    codonList.append('ATG')

    #Add HAR codons
    for codonBase1 in range(0, len(HAR), 3):
        codonList.append(HAR[codonBase1:codonBase1+3])
    
    #Define gene strand:
    if winnerdf["Strand"] == '-': #If on the minus strand, take revComp per codon
        for ind, codon in enumerate(codonList):
            codonList[ind] = revComp(codon)

    return codonList

def codonReverseFragmenter(codonsList, winnerdf):
    """
    Will take a fragmented list of codons of the start/stop site region and replace the mutated HAL/HAR arms in tfsDF.
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
    HAL = ""
    HAR = ""

    #return to + strand if the gene is on -
    if winnerdf["Strand"] == '-':
        for ind, codon in enumerate(codonsList):
            codonsList[ind] = revComp(codon)
    
    #replace HAL
    for codon in range(0, 7):
        HAL += codonsList[codon]

    #replace HAR
    for codon in range(8, 15):
        HAR += codonsList[codon]

    #replace mutated HAL and HAR
    winnerdf.at["upstreamHA"] = winnerdf["upstreamHA"][:-21] + HAL
    winnerdf.at["downstreamHA"] = HAR + winnerdf["downstreamHA"][21:]

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

def sgRNArunner():
    """
    Check functions so far work with the proposed changes.
    """
    import pandas as pd
    #set up reference sequence Bio SeqIO element
    refSeqPerChromosome = refSeq()

    #Make the dataframe for all transcription factor start/stop site infp
    TFsdf = make_dataframe_from_TFs_list("inputfiles/TFs.xlsx", refSeqPerChromosome)
    
    #set up the output DF - this is the winning sgRNA per site in TFsdf
    TFsdfWinnersandMutated = pd.DataFrame(columns=["fmin", "fmax", "#chr", "strand", "sgRNA_sequence", "Gene_ID",
                                                   "Transcript_ID", "Chromosome", "Gene_Region", "Start", "Stop",
                                                   "Strand", "Reference_Seq", "upstreamHA", "downstreamHA", "positionScore",
                                                    "PAM_in_start/stop", "<15_bp3’_overhang", "PAM_in_CDS", "PAM_outside_CDS",
                                                    "CDS_boundary", "lastG", "cutSite", "mutated"])
    
    #Per row of this dataframe, select a guideRNA
    for ind, row in TFsdf.iterrows():
        print(f"Selecting guide RNA for TFsdf row {ind}")

        filtered_sgRNA = gRNA_stringencyIterator(row, refSeqPerChromosome)

        #Score the sgRNAs for this site
        sgRNAdf = sgRNApositionCheck(filtered_sgRNA)

        #Add a column to establish whether mutation has occurred. This will be set to 'True' in the mutator function and starts as False by default.
        sgRNAdf['mutated'] = False

        #If no sgRNAs were found at any stringency
        if len(sgRNAdf) == 0:
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
                winnerdf = mutator(winnerdf) #this will mutate HAL or HAR as needed and return the original DF with mutated sequences, and the mutated column set to True
        
            #Add the winning sgRNA into the output df for this TF start/stop site
            TFsdfWinnersandMutated.loc[ind] = winnerdf
    
        #Return dataframe as an excel file (This should be unindented one, but while testing I'd like to see the output file updated after each row of the TFsdf)
        TFsdfWinnersandMutated.to_excel("outputFiles/winningsgRNAs.xlsx")

    return TFsdfWinnersandMutated

sgRNArunner()