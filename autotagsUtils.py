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

#Primer functions
def make_primers_dataframe_from_TFs_list(TF_list, ref_genome = "inputfiles/dmel-all-chromosome-r6.48.fasta", annotation = "inputfiles/dmel-all-r6.48.gtf"):
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
    
    #This is to reformat the "Attribute" category in refGenomeAnnotation, to extract Gene_ID, Gene_Symbol, Transcript ID, and Transcript Symbol
    index = 0

    #Add new categories to the dataframe
    refGenomeAnnotation = refGenomeAnnotation.assign(Gene_ID = "", Gene_Symbol = "", Transcript_ID = "", Transcript_Symbol = "")

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

    TFsdf = refGenomeAnnotation[["Gene_ID", "Gene_Symbol", "Transcript_ID", "Transcript_Symbol", "Chromosome", "Gene_Region", "Start", "Stop", "Strand"]].loc[refGenomeAnnotation["Gene_ID"].isin(queryTFsdf["Flybase_ID"])]

    #Add reference genome sequence per gene region
    #This will correspond to 1.6kb upstream and downstream of ATG/stop codon 
    TFsdf = TFsdf.assign(Reference_Seq = "")

    for index, rowcontents in TFsdf.iterrows():
        if rowcontents["Strand"] == "+":

            #Define 2.6kb gene region
            regionStart = rowcontents["Start"] - 1601
            regionStop = rowcontents["Stop"] + 1600

            #Add reference sequence
            TFsdf.at[index,"Reference_Seq"] = str(refSeqPerChromosome[rowcontents["Chromosome"]][regionStart:regionStop])

        if rowcontents["Strand"] == "-":

            #Define 2.6kb gene region
            regionStart = rowcontents["Start"] - 1601
            regionStop = rowcontents["Stop"] + 1600

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

def designPrimers(stringency, primer_type, template, useExtendedRegion = False):
 
  '''
  This function calculates a primers given the desired stringency and primer type.

  params:
    stringency: value from 0-2 for appropriate stringency conditions.
    primer_type: one of 'VAL-F', 'HAL-F', 'HAL-R', 'HAR-F', 'HAR-R', 'VAL-R'
    template: gene region of the start/stop codon (1700bp upstream and downstream)
    useExtendedRegion: if True, uses an extended region of the template for each primer.

  output:
    potential_primers = list of string(s) of the primer sequences.
    fullprimer3results = dataframe of results from primer3 with quality metrics of primer produced.
  
  '''
  import pandas as pd
  import primer3 as p3
  
  #Import the relevant specification files
  #Stringency Rules (stringency levels/index is 0-2)
  stringencyRulesdf = pd.read_excel("inputfiles/Primer_Stringency_Rules.xlsx")
  #Primer specifications
  primerNameandRegiondf = pd.read_excel("inputfiles/Primer_Name_and_Regions.xlsx")  
  #Reset index to primer_type
  primerNameandRegiondf = primerNameandRegiondf.set_index("primer_type", drop = True)

  #Extract appropriate stringency rules as dict
  stringencyRulesdf = stringencyRulesdf.loc[stringency]
  s = stringencyRulesdf.to_dict()

  #Extract necessary values from primerNameandRegiondf based on primer_type as dict
  primerNameandRegiondf = primerNameandRegiondf.loc[primer_type]
  p = primerNameandRegiondf.to_dict()

  #Check for shorter templates (if this gene region is near the end of the chromosome, might be missing a few bases)
  if len(template) == 3203 or primer_type in ['HAL-F','HAL-R','HAR-F','HAR-R']:
     extendedLength = p["extended_region_length"]
  else:
     missingRefSeqbases = 3203 - len(template)
     extendedLength = p["extended_region_length"] - missingRefSeqbases

  #Define whether to use the intial search region or extended search region
  if useExtendedRegion == False:
    gene_region = [p["initial_region_start"], p["initial_region_length"]]
  else:
    gene_region = [p["extended_region_start"], extendedLength]

  primer3Results = p3.bindings.designPrimers(
    {
        'SEQUENCE_TEMPLATE': template,
        'SEQUENCE_INCLUDED_REGION': gene_region, 
    },
    {
        'PRIMER_NUM_RETURN':p["primers_needed"],

        'PRIMER_TASK': "generic",
        'PRIMER_PICK_LEFT_PRIMER': p["left"],
        'PRIMER_PICK_INTERNAL_OLIGO': 0,
        'PRIMER_PICK_RIGHT_PRIMER': p["right"],

        'PRIMER_MIN_GC': s["GC_min"],
        'PRIMER_OPT_GC_PERCENT': s["GC_opt"],
        'PRIMER_MAX_GC': s["GC_max"],
     
        'PRIMER_MIN_SIZE': s["size_min"],
        'PRIMER_OPT_SIZE': s["size_opt"],
        'PRIMER_MAX_SIZE': s["size_max"],
     
        'PRIMER_MAX_END_GC': s["max_end_GC"],
     
        'PRIMER_GC_CLAMP': s["GC_clamp"],    
             
        'PRIMER_MAX_HAIRPIN_TH': s["TH_max_hairpin"],
    
        'PRIMER_MAX_POLY_X': s["max_polyx"],
    })
  
  #Filter the results dictionary for the keys explaining the number of primers returned
  primersReturnedCount = {k:v for (k,v) in primer3Results.items() if "NUM_RETURNED" in k}
  
  #return success = False (no primer found), or success = True, and the primers found
  if all(x == 0 for x in primersReturnedCount.values()) == True: #this is the output format if no primers are found
    success = False
    potential_primers = ["NA"]
  else: #if primers were found, compile these sequences into a list
    success = True
    potential_primers = []
    for columnName in primer3Results.keys():
      if columnName.endswith("SEQUENCE"):
        potential_primers.append(primer3Results[columnName])

  return success, potential_primers

def stringencyIterate(template, primer_type, useExtRegion = False, stringencyStart = 0):
  """
  Iterate through stringency levels to find primers. Primer mounting will be checked for "HAL-R" and "HAL-F".

  params:
    template: sequence of the full gene region as a string
    primer_type: type of primer being designed
    useExtRegion: if True, will use expanded search regions for primers.
    stringencyStart: what stringency to start the search at.
  
  output:
    stringency: stringency value for which a primer was found.
    potential_primers: sequence of the primer found. If none are found, returns NA.
    
  """
  success = False

  bestPrimerNotMounted = []

    #iterate through stringencies
  for stringency in range(stringencyStart, 3):
    success, potential_primers = designPrimers(stringency, primer_type, template, useExtendedRegion = useExtRegion)
    
    if primer_type == "HAL-R" or primer_type =="HAR-F": #check mounting for special case
      if success == True:
        bestPrimerNotMounted.append(potential_primers[0]) #keep the best primer of the first iteration

        #Check mounting of the returned primers
        potential_primers = mountedPrimers(primer_type, potential_primers, template)

        #If this is the last stringency, use makeMount to generate a mounted primer (first, try to extend best primer, otherwise mount to 30bp next to start/stop)
        if potential_primers == [] and stringency == 2:
            potential_primers = mountedPrimers(primer_type, bestPrimerNotMounted, template, extendMount = True)
            stringency = "extendMount"
            success = True 
        #Otherwise if not primer is correctly mounted, return false to move to the next stringency
        elif potential_primers == [] and stringency < 2:
            success = False

      #If no primer is found even at the last stringency for HAL-R/HAR-F, makemount
      elif success == False and stringency == 2:
            potential_primers = mountedPrimers(primer_type, potential_primers, template, makeMount = True)
            stringency = "makeMount"
            success = True 

    #if no primer obtained, re-design primer (stay in loop)
    if success == True:
       break

  #If still not returning a primer (for non-mount primers), the sequence will be NA
  if success == False:
    print("Could not find a primer at any stringency levels.") #don't want this every time
    potential_primers = ["NA"]
    stringency = "none"
  
  return stringency, potential_primers

def mountedPrimers(primer_type, potential_primers, template, makeMount = False, extendMount = False):
  """
  Perform special conditions check for primers that must be mounted to the stop/start site.

  params:
    primer_type: string variale specifying primer type (e.g. "HAL_R")
    potential_primers: list of candidate sequences for this primer
    template: 3200 bp region of the reference genome for this start/stop site
    makeMount = if True, will manually assign 25 bp from the intended primer start site.
    extendMount = if True, will take the best primer that was not mounted and extend it to start at the stop/start adjacent site.

  output:
    mountedPrimer: a single primer that is mounted at the start/stop-adjacent site.

  """
  import pandas as pd

  #Import relevant files and extract values for the primer_type from the dataframe:
  #Primer specifications
  primerNameandRegiondf = pd.read_excel("inputfiles/Primer_Name_and_Regions.xlsx")  
  #Reset index to primer_type
  primerNameandRegiondf = primerNameandRegiondf.set_index("primer_type", drop = True)
 
  primerNameandRegiondf = primerNameandRegiondf.loc[primer_type]
  p = primerNameandRegiondf.to_dict()

  mountedPrimer = []

  #If no special conditions specified, just check positioning of all potential primers and keep if mounted correctly.
  if makeMount == False and extendMount == False:
    for primer in potential_primers:
      if primer_type == "HAL-R": #This is a reverse primer, so need reverse complement to match to template
        searchPrimer = revComp(primer)
        stopIndex = template.find(searchPrimer) + len(searchPrimer)
        if stopIndex == p["initial_region_start"]+p["initial_region_length"]:
          mountedPrimer.append(primer)
      else:
        searchPrimer = primer
        startIndex = template.find(searchPrimer)
        if startIndex == p["initial_region_start"]:
          mountedPrimer.append(primer)
  
  #If no primer is found for all stringency levels then manually design.
  elif makeMount == True:
    mountedPrimer = template[p["initial_region_start"]: p["initial_region_start"]+p["initial_region_length"]] #Note here that these regions will not change even if the extended region is used
    if primer_type == "HAL-R": #This is a reverse primer, so need reverse complement to match to template
        mountedPrimer = revComp(mountedPrimer)
    mountedPrimer = [mountedPrimer]

  #In this case, the input is the best primer that is not mounted to the start site - extend this until the start.
  elif extendMount == True:
    bestPrimerNotMounted = potential_primers[0]
    #For HAL-R, find template region of the reverse complement and extend until desired 3' site.
    if primer_type == "HAL-R": #This is a reverse primer, so need reverse complement to match to template
        searchPrimer = revComp(bestPrimerNotMounted)
        startIndex = template.find(searchPrimer)
        mountedPrimer = template[startIndex: p["initial_region_start"]+p["initial_region_length"]]
        mountedPrimer = revComp(mountedPrimer)
    #For HAR-F, find primer within the template, and extend at the 5' end until desired site.
    else:
       searchPrimer = bestPrimerNotMounted
       stopIndex = template.find(searchPrimer) + len(searchPrimer)
       mountedPrimer = template[p["initial_region_start"]: stopIndex]
    mountedPrimer = [mountedPrimer]

  return mountedPrimer

def BLASTprimer(primer, chromosome):
  """
  BLAST checks a primer to the appropriate transgenic reference genome.

  params:
    primer: a single primer (str)
    chromosome: the chromosome that the template for the primers is on. e.g. 'X' (str)

  output: 
    BLASTcount: number of matches to the transgenic reference genome
  """
  from Bio import SeqIO

  #Upload transgenic reference file
  #This will becas_on_2 for chromosome 3, and cas_on_3 for all other chromosomes.
  if chromosome == 3:
      transgenicRefFile = SeqIO.parse("inputfiles/dmel6-nos-Cas9_on_2.fasta", "fasta")    
  else:
      transgenicRefFile = SeqIO.parse("inputfiles/dmel6-nos-Cas9_on_3.fasta", "fasta")    

#add the full genome (both strands)
  fullTransgenicRefSequence = ""
  for fasta in transgenicRefFile:
    currentSequence = fasta.seq
    fullTransgenicRefSequence += f"{currentSequence} "
    reverseComp = currentSequence.reverse_complement()
    fullTransgenicRefSequence += f"{reverseComp} "
  
  #If a primer is NA (none has been found) - don't BLAST and set value to 1 as this will need to be the accepted primer (there are no others).
  if primer == 'NA':
    print(f"Skipping BLAST for 'NA' primer, {primer}")
    BLASTcount = 1
  #Count instances of primer against the full transgenic reference
  else:
    BLASTcount = fullTransgenicRefSequence.count(primer)

  return BLASTcount

def findsixPrimers(template, forceExtend = False):
    """
    Calculates primers of each type for one site. Automatically uses extended regions if primers are not found for any one type in intial search.
    Extended regions can also be manually set to True if needed.

    params: 
        template: template of the gene region
        forceExtend: if True, the extended region will be used to calculate primers.

    output:
        sixPrimers: a dictionary of the six primers for this site, including metrics about the primer stringency, if extended region was used, if any were BLAST skipped, and a final primer space.
    """
    sixPrimers = {"VAL-F":[[],0,'F',0,''], "HAL-F":[[],0,'F',0,''], "HAL-R":[[],0,'F',0,''], "HAR-F":[[],0,'F',0,''], "HAR-R":[[],0,'F',0,''], "VAL-R":[[],0,'F',0,'']}

    #Calculate the six primers and record their parameters, iterating through stringencies.
    for primerType in sixPrimers.keys():
        stringency, potential_primers = stringencyIterate(template, primerType)
        sixPrimers[primerType][0] = potential_primers
        sixPrimers[primerType][1] = stringency
    
    #Use extended region if there's any NA primers in the result, or if specified in function input.
    for row in sixPrimers.values():
        if row[0] == ['NA'] or row[0] =='NA':
            forceExtend = True

    if forceExtend == True:   
        print('Using extended region.')     
        for primerType in sixPrimers.keys():
            stringency, potential_primers = stringencyIterate(template, primerType, useExtRegion= forceExtend)
            sixPrimers[primerType][0] = potential_primers
            sixPrimers[primerType][1] = stringency
            sixPrimers[primerType][2] = 'T'

    return sixPrimers

def BLASTiterator(sixPrimers, Chromosome, template, useextended = False):
    """
    BLAST through potential primers for each primer type. This will skip primers that do not match to the transgenic reference genome, and will increase stringency if none in potential_primers is BLAST accepted.
    If still no primer found and the extended region hasn't been used already, it will be used.

    params:
        sixPrimers: dictionary of potential primers and metrics in format sixPrimers = {"VAL-F":[[],0,'F',0,'']...}
        Chromosome: chromosome on which the gene region is found.
        template: gene sequence for this gene region.
        useextended: default is false. If true, primers will be recalculated with extended regions.

    returns:
        sixPrimers: in the same format as the input sixPrimers, but now BLAST adjusted.
    """

    #Iterate through initial primers and count matches for each primer type.
    for primerType, row in zip(sixPrimers.keys(),sixPrimers.values()):
      potential_primers = row[0]
      stringency = row[1]
      BLASTskip = 0
      matchCount = BLASTprimer(potential_primers[0], chromosome = Chromosome) #note that NAs are accepted as 1 in this function
      #For the case where matchCount is 0, need to iterate through the next primers in the list.
      #Will keep looping through and taking next primer along in potential_primers or increasing stringency to find a new list
      while stringency in [0,1,2] and matchCount == 0:
          print("In while loop")
          if BLASTskip<len(potential_primers)-1:
              BLASTskip +=1
              matchCount = BLASTprimer(potential_primers[BLASTskip], chromosome = Chromosome)
              print(f"new matchcount is {matchCount}")
          elif BLASTskip >= len(potential_primers)-1:
              stringency += 1
              stringency, potential_primers = stringencyIterate(template, primerType, stringencyStart= stringency, useExtRegion= useextended)
              print(f"using new stringency {stringency}")

      if matchCount == 0 and useextended == True:
        print("No primers after using all BLASTskip, stringencies, and extended region.")
        sixPrimers[primerType][4] = "NA"
        sixPrimers[primerType][3] = "ALL"
    
      if matchCount == 0 and useextended == False:
        sixPrimers[primerType][0] = "useExt"
        print("Using extended region as BLAST skips have been maxed out.") 

      if matchCount > 0:
        sixPrimers[primerType][3] = BLASTskip
        sixPrimers[primerType][4] = potential_primers[BLASTskip]

    return sixPrimers

def addSixPrimers(template, Chromosome):
  """
  Find sequences for each of the six primers needed for a start/stop site.

  params:
    template: sequence of the full gene region (str)
    Chromsome: chromosome on which the gene is found. (str)

  output:
    sixPrimers: dictionary of parameters and final primer for each primer type for that gene region.

  """
  import pandas as pd

  #Calculate an initial set of six primers
  sixPrimers = findsixPrimers(template) #this will already use the extended region if needed
  #BLAST check these
  BLASTEDsixPrimers = BLASTiterator(sixPrimers, Chromosome, template, useextended = False)
  
  #In the standard case, not using the extended region
  extend = False

  #If the BLASTiterator has indicated that no current primers at any stringency have a match to the transgenic ref, use the extended region and re-BLAST iterate.
  for row in BLASTEDsixPrimers.values():
    if row[0] == "useExt":
      extend = True

  if extend == True:
    sixPrimers = findsixPrimers(template, forceExtend= True)
    BLASTEDsixPrimers = BLASTiterator(sixPrimers, Chromosome, template, useextended = True)

    #Finally, remove the potential primers (keeping only the final picked primer)
    #Format at this point is {"VAL-F": [stringency, T/F of extended region used, BLAST skips, final primer sequence]}
  for primerType in BLASTEDsixPrimers.keys():
    BLASTEDsixPrimers[primerType] = BLASTEDsixPrimers[primerType][1:]

  return BLASTEDsixPrimers

def finalPrimersdf(TFsdf, outputFile, returnParameters = False):
  """
  Compile all primers for all transcription factor start/stop sites per isoform.

  params:
    TFsdf: full dataframe contaning transcript ID, gene region type (start/stop), chromosome, reference sequence (+-3200bp).
    outputFile: desired path and name for excel output file.
    returnParameters: if False, returns only primer IDs and primers. if true, returns parameters of the designed primers (strigency, whether the extended region was used, how many primers were skipped due to BLAST failure)

  output:
    output file will be created.

  """
  import pandas as pd
  
  #create a dataframe for primer inputs
  allPrimers = pd.DataFrame(columns = ["VAL-F", "HAL-F", "HAL-R", "HAR-F", "HAR-R", "VAL-R"])

  if returnParameters == True: #Additional metrics will be provided in the dataframe if requested.
    allPrimers[["Stringency Level", "Extended Region Used", "No. of BLAST Excluded"]] = ""

  #add primers for every row in the transcription factor dataframe
  for index, row in TFsdf.iterrows():
    print(f"now calculating index {index}")
    sixPrimersdict = addSixPrimers(row["Reference_Seq"], row["Chromosome"])
    finalsixPrimers = []
    stringenciesList = []
    extended = 'F'
    BLASTskipList = []
    
    for row in sixPrimersdict.values():
       finalsixPrimers.append(row[3])
       stringenciesList.append(row[0]) 
       extended = row[1]
       BLASTskipList.append(row[2])

    if returnParameters == False:
      allPrimers.loc[len(allPrimers)] = finalsixPrimers
    else:
      sixprimersandParameters = finalsixPrimers + [stringenciesList] + [extended] + [BLASTskipList]
      print(sixprimersandParameters)
      allPrimers.loc[len(allPrimers)] = sixprimersandParameters

  #Remove no longer needed TFdfs columns
  TFsdf = TFsdf.drop(columns = ["Gene_Symbol", "Chromosome", "Start", "Stop", "Strand", "Reference_Seq"])

  #Merge the dataframes
  TFsdf = TFsdf.join(allPrimers)

  TFsdf.to_excel(outputFile)

  return TFsdf

def primersRunner(TFs = 'inputfiles/TFs.xlsx', refgenome = 'inputfiles/dmel-all-chromosome-r6.48.fasta', annotation = 'inputfiles/dmel-all-r6.48.gtf', outputFile = "inputfiles/mockMaterials/TFsfullOutput.xlsx"):
    import pandas as pd

    TFsdf, TFsdict_of_dict = make_primers_dataframe_from_TFs_list(TFs, refgenome, annotation)

    TFsdf = finalPrimersdf(TFsdf, outputFile, returnParameters= True)

#Homology Arm Fragments

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

#sgRNA functions

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
            "sgRNA_strand" = '-'
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

            # if the gRNA lies on the minus strand than provide reverse complement sequence

            if GenomeCoordinates['strand'].iloc[gRNA] == '-':

                gRNA_string = revComp(gRNA_string)

            sgRNA_list_values.append(gRNA_string)

    TF_dict["sgRNA_list_positions"] = sgRNA_list_positions
    TF_dict["sgRNA_list_values"] = sgRNA_list_values
    TF_dict["sgRNA_strand"] = GenomeCoordinates['strand'].iloc[gRNA]


    return(TF_dict)

def positionScore(df):
    """
    Boolean evaluations of whether sgRNAs for a site meets required best-pick conditions and positional information about these sgRNAs.

    params:
        df: a dataframe for one start/stop site of the format: 
            df= {"start/stop":"start", #is it N or C termini -> do we need to look at start or stop codon for teh cut 
            'genome_start_codon_pos':403, 
            'genome_stop_codon_pos':406, # or only 1 of those 
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
    sgRNACatalogue["positionScore"] = sgRNACatalogue["sgRNA_stop"] - df["genome_stop_codon_pos"] #This is fmax - stop

    #Dataframe containing parameter ranges to interpret the positon score, based on gene strand, sgRNA strand, and start/stop
    positionScoreParameters = pd.read_excel("inputfiles/fmaxStopScore.xlsx")

    #Add the boolean columns to our output dataframe
    booleanColumns = ["P6inCDS", "Mutate2", "Mutate1", "OutsideCDS", "CutSiteInCDS", "PAMinCDS", "15bpOverhangL", "15bpOverhangR"]
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

def checkCDSCutandOrder(sgRNACatalogue):
    """
    Given an sgRNACatalogue or subset, will calculate firstly rows where the sgRNA cuts inside CDS. If multiple, selects that which cuts closest to start/stop.
    If none, selects closest cut sgRNA that is outside CDS.
    
    """
    #C. Check cut site in CDS
    conditionC = sgRNACatalogue[sgRNACatalogue["CutSiteInCDS"] == True]
    conditionCclosestCut = conditionC.sort_values('PAMRelativeEnd') #This value is the absolute value of positionScore, so indicates those that cut closest to start/stop

    if len(conditionCclosestCut) >1: #cuts in CDS, closest cut (C1, C2)
        conditionCclosestCut = conditionCclosestCut.reset_index(drop = True) #reset index
        winnersgRNA = conditionCclosestCut.at[0, "sgRNA_sequence"]
        winnerFound = True
    else: #not sgRNAs cut in CDS, select closest that still met condition B (C3)
        nonCDSclosestCut = sgRNACatalogue.sort_values('PAMRelativeEnd')
        winnersgRNA = nonCDSclosestCut.at[0, "sgRNA_sequence"]
        winnerFound = True

    return winnersgRNA, winnerFound

def find_best_gRNA(df):
    """
    params: df={
            "start/stop":"start",#is it N or C termini -> do we need to look at start or stop codon for teh cut 
            'genome_start_codon_pos':400, 
            'genome_stop_codon_pos':700, # or only 1 of those 
            'strand_type':'+',
            "sgRNA_list_positions":[[401,425],[456,467],[478,489],[395,415]],#those wil be as genome positions -assumptions - the coordinates correspond to the 1st and last bp of the strand to which the gsRNA will be complementary to
            "sgRNA_list_values":["AAGCGACTA","AAAAAAAATAAAAA","ATATATTTTTTTTTTAAAAA","AGCGCGAAATAATA"]
        }
    
    """
    #Score the sgRNAs for this site
    sgRNACatalogue = positionScore(df)

    #Set up the winning guide RNA
    winnerFound = False
    mutationNeeded = False
    winnersgRNA = "" #Could update this to be more information about the guide e.g. strand, position

    #If there is no guide RNAs at this site, flag this. The resulting sgRNA will be an empty string.
    if len(sgRNACatalogue) == 0:
        print(f"No guide RNAs were found for the start/stop site beginning at position {df['genome_start_codon_pos']}")
        winnerFound = True #To avoid looping through other conditions if there are no sgRNAs

    #A. Ideal condition - PAM in start/stop
    if winnerFound == False:
        conditionA = sgRNACatalogue[sgRNACatalogue["PAMinStartStop"] == True] #This is the subset df for which PAM is in the start/stop
        if len(conditionA) > 0: #if there is one or more sgRNAs for this condition, select the first as the winner
            conditionA = conditionA.reset_index(drop = True) #reset index
            winnersgRNA = conditionA.at[0, "sgRNA_sequence"]
            winnerFound = True

    #B. sgRNA overhang is less than 15bp
    if winnerFound == False:
        conditionB = sgRNACatalogue[sgRNACatalogue["15bpOverhangL"] == False & sgRNACatalogue["15bpOverhangR"] == False]
        if len(conditionB) == 1:
            winnersgRNA = conditionB
            winnerFound = True
        elif len(conditionB) > 1: #more than one, select in CDS preferencially, and closest cut
            winnersgRNA, winnerFound = checkCDSCutandOrder(conditionB)    

    #D. sgRNA overhang is more than 15bp, need to mutate
    if winnerFound == False:
        if len(sgRNACatalogue) == 1:
            winnersgRNA = sgRNACatalogue.at[0, "sgRNA_sequence"]
        else: #more than 1 sgRNA
            winnersgRNA, winnerFound = checkCDSCutandOrder(sgRNACatalogue)            

        mutationNeeded = True


    return winnersgRNA, mutationNeeded, sgRNACatalogue[sgRNACatalogue["sgRNA_sequence"] == winnersgRNA] #could adapt return to be just the true/false for mutation and the catalogue

##incomplete
def stringencyFilesIterator(TFsdict_of_dic):
    """"
    Loop throgh stringency files to select sgRNA for one start/stop position. 
    """
    sgRNAfiles = [""] #sgRNA files in stringency order (low to high)
    stringency = 0
    sgRNAselected = False

    #Loop through the different files and stop if an sgRNA is found
    while sgRNAselected == False:
        #From this, filter guideRNAs nearby
        df = filter_gRNA(sgRNAfiles[stringency], TFsdict_of_dic)

        #Extract necessary boolean and position values
        sgRNACatalogue = positionScore(df)

        #find best guideRNA
        winnersgRNA, mutationNeeded, winnerCatalogue = find_best_gRNA(sgRNACatalogue)

        #Check if a winner was found, go to next stringency level
        if winnersgRNA != "":
            sgRNAselected == True
        else:
            stringency +=1
            if stringency > 3:
                print("no sgRNAs could be found for any stringency file.")
                sgRNAselected = True #To get out of the while loop
    
    return "something"

def codonFragmenter(sequenceString, type = "homologyArm", direction = 'HAL', geneStrand = '+'):
    """
    Fragments codons in order moving away from the start or stop site. If gene on - strand, will take the reverse complement codon.

    params:
        sequence: sequence to fragment into codons as string.
        type: one of "homologyArm" or "primer". If primer, will truncate at far end into a fragment divisible by 3 before fragmenting into codons.
        direction: 'HAL' or 'HAR'. If on left homology arm, codon order is from right to left. If on right, codon order is from left to right.
        geneStrand: if +, codons are the original + strand sequence. If -, each individual codon has been reverse complemented.
    """

    #If sequence is for primer HAL-R, will take the reverse complement so that all fragments/primers used in code below are in the same orientation.
    if type =="primer" & direction == 'HAL':
        sequenceDirected = revComp(sequenceString)
    else:
        sequenceDirected = sequenceString

    #Define the sequence string to use/mutate. If primer, truncate until divisible by 3.
    if type == "primer":
        remainder = len(sequenceDirected) % 3
        if remainder == 1:
            sequence = sequenceDirected[:-1]
        elif remainder == 2:
            sequence = sequenceDirected[:-2]
        else:
            sequence = sequenceDirected #if already divisible by 3
    elif type == "homologyArm":
        sequence = sequenceDirected #no modification needed, always 225bp so divisible by 3
    else:
        print("Sequence type invalid.")

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

##incomplete
def mutator(UpstreamFrag, downstreamFrag, df, sequenceType = "homologyArm"):
    """
    In the case where a fragment or primer needs to be mutated, will mutate in CDS (preferably PAM, if not in the sgRNA). If not possible, will mutate PAM outside of CDS to NGT.

    params:
        sequenceToMutate: the homology arm or primer fragment that will be mutated.
        winner_sgRNACatalogue: sgRNACatalogue in format as above, with only the row for the winner sgRNA selected.
        sequenceType: one of 'homologyArm' or 'primer'.
    """

    positionScore(df)

    #sequenceToMutate, winner_sgRNACatalogue, sequenceType = "homologyArm"

    #Check mutable conditions and posiions
    if positionBoolean.at[]

    #1. Mutate PAM in CDS
    if positionBoolean.at[1,"PAMinCDS"] is True:
        
    #after mutating, revComp back

##Final runner

##incomplete
def runnersgRNAandPrimer(TF_list, type = "homologyArms", sgRNAdf = "none"):

    """
    TF_list: .xlsx file of TFs. Either the full TFs list in the first instance, or a shorter list for which to calculate primers.
    type: One of "HomologyArm" or "primers" - homology arms will be returned in the normal instance.
    sgRNAdf = if none, will run the sgRNA search function per TF in TF_list. If provided, will use the provided sgRNAs.

    """
    
    #make TFsdf
    TFsdf, TFsdict_of_dic = make_dataframe_from_TFs_list(TF_list)

    #sgRNAs
    if sgRNAdf = none
        for dict in TFsdict_ofdic = sgRNA:
            #files iterator
    #Find homology arms or primers

    #Find homology arms
    if type == "homologyArms":
        tfsDf = make_homology_arm_fragments(tfsDF)

    #Find primers
    elif type == "primers":
        TFsdf = primersRunner(TFs = 'inputFiles/truncatedTFsinput.xlsx') #replace this with the TFs list of only those where you need primers

    else:
        print("invalid type. Should be one of 'homologyArms' or 'primers')
              
    #Mutate as needed
    if mutationNeeded == True:
        mutator(sequenceToMutate, winnerCatalogue, sequenceType = "homologyArm")

    #return sgRNA selected, HAR fragments or primers, and I assume metrics
