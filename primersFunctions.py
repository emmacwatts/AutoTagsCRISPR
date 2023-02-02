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

def designPrimers(stringency, primer_type, template, useExtendedRegion = False):
 
  '''
  This function calculates a primers given the desired stringency and primer type.

  params:
    stringency: value from 0-2 for appropriate stringency conditions.
    primer_type: one of 'VAL-F', 'HAL-F', 'HAL-R', 'HAR-F', 'HAR-R', 'VAL-R'
    template: gene region of the start/stop codon (1700bp upstream and downstream)
    useExtendedRegion: if True, uses an extended region of the template for each primer.
    blastSkip: standard is 0, meaning no primers to exclude. If <1, will exclude that many potential primers from the beginning of the primer3 result. 
    stringencyRulesdf: primer specification for descending order of stringencies.
    primerNameandRegiondf: dataframe of specifications for each primer type

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

  #Define whether to use the intial search region or extended search region
  if useExtendedRegion == False:
    gene_region = [p["initial_region_start"], p["initial_region_length"]]
  else:
    gene_region = [p["extended_region_start"], p["extended_region_stop"]]

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
  
  #return success = False (no primer found), or success = True, and the primers found
  if all(x == 0 for x in primer3Results.values()) == True: #this is the output format if no primers are found
    success = False
    potential_primers = "NA"
  else: #if primers were found, compile these sequences into a list
    success = True
    potential_primers = []
    for columnName in primer3Results.keys():
      if columnName.endswith("SEQUENCE"):
        potential_primers.append(primer3Results[columnName])

    return success, potential_primers

def stringencyIterate(template, primer_type, useExtRegion = False):
  """
  Iterate through stringency levels to find primers. Primer mounting will be checked for "HAL-R" and "HAL-F".

  params:
    template: sequence of the full gene region as a string
    primer_type: type of primer being designed
    useExtRegion: if True, will use expanded search regions for primers.
  
  output:
    potential_primers: sequence of the primer found. If none are found, returns NA
    
  """
  success = False

    #iterate through stringencies
  for stringency in range(0,2):
    success, potential_primers = designPrimers(stringency, primer_type, template, useExtendedRegion = useExtRegion)
    if success == True and primer_type == "HAL-R" or primer_type =="HAR-F": #check mounting for special case
      potential_primers = mountedPrimers(primer_type, potential_primers, template)
      #If this is the last stringency, use makeMount to generate a mounted primer
      if potential_primers == [] and stringency == 2:
          potential_primers = mountedPrimers(primer_type, potential_primers, template, makeMount = True)
          success = True 
      #Otherwise if not primer is correctly mounted, return false to move to the next stringency
      elif potential_primers == [] and stringency < 2:
          success = False
    #if no primer obtained, re-design primer (stay in loop)
    #else, keep this potential_primers value
    if success == True:
      break

  #If still not returning a primer (for non-mount primers), the sequence will be NA
  if success == False:

    print("Could not find a primer at any stringency levels.")
    potential_primers = "NA"
    stringency = "none"
  
  return stringency, potential_primers

def mountedPrimers(primer_type, potential_primers, template, makeMount = False):
  """
  Perform special conditions check for primers that must be mounted to a site in the template, without a force-mount to maintain primer quality.

  params:
    primer_type: string variale specifying primer type (e.g. "HAL_R")
    potential_primers: list of candidate sequences for this primer
    template: 3400 bp region of the reference genome for this start/stop site
    makeMount = if True, will manually assign 25 bp from the intended primer start site.

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

  if makeMount == False:

    for primer in potential_primers:
      if primer_type == "HAL-R": #This is a reverse primer, so need reverse complement to match to template
        primer = revComp(primer)
      startIndex = template.find(primer)
      if startIndex == p["initial_region_start"]:
        mountedPrimer.append(primer)
  
  #If no primer is found for all stringency levels then manually design.
  else:
    mountedPrimer = template[p["initial_region_start"]: p["initial_region_start"]+p["initial_region_length"]]
    if primer_type == "HAL-R": #This is a reverse primer, so need reverse complement to match to template
        mountedPrimer = [revComp(mountedPrimer)]
  
  return mountedPrimer

def addSixPrimers(template, useExtendedRegion = False, blastSkip = {}):
  """
  Find sequences for each of the six primers needed per transcript per start/stop region.

  params:
    template: sequence of the full gene region as a string
    useExtendedRegion: if true, the primers will be designed using a larger region around the start/stop codon.
    blastSkip: dictionary of primer types: number of BLAST-failed primers to skip for that type.

  output:
    stringencyList: the stringency used per primer type.
    sixPrimers: the six primer sequences as a list of strings.

  """
  primersNeeded = ["VAL-F", "HAL-F", "HAL-R", "HAR-F", "HAR-R", "VAL-R"]

  #setting up lists for desired output
  sixPrimers = []
  stringencyList = []

  #Obtain each of the 6 primers one by one
  for primer_type in primersNeeded:
    #Note here that potential_primers could be NA, many primers, or a single primer in a list
    stringency, potential_primers = stringencyIterate(template, primer_type, useExtRegion = useExtendedRegion)

    #Record the stringency used per primer type
    stringencyList.append(stringency)

    #If BLAST has indicated that we need to skip primers, blastSkip will increase in value for that primer type and this will call the next primer along in potential_primers.
    if blastSkip[primer_type] > len(primer_type):
      sixPrimers.append("NA")
      print("All potential primers have failed the BLAST check.")
    else:
      sixPrimers.append(potential_primers[blastSkip[primer_type]])

  return stringencyList, sixPrimers

def checkSixPrimers(template, chromosome):
  """
  Quality control and redesign if needed for the six primers per template.

  params:
    template: string of the 3400 bp region surrounding the start/stop site.
    chromosome: chromosome number of the template strand.
  
  output:
    blastSkip: a dictionary of primer type : how many primers failed BLAST check for that primer.
    stringencyList: list of stringencies for each primer in the set (from VAL-F -> VAL-R)
    extended: whether or not the extended region was used to design these primers.
    sixPrimers: the final set of six primers designed (from VAL-F -> VAL-R)

  """
  #Standard case - no extended region needed, and BLAST hasn't yet indicated any primers that need to be skipped
  extended = False
  blastSkip = {"VAL-F":0, "HAL-F":0, "HAL-R":0, "HAR-F":0, "HAR-R":0, "VAL-R":0}

  #run addSixPrimers
  stringencyList, sixPrimers = addSixPrimers(template)

  #If any of the 6 primers are 'NA', redesign the set of 6 with extended regions.
  if 'NA' in sixPrimers:
    extended = True
    stringencyList, sixPrimers = addSixPrimers(template, useExtendedRegion = extended)

  #BLAST check the resulting primers    
  faultyPrimers = BLASTSixPrimers(sixPrimers, chromosome)
  while faultyPrimers != []: #This indicates that there is faulty primers present
    for primer in faultyPrimers:
      blastSkip[primer] += 1
    stringencyList, sixPrimers = addSixPrimers(template, useExtendedRegion = extended, blastSkip = blastSkip)
    faultyPrimers = BLASTSixPrimers(sixPrimers, chromosome)

  return blastSkip, stringencyList, extended, sixPrimers

def BLASTSixPrimers(sixPrimers, chromosome):
  """
  BLAST check the primer set to the appropriate transgenic reference genome.

  params:
    sixPrimers: the set of 6 primers in a list of strings.
    chromosome: the chromosome that the template for the primers is on.

  output: 
    faultyPrimers: list of the primer types of any primers that match more or less than once to the transgenic ref genome.
  """

  import pandas as pd
  from Bio import SeqIO

  #This corresponds to the order of the list, sixPrimers
  primerTypes = ["VAL-F", "HAL-F", "HAL-R", "HAR-F", "HAR-R", "VAL-R"]
  faultyPrimers = []

  #Upload transgenic reference file
  #This will be chromosome 3 of cas_on_2, and all other chromosomes from cas_on_3.
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
          faultyPrimers.append(primerTypes[sixPrimers.index(primer)])
      else:
          continue

  return faultyPrimers

def finalPrimersdf(TFsdf, returnParameters = False):
  """
  Compile all primers for transcription factors.

  params:
    TFsdf: full dataframe contaning transcript ID and gene region type (start/stop)
    returnParameters: if False, returns only primer IDs and primers. if true, returns parameters of the designed primers (strigency, whether the extended region was used, how many primers were skipped due to BLAST removal)

  output:
    TFsdf: a dataframe of primer information.

  """
  #Start blank columns for the primers we need to add per transcript per start/stop (so for every row in TFsdf)
  import pandas as pd

  #This will contain all primers identified by primer type (6 columns)
  allPrimers = pd.DataFrame(columns = ["VAL-F", "HAL-F", "HAL-R", "HAR-F", "HAR-R", "VAL-R"])
  if returnParameters == True: #Additional metrics will be provided in the dataframe if the user requests this.
    allPrimers[["Stringency Level", "Extended Region Used", "Number of Primers BLAST Excluded"]]
  else: #Otherwise, a truncated TFs information list is returned (Transcript_ID, Gene_Region)
    TFsdf = TFsdf.drop(["Gene_ID", "Gene_Symbol", "Transcript_ID", "Chromosome", "Start", "Stop", "Strand", "Reference_Seq"])

  for index, row in TFsdf.iterrows():
    blastSkip, stringencyList, extended, sixPrimers = checkSixPrimers(row["Reference_Seq"], row.at[index, "Chromosome"])
    if returnParameters == False:    
      allPrimers.loc[len(allPrimers)] = sixPrimers
    else:
      fullResults = sixPrimers + stringencyList + extended + blastSkip
      allPrimers.loc[len(allPrimers)] = fullResults

  #Add primer sequence df to original (or truncated) TFsdf
  TFsdf = TFsdf.join(allPrimers)

  return TFsdf

