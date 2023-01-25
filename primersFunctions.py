def designPrimers(stringency, primer_type, template, extend = False, nextPrimer = 0, printFullResults = False):
 
  '''
  This function calculates a primers given the desired stringency and primer type.

  params:
    stringency: value from 0-2 for appropriate stringency conditions.
    primer_type: one of 'VAL-F', 'HAL-F', 'HAL-R', 'HAR-F', 'HAR-R', 'VAL-R'
    template: gene region of the start/stop codon (1700bp upstream and downstream)
    extend: if True, uses an extended region of the template for each primer.
    nextPrimer: standard is 0, meaning no primers to exclude. If <1, will exclude that many potential primers from the beginning of the primer3 result. 
    stringencyRulesdf: primer specification for descending order of stringencies.
    primerNameandRegiondf: dataframe of specifications for each primer type

  output:
    primer_sequence = list of string(s) of the primer sequences.
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
  if extend == False:
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
  
  if printFullResults == True:
    return primer3Results
  
  else: #return success = False (no primer found), or success = True, and the primer(s) found
    if all(x == 0 for x in primer3Results.values()) == True: #this is the output format is no primers are found
      success = False
      primer_sequence = "NA"
    else: #if primers were found, compile the sequences of these into a list
      success = True
      primer_sequence = []
    # fullprimer3results = pd.DataFrame.from_dict(primer) #can modify the code to also return full dataframe information is desired
      for columnName in primer3Results.keys():
        if columnName.endswith("SEQUENCE"):
          primer_sequence.append(primer3Results[columnName])

      #nextPrimers variable allows removal of the top hit where this primer needs to be excluded.
      if nextPrimer == 0: 
        primer_sequence = primer_sequence
      elif nextPrimer <= len(primer_sequence):
    #If, after BLAST we need to exclude the first x primers:
        primer_sequence = primer_sequence[nextPrimer:] #remove the top hit(s)
      else:
        print("After excluding BLAST-searched primers, no remaining primers are compatible.")    
        success = False
        primer_sequence = "NA"
    
    return success, primer_sequence

def stringencyAndExtension(template, primer_type, nextprimer = 0):
  """
  Iterate through stringency levels and extended sequence regions to find primers.

  params:
    template: sequence of the full gene region as a string
    primer_type: type of primer being designed
  
  output:
    primer_sequence: sequence of the primer found. If none are found, returns NA
    
  """
    #start at stringency 0
  for stringency in range(0,2):
    success, primer_sequence = designPrimers(stringency, primer_type, template, nextPrimer = nextprimer)
    #if no primer obtained, re-design primer (stay in loop)
    #else, keep this primer_sequence value
    if success == True:
      break
  
  #After iterating through stringencies, if success is still False, check extended regions
  if success == False:
    #Again, interate through stringencies
    for stringency in range(0,2):
      success, primer_sequence = designPrimers(stringency, primer_type, template, extend = True, nextPrimer = nextprimer)
      if success == True:
        break

  #If still not returning a primer, the sequence will be NA
  if success == False:
    print("Could not find a primer within the defined region at these stringency levels.")
    primer_sequence = "NA"
  
  return primer_sequence

def addSixPrimers(template):
  """
  Find sequences for each of the six primers needed per transcript per start/stop region.

  params:
    template: sequence of the full gene region as a string

  output:
    sixPrimers: the six primer sequences as a list of strings.
  """
  primersNeeded = ["VAL-F", "HAL-F", "HAL-R", "HAR-F", "HAR-R", "VAL-R"]
  #sixPrimers is where the primer sequences will go, and will correspond in order with 'primersNeeded' above
  sixPrimers = []

  #Obtain each of the 6 primers one by one
  for primer_type in primersNeeded:
    #Note here that primer_sequuence could be NA, many primers, or a single primer in a list
    primer_sequence = stringencyAndExtension(template, primer_type)

    if primer_type == "HAL-R" or primer_type =="HAR-F":
      #In these two cases, primer_sequence should return about 1000 primers.
      primer_sequence = mountedPrimers(primer_type, primer_sequence, template)
    
    if primer_sequence == []:
      primer_sequence = ['NA']
      
    #If there's many primers (and we've already filtered for special condition checks, just take the best/first one)
    sixPrimers.append(primer_sequence[0])    

  return sixPrimers

def mountedPrimers(primer_type, primer_sequence, template):
  """
  Perform special conditions check for primers that must be mounted to a site in the gene region, without a force-mount to maintain primer quality.

  params:
    primer_type: string variale specifying primer type of the six types (e.g. "HAL_R")
    primer_sequence: list of candidate sequences for this primer
    template: region of the reference genome for this start/stop site

  output:
    mountedPrimer: a single primer that is mounted at the start/stop-adjacent site.

  """
  import pandas as pd

  #Import relevant files and extract  values for the primer_type from the dataframe:
  #Primer specifications
  primerNameandRegiondf = pd.read_excel("inputfiles/Primer_Name_and_Regions.xlsx")  
  #Reset index to primer_type
  primerNameandRegiondf = primerNameandRegiondf.set_index("primer_type", drop = True)
 
  primerNameandRegiondf = primerNameandRegiondf.loc[primer_type]
  p = primerNameandRegiondf.to_dict()

  foundMountedPrimer = False

  for primer in primer_sequence:
    startIndex = template.find(primer)
    if startIndex == p["initial_region_start"]:
      foundMountedPrimer = True
      mountedPrimer = primer
      return [mountedPrimer]

  #If no primer is found in the up to 1000 primers requested from primer3, then manually design.
  if foundMountedPrimer == False:
    mountedPrimer = template[p["initial_region_start"]: p["initial_region_start"]+p["initial_region_length"]]
    return [mountedPrimer]


def finalPrimersdf(TFsdf):
  """
  Compile all primers for transcription factors.

  params:
    TFsdf: full dataframe contaning transcript ID and gene region type (start/stop)

  output:
    allPrimers: a dataframe of 6 primers per site, with indexes matching TFsdf.
    TFsdf: a dataframe of allPrimers joined to the input TFsdf.

  """
  #Start blank columns for the primers we need to add per transcript per start/stop (so for every row in TFsdf)
  import pandas as pd

  #This will contain all primers identified by primer type (6 columns)
  allPrimers = pd.DataFrame(columns = ["VAL-F", "HAL-F", "HAL-R", "HAR-F", "HAR-R", "VAL-R"])

  for index, row in TFsdf.iterrows():
    sixPrimersList = addSixPrimers(row["Reference_Seq"])
    #Add this list to the end of the full primers dataframe
    allPrimers.loc[len(allPrimers)] = sixPrimersList

  #Add primer sequences to the TFsdf full dataframe, which contains transcript and gene region IDs
  TFsdf = TFsdf.join(allPrimers)

  #This output gives the option to see all primer sequences only, or the primer sequences added to TFsdf
  return allPrimers, TFsdf

def primersBLASTrecalc(problemPrimerIndex, problemPrimerType, TFsdf, allPrimers):
  """
  Recalculates a new primer, excluding the top hit(s) where BLAST results indicate that they should be ruled out.
  params:
    problemPrimerIndex: index number corresponding to TFsdf and allPrimers for the erroneous primer.
    problemPrimerType: primer type (in format e.g. "VAL-F") for the erroneous primer.
    TFsdf: dataframe of full information about each TF, including designated primers.
    allPrimers: dataframe of all primers only, with index corresonding to TFsdf.
  output:
  """
  import pandas as pd

  #Extract the template again for the problem primer
  template = TFsdf.at[problemPrimerIndex, "Reference_Seq"]
  nextprimer = 1
  primer_sequence = stringencyAndExtension(template, primer_type = problemPrimerType, nextprimer = nextprimer)
  #Fix the primer sequence
  TFsdf.at[problemPrimerIndex, "problemPrimerType"] = primer_sequence[0]
  allPrimers.at[problemPrimerIndex, "problemPrimerType"] = primer_sequence[0]

  print(f"Problem primer {problemPrimerType} at index {problemPrimerIndex} has been fixed in all dataframes.")
  
  return TFsdf, allPrimers