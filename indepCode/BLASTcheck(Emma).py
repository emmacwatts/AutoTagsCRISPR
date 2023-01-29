#Primers output will be a dataframe, for BLAST we need FASTA:
def dftoFASTA(finalTFsdfwithPrimers, transgenicFileon3, transgenicFileonOther):
  """
  Reformats a primers dataframe to the appropriate format for a BLAST search. 

  params: 
    primersdf: a dataframe of primers with sequence in the "primer_sequence" column and chromosome number in the "Chromosome" column.
    transgenicFileon3: path and desired filename for fasta file to be created for primers on chromosome 3.
    transgenicFileonOther: path and desired filename for fasta file to be created for primers on all other chromosomes.

  returns: the two fasta files defined in params, in format:
      >"primer{index from primersdf}"
      primer sequence

  """
  import pandas as pd
  from Bio import SeqIO

  #write FASTA file for primers of chromosome3
  primersFASTA = open(transgenicFileon3, "w")

  #Write this into the fasta file as>header \n sequence
  for index, rowcontents in primersdf.iterrows():
    if rowcontents["Chromosome"] == 3 or '3':
      header = (">" + "primer" + str(index))
      sequence = rowcontents["primer_sequence"]

  primersFASTA.write(str(header) + "\n" + str(sequence) + "\n")

  primersFASTA.close()

  #Write FASTA file for primers on all other chromosomes
  primersFASTA = open(transgenicFileonOther, "w")

  #Per row in the dataframe, concatenate all information as a header, and primer sequence as sequence
  #Write this into the fasta file as>header \n sequence
  for index, rowcontents in primersdf.iterrows():
    if rowcontents["Chromosome"] != 3 or '3':
      header = (">" + "primer" + str(index))
      sequence = rowcontents["primer_sequence"]

  primersFASTA.write(str(header) + "\n" + str(sequence) + "\n")

  primersFASTA.close()

#In utils, run this twice, one for each transgenic reference    
def runBLAST(queryPrimersfasta, transgenicreffasta):
  """
  A subprocess shell command to run BLAST. BLAST must be installed.

  params:
    queryPrimersfasta: fasta file of primers to query
    transgeicreffasta: fasta file of the reference genome

  return:
    BLASToutput: BLAST output file in output format 6 as a .txt file

  """ 

  #still to be completed 
  import subprocess

  subprocess.run["BLAST.sh", f"--input_file = {queryPrimersfasta}", 
  f"--output_file = {transgenicreffasta}"]

def primerCountBLAST(BLASTresults, primersdf):
  '''
  Checks for the number of exact matches of a primer against a reference genome.

  Input:
  BLASTresults = txt file of BLAST results in -outfmt 6
  primersdf = p

  Output:
  faultyPrimers: the primers from primersdf that match zero times or more than once to the reference genome. Given in dataframe format.
  '''  
  import pandas as pd
  from Bio import SeqIO

  #Append the primers df with a new column counting the number of instances of that primer against the reference genome
  primersdf["Match Count"] = "0"

  #Open BLAST results
  primerBLAST = open(BLASTresults, "r")
  primerBLASTlines = primerBLAST.readlines()

  #Each 100% match will be labelled with the primer name, so count the instances of each 
  #primer name in the file and append this to "Match Count" in primersdf
  for line in primerBLASTlines:
    lineattributes = (line.split("\t"))
    primerindex = int(lineattributes[0].replace("primer", ""))
    currentcount = primersdf.at[primerindex, "Match Count"]
    primersdf.at[primerindex, "Match Count"] = int(currentcount) + 1
  
  #Keep only faulty primers
  faultyPrimers = primersdf
  for index, rowcontents in faultyPrimers.iterrows():
    if rowcontents["Match Count"] == 1:
      faultyPrimers.drop([index], axis = 0, inplace = True)

  return faultyPrimers

  ##This is the primer BLAST code for per chromosome check
  def BLASTSixPrimers(sixPrimers, chromosome):
    
    import pandas as pd
    from Bio import SeqIO

    #This corresponds to the order of the list, sixPrimers
    primerTypes = ["VAL-F", "HAL-F", "HAL-R", "HAR-F", "HAR-R", "VAL-R"]

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
            faultyPrimer = primerTypes[sixPrimers.index(primer)]
        else:
            faultyPrimer = "none"

    return faultyPrimer

##Just storing this here in case I need it
        #nextPrimers variable allows removal of the top hit where BLAST has indicated primer needs to be excluded.
      if nextPrimer == 0: 
        potential_primers = potential_primers
      elif nextPrimer <= len(potential_primers):
    #If, after BLAST we need to exclude the first x primers:
        potential_primers = potential_primers[nextPrimer:] #remove the top hit(s)
      else:
        print("After excluding BLAST-searched primers, no remaining primers are compatible.")    
        success = False
        potential_primers = "NA"
    
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
    from gtfparse import read_gtf

    #This is the input file containing the TFs we want to query
    #Imported as a pandas dataframe
    queryTFsdf = pd.read_excel(TF_list)

    #This is the .gtf file with annotations for each gene on the reference genome
    refGenomeAnnotations = read_gtf("gene_annotations.gtf")
    
    #This is the FASTA file of the reference genome sequence
    refSeqPerChromosome = {}
    for seq in SeqIO.parse(open(ref_genome), 'fasta'):
        refSeqPerChromosome[seq.id] = seq.seq

    
    refGenomeAnnotation = refGenomeAnnotation.loc[refGenomeAnnotation["feature"].isin(["start_codon", "stop_codon"])]

    TFsdf = refGenomeAnnotation[["Gene_ID", "Gene_Symbol", "Transcript_ID", "Transcript_Symbol", "Chromosome", "Gene_Region", "Start", "Stop", "Strand"]].loc[refGenomeAnnotation["gene_id"].isin(queryTFsdf["Flybase_ID"])]
    
