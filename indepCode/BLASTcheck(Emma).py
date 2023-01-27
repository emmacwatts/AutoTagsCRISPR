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