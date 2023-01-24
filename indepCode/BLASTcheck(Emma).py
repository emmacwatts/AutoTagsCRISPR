#Primers output will be a dataframe, for BLAST we need FASTA:
def dftoFASTA(primersdf, transgenicFile):
  """
  Reformats a primers dataframe to the appropriate format for a BLAST search. 

  params: 
    primersdf: a dataframe of primers
    transgenicFile: filename and path for fasta file to be created

  returns: transgenicFile - this is a fasta file of the primers identified by all other rows in the dataframe.

  """
  import pandas as pd
  from Bio import SeqIO

  primersFASTA = open(transgenicFile, "w")

  #Per row in the dataframe, concatenate all information as a header, and primer sequence as sequence
  #Write this into the fasta file as>header \n sequence
  for index, rowcontents in primersdf.iterrows():
    split = (str(rowcontents).split())
    header = (">" + "primer" + str(index))
    sequence = split[-7]

    primersFASTA.write(str(header) + "\n" + str(sequence) + "\n")

  primersFASTA.close()

#use shell to run BLAST on primers against custom transgenic reference file
%%shell

#I've made a custom transgenic reference FASTA which uses the appropriate transgenic strain for each chromosome
#(has chromosome 3 from nos-Cas9_on_2 and all other chromosomes from nos-Cas9_on_3)
#Use this to BLAST against - "transgenicRef.fasta"

#install BLAST (if necessary)
gunzip 'ncbi-blast-2.13.0+-x64-linux.tar.gz'
tar -xvf 'ncbi-blast-2.13.0+-x64-linux.tar'
#change directory to ncbi-blast-2.13.0+/

#run BLAST
makeblastdb -in transgenicRef.fasta -dbtype nucl -parse_seqids 
blastn -query queryprimers.fasta -subject transgenicRef.fasta -evalue 1e-23 -outfmt 6 -out primerblast.txt -max_target_seqs 2
  #This output format is as simple as possible, so that we get as small a file when querying many primers at once

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