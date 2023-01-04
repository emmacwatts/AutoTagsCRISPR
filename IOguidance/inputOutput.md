
## For initial data formatting (Emma):
    
   Input: 
     TFs_list: excel file of query sequences with Gene_ID and Transcript_ID
     ref_genome: fasta file for reference genome
     annotation: .gtf file for ref_genome
    
   Output
     TFsdf: dataframe of TF information and sequences
     TFsdict_of_dict: TFsdf as a dictionary of dictionaries per index
  
![image](https://user-images.githubusercontent.com/120821707/210603565-492daac8-42fc-467c-abb8-ff107856d6b4.png)

    The dataframe format can be viewed with display() to show it more clearly:
    
![image](https://user-images.githubusercontent.com/120821707/210604459-6cfd0fa3-d152-4585-984c-5484b05309a1.png)

   
## For primers BLAST (Emma):

dftoFASTA function:

  Input:
    primersdf: a dataframe of primers
    transgenicFile: filename and path for fasta file to be created

  Output: transgenicFile - this is a fasta file of the primers identified by all other rows in the dataframe.

![image](https://user-images.githubusercontent.com/120821707/210605340-f6d6edc9-2176-4031-b0e8-5e53b41c3c66.png)

 shell BLAST section:
  Input:
    queryprimers.fasta - the fasta file of your query primers created with the dftoFASTA function
    transgenicRef.fasta - the full genome for the merged transgenic reference
    
  Output:
    primerblast.txt - a truncated BLAST search result for all primers
 
![image](https://user-images.githubusercontent.com/120821707/210605999-e1614a9c-ec02-4c3f-9592-5476f9168dd2.png)

primersCountBLAST function:

  Input:
    BLASTresults = txt file of BLAST results (primerblast.txt)
    primersdf = dataframe of primers to be appended to

  Output:
    faultyPrimers: the primers from primersdf that match zero times or more than once to the reference genome. Given in dataframe format.
    This will look roughly like the image below, although I've not yet run this with Giulia's output.
    
    ![image](https://user-images.githubusercontent.com/120821707/210606850-fd4b6683-6ac2-493e-b3bf-840445143cc0.png)
