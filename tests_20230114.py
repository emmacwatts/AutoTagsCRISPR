from utils_2023_01_142 import revComp
from utils_2023_01_142 import make_dataframe_from_TFs_list

def test_revComp():

    sequence = "ATTGC"
    reverse_complement = "GCAAT"

    if revComp(sequence) == reverse_complement:

        print("The function revComp is working!")

    else: print("The function revComp is broken. Do not run the whole pipeline until this bug is fixed!")

    return

def test_make_dataframe_from_TFs_list():

    TF_name = "TFs_mock.xlsx"
    ref_gene = "dmel-all-chromosome-r6.48.fasta"
    annotation = "dmel-all-r6.48.gtf"
    ref_genome = "dmel-all-chromosome-r6.48.fasta"
    transgenic_genome_chr2 = "dmel6-nos-Cas9_on_2.fasta"
    transgenic_genome_chr3 = "dmel6-nos-Cas9_on_3.fasta"

    dataframe = make_dataframe_from_TFs_list(TF_name, ref_genome, annotation, transgenic_genome_chr2, transgenic_genome_chr3)
    
    print(len(dataframe))
    
    if len(dataframe) == 30: # this TF has 15 isoforms and the dataframe outputs two dictionaries per isoform for each start and stop

        if dataframe[0]["Gene_ID"] == dataframe[1]["Gene_ID"]:

            if dataframe[0]["Transcript_ID"] == dataframe[1]["Transcript_ID"]:

                if dataframe[0]["Gene_Region"] == "start_codon" and dataframe[1]["Gene_Region"] == "stop_codon":

                    if dataframe[0]["Start"] == 18414273 and dataframe[1]["Stop"] == 18545586: 

                        if dataframe[0]["Strand"] == "-":

                            if dataframe[2]["Transcript_ID"] == dataframe[0]["Transcript_ID"]
        ["'Gene_ID': 'FBgn0004652',
                    'Transcript_ID': 'FBtr0083651', 
                    'Chromosome': '3R', 
                    'Gene_Region': 'stop_codon', 
                    'Start': 18426145, 
                    'Stop': 18426147, 
                    'Strand': '-'"]
        print("The function make_dataframe_from_TFs_list is working!")

    else: print("The function make_dataframe_from_TFs_list is broken. Do not run the whole program until this error is fixes!")
    
    return