from utils_2023_01_142 import filter_gRNA
from utils_2023_01_142 import make_dataframe_from_TFs_list
from utils_2023_01_142 import find_best_gRNA
from tests_2023_01_14 import test_revComp
from tests_2023_01_14 import test_make_dataframe_from_TFs_list
import pandas as pd

gRNA_file = "NoOffTarget_high_stringency.gff"
TF_names_list = "TFs_mock.xlsx"
TF_genes_list = "TF_log_mock.csv"
annotation = "dmel-all-r6.48.gtf"
ref_genome = "dmel-all-chromosome-r6.48.fasta"
transgenic_genome_chr2 = "dmel6-nos-Cas9_on_2.fasta"
transgenic_genome_chr3 = "dmel6-nos-Cas9_on_3.fasta"

# run initial testing before running the main program

test_revComp()
test_make_dataframe_from_TFs_list()

'''
TF_dict_of_dict = make_dataframe_from_TFs_list(TF_names_list, ref_genome, annotation, transgenic_genome_chr2, transgenic_genome_chr3)

for TF_dict in TF_dict_of_dict:

    # dictionary that translates the dictionary from filter_gRNA to find_best_gRNA function
    gRNA_dict = {}
    gRNA_dict["start/stop"] = TF_dict_of_dict[TF_dict]["Gene_Region"]

    if gRNA_dict["start/stop"] == "start_codon":

        gRNA_dict["genome_start_codon_pos"] = TF_dict_of_dict[TF_dict]["Start"]
        gRNA_dict["genome_stop_codon_pos"] = 0

    else:
        gRNA_dict["genome_start_codon_pos"] = 0
        gRNA_dict["genome_stop_codon_pos"] = TF_dict_of_dict[TF_dict]["Stop"]

    gRNA_dict["strand_type"] = TF_dict_of_dict[TF_dict]["Strand"]
    
    single_TF_dict= TF_dict_of_dict[TF_dict]

    gRNA_hits = filter_gRNA(gRNA_file,single_TF_dict)

    print(gRNA_hits)

    gRNA_dict["sgRNA_list_pos"] = gRNA_hits["sgRNA_list_pos"]
    gRNA_dict["sgRNA_list_values"] = gRNA_hits["sgRNA_list_values"]
    
    if gRNA_dict["sgRNA_list_pos"]:

        winner_gRNA = find_best_gRNA(gRNA_dict)

    else: print("not a winner")

print(winner_gRNA)
'''