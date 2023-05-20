from utils_20230401 import filter_gRNA
from utils_20230401 import make_dataframe_from_TFs_list
from utils_20230401 import find_best_gRNA
from tests_20230401 import test_revComp
from tests_20230401 import test_make_dataframe_from_TFs_list
from tests_20230401 import test_filter_gRNA
from tests_20230401 import test_check_start_stop_NGG
from tests_20230401 import test_check_over_15
from tests_20230401 import test_select_closest
from tests_20230401 import test_check_cutting_site_inside_CDS
from tests_20230401 import test_find_synonymous_codons
from tests_20230401 import test_mutate_PAM_in_codon
from tests_20230401 import test_make_synonymous_mutation
from tests_20230401 import test_translate_nucleotide_position_into_codon_position
from tests_20230401 import test_mutate_PAM_in_HDR_primer
import time
import pandas as pd
import os

# input file paths - adjust this accordingly

sgRNA_files = ["inputfiles/NoOffTarget_high_stringency.gff", "inputfiles/NoOffTarget_med_stringency.gff", "inputfiles/NoOffTarget_low_stringency.gff", \
    "inputfiles/1to3NonCdsOffTarget_low_stringency.gff", "inputfiles/ManyOffTarget_low_stringency.gff"]
TF_names_list = "inputfiles/TFs_mock.xlsx"
TF_genes_list = "inputfiles/TF_log_mock.csv"
annotation = "inputfiles/dmel-all-r6.48.gtf"
ref_genome = "inputfiles/dmel-all-chromosome-r6.48.fasta"
transgenic_genome_chr2 = "inputfiles/dmel6-nos-Cas9_on_2.fasta"
transgenic_genome_chr3 = "inputfiles/dmel6-nos-Cas9_on_3.fasta"
df_sgRNA = pd.DataFrame(columns=['gene_ID','transcript_ID', 'chromosome', 'strand_type', 'start/stop', 'genome_start_codon_pos', 'genome_stop_codon_pos', 'must_PAM_be_mutated_in_HDR_plasmid?', 'cut_outside_of_CDS', 'sgRNA_type', 'sgRNA_seq', 'sgRNA_coordinates'])
timestr = time.strftime("%Y%m%d-%H%M%S")
start = time.time()


# run initial testing before running the main program

test_revComp()
#test_make_dataframe_from_TFs_list() # just commented out because it takes long to run
test_filter_gRNA()
test_check_start_stop_NGG()
test_check_over_15()
test_select_closest()
test_check_cutting_site_inside_CDS()
test_find_synonymous_codons()
test_mutate_PAM_in_codon()
test_make_synonymous_mutation()
test_translate_nucleotide_position_into_codon_position()
test_mutate_PAM_in_HDR_primer()

#@TODO: write a function for this - does not look pretty

TFsdf, TF_dict_of_dict = make_dataframe_from_TFs_list(TF_names_list, ref_genome, annotation, transgenic_genome_chr2, transgenic_genome_chr3)

for TF_dict in TF_dict_of_dict:

    print("TF_dict:",TF_dict)

    # dictionary that translates the dictionary from filter_gRNA to find_best_gRNA function
    gRNA_dict = {}
    gRNA_dict["gene_ID"] = TF_dict_of_dict[TF_dict]["Gene_ID"]
    gRNA_dict["transcript_ID"] = TF_dict_of_dict[TF_dict]["Transcript_ID"]
    gRNA_dict["chromosome"] = TF_dict_of_dict[TF_dict]["Chromosome"]
    gRNA_dict["start/stop"] = TF_dict_of_dict[TF_dict]["Gene_Region"]
    gRNA_dict["strand_type"] = TF_dict_of_dict[TF_dict]["Strand"]

    if gRNA_dict["start/stop"] == "start_codon":

        gRNA_dict["genome_start_codon_pos"] = TF_dict_of_dict[TF_dict]["Start"]
        gRNA_dict["genome_stop_codon_pos"] = "N/A"

    else:
        gRNA_dict["genome_start_codon_pos"] = "N/A"
        gRNA_dict["genome_stop_codon_pos"] = TF_dict_of_dict[TF_dict]["Start"]
    
    single_TF_dict= TF_dict_of_dict[TF_dict]

    for sgRNA_file in sgRNA_files:

        gRNA_hits = filter_gRNA(sgRNA_file,single_TF_dict)

        gRNA_dict["sgRNA_list_positions"] = gRNA_hits["sgRNA_list_positions"]
        gRNA_dict["sgRNA_list_values"] = gRNA_hits["sgRNA_list_values"]
        
        # In case the filter_gRNA() function finds more than one suitable sgRNAs then the funtion find_best_gRNA()
        # is executed to compare which of the sgRNAs is the most suitable one

        if len(gRNA_dict["sgRNA_list_values"]) > 1:

            winner_gRNA = find_best_gRNA(gRNA_dict)
            print(winner_gRNA)
            winner_gRNA["sgRNA_type"] = sgRNA_file

            # The sgRNA_list_positions and sgRNA_list_values have to be renamed because the otherwise the dataframe cannot be 
            # printed as a table in Excel since the arrays have different length

            winner_gRNA["sgRNA_seq"] = winner_gRNA["sgRNA_list_values"]
            winner_gRNA["sgRNA_start"] = winner_gRNA["sgRNA_list_positions"][0]
            winner_gRNA["sgRNA_end"] = winner_gRNA["sgRNA_list_positions"][1]
            del winner_gRNA["sgRNA_list_values"]
            del winner_gRNA["sgRNA_list_positions"]
            df_sgRNA = df_sgRNA.append(winner_gRNA, ignore_index = True)
            print("we have found a winner for this transcript of this transcription factor. No more cycling")
        
            break

        elif len(gRNA_dict["sgRNA_list_values"]) == 1:

            winner_gRNA = gRNA_dict
            winner_gRNA["sgRNA_type"] = sgRNA_file
            winner_gRNA["sgRNA_seq"] = winner_gRNA["sgRNA_list_values"][0]
            winner_gRNA["sgRNA_start"] = winner_gRNA["sgRNA_list_positions"][0][0]
            winner_gRNA["sgRNA_end"] = winner_gRNA["sgRNA_list_positions"][0][1]
            del winner_gRNA["sgRNA_list_values"]
            del winner_gRNA["sgRNA_list_positions"]
            print(winner_gRNA)
            df_sgRNA = df_sgRNA.append(winner_gRNA, ignore_index = True)
            print("we have found a winner for this transcript of this transcription factor. No more cycling")
            

            break

        elif len(gRNA_dict["sgRNA_list_values"]) < 1: 
            
            print("not a winner")
            print("cycling through next sgRNA file")
    
    print("done cycling through sgRNA files")

df_sgRNA.to_excel(f'output/optimal_sgRNAs_{timestr}.xlsx', index=False, header=True)
end = time.time()
print("", end - start)
print("Program finished! The program took", time.strftime("%H:%M:%S", time.gmtime(end - start)) , "to run!")
