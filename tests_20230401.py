from utils_20230401 import revComp
from utils_20230401 import make_dataframe_from_TFs_list
from utils_20230401 import filter_gRNA
from utils_20230401 import check_start_stop_NGG
from utils_20230401 import check_over_15
from utils_20230401 import select_closest
from utils_20230401 import check_cutting_site_inside_CDS
from utils_20230401 import find_synonymous_codons
from utils_20230401 import mutate_PAM_in_codon
from utils_20230401 import make_synonymous_mutation
from utils_20230401 import translate_nucleotide_position_into_codon_position
from utils_20230401 import mutate_PAM_in_HDR_plasmid
from os import sys

def test_revComp():

    sequence = "ATTGC"
    reverse_complement = "GCAAT"

    if revComp(sequence) == reverse_complement:

            print("The function revComp() is working!")

    else: 
        
        print("The function revComp is broken. Do not run the whole pipeline until this bug is fixed!")

        sys.exit()

    return

def test_make_dataframe_from_TFs_list():

    TF_name = "inputfiles/TFs_mock.xlsx"
    ref_gene = "inputfiles/dmel-all-chromosome-r6.48.fasta"
    annotation = "inputfiles/dmel-all-r6.48.gtf"
    ref_genome = "inputfiles/dmel-all-chromosome-r6.48.fasta"
    transgenic_genome_chr2 = "inputfiles/dmel6-nos-Cas9_on_2.fasta"
    transgenic_genome_chr3 = "inputfiles/dmel6-nos-Cas9_on_3.fasta"

    TFsdf, dataframe = make_dataframe_from_TFs_list(TF_name, ref_genome, annotation, transgenic_genome_chr2, transgenic_genome_chr3)
    
    if not len(dataframe) == 30: # this TF has 15 isoforms and the dataframe outputs two dictionaries per isoform for each start and stop

        print("Length of the dataframe generated by transcription factor isoforms is incorrect. Do not run the whole pipeline until this bug in function make_dataframe_from_TFs_list() is fixed!")

        sys.exit()
    
    elif not dataframe[list(dataframe.keys())[0]]["Gene_ID"] == dataframe[list(dataframe.keys())[1]]["Gene_ID"]:

        print("Gene_ID is not matching. Do not run the whole pipeline until this bug in function make_dataframe_from_TFs_list() is fixed!")

        sys.exit()

    elif not dataframe[list(dataframe.keys())[0]]["Transcript_ID"] == dataframe[list(dataframe.keys())[1]]["Transcript_ID"]:

        print("Transcript_IDs for start and stop are not matching. Do not run the whole pipeline until this bug in function make_dataframe_from_TFs_list() is fixed!")

        sys.exit()

    elif not dataframe[list(dataframe.keys())[0]]["Gene_Region"] == "start_codon" and not dataframe[list(dataframe.keys())[1]]["Gene_Region"] == "stop_codon":

        print("Transcript does not contain start and stop codon. Do not run the whole pipeline until this bug in function make_dataframe_from_TFs_list() is fixed!")

        sys.exit()

    elif not dataframe[list(dataframe.keys())[0]]["Start"] <= 18545586 and not dataframe[list(dataframe.keys())[1]]["Stop"] >= 18414273: 
        
        print("Start and Stop of Transcript have wrong genomic coordinates. Do not run the whole pipeline until this bug in function make_dataframe_from_TFs_list() is fixed!")

        sys.exit()

    elif not dataframe[list(dataframe.keys())[0]]["Strand"] == "-":

        print("Strand is incorrect. Do not run the whole pipeline until this bug in function make_dataframe_from_TFs_list() is fixed!")

        sys.exit()

    elif not int(dataframe[list(dataframe.keys())[2]]["Transcript_ID"][-7:]) == int(dataframe[list(dataframe.keys())[0]]["Transcript_ID"][-7:]) + 1:

        print("Transcript_ID´s are not consecutive. Do not run the whole pipeline until this bug in function make_dataframe_from_TFs_list() is fixed!")

        sys.exit()

    elif not dataframe[list(dataframe.keys())[0]]["Chromosome"] == "3R":

        print("Incorrect chromosome. Do not run the whole pipeline until this bug in function make_dataframe_from_TFs_list() is fixed!")

        sys.exit()

    elif not len(dataframe[list(dataframe.keys())[0]]["Reference_Seq"]) == 3203:

        print("The Referense_Seq is only", len(dataframe[list(dataframe.keys())[0]["Reference_Seq"]]), "long. Do not run the whole pipeline until this bug in function make_dataframe_from_TFs_list() is fixed!")

        sys.exit()

    elif not dataframe[list(dataframe.keys())[0]]["Reference_Seq"][1600 : 1600 + 3] == "ATG":

        print("The position of the start codon does not add up. Do not run the whole pipeline until this bug in function make_dataframe_from_TFs_list() is fixed!")

        sys.exit()

    else: 
        
        print("The function make_dataframe_from_TFs_list() is working!")
    
    return

def test_filter_gRNA():

    TF_dict1 = {
        'Gene_ID': 'FBgn0004652',
        'Transcript_ID': 'FBtr0083651', 
        'Chromosome': 'Y', 
        'Gene_Region': 'stop_codon', 
        'Start': 29027+19, 
        'Stop': 29027+19+2, 
        'Strand': '+', 
        'Reference_Seq': 1600 * "C" + "TAG" + "G" + 1599 * "C", 
    }

    sgRNA_file = "inputfiles/sgRNA_mock.gff"

    solution1 = {
        'Gene_ID': 'FBgn0004652',
        'Transcript_ID': 'FBtr0083651', 
        'Chromosome': 'Y', 
        'Gene_Region': 'stop_codon', 
        'Start': 29027+19, 
        'Stop': 29027+19+2, 
        'Strand': '+', 
        'Reference_Seq': 1600 * "C" + "TAG" + "G" + 1599 * "C", 
        'sgRNA_list_positions': [[29027, 29049]],
        'sgRNA_list_values': [19*'C'+'TAGG']
    }

    solution2 = {
        'Gene_ID': 'FBgn0004652',
        'Transcript_ID': 'FBtr0083651', 
        'Chromosome': '3R', 
        'Gene_Region': 'stop_codon', 
        'Start': 29027+19, 
        'Stop': 29027+19+2, 
        'Strand': '+', 
        'Reference_Seq': 1600 * "C" + "TAG" + "G" + 1599 * "C", 
        'sgRNA_list_positions': [[29027, 29049]],
        'sgRNA_list_values': [19*'C'+'TAGG']
    }

    outcome1 = filter_gRNA(sgRNA_file, TF_dict1)

    if  outcome1 == solution1:

        if outcome1 != solution2:

            print("The function filter_sgRNA() is working!")
    
        else:

            print("Mismatching chromosomes are not discovered. Do not run the whole pipeline until this bug in function filter_sgRNA() is fixed!")

            sys.exit()

    else:
        
        print("Do not run the whole pipeline until this bug in function filter_sgRNA() is fixed!")

        sys.exit()

    return

def test_check_start_stop_NGG():
    
    TF_dict1 = {
        'Gene_ID': 'FBgn0004652',
        'Transcript_ID': 'FBtr0083651', 
        'Chromosome': 'Y', 
        'start/stop': 'stop_codon', 
        'genome_stop_codon_pos': 29027+19, 
        'genome_start_codon_pos': 0, 
        'strand_type': '+', 
        'Reference_Seq': 1600 * "C" + "TAG" + "G" + 1599 * "C", 
        'sgRNA_list_positions': [[29027,29049]],
        'sgRNA_list_values': [19*'C'+'TAGG']
    }

    solution1 = {
        'Gene_ID': 'FBgn0004652',
        'Transcript_ID': 'FBtr0083651', 
        'Chromosome': 'Y', 
        'start/stop': 'stop_codon', 
        'genome_stop_codon_pos': 29027+19, 
        'genome_start_codon_pos': 0, 
        'strand_type': '+', 
        'Reference_Seq': 1600 * "C" + "TAG" + "G" + 1599 * "C", 
        'sgRNA_list_positions': [29027,29049],
        'sgRNA_list_values': 19*'C'+'TAGG',
        'must_PAM_be_mutated_in_HDR_plasmid?' : 'no'
    }

    TF_dict2 = {
        'Gene_ID': 'FBgn0004652',
        'Transcript_ID': 'FBtr0083651', 
        'Chromosome': 'Y', 
        'start/stop': 'stop_codon', 
        'genome_stop_codon_pos': 29027+19, 
        'genome_start_codon_pos': 0, 
        'strand_type': '+', 
        'Reference_Seq': 1600 * "C" + "TAG" + "G" + 1599 * "C", 
        'sgRNA_list_positions': [[29027,29049]],
        'sgRNA_list_values': [18*'C'+'TAGG'+'C']
    }
    
    TF_dict3 = {
        'Gene_ID': 'FBgn0004652',
        'Transcript_ID': 'FBtr0083651', 
        'Chromosome': 'Y', 
        'start/stop': 'start_codon', 
        'genome_stop_codon_pos': 0, 
        'genome_start_codon_pos': 29027+19, 
        'strand_type': '+', 
        'Reference_Seq': 1600 * "C" + "ATG" + "G" + 1599 * "C", 
        'sgRNA_list_positions': [[29027,29049]],
        'sgRNA_list_values': [19*'C'+'ATGG']
    }

    solution3 = {
        'Gene_ID': 'FBgn0004652',
        'Transcript_ID': 'FBtr0083651', 
        'Chromosome': 'Y', 
        'start/stop': 'start_codon', 
        'genome_stop_codon_pos': 0, 
        'genome_start_codon_pos': 29027+19, 
        'strand_type': '+', 
        'Reference_Seq': 1600 * "C" + "ATG" + "G" + 1599 * "C", 
        'sgRNA_list_positions': [29027,29049],
        'sgRNA_list_values': 19*'C'+'ATGG',
        'must_PAM_be_mutated_in_HDR_plasmid?' : 'no'
    }

    outcome1 = check_start_stop_NGG(TF_dict1)
    outcome2 = check_start_stop_NGG(TF_dict2)
    outcome3 = check_start_stop_NGG(TF_dict3)

    if outcome1 == solution1 and outcome3 == solution3:

        if not outcome2:

            print("The function test_check_start_stop_NGG() is working!")
        
        else:   
            
            print("False positive recognition of best case scenario. Do not run the whole pipeline until this bug in function check_start_stop_NGG() is fixed!")

            sys.exit()
    else: 

        print ("False negative recognition of best case scenario. Do not run the wole pipeline until this bug in check_start_stop_NGG() is fixed!")

        sys.exit()

def test_check_over_15():

    # the PAM in the HDR plasmid for this sgRNA does not need to be mutated

    TF_dict1 = {
        'Gene_ID': 'FBgn0004652',
        'Transcript_ID': 'FBtr0083651', 
        'Chromosome': 'Y', 
        'start/stop': 'stop_codon', 
        'genome_stop_codon_pos': 29049-17, 
        'genome_start_codon_pos': 0, 
        'strand_type': '+', 
        'Reference_Seq': 1600 * "C" + "TAG" + 12 * "C"+ "TGG" + 1585 * "C", 
        'sgRNA_list_positions': [[29027,29049]],
        'sgRNA_list_values': [5 * "C" + "TAG" + 12 * "C"+ "TGG"]
    }

    # the PAM in the HDR plasmid for this sgRNA does (!) need to be mutated

    TF_dict2 = {
        'Gene_ID': 'FBgn0004652',
        'Transcript_ID': 'FBtr0083651', 
        'Chromosome': 'Y', 
        'start/stop': 'stop_codon', 
        'genome_stop_codon_pos': 29049-18, 
        'genome_start_codon_pos': 0, 
        'strand_type': '+', 
        'Reference_Seq': 1600 * "C" + "TAG" + 13 * "C"+ "TGG" + 1584 * "C", 
        'sgRNA_list_positions': [[29027,29049]],
        'sgRNA_list_values': [4 * "C" + "TAG" + 13 * "C"+ "TGG"]
    }
    
    # the PAM in the HDR plasmid for this sgRNA does not need to be mutated

    TF_dict3 = {
        'Gene_ID': 'FBgn0004652',
        'Transcript_ID': 'FBtr0083651', 
        'Chromosome': 'Y', 
        'start/stop': 'start_codon', 
        'genome_stop_codon_pos': 0, 
        'genome_start_codon_pos': 29049-17, 
        'strand_type': '+', 
        'Reference_Seq': 1600 * "C" + "ATG" + 12 * "C"+ "TGG" + 1585 * "C", 
        'sgRNA_list_positions': [[29027,29049]],
        'sgRNA_list_values': [5 * "C" + "ATG" + 12 * "C"+ "TGG"]
    }
    
    solution1 = check_over_15(TF_dict1)
    solution2 = check_over_15(TF_dict2)
    solution3 = check_over_15(TF_dict3)

    if solution1["must_PAM_be_mutated_in_HDR_plasmid?"] == "no":

        if solution3["must_PAM_be_mutated_in_HDR_plasmid?"] == "no":

            if not solution2["must_PAM_be_mutated_in_HDR_plasmid?"] == "no":

                print("The function test_check_over_15() is working!")

            else: 

                print("sgRNA that is too far away from stop/start was wrongfully classified. Do not run the wole pipeline until this bug in check_over_15() is fixed!")

                sys.exit()

        else: 

            print("sgRNA that is close enough to start codon was wrongfully classified. Do not run the wole pipeline until this bug in check_over_15() is fixed!")

            sys.exit()

    else:

        print("sgRNA that is close enough to stop codon was wrongfully classified. Do not run the wole pipeline until this bug in check_over_15() is fixed!")

        sys.exit()

    return

def test_select_closest():

    TF_dict1 = {
        'Gene_ID': 'FBgn0004652',
        'Transcript_ID': 'FBtr0083651', 
        'Chromosome': 'Y', 
        'start/stop': 'stop_codon', 
        'genome_stop_codon_pos': 29049-17, 
        'genome_start_codon_pos': 0, 
        'strand_type': '+', 
        'Reference_Seq': 1600 * "C" + "TAG" + 12 * "C"+ "TGG" + 1585 * "C", 
        'sgRNA_list_positions': [[29026,29049], [29027, 29050]],
        'sgRNA_list_values': [6 * "C" + "TAG" + 11 * "C"+ "TGG", 5 * "C" + "TAG" + 12 * "C"+ "TGG"]
    }

    TF_dict2 = {
        'Gene_ID': 'FBgn0004652',
        'Transcript_ID': 'FBtr0083651', 
        'Chromosome': 'Y', 
        'start/stop': 'stop_codon', 
        'genome_stop_codon_pos': 29049-17, 
        'genome_start_codon_pos': 0, 
        'strand_type': '+', 
        'Reference_Seq': 1600 * "C" + "TAG" + 12 * "C"+ "TGG" + 1585 * "C", 
        'sgRNA_list_positions': [[29004, 29027], [29002, 29025]],
        'sgRNA_list_values': [23 * "C", 6 * "C" + "TAG" + 11 * "C"+ "TGG"]
    }

    solution1 = select_closest(TF_dict1)
    solution2 = select_closest(TF_dict2)

    if solution1["sgRNA_list_values"] == 6 * "C" + "TAG" + 11 * "C"+ "TGG":

        if solution2["sgRNA_list_values"] == 23 * "C":

            print("The function select_closest() is working!")

        else:

            print("The wrong sgRNA was selected from two sgRNAs with different distances from start/stop codon. \
            Do not run the wole pipeline until this bug in select_closest() is fixed!")
            
            sys.exit()

    else: 

        print("The wrong sgRNA was selected from two sgRNAs with different distances from start/stop codon. \
            Do not run the wole pipeline until this bug in select_closest() is fixed!")

        sys.exit()
        
    return

def test_check_cutting_site_inside_CDS():

    TF_dict1 = {
        'Gene_ID': 'FBgn0004652',
        'Transcript_ID': 'FBtr0083651', 
        'Chromosome': 'Y', 
        'start/stop': 'stop_codon', 
        'genome_stop_codon_pos': 29049-17, 
        'genome_start_codon_pos': 0, 
        'strand_type': '+', 
        'Reference_Seq': 1600 * "C" + "TAG" + 12 * "C"+ "TGG" + 1585 * "C", 
        'sgRNA_list_positions': [[29004, 29027], [29026, 29049]],
        'sgRNA_list_values': [23 * "C", 6 * "C" + "TAG" + 12 * "C"+ "TGG"]
    }

    TF_dict2 = {
        'Gene_ID': 'FBgn0004652',
        'Transcript_ID': 'FBtr0083651', 
        'Chromosome': 'Y', 
        'start/stop': 'start_codon', 
        'genome_stop_codon_pos': 0, 
        'genome_start_codon_pos': 29049-17, 
        'strand_type': '+', 
        'Reference_Seq': 1600 * "C" + "ATG" + 12 * "C"+ "TGG" + 1585 * "C", 
        'sgRNA_list_positions': [[29004, 29027], [29026, 29049]],
        'sgRNA_list_values': [23 * "C", 6 * "C" + "ATG" + 12 * "C"+ "TGG"]
    }

    
    TF_dict3 = {
        'Gene_ID': 'FBgn0004652',
        'Transcript_ID': 'FBtr0083651', 
        'Chromosome': 'Y', 
        'start/stop': 'start_codon', 
        'genome_stop_codon_pos': 0, 
        'genome_start_codon_pos': 29049-17, 
        'strand_type': '+', 
        'Reference_Seq': 1600 * "C" + "ATG" + 12 * "C"+ "TGG" + 1585 * "C", 
        'sgRNA_list_positions': [[29004, 29027]],
        'sgRNA_list_values': [23 * "C"]
    }


    solution1 = check_cutting_site_inside_CDS(TF_dict1)
    solution2 = check_cutting_site_inside_CDS(TF_dict2)
    solution3 = check_cutting_site_inside_CDS(TF_dict3)
    
    if solution1["cut_outside_of_CDS"] == "no" and solution1["sgRNA_list_values"] == [23 * "C"]:

        if solution2["cut_outside_of_CDS"] == "no" and solution2["sgRNA_list_values"] == [6 * "C" + "ATG" + 12 * "C"+ "TGG"]:

            if solution3["cut_outside_of_CDS"] == "yes" and solution3["sgRNA_list_values"] == [23 * "C"]:

                print("The function check_cutting_site_inside_CDS() is working!")

            else: 

                print("Error! A sgRNA that cuts outside the CDS was wrongly classified. Do not run the whole pipeline until this error in check_cutting_site_inside_CDS() is fixed!")

                sys.exit()

        else:

            print("Error! A sgRNA that cuts inside the CDS was wrongly classified. Do not run the whole pipeline until this error in check_cutting_site_inside_CDS() is fixed!")

            sys.exit()
    else: 
        
        print("Error! A sgRNA that cuts inside the CDS was wrongly classified. Do not run the whole pipeline until this error in check_cutting_site_inside_CDS() is fixed!")

        sys.exit()
    
    return

def test_find_synonymous_codons():

    codon_table_excel = "inputfiles/codon_table.xlsx"
    codon_to_mutate = 'GGG'

    synonymous_codons = find_synonymous_codons(codon_to_mutate, codon_table_excel)

    if synonymous_codons == ['GGT', 'GGC', 'GGA']:

        print("The function find_synonymous_codons() is working!")
    
    else: 

        print("Error! find_synonymous_codons is not returning the correct list of codons that encode for the same amino acid as the query codon.")

        sys.exit()
    
    return

def test_mutate_PAM_in_codon():

    codon_to_mutate_1 = 'GGG'
    codon_to_mutate_2 = 'GAA'
    synonymous_codons_1 = ['GGT', 'GGC', 'GGA']
    synonymous_codons_2 = ['AAA']
    synonymous_codons_3 = ['GTT']

    selected_codon_1 = mutate_PAM_in_codon(codon_to_mutate_1, synonymous_codons_1)
    selected_codon_2 = mutate_PAM_in_codon(codon_to_mutate_2, synonymous_codons_2)
    selected_codon_3 = mutate_PAM_in_codon(codon_to_mutate_2, synonymous_codons_3)

    if selected_codon_1 == 'GGT':

        if selected_codon_2 == 'AAA':

            if not selected_codon_3:

                print("The function mutate_PAM_in_codon is working!")

            else: 

                print("Error! A codon with a G as first nucleotide was selected as synonymous codon for a codon with a G, as first letter of the PAM, as first nucleotide.")
                sys.exit()

        else: 

            print("Error! A suitable codon was not selected as synonymous codon for the query codon, which had a G as part of the PAM as the first nucleotide of the codon.")
            sys.exit()

    else: 

        print("Error! A suitable codon was not selected as synonymous codon for the query codon, which had a G as part of the PAM as the last nucleotide of the codon.")
        sys.exit()

    return

def test_make_synonymous_mutation():

    sequence_1 = "AAAAAAAGG"
    position_of_mutation_1 = 8

    sequence_2 = "AAAAAAGGA"
    position_of_mutation_2 = 6


    mutated_sequence_1 = "AAAAAACGT"
    mutated_sequence_2 = ""
    codon_table_excel = "inputfiles/codon_table.xlsx"

    query_2  = make_synonymous_mutation(sequence_2, position_of_mutation_2, codon_table_excel)

    if make_synonymous_mutation(sequence_1, position_of_mutation_1, codon_table_excel) == mutated_sequence_1:

        if query_2 == mutated_sequence_2:

            print("The function make_synonymous_mutation() is working!")

        else:
            
            print("Error! The function make_synonymous_mutation() has returned a synonymous mutation even though none should exist")
            
            sys.exit()


    else:
    
        print("Error! The function make_synonymous_mutations() has returned the wrong synonymous mutation for this sequence")

        sys.exit()

    return

def test_translate_nucleotide_position_into_codon_position():

    sequence = "qwertzuiop"
    
    if translate_nucleotide_position_into_codon_position(sequence, 0) == 0:

        if translate_nucleotide_position_into_codon_position(sequence,2) == 0:

            if translate_nucleotide_position_into_codon_position(sequence,3) ==1:

                print("The function translate_nucleotide_position_into_codon_position() is working!")

            else:

                print("Error! The function translate_nucleotide_position_into_codon_position() has returned the wrong codon position")

                sys.exit()
        
        else:

                print("Error! The function translate_nucleotide_position_into_codon_position() has returned the wrong codon position")

                sys.exit()
    
    else:

                print("Error! The function translate_nucleotide_position_into_codon_position() has returned the wrong codon position")

                sys.exit()

    return

def test_mutate_PAM_in_HDR_plasmid():

    df_1 = {
        'Gene_ID': 'FBgn0004652',
        'Transcript_ID': 'FBtr0083651', 
        'Chromosome': 'Y', 
        'start/stop': 'stop_codon', 
        'genome_stop_codon_pos': 29026, 
        'genome_start_codon_pos': 0, 
        'strand_type': '+', 
        'Reference_Seq': 1600 * "C" + "TAG" + 12 * "C"+ "TGG" + 1585 * "C", 
        'sgRNA_list_positions': [29004, 29027],
        'sgRNA_list_values': 23 * "C"
        }
    
    df_2 = {
        'Gene_ID': 'FBgn0004652',
        'Transcript_ID': 'FBtr0083651', 
        'Chromosome': 'Y', 
        'start/stop': 'stop_codon', 
        'genome_stop_codon_pos': 29025-16, 
        'genome_start_codon_pos': 0, 
        'strand_type': '+', 
        'Reference_Seq': 1600 * "C" + "TAG" + 12 * "C"+ "TGG" + 1585 * "C", 
        'sgRNA_list_positions': [29004, 29027],
        'sgRNA_list_values': 23 * "C"
        }

    df_3 = {
        'Gene_ID': 'FBgn0004652',
        'Transcript_ID': 'FBtr0083651', 
        'Chromosome': 'Y', 
        'start/stop': 'stop_codon', 
        'genome_stop_codon_pos': 29026, 
        'genome_start_codon_pos': 0, 
        'strand_type': '+', 
        'Reference_Seq': 1600 * "C" + "TAG" + 12 * "C"+ "TGG" + 1585 * "C", 
        'sgRNA_list_positions': [29004, 29027],
        'sgRNA_list_values': 23 * "C"
        }

    HAL_R_1 = 12 * "C"+ "CGG"

    HAR_F_1 = 15 * "C"

    HAL_R_2 = 15 * "C"

    HAR_F_2 = 12 * "C"+ "CGG"

    HAL_R_3 = 12 * "C"+ "TGG"

    HAR_F_3 = 15 * "C"

    df_out_1 = {
        'Gene_ID': 'FBgn0004652',
        'Transcript_ID': 'FBtr0083651', 
        'Chromosome': 'Y', 
        'start/stop': 'stop_codon', 
        'genome_stop_codon_pos': 29026, 
        'genome_start_codon_pos': 0, 
        'strand_type': '+', 
        'Reference_Seq': 1600 * "C" + "TAG" + 12 * "C"+ "TGG" + 1585 * "C", 
        'sgRNA_list_positions': [29004, 29027],
        'sgRNA_list_values': 23 * "C",
        'HAL-R': 12 * "C"+ "CGT",
        'HAR-F': 15 * "C"
        }

    df_out_2 = {
        'Gene_ID': 'FBgn0004652',
        'Transcript_ID': 'FBtr0083651', 
        'Chromosome': 'Y', 
        'start/stop': 'stop_codon', 
        'genome_stop_codon_pos': 29025-16, 
        'genome_start_codon_pos': 0, 
        'strand_type': '+', 
        'Reference_Seq': 1600 * "C" + "TAG" + 12 * "C"+ "TGG" + 1585 * "C", 
        'sgRNA_list_positions': [29004, 29027],
        'sgRNA_list_values': 23 * "C",
        'HAL-R': 15 * "C",
        'HAR-F': 12 * "C"+ "CGT"
        }

    df_out_3 = {
        'Gene_ID': 'FBgn0004652',
        'Transcript_ID': 'FBtr0083651', 
        'Chromosome': 'Y', 
        'start/stop': 'stop_codon', 
        'genome_stop_codon_pos': 29026, 
        'genome_start_codon_pos': 0, 
        'strand_type': '+', 
        'Reference_Seq': 1600 * "C" + "TAG" + 12 * "C"+ "TGG" + 1585 * "C", 
        'sgRNA_list_positions': [29004, 29027],
        'sgRNA_list_values': 23 * "C"
    }

    mutated_primer_1 = mutate_PAM_in_HDR_plasmid(HAL_R_1, HAR_F_1, df_1)
    mutated_primer_2 = mutate_PAM_in_HDR_plasmid(HAL_R_2, HAR_F_2, df_2)
    mutated_primer_3 = mutate_PAM_in_HDR_plasmid(HAL_R_3, HAR_F_3, df_3)  
    
    if mutated_primer_1 == df_out_1:

        if mutated_primer_2 == df_out_2:
        
            if mutated_primer_3 == df_out_3:

                print("The function mutate_PAM_in_HDR_plasmid() is working!")

            else: 

                print("The function mutate_PAM_in_HDR_plasmid() is returning mutated primer even though there is no synonymous codon for the PAM codon")

                sys.exit()
            
        else:
            
            print("The function mutate_PAM_in_HDR_plasmid() is returning the wrong mutated primer for PAM in the HAR_F")

            sys.exit()

    else:
    
        print("The function mutate_PAM_in_HDR_primer() is returning the wrong mutated primer for PAM in the HAL_R")

        sys.exit()

    