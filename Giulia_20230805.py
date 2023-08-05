def concatenate_for_primer_mutation(primers_output_file, optimal_sgRNA_output_file):
    """"
    this function concatenates the output file of the primers functions with the optimal sgRNA files
    
    input:
    primers_output_file: excel file with primers. This is the output of the functions that design primers.
    optimal_sgRNA_output_file: excel file with sgRNA. Output if sgRNA design function.

    output:
    concat_primers: concatenated excel with the inputs
    
    """""
    import pandas as pd

    sgRNA_doc = pd.read_excel(optimal_sgRNA_output_file)
    primer_doc = pd.read_excel(primers_output_file)
    sgRNA_doc.rename(columns={'gene_ID': 'Gene_ID','transcript_ID':'Transcript_ID','start/stop':'Gene_Region'}, inplace=True)

    #concatenate excel files
    concat_primers = pd.merge(sgRNA_doc, primer_doc, on=['Gene_ID','Transcript_ID','Gene_Region'])
    concat_primers = concat_primers.drop(columns=["Unnamed: 0"])
    concat_primers = concat_primers.reset_index()

    concat_primers.to_excel("outputFiles\output_concat_primers_and_sgRNA_prior_to_mutation.xlsx")
    return concat_primers


def revComp(inputSeq):
    """
    This function takes an input sequence and returns the reverse complement.

    Input: inputSeq in str format
    Output: revComp in str format

    """
    complement = {'A': 'T', 'C': 'G', 'G': 'C', 'T': 'A'}

    revComp = ""
    for base in inputSeq[::-1]:
        revComp += complement[(base.upper())]

    return revComp


def primer_parser(primer_sequence, primer_type, strand_type):

    primer_length = len(primer_sequence)
    leftover_size = primer_length%3
    splitat = primer_length - leftover_size

    if strand_type == '-':

        if primer_type == 'hal-r':
            shortened_primer, leftover = primer_sequence[:splitat], primer_sequence[splitat:]
        if primer_type == 'har-f':
            primer_sequence_r = revComp (primer_sequence)
            shortened_primer, leftover = primer_sequence_r[:leftover_size], primer_sequence_r[leftover_size:]

    if strand_type == '+':
        if primer_type == 'hal-r':
            primer_sequence_r = revComp(primer_sequence)
            shortened_primer, leftover = primer_sequence_r[:leftover_size], primer_sequence_r[leftover_size:]
        if primer_type == 'har-f':
            shortened_primer, leftover = primer_sequence[:splitat], primer_sequence[splitat:]

    return shortened_primer, leftover


def mutate_PAM_in_primers(concat_primers):

    concat_primers['mutated_primers_?'] = 'mutation was not necessary'

    for index, row in concat_primers.iterrows():

        if row['must_PAM_be_mutated_in_HDR_plasmid?'] != 'no':

            #@TODO: must make sure that df is being passed in the right df format and that nothing is being left out
            df = {
                "start/stop": row['Gene_Region'],
                'genome_start_codon_pos': row['genome_start_codon_pos'],
                'genome_stop_codon_pos': row['genome_stop_codon_pos'],
                'strand_type': row['strand_type'],
                "sgRNA_list_positions": [row['sgRNA_start'], row['sgRNA_start']],
                "sgRNA_list_values": row['sgRNA_seq']
            }

            hal_r_parsed, hal_r_leftover = primer_parser(row['HAL-R'],'hal-r',row['strand_type'])
            har_f_parsed, har_f_leftover = primer_parser(row['HAR-F'],'har-f',row['strand_type'])

            #@TODO - uncomment this - df_updated = mutate_PAM_in_HDR_plasmid(hal_r_parsed,har_f_parsed,df)

            df_updated = {
                "start/stop": row['Gene_Region'],
                'genome_start_codon_pos': row['genome_start_codon_pos'],
                'genome_stop_codon_pos': row['genome_stop_codon_pos'],
                'strand_type': row['strand_type'],
                "sgRNA_list_positions": [row['sgRNA_start'], row['sgRNA_start']],
                "sgRNA_list_values": row['sgRNA_seq'],
                "HAR-F": "ParsedAndMutated",
                "HAL-R": "ParsedAndMutated",
                "mutated?": "MutationPlaceholder"
            }

            #1. check what was mutated - add this info to the extra column
            #add column for the Pam in concatenate_for_primer_mutation()
            #fill in

            #row['HAR-F']=df_updated['HAR-F']
            #row['HAL-R'] = df_updated['HAL-R']


    #concat_primers.to_excel("outputFiles\output_concat.xlsx")
    print(concat_primers)
    return concat_primers


optimal_sgRNA_output_file = "inputfiles\mockMaterials\optimal_sgRNAs_mock.xlsx"
primers_outputfile_nomutation = "inputfiles\mockMaterials\mockTFsdfwithPrimers.xlsx"

concatenated_primers_sgRNA = concatenate_for_primer_mutation(primers_outputfile_nomutation, optimal_sgRNA_output_file)
mutated_primers = mutate_PAM_in_primers(concatenated_primers_sgRNA)

#check if the PAM is