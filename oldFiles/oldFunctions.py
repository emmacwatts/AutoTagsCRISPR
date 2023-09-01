
def pamPosition(HAL_R, HAR_F, df):

    from os import sys

    codon_table_excel = "inputfiles/codon_table.xlsx"

    if df["start/stop"] == "start_codon":

        pos_of_interest = df["genome_start_codon_pos"]

    if df["start/stop"] == "stop_codon":

        pos_of_interest = df["genome_stop_codon_pos"]

    if df['sgRNA_strand'] == df['strand_type']: #same strand

        PAM_pos = df["sgRNA_list_positions"][1] #use end of sgRNA

        # check whether PAM in HAL-R

        if PAM_pos < pos_of_interest:

            # calculate distance between end of HAL-R and end of PAM

            PAM_pos_in_HAL_R = PAM_pos - pos_of_interest

            # try mutating the PAM in the HDR fragment
            mutated_HAL_R = make_synonymous_mutation(HAL_R, PAM_pos, codon_table_excel)

        elif PAM_pos-2 > pos_of_interest + 2:


    # check whether PAM in HAR-F
    
    elif df["sgRNA_list_positions"][1] - 2 - pos_of_interest >= 16:

                
    return df

