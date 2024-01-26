def mutate_PAM_in_HDR_plasmid(HAL_R, HAR_F, df):
    from os import sys

    if df["start/stop"] == "start_codon":
        pos_of_interest = df["genome_start_codon_pos"]
    elif df["start/stop"] == "stop_codon":
        pos_of_interest = df["genome_stop_codon_pos"]

    codon_table_excel = "inputfiles/codon_table.xlsx"

    # check whether PAM in HAL-R
    if 0 >= df["sgRNA_list_positions"][1] - 2 - pos_of_interest:
        PAM_pos = HAL_R.find("GG")

        # sanity check since function x.find() returns -1 if string is not found
        if PAM_pos == -1:
            #print(HAL_R, df)
            #print("Error! No GG found in HAL-R even though PAM should be located in HAL-R")
            df["HAL-R"] = HAL_R
            df["HAR-F"] = HAR_F
            df["mutated?"] = "no"
            df["PAM in primer?"] = "not found in HAL-R"


        else:
            # try mutating the PAM in the primer/HDR fragment
            mutated_HAL_R = make_synonymous_mutation(HAL_R, PAM_pos, codon_table_excel)
            # if you were able to mutate the PAM in the primer/ fragment sequence, add mutated primer
            # and information about which primer was mutated to the dataframe
            if mutated_HAL_R:
                df["HAL-R"] = mutated_HAL_R
                df["HAR-F"] = HAR_F
                df["mutated?"] = "HAL_R"
                df["PAM in primer?"] = "yes"
            # if you were not successful in mutating the PAM in the primer/ fragment sequence,
            # try mutating the sgRNA recognition sequence in the primer/ fragement sequence
            else:
                # calculate the distance between the start/ stop condon and the end of the PAM
                # to determine where the sgRNA recognition site is located in the primer/ fragement sequence
                distance = pos_of_interest - df["sgRNA_list_positions"][1] - 2
                mutated_HAL_R = mutate_sgRNA_recognition_site_in_HDR_plasmid('HAL_R', HAL_R, distance,)

                if mutated_HAL_R:
                    # if mutation was successful,  update the sgRNA dictionary with mutated primers
                    # and add information about mutation status
                    df["HAR-F"] = HAR_F
                    df["HAL-R"] = mutated_HAL_R
                    df["mutated?"] = "HAL_R"
                    df["PAM in primer?"] = "yes"
                else:
                    # if mutation was not successful, add information about mutation status
                    # and keep old unmutated primer
                    df["HAR-F"] = HAR_F
                    df["HAL-R"] = HAL_R
                    df["mutated?"] = "no"
                    df["PAM in primer?"] = "yes"

    # check whether PAM in HAR-F
    elif df["sgRNA_list_positions"][1] - 2 - pos_of_interest >= 16:
        # sanity check since funtion x.find() returns -1 if string is not found
        PAM_pos = HAR_F.find("GG")

        if PAM_pos == -1:
            #print("Error! No GG found in HAR-F even though PAM should be located in HAR-F")
            #sys.exit()
            df["HAL-R"] = HAL_R
            df["HAR-F"] = HAR_F
            df["mutated?"] = "no"
            df["PAM in primer?"] = "not found in HAR-F"
        else:
            # try mutating the PAM in the primer/ fragment sequence
            mutated_HAR_F = make_synonymous_mutation(HAR_F, PAM_pos, codon_table_excel)

            # if you were able to mutate the PAM in the primer sequence, add mutated primer
            # and information about which primer was mutated to the dataframe
            if mutated_HAR_F:
                df["HAR-F"] = mutated_HAR_F
                df["HAL-R"] = HAL_R
                df["mutated?"] = "HAR_F"
                df["PAM in primer?"] = "yes"
            # if there was no synonymous mutation for the PAM, mutate the sgRNA recognition site

            else:
                # calculate distance from the beginning of the PAM to the beginning of the start/ stop codon
                distance = df["sgRNA_list_positions"][1] - 2 - pos_of_interest - 3
                # try to mutate sgRNA recognition site in primer
                mutated_HAR_F = mutate_sgRNA_recognition_site_in_HDR_plasmid('HAR_F', HAR_F, distance)

                if df:
                    # if mutation was successful,  update the sgRNA dictionary with mutated primers
                    # and add information about mutation status
                    df["HAR-F"] = mutated_HAR_F
                    df["HAL-R"] = HAL_R
                    df["mutated?"] = "HAR_F"
                    df["PAM in primer?"] = "yes"
                else:
                    # if mutation was not successful, add information about mutation status
                    # and keep old unmutated primer
                    df["HAR-F"] = HAR_F
                    df["HAL-R"] = HAL_R
                    df["mutated?"] = "no"
                    df["PAM in primer?"] = "yes"
    return df



def make_synonymous_mutation(sequence, position_of_mutation, codon_table_excel):
    '''
    Mutates a G from the PAM in a sequence in a way, that the codon containing the G still encodes for the same amino acid.

    params:
        sequence: string, containing the sequence that is supposed to be mutated.
        positions_of_mutation: int, position of nucleotide that is supposed to be mutated in sequence starting with 0.
    returns: string, mutated sequence.
    '''

    # a codon is 3 nucleotides long
    x = 3
    list_of_codons = [sequence[y - x:y] for y in range(x, len(sequence) + x, x)]
    codon_position = translate_nucleotide_position_into_codon_position(sequence, position_of_mutation)
    codon_to_mutate = list_of_codons[codon_position]
    synonymous_codons = find_synonymous_codons(codon_to_mutate, codon_table_excel)
    selected_codon = mutate_PAM_in_codon(codon_to_mutate, synonymous_codons)

    if selected_codon:
        list_of_codons[codon_position] = selected_codon
        mutated_sequence = "".join(list_of_codons)

    else:
        mutated_sequence = ""

    return mutated_sequence


def translate_nucleotide_position_into_codon_position(sequence, nucleotide_position):

    '''
    params:
        sequence: string, nucleotide sequence
        nucleotide_position: integer, position of nucleotide that you want to translate into the codon position
    returns:
        count: integer, position of codon containing nucleotide from nucleotide position
    '''

    count = 0
    # since position of first letter in string is zero, function has to be adapted accordingly
    for n in range(0,len(sequence)):

        if nucleotide_position - 3 >= 0:
            count = count +1
            nucleotide_position = nucleotide_position - 3
        else:
            break

    return(count)


def find_synonymous_codons(query_codon, codon_table_excel):
    '''
    Uses the amino acids table to select codons that encode for the same amino acid as the query codon.

    Params:
        codon: string, codon to select synonymous codons for
        codon_table_excel: string, path to an excel file that lists per codon which amino acid that codon encodes.
    Returns:
        synonymous_codons: list of strings, each string is a codon that encodes for the same amino acid as the query codon.

    '''
    import pandas as pd
    synonymous_codons = []
    codon_table = pd.read_excel(codon_table_excel)
    amino_acid = codon_table.query(f"codon =='{query_codon}'")
    amino_acid_query = amino_acid.iloc[0]['amino_acid']
    codons = codon_table.query(f"amino_acid=='{amino_acid_query}' & codon !='{query_codon}'")

    for codon in range(0, len(codons.index)):
        synonymous_codon = codons.iloc[codon]['codon']
        synonymous_codons.append(synonymous_codon)

    return synonymous_codons


def mutate_PAM_in_codon(query_codon, synonymous_codons):
    '''
    Mutates one of the G`s from the PAM in the query codon to another nucleotide based on the codon list provided.
    Returns an empty list if mutation was not possible.

    params:
        query_codon: string, codon containing the first G of the PAM.
        synonymous_codons: list of strings, codons that incode for the same amino acid as the query codon.
    returns:
        selected_codon: string, codon that is the same as the query codon apart from one G, which also encodes for the same amino acid as the query codon.
    '''

    if synonymous_codons:
        list_query_codon = list(query_codon)

        if list_query_codon[2] == 'G':

            for synonymous_codon in synonymous_codons:
                list_synonymous_codon = list(synonymous_codon)

                if list_synonymous_codon[2] != 'G':
                    selected_codon = synonymous_codon
                    break
                else:
                    selected_codon = ''

        elif list_query_codon[0] == 'G':

            for synonymous_codon in synonymous_codons:
                list_synonymous_codon = list(synonymous_codon)

                if list_synonymous_codon[0] != 'G':
                    selected_codon = synonymous_codon
                    break
                else:
                    selected_codon = ''
        else:
            print('There is no synonymous codon that can be used to mutate the PAM. Searching for other sgRNAs.')
            selected_codon = ''
    else:
        print('There is no synonymous codon that can be used to mutate the PAM. Searching for other sgRNAs.')
        selected_codon = ''

    return selected_codon


def mutate_sgRNA_recognition_site_in_HDR_plasmid(sequenceType, sequenceToMutate, positionCount):
    """
    Mutates HDR arm fragments in the sgRNA recognition site in the 2-3 codons closest to PAM in the case where PAM itself cannot be mutated. 1-2 mutations are made in separate amino acids.

    params:
        sequenceType: 'HAL-R' or 'HAR-F' (corresponding to left or right homology arm)
        sequenceToMutate: the primer/fragment sequence. If reverse primer, this should be the + strand (revComp)
        positionCount: Distance of start/stop from PAM (excluding NGG).

    returns:
        newSequence: The mutated HDR arm. If no mutations were not possible, an empty string will be returned.

    """
    # Define all codons in the homology arm fragment
    allCodons = []
    for codonStart in range(0, len(sequenceToMutate), 3):
        allCodons.append(sequenceToMutate[codonStart:codonStart + 3])

    if sequenceType == 'HAL-R':
        allCodons.reverse()  # For HAL-R, order codons moving away from te STOP site

    # Take the 1-3 codons closest to PAM in order
    if positionCount <= 3:  # Can only attempt to mutate one codon if positionCount is less than three
        orderedMutableCodons = {0: allCodons[0]}
    elif positionCount <= 6:  # Can only attempt to mutate two codon if positionCount is less than six
        orderedMutableCodons = {0: allCodons[0], 1: allCodons[1]}
    else:
        PAMcodon = int((positionCount - 1) / 3)  # This is the codon number that PAM is in
        if positionCount % 3 != 0:  # If PAM is not in frame, take the codon that it is in and the two after that
            orderedMutableCodons = {PAMcodon: allCodons[PAMcodon], PAMcodon - 1: allCodons[PAMcodon - 1],
                                    PAMcodon - 2: allCodons[PAMcodon - 2]}
        else:  # If in frame, just take the 2 codons after PAM
            orderedMutableCodons = {PAMcodon - 1: allCodons[PAMcodon - 1], PAMcodon - 2: allCodons[PAMcodon - 2]}

            # This is the count of mutations completed in the sequence.
    mutatedCount = 0
    newSequence = sequenceToMutate  # Make copy to avoid mutating original sequence

    # loop through the codons and mutate if possible to a maximum of two mutations.
    for key, codon in orderedMutableCodons.items():
        # find synonymous codons - will provide list of synonymous codons
        codonOptions = find_synonymous_codons(codon)

        # if the list is not empty, select the first synonymous codon and replace this within the mutable region. Add 1 to mutatedCount.
        if codonOptions != []:
            newcodon = codonOptions[0]
            newSequence = newSequence[:(len(allCodons) - key - 1) * 3] + newcodon + newSequence[(len(
                allCodons) - key) * 3:]  # replace the appropriate codon in sequenceToMutate
            mutatedCount += 1
        # Once you reach 2, stop mutating
        if mutatedCount == 2:
            print(mutatedCount)
            break

    if mutatedCount == 0:  # If nothing was mutated, return an empty string
        return ""
    else:
        return newSequence  # Otherwise, return the mutated sequence


def retrieve_HDR_arm(df):
    """
    Accesses the HDR arm fragment from the dataframe with the start/ stop information about the transcription factors

    params:
            df: dataframe with start/stop information about the transcription factors

    returns:
            fragment: fragment to synthesize for obtaining the HDR arm
    """

    HAL = df["upstreamHA"]
    HAR = df["downstreamHA"]

    return HAL, HAR