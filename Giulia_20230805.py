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

def mutate_PAM_in_primers(concat_primers):

    #for index, row in concat_primers.iterrows():

        #body of the function, or call out to fucntion that actually mutates
        # if column says yes then pass to function
        #df_updated = fct.mutate_PAM_in_HDR_plasmid(row['HAL-R'],row['HAR-F'],df)
        #1. check what was mutated - add this info to the extra column
        #add column for the Pam in concatenate_for_primer_mutation()
        #fill in

        #row['HAR-F']=df_updated['HAR-F']
        #row['HAL-R'] = df_updated['HAL-R']


    #concat_primers.to_excel("outputFiles\output_concat.xlsx")
    return concat_primers


optimal_sgRNA_output_file = "inputfiles\mockMaterials\optimal_sgRNAs_mock.xlsx"
primers_outputfile_nomutation = "inputfiles\mockMaterials\mockTFsdfwithPrimers.xlsx"

concatenated_primers_sgRNA = concatenate_for_primer_mutation(primers_outputfile_nomutation, optimal_sgRNA_output_file)
print(concatenated_primers_sgRNA)

