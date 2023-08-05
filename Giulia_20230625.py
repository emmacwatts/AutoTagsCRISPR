def Check_and_mutate_HDR_plasmid(primers_output_file, optimal_sgRNA_output_file):
    import pandas as pd
    import numpy as np
    import utils_20230401 as fct

    sgRNA_doc = pd.read_excel(optimal_sgRNA_output_file)
    primer_doc = pd.read_excel(primers_output_file)
    sgRNA_doc.rename(columns={'gene_ID': 'Gene_ID','transcript_ID':'Transcript_ID','start/stop':'Gene_Region'}, inplace=True)

    #concatenate excel files
    concat_primers = pd.merge(sgRNA_doc, primer_doc, on=['Gene_ID','Transcript_ID','Gene_Region'])
    concat_primers = concat_primers.drop(columns=["Unnamed: 0"])
    concat_primers = concat_primers.reset_index()

    print(concat_primers)

    for index, row in concat_primers.iterrows():
        df = {
        "start/stop": row['Gene_Region'],
        'genome_start_codon_pos': row['genome_start_codon_pos'],
        'genome_stop_codon_pos': row['genome_stop_codon_pos'],
        'strand_type': row['strand_type'],
        "sgRNA_list_positions": [row['sgRNA_start'], row['sgRNA_start']],
        "sgRNA_list_values": row['sgRNA_seq']
        }

        print(df)

        #df_updated = fct.mutate_PAM_in_HDR_plasmid(row['HAL-R'],row['HAR-F'],df)
        df_updated = {
        "start/stop": row['Gene_Region'],
        'genome_start_codon_pos': row['genome_start_codon_pos'],
        'genome_stop_codon_pos': row['genome_stop_codon_pos'],
        'strand_type': row['strand_type'],
        "sgRNA_list_positions": [row['sgRNA_start'], row['sgRNA_start']],
        "sgRNA_list_values": row['sgRNA_seq'],
        "HAR-F" : "AGCCCTT",
        "HAL-R" : "AGGCCTTGG"
        }

        row['HAR-F']=df_updated['HAR-F']
        print(row['HAR-F'])
        row['HAL-R'] = df_updated['HAL-R']


    concat_primers.to_excel("outputFiles\output_concat.xlsx")
    return concat_primers


optimal_sgRNA_output_file = "inputfiles\mockMaterials\optimal_sgRNAs_mock.xlsx"
outputfile = "inputfiles\mockMaterials\mockTFsdfwithPrimers.xlsx"

df = Check_and_mutate_HDR_plasmid(outputfile, optimal_sgRNA_output_file)








