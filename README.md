# Introducing AutoTagsCRISPR

AutoTagsCRISPR is a pipeline for automated CRISPR construct design to tag genes at all annotated termini (i.e., at every start and stop codon). The CRISPR/Cas9 system allows to insert tags (e.g. fluorescent tags) into the genome to investigate the function of previously unstudied genes or transcripts. When binding to its DNA target, the CRISPR/Cas9 system cuts at a defined position. A donor DNA fragment with flanking regions homologous to the flanking regions of the cut site can then be integrated into the genome via a process referred to as Homology Directed Repair (HDR). To target the CRISPR/Cas9 system to the correct locus within the genome, a suitable sgRNA needs to be designed. To insert the donor DNA fragment, a Homology Arms (HA) need to be designed. Part of the design process is to ensure that the designed sgRNA does not cut inside the HA since this would inhibit the HDR process required for inserting the tag. AutoTagsCRISPR automates both the sgRNA and HA design and ensures the HAs are cut-proof.

Here, we will demonstrate the usability and logic of AutoTagsCRISPR using a Drosophila Melanogaster project as an example. :fly:
To understand the function of a transcription factor (TF), it is necessary to study its tissue distribution, binding characteristics at physiological concentration and chromatin accessibility state, effect on transcription, and protein interactions. To investigate these properties, S. Kittelmann et al.,  propose to generate a biological resource to enable the in-depth study of TF function in Drosophila. This resource will consist of three parts: a set of plasmids for tagging all Drosophila TFs with an exchangeable epitope, fly lines in which TFs havebeen tagged, and a database with expression and binding information for a subset of previously unstudied TFs. S. Kittelmann et al., will insert a superfolder-GFP (sfGFP) to 1) tag specifically TFs; 2) tag at the endogenous genomic location, capturing all regulatory information; 3) tag all isoforms with different N and C-termini; 4) allow easy tag exchange; 5) support easy removal of transgenic markers, allowing virtually scarless gene editing. There are on average 2.56 termini annotated per TF gene, amounting to 1,915 CRISPR constructs in total for the 753 TF genes in the Drosophila genome. AutoTagsCRISPR will allow for an automated CRISPR construct design to speed up the design process and minimise costs.

![image](https://user-images.githubusercontent.com/120821707/210607784-b8ccab0c-a99f-46fc-afdc-5f6f702fe3a1.png)

---

## Table of Contents

### How to set up the pipeline and run the pipeline for our example use case
```bash
# clone GitHub repository and move to workspace
git clone https://github.com/emmacwatts/AutoTagsCRISPR.git
cd AutoTagsCRISPR

# install dependencies
conda env create -f environment.yml
conda activate AutoTagsEnv
```
Go to our GoogleDrive repositiory and 
Now you need to save the following files in the folder /AutoTagsCRISPR/inputfiles:
- dmel-all-r6.48.gtf
- dmel-all-chromosome-r6.48.fasta

And the following files in the folder /AutoTagsCRISPR/inputfiles/sgRNAFiles:
- NoOffTarget_high_stringency.gff
- NoOffTarget_med_stringency.gff
- NoOffTarget_low_stringency.gff
- 1to3NonCdsOffTarget_low_stringency.gff
- ManyOffTarget_low_stringency.gff
  
```bash
# run jupyter notebook tests
jupyter notebook sgRNA_tests_window_21_pb.ipynb
jupyter notebook sgRNA_tests_window_42_pb.ipynb

# run tests with shorter mock files for a window around start/stop codon of 21 bp and 42 bp
python sgRNArunner.py -window=21 -inputfile="inputfiles/mockMaterials/TFsTruncatedLong.xlsx"
python sgRNArunner.py -window=42 -inputfile="inputfiles/mockMaterials/TFsTruncatedLong.xlsx"

# run whole pipeline
python sgRNArunner.py -window=21 -inputfile="inputfiles/TFs.xlsx"
python sgRNArunner.py -window=42 -inputfile="inputfiles/TFs.xlsx"
```
### [IOguidance](https://github.com/emmacwatts/AutoTagsCRISPR/tree/main/IOguidance)

Demonstrate input and output formats for your functions here with screenshots/descriptions. This will ensure our different code sections match up.

### [finalCode](https://github.com/emmacwatts/AutoTagsCRISPR/tree/main/finalCode)

This is where final, compiled scripts can go. Ideally we want everything we want to package in this directory only.

### [indepCode](https://github.com/emmacwatts/AutoTagsCRISPR/tree/main/indepCode)

This directory contains indpendent working spaces - this is not the final code.
Might be a good idea to just put a copy of your independent code here so others can see it/move snippets to our final code, etc.

### [inputFiles](https://github.com/emmacwatts/AutoTagsCRISPR/tree/main/inputfiles)

This contains the input and advice files we were provided with. Note that .fasta and.gtf files were too large to save here.


## Dependencies

Latest versions of biopython, primer3-py, pandas, and tqdm.

## Primers function update - key information sheet

![image](https://github.com/emmacwatts/AutoTagsCRISPR/files/10706002/CRISPR.primers.workflow.pdf)
![image](https://github.com/emmacwatts/AutoTagsCRISPR/files/10706079/CRISPR.primers.workflow.pdf)

### Note on position count

positionCount in utils20230401 (mutate_sgRNA_recognition_site_in_HDR_plasmid function) indicates the number of bases in the segment highlighted in red in the two possible scenarios.

![image]("https://github.com/emmacwatts/AutoTagsCRISPR/assets/120821707/09faf6d0-a22b-4d64-82d9-071158f65e63">)

##Notes for meeting with S/K
* Problem - PAM not in sgRNA sequence given - affects downstream mutation
* Mutating PAM - would need to identify outside? or could mutate 3 codons in every case, but then primers might not work
