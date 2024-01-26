# AutoTagsCRISPR

Building a pipeline for automated CRISPR construct design in Drosophila Melanogaster. :fly:

Data is organised from input files with gene IDs, reference genome sequences and annotations, and all potential guideRNA sequences.
The best guide RNA is chosen for the given parameters, and primers are designed for the left and right homology arms of this region
as well as two verification primers outside of this region.

![image](https://user-images.githubusercontent.com/120821707/210607784-b8ccab0c-a99f-46fc-afdc-5f6f702fe3a1.png)

---

## Table of Contents

### How to set up the pipeline
```bash
# clone GitHub repository and move to workspace
git clone https://github.com/emmacwatts/AutoTagsCRISPR.git
cd AutoTagsCRISPR

# install dependencies
conda env create -f environment.yml
conda activate AutoTagsEnv

# run jupyter notebook tests
jupyter notebook blob/main/sgRNA_tests_window_21_pb.ipynb
jupyter notebook blob/main/sgRNA_tests_window_42_pb.ipynb

# run whole pipeline
python sgRNArunner.py
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
