# AutoTagsCRISPR

Building a pipeline for automated CRISPR construct design in Drosophila Melanogaster. :fly:

Data is organised from input files with gene IDs, reference genome sequences and annotations, and all potential guideRNA sequences.
The best guide RNA is chosen for the given parameters, and primers are designed for the left and right homology arms of this region
as well as two verification primers outside of this region.

![image](https://user-images.githubusercontent.com/120821707/210607784-b8ccab0c-a99f-46fc-afdc-5f6f702fe3a1.png)

---

## Table of Contents

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
