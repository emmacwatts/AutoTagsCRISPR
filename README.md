# Introducing AutoTagsCRISPR

AutoTagsCRISPR is a pipeline for automated CRISPR construct design to tag genes at all annotated termini (i.e., at every start and stop codon). The CRISPR/Cas9 system allows insertion of tags (e.g. fluorescent tags) into the genome to investigate the function of previously unstudied genes or transcripts. When binding to its DNA target, the CRISPR/Cas9 system cuts at a defined position. A donor DNA fragment with flanking regions homologous to the flanking regions of the cut site can then be integrated into the genome via a process referred to as Homology Directed Repair (HDR). To target the CRISPR/Cas9 system to the correct locus within the genome, a suitable sgRNA needs to be designed. For the donor DNA fragment, Homology Arms (HA) need to be designed. Part of the design process is to ensure that the designed sgRNA does not cut inside the HA since this would inhibit the HDR process. AutoTagsCRISPR automates both the sgRNA and HA design and ensures that the HAs are cut-proof.

Here, we will demonstrate the usability and logic of AutoTagsCRISPR using a Drosophila Melanogaster project as an example. :fly:
To understand the function of a transcription factor (TF), it is necessary to study its tissue distribution, binding characteristics at physiological concentration, chromatin accessibility state, effect on transcription, and protein interactions. To investigate these properties, S. Kittelmann et al.,  propose to generate a biological resource to enable the in-depth study of TF function in Drosophila. This resource will consist of three parts: a set of plasmids for tagging all Drosophila TFs with an exchangeable cargo, fly lines in which TFs have been tagged, and a database with expression and binding information for a subset of previously unstudied TFs. S. Kittelmann et al., will insert a superfolder-GFP (sfGFP) to 1) tag TFs specifically; 2) tag at the endogenous genomic location, capturing all regulatory information; 3) tag all isoforms with different N and C-termini; 4) allow easy tag exchange; 5) support easy removal of transgenic markers, allowing virtually scarless gene editing. There are on average 2.56 termini annotated per TF gene, amounting to 1,915 CRISPR constructs in total for the 753 TF genes in the Drosophila genome. AutoTagsCRISPR will allow for an automated CRISPR construct design to speed up the design process and minimise costs.

For further information read [PipelineLogic.pptx](./PipelineLogic.pptx) in Slide Show mode.

![image](https://user-images.githubusercontent.com/120821707/210607784-b8ccab0c-a99f-46fc-afdc-5f6f702fe3a1.png)

---
# Introducing the folder structure
## [inputfiles](./inputfiles)
Contains genome sequence files, genome annotation files, codon table files and files describing the realtive position of sgRNA and the start/stop codon.

## [outputFiles](./outputFiles)
Location where the pipeline output is saved. Output consists of excel sheet listening information about the TF transcripts and their respective sgRNA and HA.

## [oldFiles](./oldFiles)
Files that are NOT necessary to run the pipeline but are relicts of the development process. Should be deleted if possible.

## [primerScripts](./primerScripts)
Files that are NOT necessary to run the pipeline but are relicts of the development process to design primers. These scripts are relevant for the case that HAs cannot be synthesized and have to be cloned from the Drosophila genome. The scripts allow to design suitable primers to 1) clone the Homology arms and 2) verify the successful tagging of the start/stop codon.

## [PipelineLogic.pptx](./PipelineLogic.pptx)
Guide describing the logic and usability of the AutoTagsCRISPR pipeline. Before building on this pipeline it would be advisable to carefully read this document. Read the document in Slide Show mode as there are animations that will help you understand what is going on.

---
# How to set up the pipeline and run the pipeline for our example use case

1. Set up a local version of this GitHub repository

   ```bash
   # clone GitHub repository and move to workspace
   git clone https://github.com/emmacwatts/AutoTagsCRISPR.git
   cd AutoTagsCRISPR

   # install dependencies
   conda env create -f environment.yml
   conda activate AutoTagsEnv
   ```

2. Go to our [Google Drive](https://drive.google.com/drive/u/2/folders/162uoqatwaCwedfsUMTlAHNHkJTP7n6ki)
   
3. Download and save the following files in [inputfiles](inputfiles):
   - dmel-all-r6.48.gtf
   - dmel-all-chromosome-r6.48.fasta

4. Create a new folder in [inputfiles](inputfiles) called sgRNAFiles:

   ```bash
   # make sgRNAFiles folder
   mkdir inputfiles/sgRNAFiles
   ```
    
5. Download and save the following files in the freshly created inputfiles/sgRNAFiles folder:
   - NoOffTarget_high_stringency.gff
   - NoOffTarget_med_stringency.gff
   - NoOffTarget_low_stringency.gff
   - 1to3NonCdsOffTarget_low_stringency.gff
   - ManyOffTarget_low_stringency.gff

6. Test for successful implementation
   
   ```bash
   # run jupyter notebook tests
   jupyter notebook sgRNA_tests_window_21_pb.ipynb
   jupyter notebook sgRNA_tests_window_42_pb.ipynb
   
   # run tests with shorter mock files for a window left and right of the annotated termini of 21 bp and 42 bp
   # print statements will get stored in text files in outputFiles
   python sgRNArunner.py "inputfiles/mockMaterials/TFsTruncatedLong.xlsx" "inputfiles/dmel-all-chromosome-r6.48.fasta" "inputfiles/dmel-all-r6.48.gtf" "inputfiles/sgRNAFiles" 21 >> outputFiles/outMock21.txt
   python sgRNArunner.py "inputfiles/mockMaterials/TFsTruncatedLong.xlsx" "inputfiles/dmel-all-chromosome-r6.48.fasta" "inputfiles/dmel-all-r6.48.gtf" "inputfiles/sgRNAFiles" 42 >> outputFiles/outMock42.txt
   ```

8. Run AutoTagsCRISPR for all annotated termini of the 753 Drosophila TFs.
   **_⚠️ Warning:_** This takes about 3 days to run. Do not open excel sheets from this repository while the pipeline is running. Opening excel sheets would cause the run to crash.
     
   ```bash
   # run whole pipeline
   # of note, you can specify the window (number of bp) left and right of the annotated termini in which you would like to design the sgRNA
   # the number of bp has to be divisible by 3 and can be maximum 42
   # print statements will get stored in text files in outputFiles
   python sgRNArunner.py "inputfiles/TFs.xlsx" "inputfiles/dmel-all-chromosome-r6.48.fasta" "inputfiles/dmel-all-r6.48.gtf" "inputfiles/sgRNAFiles" 21 >> outputFiles/outFull21.txt
   python sgRNArunner.py "inputfiles/TFs.xlsx" "inputfiles/dmel-all-chromosome-r6.48.fasta" "inputfiles/dmel-all-r6.48.gtf" "inputfiles/sgRNAFiles" 42 >> outputFiles/outFull42.txt
   ```

   ---
# Future project ideas
- Improve code run time using multiprocessing
- Change HA and sgRNA output into the format used for ordering the fragments
