from sgRNAutils import sgRNArunner
import sys

# If you choose a window size larger than 42pb, you will need to extend the mutation space in fmaxStopScoreML.xlsx
# the input for the function sgRNArunner should be geneFile, fastaFile, annotationFile, sgRNAFolder and window size
#sgRNArunner(sys.argv[1], sys.argv[2], sys.argv[3], sys.argv[4], int(sys.argv[5]))
sgRNArunner("inputfiles/mockMaterials/TFsTruncatedLong.xlsx", "inputfiles/dmel-all-chromosome-r6.48.fasta",
            "inputfiles/dmel-all-r6.48.gtf", "inputfiles/sgRNAFiles", 21)
