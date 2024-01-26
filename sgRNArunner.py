from sgRNAutils import sgRNArunner

# If you choose a window size larger than 42pb, you will need to extend the mutation space in fmaxStopScoreML.xlsx
sgRNArunner('inputfiles/TFs.xlsx', window = 21)