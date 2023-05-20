def homologyArmsFragments(tfsDF):
    """
    Takes in the dataframe of information about start/stop codon regions, 
    and appends with columns for the 225 bp upstream and downstream of
    this site (the homlogy arm fragments).

    input: tfsDF - as defined in previous function and on GitHub README, with Reference_Seq 1700bp either side of the start/stop.
    output: tfsDF appended with upstreamHA and downstreamHA

    """

    tfsDF["upstreamHA"] = tfsDF.Reference_Seq.str[1475:1700]
    tfsDF["downstreamHA"] = tfsDF.Reference_Seq.str[1704:1929]

    return tfsDF

#Tested and returns right fragments.