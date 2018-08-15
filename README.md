# ChIPSeq data analysis and visualisation

## motif_location_on_sequence 
This function helps you to locate motifs on the reference fasta sequences provided by the user. 
As a prerequisite user have to provide the list of motifs and reference fasta sequence. 
User have to mention if he wants to display percentile location i.e. in case of ChIPSeq peaks when the binding site have variable lengths displaying percent position makes more sense. Or you can just plot the start of the motif on the sequence.
This function takes into account the reverse complement of the motif while finding the motif so no need to worry about the motifs on complementary strands.

Input files
```
1. fastaSequence file
>GeneA
GTGCAACGCGTATAAACCTTTGGGCCTCTCTCACGATGAATAGGGGCAGTAGCATTACACTCTGGCATAACAATAGGGCCCTGGCAATAGGGCTCTGGCAATAGGGCTCTGGCAATAGACTGTGGGAGGGCCATACTG
>GeneB
TAGAAGAAAAAAAGACAAGCAAAACAAAATCACACAACACATAATCAAAAAGCATTACACTCTGGCATAACAATAGGGCCCTGGCAATAGGGCTCTGGCAATAGGGCTCTGGCAATAGACTGTGGGAGGGCCATACTG

2. mymotifs (can be copied from clipboard)
TCAATAGGGGATT
TGAATAGGGGC
ACAATAGGGCC
GCAATAGGGCT
GCAATAGGGCT
```

Run the function as:
```
motif_location_on_sequence(fastaSequence, mymotifs,"outfile",x = "start")
```
## Output
1. outfile_GenesWithMotifs.txt
2. Motif location visualisation
![](motif_location_on_sequence.png)
