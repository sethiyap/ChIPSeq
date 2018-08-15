###############
## Author : Pooja Sethiya
## Institute : Chris Lab Faculty of Health Science / University of Macau. 
## Email : yb57662@umac.mo
## Date : 15 Aug 2018
###############
#--- This function also considers the reverse complement of the motif while finding the motif
#--- 1. Provide fastaSequence in which the motif is to be searched
#--- 2. List of motifs or combination of motifs 
#         AAAATAGGGAT
#         TCAATAGGGGA
#         TAAATAGGGAT
#--- 3. Provide output file name, this will contain the location of motif and which motif is found in which gene
#--- 4. Provide x-axis to be plotted i.e. "percent" position or just the "start" of the motif on the reference sequence
#--- Run the function as: motif_location_on_sequence(fastaSequence, mymotifs,"PPP_Msn4",y = "start")

#-------------------------- INPUT ------------------------------------

fastaSequence = readDNAStringSet("PPP_Promoter_500bpFromATG.fa",format = "fasta")
mymotifs = read.clipboard(header=FALSE)

#------------------------- FUNCTION ----------------------------------
motif_location_on_sequence = function(fastaSequence, mymotifs, output_name,x="Percent"){

require("tidyverse")
require("ggpubr")
require("Biostrings")

DNA_mymotifs <- DNAStringSet(as.matrix(mymotifs))

#--- Compute reverse complement of the given motifs
revComplement <- reverseComplement(DNA_mymotifs)

revComplement <- as.data.frame(revComplement)

#--- Combine the motifs and their reverse complements
all_motifs <- cbind(mymotifs$V1, revComplement$x)
all_motifs <- as.matrix(all_motifs[,2])

#--- compute the location of the motifs on the given sequences
ll <- list()
tt  <-  list()
for(i in seq_along(all_motifs)){
          #i=1
          
          mi0 <- vmatchPattern(all_motifs[i], fastaSequence)     
          
          coords = as.data.frame(mi0)
          
          nmatch_per_seq <- elementNROWS(mi0)
          
          pos = which(nmatch_per_seq>0)
          
          tt[[i]] = table(nmatch_per_seq)
          
          Freq = nmatch_per_seq[coords$group]
          
          genes = names(mi0)[coords$group]
          
          start=coords[,3]
          end=coords[,4]
          
          ll[[i]] = as.data.frame(cbind(genes,Freq,start,end))
          #print(ll[[i]])
}


names(ll)=mymotifs[,1]
names(tt)=mymotifs[,1]


#---- Unique genes associated with the motifs
genes_with_motifs <- as.tibble(do.call("rbind", ll))
genes_with_motifs$gene_width = width(fastaSequence[genes_with_motifs$genes])

message("Total Genes with motif ", length(unique(genes_with_motifs$genes)))

write.table(genes_with_motifs, paste(output_name,"_GenesWithMotifs.txt", sep=""), sep="\t",col.names = TRUE, row.names = TRUE, quote = FALSE)

#---- If you already have motifs mapped to the sequences in the following format upload them on clipboard:
#  motif_data <-  read_delim(pipe("pbpaste"), delim = "\t", col_names =TRUE)

motif_data <- genes_with_motifs  %>% dplyr::select(genes, start,gene_width)

# genes                                                   start gene_width
#CAGL0I02200g::ChrI_C_glabrata_CBS138:188657-189157(+)   139          500
#CAGL0I02200g::ChrI_C_glabrata_CBS138:188657-189157(+)   312          500
#CAGL0L03740g::ChrL_C_glabrata_CBS138:435271-435771(-)   143          500
#CAGL0L05478g::ChrL_C_glabrata_CBS138:601708-602208(+)   276          500

motif_data$start <- as.numeric(levels(motif_data$start))[motif_data$start]

#-- Compute percent location for the fasta sequences with variable lengths
motif_data <- motif_data %>% dplyr::mutate(Percent=100*(start/gene_width))

motif_data <- motif_data[order(motif_data$start, decreasing = TRUE),]

ggdotchart(motif_data,x="genes", y=x,
           color = "genes",                                # Color by groups
           #palette = rainbow(nrow(motif_data)) , # Custom color palette
           sorting = "descending",                       # Sort value in descending order
           rotate = FALSE,                                # Rotate vertically
           dot.size = 4,                                 # Large dot size
           y.text.col = FALSE,                            # Color y text by groups
           ggtheme = theme_pubr())+
          scale_y_reverse()+ 
          ylab("")+
          coord_flip()+
          theme(legend.text = element_blank(),legend.position = "none",
                axis.text.x = element_text(face="bold", colour="black", size=16,angle=0),
                axis.text.y = element_text(face="bold", colour="black", size=12,angle=0))+
          guides(fill = guide_legend(title = ""),color = guide_legend(title = ""))+
          theme_cleveland()

}
