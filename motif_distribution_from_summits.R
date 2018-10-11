###############
## Author : Pooja Sethiya
## Institute : Chris Lab Faculty of Health Science / University of Macau. 
## Email : yb57662@umac.mo
## Date : 08 Oct 2018
###############
#--- Input
# 1. Genome file
# 2. summits file in bed format (preferably macs output)
# 3. copy paste motifs list 
# 4. flankong base pair from summit (e.g. 100)
# 5. want to compute reverse complement for your motifs list? Work only for DNA (A,T,G,C,N) sequence
genome_fasta <- "ABC_chromosomes.fasta"
macs2_summit_file <- "ABC_summits.bed"
mymotifs = read.clipboard(header=FALSE)
# SYGGRG	CTGGAG				    
# SYGGRG	CTGGGG
# SYGGRG	CCGGAG
# SYGGRG	CCGGGG
# SYGGRG	GTGGAG
# SYGGRG	GTGGGG
# SYGGRG	GCGGAG
# SYGGRG	GCGGGG
# GGCSS	GGCCC
# GGCSS	GGCCG
# GGCSS	GGCGC
# GGCSS	GGCGG

#--- run as
motif_distirbution_from_summits(macs2_summit_file,genome_fasta,mymotifs,flank_from_summit = 100, motif_revComplement="FALSE")

#-- load function first
motif_distirbution_from_summits <- function(macs2_summit_file,genome_fasta, mymotifs, flank_from_summit, motif_revComplement="FALSE"){
          
          #--- Load packages
          library(GenomicRanges)
          library(IRanges)
          library(tidyverse)
          library(seqinr)
          library(rtracklayer)
          library(BSgenome)
          library(psych)
          library(reshape2)
          
          macs2_summits <- import.bed(macs2_summit_file)
          print(head(macs2_summits))
          
          #--- get 100bp from summit
          summit_100bp <- macs2_summits+flank_from_summit
          print(head(summit_100bp))
          #---- Get the sequence
          
          dna <- readDNAStringSet(genome_fasta)
          names(dna) <- gsub(' .*', '',names(dna))
          
          #---- Check whether the region boundaries are within genome, remove if out of range
          dd = data.frame(cbind(names(dna), width(dna))) %>% mutate(Start=rep(1,length(names(dna)))) %>% dplyr::select(c("X1","Start","X2"))
          colnames(dd)=c("Chr","Start","End")    
          dd$Chr <- gsub(' .*', '',dd$Chr)
          dd = makeGRangesFromDataFrame(dd)
          flank_region_within_bound <- subsetByOverlaps(summit_100bp,dd,type = "within")
          
          message("Binding sites within genomic range: ",length(flank_region_within_bound))
          #---- Get Sequence of within range regions
          flank_seq <- getSeq(dna, flank_region_within_bound)
          names(flank_seq) = flank_region_within_bound$name
          
          #--- Compute reverse complement of the given motifs
          DNA_mymotifs <- DNAStringSet(as.matrix(mymotifs))
          
          if(motif_revComplement=="TRUE"){
                    revComplement <- reverseComplement(DNA_mymotifs)
                    revComplement <- data.frame(revComplement)
                    
                    #--- Combine the motifs and their reverse complements
                    all_motifs <- rbind.DataFrame(mymotifs$V1,revComplement$revComplement)
                    all_motifs <- as.matrix(all_motifs$X)   
          }
        
          else{
                    all_motifs <- as.matrix(mymotifs)
          }
          
          #--- compute the location of the motifs on the given sequences
          ll <- list()
          tt  <-  list()
          for(i in seq_along(all_motifs[,1])){
                    #i=1
                    mi0 <- vmatchPattern(all_motifs[i,2], flank_seq,fixed="subject")     
                    
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

          names(ll)=all_motifs[,1]
          
          
          
          #---- Unique genes associated with the motifs
          genes_with_motifs <- as.tibble(do.call("rbind", ll)) %>% mutate(motif=rownames(.)) 
          genes_with_motifs$motif <- str_replace(genes_with_motifs$motif,"\\..*","")
          
          write_delim(genes_with_motifs, paste(basename(macs2_summit_file),"_peaks_with_motif.tab", sep=""), delim ="\t",col_names = TRUE )
          
          genes_with_motifs$gene_width = width(flank_seq[genes_with_motifs$genes])
          
          message("motif occurrences: ", nrow(genes_with_motifs))
          print(table(genes_with_motifs$motif))
          
          
          genes_with_motifs$start <- as.numeric(levels(genes_with_motifs$start))[genes_with_motifs$start]
          
          #-- plot the motifs on the binding site
         gg <- ggplot(genes_with_motifs,aes(x=start, y=genes, color=motif,shape=motif))+
                    geom_point(alpha=0.8, size=1.8)+
                    geom_vline(data = genes_with_motifs, aes(xintercept=gene_width/2),color="blue",size=2)+
                    ylab("binding sites")+
                    xlab("")+
                    theme_classic()+
                    scale_x_continuous(limits = c(0, genes_with_motifs$gene_width[1]),breaks=c(0,51,101,151,201), labels=c("-100bp", "-50bp", "summit", "50bp", "100bp"))+
                    theme(legend.position = "top",
                          axis.ticks.y = element_blank(),
                          legend.text = element_text(face="bold", colour="black", size=12,angle=0),
                          axis.text.y = element_blank(),
                          axis.text.x = element_text(face="bold", colour="black", size=12,angle=0))
          
         print(gg)
         #--- save output file
         ggsave(paste(basename(macs2_summit_file),"_motifdistribution.pdf", sep=""),device = "pdf", width=5, height = 7)
          
}
