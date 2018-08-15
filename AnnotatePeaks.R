###############
## Author : Pooja Sethiya
## Institute : Chris Lab Faculty of Health Science / University of Macau. 
## Email : yb57662@umac.mo
## Date : 02-04-2018
###############
# subject.bed = Path/To/Bed/File
######## subject file format##########
# Chr1	3275	4063	Gene1	0	+
# Chr2	4131	4571	Gene2	0	+
# Chr2	4571	5524	Gene3	0	+
# Chr2	9929	10738	Gene4	0	-
# chr3	12469	12699	Gene5	0	-
############################
# query <- "path/to/Query/file"
########## Query input file format######
# Chr	Start	End	PeakID
# Chr1	1791855	1792151	Peak-1
# Chr1	57460	57756	Peak-2
# Chr2	892611	892907	Peak-3
# chr3	1297088	1297384	Peak-4
# chr4	1280919	1281215	Peak-5
###### Run function as
#AnnotatePeaksInR(subject.bed,query,"out")
####################

subject= "/Users/Pooja/Documents/Data-Analysis/PhD/Hsf1/ChIP-Seq/Calbicans/CaHsf1_ChIP_Set1/Reference/Ca_orf_6221_forGRanges.bed"
#subject="/Users/Pooja/Documents/Data-Analysis/PhD/Pol2-RNASeq/CG-RNASeq-PolI/Cg-Pol2-All-Time-RNASeq2hr-rawdata/Genome-Reference-Files/Cg_s02_m07_r06_Orf_withOutHeader.bed"
subject= rtracklayer::import(subject,format = "bed")
query = "Glucose_42C_20min_B_narrow_TotalPeaks.tab"
query.input <- read.table(query,sep="\t", header=T)
outfileName <- "Glucose_42C_20min_B_ReSeq_Combined_TotalPeaks"

AnnotatePeaksInR=function(subject,query.input,outfileName,MatchMotif="FALSE"){
          library("IRanges")
          library("GenomicRanges")
          library("rtracklayer")
          library(dplyr)
          
          ## Read Input files
          query.input <- read.table(query.input,sep="\t", header=T)
          query <- makeGRangesFromDataFrame(query.input, keep.extra.columns = T)
          subject= rtracklayer::import(subject,format = "bed")
          #### preceeding gene to query
          pp = GenomicRanges::precede(query,subject,ignore.strand=F,select="all")
          df.p = data.frame(query[queryHits(pp),], subject[subjectHits(pp),])
          nrow(df.p)
          df.p <- as_tibble(df.p)
          
          ####  CDS to query All (i.e. include complete and partial overlap)
          
          oo.1 = findOverlaps(query,subject,ignore.strand=F)
          df.o.1 = data.frame(query[queryHits(oo.1),], subject[subjectHits(oo.1),])
          nrow(df.o.1)
          
          ####overlapping CDS to query within gene body
          
          oo.2 = findOverlaps(query,subject,ignore.strand=T,select="all",type="within")
          df.o.2 = data.frame(query[queryHits(oo.2),], subject[subjectHits(oo.2),])
          
          message("Number of peaks with complete overlap with the coding region")
          print(nrow(df.o.2))
          
          ##########################
          
          partial.overlap <- subset(df.o.1, !(df.o.1[,6] %in% df.o.2[,6]))
          message("Number of peaks with Partial Overlap")
          print(nrow(partial.overlap))
          
          
          ### Filter Peak to get nearest gene
          ### Assign two genes to the partial overlapping peaks
          ### Filter targets having partial overlap at 3'end
          ### Depending upon the distance denote the peaks as PartialOverlap_Gene_1,2,3
          
          partial.overlap <- as_tibble(partial.overlap)
          pp <-  partial.overlap %>% 
                    mutate(distance =if_else(strand.1=="-", start-end.1, start.1-end)) %>% 
                    filter(abs(distance)< width) %>%
                    dplyr::group_by(V4)  %>% 
                    dplyr::mutate(min = order(distance)) %>% 
                    mutate(Category = paste("PartialOverlap_Gene",min , sep= "_")) %>% dplyr::select(-c(min))
                    
          ### Find the peaks within gene body
          ### Assign gene terminal if the ratio of distance > 0.75 i.e. if the peak is lying at extreme 3' end
          df.o.2 = as_tibble(df.o.2)
          pGeneBody <- df.o.2 %>% 
                    mutate(distance =if_else(strand.1=="-",  end.1-end,start-start.1), Ratio = (abs(distance)/width.1), Category=if_else((abs(distance)/width.1)>0.75,"AtGeneEnd","WithInGeneBody")) %>% 
                    dplyr::select(-c(Ratio))
          
          ### Find peaks in promoter region
          ### filter the peaks based on distance
          pPrecede <- mutate(df.p, distance =if_else(strand.1=="-", start-end.1, start.1-end),Category="InPromoter")  %>% 
                      group_by(V4) %>% 
                      filter(distance==min(distance))  
           pPromoter <- pPrecede
          ### Remove peaks from  Promoter if they are already showing partial overlap       
          pPrecede <- anti_join(pPrecede, pp, by = "V4") 
          
          ### Remove peaks from  Promoter if they are already showing overlap with GeneBody   
          pPrecede <- anti_join(pPrecede, pGeneBody, by = "V4")
          
          
          ### All possible targets for each peak          
          df.All <- bind_rows(pp,pGeneBody, pPrecede)    
          
          ### Peaks with distance > 2000 in promoter, assign them to nearest gene irrespective of strand (ie --> Peak <-- <-- -->Target(according To Strand))
          df.LongPromoters <- df.All %>% mutate(Condition=if_else((distance>2000 & Category=="InPromoter"), "1","0")) %>% filter(Condition=="1")  
          
          qq.subsetLongPromoters <-  subset(query, query$V4 %in% df.LongPromoters$V4)
           
          nn <- nearest(qq.subsetLongPromoters,subject,ignore.strand=T,select="all")
          df.near <- data.frame(qq.subsetLongPromoters[queryHits(nn),], subject[subjectHits(nn),])
          nrow(df.near)
          df.near <-  as_tibble(df.near)
          df.near.1 <- mutate(df.near, distance =if_else(start > end.1, start-end.1, start.1-end),Category="NearestTarget") 
          
          ### Filter peaks in Promoter if distance > 2000bp from the all targets
          df.All <- anti_join(df.All, df.near.1, by="V4")
          
          df.All <- bind_rows(df.All, df.near.1)
          
         
          ### Filter the peaks AtGeneEnd if they have promoter of next gene nearby
          df.filterGeneTerminal <- bind_rows(pGeneBody, pPromoter)  %>% filter(Category=="InPromoter"|Category=="AtGeneEnd") %>% group_by(V4) %>% filter(n()>1) %>% filter(distance==min(distance))
          
          ### All unique targets
          pPrecede <- anti_join( pPrecede,df.near.1, by="V4")
          pGeneBody <- anti_join( pGeneBody,df.filterGeneTerminal, by="V4")
          df.Unique <- bind_rows(df.filterGeneTerminal,pPrecede, pp, pGeneBody, df.near.1)
          df.Unique <- unique(df.Unique)
          #df.Unique=df.Unique[,c(1,2,3,6,7,8,9,10,11,12,13,14,15,16,19,20,21)]
          
          message("Number of Unique peaks")
          print(nrow(df.Unique))
          #colnames(df.out1) <- c("Chr","PeakStart","PeakEnd","PeakWidth","PeakStrand", "PeakId","Chr","GeneStart","GeneEnd","GeneLength","GeneStrand","GeneName", "Score","Distance","Catergory")
          
          if(MatchMotif=="TRUE"){
                    ### Match the motif with the peaks
                    hsf1_motif <- rtracklayer::import("/Users/Pooja/Documents/Data-Analysis/PhD/Hsf1/ChIP-Seq/Calbicans/CaHsf1_ChIP_Set1/Reference/Ca_GenomeWideMotifMapping_Hsf1_CombinedMotifs.bed",format = "bed")
                    hsf1.oo = findOverlaps(hsf1_motif,query,ignore.strand=T,select="all",type="within")
                    Hsf1.df = data.frame( query[subjectHits(hsf1.oo),],hsf1_motif[queryHits(hsf1.oo),])
                   # Hsf1.df = Hsf1.df[,c(1,2,3,6,8,13,14,15,18)]
                    
                    clip <- pipe("pbcopy", "w")
                    write.table(Hsf1.df, file=clip,sep="\t",quote=F,row.names = F)
                    close(clip)
                    
                    write.xlsx(x=as.data.frame(Hsf1.df), file = paste(outfileName,"AnnotatedPeaksbyR.xlsx", sep="_"), sheetName = "Hsf1Motifs", append = TRUE) 
                    
          }
          #write.xlsx(x=as.data.frame(df.All), file = paste(outfileName,"AnnotatedPeaksbyR.xlsx", sep="_"),sheetName = "AllPossibleTargets",append = TRUE)    
          write.xlsx(x=as.data.frame(df.Unique), file = paste(outfileName,"AnnotatedPeaksbyR.xlsx", sep="_"), sheetName = "UniqueTargets", append = TRUE) 
          
          tt = df.Unique %>% mutate(Category=replace(Category,str_detect(Category, "PartialOverlap_*"), "PartialOverlap"))
          tt = table(tt$Category)
          ttm=melt(tt)
          
          
         
          
          # ######### For tRNA overlap
          # 
          # query.tRNA <- rtracklayer::import("../../tRNA_A_nidulans_FGSC_A4_version_s10-m04-r07.bed",format = "bed")
          # 
          # tRNA.1 = findOverlaps(query.tRNA,query,ignore.strand=F)
          # tRNA.o.1 = data.frame( query[subjectHits(tRNA.1),],query.tRNA[queryHits(tRNA.1),])
          # nrow(tRNA.o.1)
          # tRNA.o.1=tRNA.o.1[,c(1,2,3,6,7,8,9,10,13)]
          # write.xlsx(x=tRNA.o.1, file = paste(outfileName,"AnnotatedPeaksbyR.xlsx", sep="_"), sheetName = "tRNA_Targets", append = TRUE) 
          # 
          ##### Plot the distribution of targets 
          gg = ggplot(ttm, aes("",y=value, fill=Var1))+geom_bar(width = 1, stat = "identity")+geom_text(aes(label = value), position = position_stack(vjust = 0.5),size=5)+coord_polar(theta = "y")+ylab("No. of Peaks")+xlab("")+theme_classic()
           pdf(paste(outfileName,"AnnotatedPeaksbyR.pdf", sep="_"), height = 5, width = 10)
           print(gg)
           dev.off()
          gg
          
          
          
}

AnnotatePeaksInR(subject, query.input,outfileName,MatchMotif = "TRUE")

          
