library("sangerseqR")
require("BSgenome.Hsapiens.UCSC.hg19");
require("TxDb.Hsapiens.UCSC.hg19.knownGene");
require("org.Hs.eg.db");

GeneSymb="CYP21A2";                                          # <==#####

### start and end for trimming
TRIM_5 = matrix( 0,nrow=Exon_cnt ) # bases omit from start   # Automated
TRIM_3 = matrix( 0,nrow=Exon_cnt ); # bases omit from end    # Automated

# Other Parameters
Heteratio = 0.33;                                            # <==#####
minQuality = 100;                                            # <==#####
cnt = 2;
for ( i in 1:cnt ){
  ab1_Ex_list[i]  <- read.abif( file.choose() )
  sang_Ex_list[i] <- sangerseq( ab1_Ex_list[[i]] )
  #sang_Ex_list[i] <- makeBaseCalls( sang_Ex_list[[i]] , ratio = Heteratio )
}

for ( i in 1:cnt ){
  #ab1_Ex_list[i]  <- read.abif( file.choose() )
  #sang_Ex_list[i] <- sangerseq( ab1_Ex_list[[i]] )
  sang_Ex_list[i] <- makeBaseCalls( sang_Ex_list[[i]] , ratio = Heteratio )
}

# Updates a sangerseq class object to contain primary and secondary peak calls

chromatogram( sang_Ex_list[[1]] , showhets = TRUE )
chromatogram( sang_Ex_list[[2]] , showhets = TRUE )
chromatogram( sang_Ex_list[[3]] , showhets = TRUE )
chromatogram( sang_Ex_list[[4]] , showhets = TRUE )
chromatogram( sang_Ex_list[[5]] , showhets = TRUE )

sanger2 = sang_Ex_list[[2]]
sanger2p = sanger2
sanger2p@primarySeq = sanger2p@primarySeq[102:200]
sanger2p@secondarySeq = sanger2p@secondarySeq[102:200]
sanger2p@peakPosMatrix = sanger2p@peakPosMatrix[102:200,]
sanger2p@peakAmpMatrix = sanger2p@peakAmpMatrix[102:200,]
chromatogram( showcalls ="both",sanger2p, showhets = TRUE , showtrim = TRUE ,width = 50)


a = sanger2p@peakAmpMatrix
a_1sh = rbind( a[ 2:dim(a)[1], ] , t(as.matrix(c(0,0,0,0))) )
a_2sh = rbind( a_1sh[ 2:dim(a_1sh)[1], ] , t(as.matrix(c(0,0,0,0))) )
a_3sh = rbind( a_2sh[ 2:dim(a_2sh)[1], ] , t(as.matrix(c(0,0,0,0))) )

a1 = sang_Ex_list[[1]]@traceMatrix
a2 = sang_Ex_list[[1]]@peakPosMatrix
a3 = sang_Ex_list[[1]]@peakAmpMatrix

write.table(a1,file = "Desktop/tracematrix")
write.table(a2,file = "Desktop/peakPosMatrix")
write.table(a3,file = "Desktop/peakAmpMAtrix")


indx = as.character(as.vector(sanger2p@primarySeq))
A_indx = (indx == "A") | (indx == "a")
C_indx = (indx == "C") | (indx == "c")
G_indx = (indx == "G") | (indx == "g")
T_indx = (indx == "T") | (indx == "t")
b = sanger2p@peakAmpMatrix
b[ indx == , 1]
