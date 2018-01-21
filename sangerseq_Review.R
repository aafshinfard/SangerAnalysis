# Library:
library("sangerseqR")

# importing files into R: read.abif()
# reads an abif file and create an abif object
# examples : 
abif_object1 <- read.abif( system.file("extdata", "heterozygous.ab1", package = "sangerseqR") )     # a seq from the library itself
abif_object2 <- read.abif( file.choose() )                                                          # from file

# try checking usefull fields of an abif file, for example:
abif_object1@data$DATA.1

# converting abif object into sangerseq object: sangerseq()
sanger_obj1 <- sangerseq( abif_object1 )

# try checking all fields of a sangerseq object, for example:
sanger_obj1@primarySeq
primarySeq(sanger_obj1,string = TRUE)
# and others.
PolyPeakParser()
# Basecalling for and sangerseq object :
basecalled <- makeBaseCalls( sanger_obj1 )

# chromatogram: 
chromatogram( basecalled )
