a = array(c(0,1,0,0,   0,0,0,1,   0,0,1,0,    1,0,0,0,    0,0,1,1,    1,0,0,1,   0,1,1,0,  0,0,0,1,  1,1,0,0,  0,1,0,1, 1,0,0,0,   0,1,1,0,  1,0,0,0,  1,0,0,1) , dim= c(4,15))
a = array(c(0,1,0,0, 0,0,0,1, 0,0,1,0, 1,0,0,0, 0,0,1,1, 1,0,0,1, 0,1,1,0, 0,0,0,1, 1,1,0,0, 0,1,0,1, 1,0,0,0, 0,1,1,0,  1,0,0,0, 0,0,1,0, 1,0,0,1) , dim= c(4,15))
cc =array( c(0,0,0,1, 0,1,0,0, 1,0,1,0, 1,1,0,0, 0,0,0,1, 0,1,1,0 ),dim=c(4,6))

HShift <- function( Mat, n ){
  
  if(abs(n) >= dim(Mat)[2] )
    Res = array(0 , dim = dim(Mat) )
  else if( n >= 0)
    Res = cbind( array(0 , dim = c( dim(Mat)[1] , abs(n) ) ) , t(head(t(Mat),dim(Mat)[2]-abs(n)))  )
  else {
    Res = cbind(  t(tail(t(Mat),dim(Mat)[2]-abs(n))), array(0 , dim = c( dim(Mat)[1] , abs(n) ) )  )
    colnames(Res) <- NULL
    rownames(Res) <- NULL
  }
  Res
}

HConvol <- function(Mat){
  Res = rep( 0 , (dim(Mat)[2]*2)-1 )
  j = 1
  for ( i in ((dim(Mat)[2]-1)*-1) : (dim(Mat)[2]-1)){
    Res[j] = sum(HShift( Mat , i ) * Mat)
    j = j+1
  }
  Res
}

HConvol14 <- function( Mat , n){
  Res = HShift( Mat , n ) * Mat
  Res
}

HConvol1 <- function( Mat , n){
  Res = colSums(HShift( Mat , n ) * Mat )
  Res
}


############################################
########### Sequence to Signal Simulator :

seq1 = c("A","C","G","C","T","A","G","A","T","C","G","C","G","G","A","C","A","T","T","G","A","G","T","C","C","G","T","A")
seq2 = c("A","C","G","C","T","A","A","G","A","T","C","G","G","A","C","A","T","T","G","A","G","T","C","C","G","T","A")

seq1 = c("A","C","G","C","T","A","G","A","T","C","G","G","A","C","A","T","T","G","A","G","T","C","C","G","T","A")
seq2 = c("A","C","G","C","T","A","A","G","A","T","C","G","C","G","G","A","C","A","T","T","G","A","G","T","C","C","G","T","A")

trace = array(rep(0,4*max(length(seq1),length(seq2))),dim = c(4,max(length(seq1),length(seq2))))
rownames(trace) = c("A","C","G","T")

for( i in 1:dim(trace)[2]){
  if(length(seq1) >= i)
    trace[seq1[i],i] = 1
  if(length(seq2) >= i)
    trace[seq2[i],i] = 1
}

