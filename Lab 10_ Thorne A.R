##### Lab 10 #####
## Andrew Thorne ##
# Due: 11/16/21


#A. 
##### input 

x_str <- "ATAC"   # side sequence 

y_str <- "GTGTAC" # top sequence 

match_score <- 3 

mismatch_score <- -1 

gap_penalty <- -4 


dna.letters<-c("A","C","G","T") 
num.letters <- length(dna.letters) 
S<-data.frame(matrix(0,nrow=num.letters,ncol=num.letters))  # data frame 
rownames(S)<-dna.letters; colnames(S)<-dna.letters 
for (i in 1:4){ 
  for (j in 1:4){ 
    if(dna.letters[i]==dna.letters[j]){  
      
      S[i,j]<- match_score  
    } 
    else{  
      S[i,j]<- mismatch_score
    } 
  } 
} 

#1.
# This code displays Substitution (S) matrix, which, as stated in the lab,
# Does not involve the recursion algorithm or gap penalty. 
# This only contains the match and mismatch score for each pair of letters.
# The if statements shows the match, and if not matching, the mismatch score. 
#2.
# Size of S matrix is 4X4 Matrix.
#3.
# The if statement checks if the i and j in dna.letters are the same.
# If there are, you receive the match score, if not, the mismatch score is produced.
# The blanks are "match_score" and "mismatch_score" 

#4.
S
S["A","T"]


#B.
x <- unlist(strsplit(x_str, "")) 

y <- unlist(strsplit(y_str, "")) 

x.len <- length(x)  

y.len <- length(y)  



Fmat<-matrix(0,nrow=x.len+1,ncol=y.len+1)  

Tmat<-Fmat  #0's to start




rownames(Fmat)<-c("-",x); colnames(Fmat)<-c("-",y) 

rownames(Tmat)<-c("-",x); colnames(Tmat)<-c("-",y) 



# create first row and column 

Fmat[,1]<- seq(from=0,len=x.len+1,by=-abs(gap_penalty))   

Fmat[1,]<- seq(from=0,len=y.len+1,by=-abs(gap_penalty)) 

Tmat[,1]<- rep(2,x.len+1)  # 2 means align with a gap in the upper seq 

Tmat[1,]<- rep(3,y.len+1)  # 3 means align with a gap in the side seq 

Tmat[,1]
Tmat[1,]

Tmat
Fmat
#1. This code creates two matrices, one is the Fmat 
# Matrix with the "Fmat[,1]" representing the first column
# And the "Fmat[1,] representing the first row.
# This principle applies to the Tmat matrix as well.
#2. The size of Fmat and Tmat is larger than the size of x and y by 1 because
# They must equal to the array extent. If the values were given as they were,
# The code would not run, and the observation of the matrix
# Could not even be constructed. The Smith-Waterman algorithm
# Only works if the matrix between the two proteins are of size m+1 and n+1,
# Also, no values in the scoring matrix can be negative, so S>=0
#3.
# The first rows and columns of Fmat are given what they are given because
# The code states the starting point of 0, to a length of x.len+1, with the gap penalty
# Being what was stated first in the lab, -abs(-4). This is the same 
# Example that was shown in the pairwise_alignment_dynamic_programming powerpoint lecture.
# The y length is 6, so that means there will be a sequence of 7 starting from 0, decreasing by
# -4 Each time. Same principle applies to the x.len of 4. 

#C.
my.numbers <- c(7,9,-4) 

max(my.numbers) 

which.max(my.numbers) 

#1.
#Max identifies which number is the greatest 
# In the vector while which.max identifies the position in the vector
# Of the Max. 

#2.
my.numbers <- c(9,9,-4)
which.max(my.numbers)
#Which.max chooses the first position in the vector that contains the same value
# As another spot within the vector, for example, in this illustration, 
# Which.max chooses the 1st position over the second because it comes first in the vector.


# Use the recursion formula with the max function 
#3.

for (i in 2:nrow(Fmat)){ 
  
  for (j in 2:ncol(Fmat)){    # use F recursive rules  
    
    test_three_cases <- c(Fmat[i-1,j-1]+S[x[i-1],y[j-1]],   # 1 mis/match 
                                                  
                                                  Fmat[i-1,j]+gap_penalty,    # 2 up-gap 
                                                  
                                                  Fmat[i,j-1]+gap_penalty )  # 3 left-gap 
    
    Fmat[i,j]=max(test_three_cases) 
    
    Tmat[i,j]=which.max(test_three_cases) 
    
  } 
  
} 

final_score <- Fmat[nrow(Fmat),ncol(Fmat)]
final_score

#4.
Fmat
Tmat
final_score
# This Matches Lecture
#5.
# Tmat would be - - A T A C 
#               G T G T A C 


#D.
make.alignment.matrices <- function(x_str, y_str, match_score, mismatch_score, gap_penalty){ 
  dna.letters<-c("A","C","G","T") 
  num.letters <- length(dna.letters) 
  S<-data.frame(matrix(0,nrow=num.letters,ncol=num.letters))  # data frame 
  rownames(S)<-dna.letters; colnames(S)<-dna.letters 
  for (i in 1:4){ 
    for (j in 1:4){ 
      if(dna.letters[i]==dna.letters[j]){  
        
        S[i,j]<- match_score  
      } 
      else{  
        S[i,j]<- mismatch_score
      } 
    } 
  } 
  x <- unlist(strsplit(x_str, "")) 
  
  y <- unlist(strsplit(y_str, "")) 
  
  x.len <- length(x)  
  
  y.len <- length(y)  
  
  Fmat<-matrix(0,nrow=x.len+1,ncol=y.len+1)  
  
  Tmat<-Fmat  
  
  rownames(Fmat)<-c("-",x); colnames(Fmat)<-c("-",y) 
  
  rownames(Tmat)<-c("-",x); colnames(Tmat)<-c("-",y)
  
  Fmat[,1]<- seq(from=0,len=x.len+1,by=-abs(gap_penalty))   
  
  Fmat[1,]<- seq(from=0,len=y.len+1,by=-abs(gap_penalty)) 
  
  Tmat[,1]<- rep(2,x.len+1)
   
  Tmat[1,]<- rep(3,y.len+1)
  for (i in 2:nrow(Fmat)){ 
    
    for (j in 2:ncol(Fmat)){    # use F recursive rules  
      
      test_three_cases <- c(Fmat[i-1,j-1]+S[x[i-1],y[j-1]],   # 1 mis/match 
                            
                            Fmat[i-1,j]+gap_penalty,    # 2 up-gap 
                            
                            Fmat[i,j-1]+gap_penalty )  # 3 left-gap 
      
      Fmat[i,j]=max(test_three_cases) 
      
      Tmat[i,j]=which.max(test_three_cases) 
      
    } 
    
  } 
   
  final_score <- Fmat[nrow(Fmat),ncol(Fmat)] 
  
  return(list(Fmat=Fmat, Tmat=Tmat, score_out=final_score)) 
  
} 


# load new input 

x_str <- "GATTA"  # side sequence 

y_str <- "GAATTC" # top sequence 

match_score <- 2 

mismatch_score <- -1 

gap_penalty <- -2 



align.list <- make.alignment.matrices(x_str, y_str, match_score, 
                                      
                                      mismatch_score, gap_penalty) 

align.list$Fmat 

align.list$Tmat 

align.list$score_out 


#1.
Fmat
Tmat
final_score

#2.
#install.packages("gplots") 

library(gplots) 

Fmat <- align.list$Fmat 

col = c("black","blue","red","yellow","green") 

breaks = seq(min(Fmat),max(Fmat),len=length(col)+1) 



heatmap.2(Fmat[-1,-1], dendrogram='none', density.info="none",  
          
          Rowv=FALSE, Colv=FALSE, trace='none',  
          
          breaks = breaks, col = col, 
          
          sepwidth=c(0.01,0.01), 
          
          sepcolor="black", 
          
          colsep=1:ncol(Fmat), 
          
          rowsep=1:nrow(Fmat)) 

#E.
show.alignment <- function(x_str,y_str,Tmat){  
  
  ################ create the alignment 
  
  # input Tmat and the two sequences: x side seq and y is top seq 
  
  # make character vectors out of the strings 
  
  x<-unlist(strsplit(x_str,"")) 
  
  y<-unlist(strsplit(y_str,"")) 
  
  
  
  n<-nrow(Tmat) # start at bottom right of Tmat 
  
  m<-ncol(Tmat) 
  
  alignment<-character() 
  
  while( (n+m)!=2 ){ 
    
    if (Tmat[n,m]==1){ 
      
      # subtract 1 from x and y indices because they are 
      
      # one row/col smaller than Tmat 
      
      curr_align_col <- rbind(x[n-1],y[m-1])  
      
      alignment <- cbind(curr_align_col,alignment) 
      
      n=n-1; m=m-1; # move back diagonally 
      
    }else if(Tmat[n,m]==2){  
      
      curr_align_col <- rbind(x[n-1],"-") # put gap in top seq 
      
      alignment <- cbind(curr_align_col,alignment) 
      
      n=n-1 # move up 
      
    }else{ 
      
      curr_align_col <- rbind("-",y[m-1]) # put gap in side seq 
      
      alignment <- cbind(curr_align_col,alignment) 
      
      m=m-1 # move left 
      
    }            
    
  } # end while 
  
  return(alignment) 
  
} # end function 

#1.
alignment <- show.alignment(x_str,y_str,align.list$Tmat) 

alignment <- show.alignment(x_str,y_str,Tmat) 

alignment 

# the following prints nicely 

write.table(alignment,row.names=F,col.names=F,quote=F) 

#2.
x_str <- "ATCGT"  # side sequence 

y_str <- "TGGTG" # top sequence 

match_score <- 1 

mismatch_score <- -2 

gap_penalty <- -1  

align.list <- make.alignment.matrices(x_str, y_str, match_score, 
                                      
                                      mismatch_score, gap_penalty) 
align.list$Tmat
align.list$Fmat
align.list$score_out

alignment <- show.alignment(x_str,y_str,align.list$Tmat) 

alignment <- show.alignment(x_str,y_str,Tmat) 

alignment 

Fmat <- align.list$Fmat 

col = c("black","blue","red","yellow","green") 

breaks = seq(min(Fmat),max(Fmat),len=length(col)+1) 

heatmap.2(Fmat[-1,-1], dendrogram='none', density.info="none",  
          
          Rowv=FALSE, Colv=FALSE, trace='none',  
          
          breaks = breaks, col = col, 
          
          sepwidth=c(0.01,0.01), 
          
          sepcolor="black", 
          
          colsep=1:ncol(Fmat), 
          
          rowsep=1:nrow(Fmat)) 
