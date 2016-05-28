#do the same thing 3 ways but with huge differences in speed
mat <- matrix(c(rep(1:5),rep(6:10),rep(11:15),rep(16:20)),ncol=4)
mat <- matrix(runif(1e4),nrow=100)

mapplyTst <- function(mat){
  iterator <- expand.grid(1:nrow(mat),1:ncol(mat))
  new_mat <- mapply(function(i,j,mat){mat[i,j]+mat[i,j]},i=iterator$Var1,j=iterator$Var2,MoreArgs=list(mat=mat),SIMPLIFY=TRUE)
  new_mat <- matrix(data=new_mat,nrow=nrow(mat),ncol=ncol(mat))
  return(new_mat)
}

loopTst <- function(mat){
  new_mat <- matrix(NA,nrow=nrow(mat),ncol=ncol(mat))
  for(i in 1:nrow(mat)){
    for(j in 1:ncol(mat)){
      new_mat[i,j] <- mat[i,j]+mat[i,j]
    }
  }
  return(new_mat)
}

outerHelp <- function(i,j,mat){
  return(mat[i,j]+mat[i,j])
}

outerHelpV <- Vectorize(outerHelp,c("i","j"))

outerTst <- function(mat){
  outer(1:nrow(mat),1:ncol(mat),outerHelpV,mat=mat)
}

mapplyTst(mat)
loopTst(mat)
outerTst(mat)

library(microbenchmark)
microbenchmark(mapplyTst(mat),loopTst(mat),outerTst(mat))
