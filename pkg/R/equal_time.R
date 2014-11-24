## Implements the equal-time distance between a pair of trajectories




"equal.time" <- function(t1, t2, pd=euclidian) {
  
  #t1 <- matrix(c(0,0,2,2,0,2,2,0,1,2,3,4), 4)
  #t2 <- matrix(c(0,1,1,3,0,1,0,0,1,2,3,4), 4)  
  
  t1 <- matrix(c(0,1,5,0,2,3,1,2,4), 3)
  t2 <- matrix(c(0,3,5,0,1,2,1,3,4), 3)
  
  LIt1 <- equal.time.LI(t1,t2)
  LIt2 <- equal.time.LI(t2,t1)
  
  maxdis1 <- equal.time.DisSimul(LIt1,t2,pd)
  maxdis2 <- equal.time.DisSimul(t1,LIt2,pd)
  
  maxdis <- max(maxdis1,maxdis2)
  
  maxdis
  
  return(maxdis)
}


"equal.time.LI" <- function(t1, t2) {
  #t1 <- matrix(c(0,1,5,0,2,3,1,2,4), 3)
  #t2 <- matrix(c(0,3,5,0,1,2,1,3,4), 3)

  
  row1=nrow(t1)
  col=ncol(t1)
  
  time1=t1[,col]
  
  fl <- vector("list", col -1)   ## fl: function list

  for (n in 1:(col-1)){
    fl[[n]]<- approxfun(time1, t1[,n])   ## for each dimension, we generate one approximation function
  }
  
  time2 <- t2[,col]
  row2=nrow(t2)
  LIt1 <- matrix(data=NA, nrow=row2, ncol= col)
  for (n in 1:(col-1)){
    LIt1[,n]<- fl[[n]](time2)
  }
  LIt1[,col] <- time2
  
  LIt1
  
  
  return (LIt1)
}

##the input t1 and t2 have to be recorded Simultaneously, which also means t1 and t2 have the same numbers of records
"equal.time.DisSimul" <- function(t1, t2, pd=euclidian) {
  t1 <- matrix(c(0,1,5,0,2,3,1,2,4), 3)
  t2 <- matrix(c(0,3,5,0,1,2,1,3,4), 3)
  
  t1
  t2
  
  row <- nrow(t1)
  
  maxdis <- 0
  for (n in 1:row) {
    dis <- pd(t1[n], t2[n])     
    if(dis > maxdis) {
      maxdis <- max(maxdis,dis)
    }
  }
  
  maxdis
  
  return (maxdis)
}

