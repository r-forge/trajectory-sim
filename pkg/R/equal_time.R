## Implements the equal-time distance between a pair of trajectories

"equal.time" <- function(t1, t2, pd=euclidian) {
  
  #t1 <- matrix(c(0,0,2,2,0,2,2,0,1,2,3,4), 4)
  #t2 <- matrix(c(0,1,1,3,0,1,0,0,1,2,3,4), 4)
  
  row <- nrow(t1)
  maxdis <- 0
  for (n in 1:row) {
    dis <- sqrt((t1[n,1]-t2[n,1])*(t1[n,1]-t2[n,1])+(t1[n,2]-t2[n,2])*(t1[n,2]-t2[n,2]))
    if(dis > maxdis) {maxdis <- dis}
  }
  
  return(maxdis)

  ## t1 and t2 are a pair of trajectories, for example obtained using
	#	data(example_traj)
	#	example.traj[[1]]
	## pd is a distance metric for points on the trajectories,
	## which defaults to the euclidian distance.
	
	## The output should be the equal-time distance between t1 and t2.
	#0 # Placeholder
}

