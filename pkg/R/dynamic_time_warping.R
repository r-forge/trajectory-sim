## Compute trajectory similarities according to dynamic time warping
## TODO: Check details of implementation against original paper

"DTW" <- function(trajectories, pd=euclidian) {
	trajectory.similarity(trajectories, implementation=DTW.pairwise, pd=pd, symmetric=TRUE, diagonal=0)
}

"DTW.pairwise" <- function(T1, T2, pd=euclidian, ...) {
	traj.sim.dp(T1, T2, .DTW.step.fun, pd, ...)
}

## DP step function for DTW
".DTW.step.fun" <-  function(T1, T2, i, j, prev, pd) {
	if (length(prev) == 0) { 
		p = new("ts.dp.entry", value=0, pred=NULL, data=list(sum=0, len=0))
		pd = NULL
		pv = 0
	} else {
		# Select the predecessor with the smallest sum so far
		pm <- which.min(sapply(prev, function(p) { p@data$sum }))
		pd <- names(pm)
		p <- prev[[pm]]
		pv <- p@value
	}
	
	s <- p@data$sum + pd(T1[i,], T2[j,])^2 # Squared sum of distances so far
	l <- p@data$len + 1 # Length of path through matrix
	new("ts.dp.entry", 
		value=sqrt(s)/l,
		pred=pd,
		data=list(sum=s, len=l))
}

