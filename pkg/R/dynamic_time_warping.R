## Compute trajectory similarities according to dynamic time warping
## Follows definition in:
## Chen, Lei, M. Tamer Ozsu, and Vincent Oria. 
## Robust and fast similarity search for moving object trajectories. 
## Proc. 2005 ACM SIGMOD intern. conf. Management of data.

"DTW" <- function(trajectories, pd=euclidian) {
	trajectory.similarity(trajectories, implementation=DTW.pairwise, pd=pd, symmetric=TRUE, diagonal=0)
}

"DTW.pairwise" <- function(T1, T2, pd=euclidian, ...) {
	traj.sim.dp(T1, T2, .DTW.step.fun, pd=pd, ...)
}

## DP step function for DTW
".DTW.step.fun" <-  function(T1, T2, i, j, prev, pd) {
	if (length(prev) == 0) {
		p = new("ts.dp.entry", value=0, pred=NULL)
		pred = NULL
	} else {
		# Select the predecessor with the smallest sum so far
		pm <- which.min(sapply(prev, slot, "value"))
		pred <- names(pm)
		p <- prev[[pm]]
	}
	
	new("ts.dp.entry", 
		value=p@value + pd(T1[i,], T2[j,]),
		pred=pred)
}

